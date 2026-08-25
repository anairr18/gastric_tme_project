#!/usr/bin/env python3
"""Build auditable sample-level pseudobulk from GSE189926 public matrices.

GSE189926 is a CD45-selected gastric tumour dataset collected before and after
disulfiram plus nivolumab. This script deliberately creates sample-level
immune pseudobulk only: it cannot assess endothelial or fibroblast abundance.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import shutil
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse


ACCESSION = "GSE189926"
SOFT_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189926/soft/GSE189926_family.soft.gz"
RAW_TAR_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE189nnn/GSE189926/suppl/GSE189926_RAW.tar"
SOURCE_PAGE = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE189926"


def download(url: str, destination: Path, minimum_bytes: int) -> None:
    if destination.exists() and destination.stat().st_size >= minimum_bytes:
        print(f"[reuse] {destination}", flush=True)
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": "gastric-tme-validation/1.0"})
    print(f"[download] {url}", flush=True)
    with urllib.request.urlopen(request, timeout=600) as response, temporary.open("wb") as handle:
        shutil.copyfileobj(response, handle)
    if temporary.stat().st_size < minimum_bytes:
        temporary.unlink(missing_ok=True)
        raise ValueError(f"Unexpectedly small download: {destination}")
    temporary.replace(destination)


def parse_soft(path: Path) -> pd.DataFrame:
    records, current = [], None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current is not None:
                    records.append(current)
                current = {"gsm": line.split("=", 1)[1].strip()}
            elif current is None:
                continue
            elif line.startswith("!Sample_title = "):
                current["title"] = line.split("=", 1)[1].strip()
            elif line.startswith("!Sample_supplementary_file_1 = "):
                current["matrix_url"] = line.split("=", 1)[1].strip().replace("ftp://", "https://")
            elif line.startswith("!Sample_characteristics_ch1 = "):
                key, _, value = line.split("=", 1)[1].strip().partition(":")
                current[key.strip().lower()] = value.strip()
    if current is not None:
        records.append(current)
    metadata = pd.DataFrame(records)
    if metadata.empty or "title" not in metadata or "matrix_url" not in metadata:
        raise ValueError("GSE189926 SOFT metadata lacked sample titles or matrix URLs.")
    identity = metadata["title"].str.extract(r"^(?P<patient_id>DSF\d+)\s+(?P<timepoint>pre-treatment|post-treatment)", expand=True)
    metadata = pd.concat([metadata, identity], axis=1)
    metadata["sample_id"] = metadata["gsm"].astype(str)
    metadata["timepoint"] = metadata["timepoint"].map({"pre-treatment": "baseline", "post-treatment": "on_treatment"})
    metadata["response"] = metadata.get("outcome", pd.Series(index=metadata.index, dtype=str)).astype(str).str.upper().str.strip()
    metadata["assay_scope"] = "CD45_selected_immune_cells"
    required = ["sample_id", "patient_id", "timepoint", "response", "matrix_url"]
    if metadata[required].isna().any().any():
        raise ValueError("Could not extract complete patient/timepoint/response metadata from GSE189926 SOFT.")
    return metadata[required + ["title", "tissue", "cell type", "treatment", "assay_scope"]]


def sum_matrix(path: Path) -> pd.Series:
    """Read a dense gene-by-cell GEO TXT matrix in bounded row chunks."""
    sums: list[pd.Series] = []
    for chunk in pd.read_csv(path, sep="\t", compression="gzip", index_col=0, chunksize=500, low_memory=False):
        numeric = chunk.apply(pd.to_numeric, errors="coerce").fillna(0.0)
        sums.append(numeric.sum(axis=1))
    if not sums:
        raise ValueError(f"No rows found in {path}")
    result = pd.concat(sums)
    result.index = result.index.astype(str).str.upper()
    return result.groupby(level=0).sum()


def read_matrix_as_anndata(path: Path, sample: pd.Series) -> tuple[pd.Series, ad.AnnData]:
    """Read one GEO dense gene-by-cell matrix and retain raw integer counts."""
    table = pd.read_csv(path, sep="\t", compression="gzip", index_col=0, low_memory=False)
    table.index = table.index.astype(str).str.upper()
    table = table.apply(pd.to_numeric, errors="coerce").fillna(0)
    if table.index.has_duplicates:
        table = table.groupby(level=0, sort=False).sum()
    counts = table.sum(axis=1)
    barcodes = table.columns.astype(str)
    matrix = sparse.csr_matrix(table.to_numpy(dtype=np.int32, copy=False).T)
    obs = pd.DataFrame({
        "cell_barcode": barcodes,
        "sample_id": str(sample["sample_id"]),
        "patient_id": str(sample["patient_id"]),
        "timepoint": str(sample["timepoint"]),
        "response": str(sample["response"]),
        "tissue": "tumor",
        "assay_scope": str(sample["assay_scope"]),
        "dataset_id": ACCESSION,
    }, index=pd.Index(str(sample["sample_id"]) + "::" + barcodes, name="cell_id"))
    return counts, ad.AnnData(X=matrix, obs=obs, var=pd.DataFrame(index=table.index))


def matrix_path(root: Path, sample: pd.Series) -> Path:
    filename = Path(str(sample["matrix_url"])).name
    return root / "matrices" / filename


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("/content/drive/MyDrive/data/external/validation_extensions/GSE189926"),
    )
    parser.add_argument("--download-matrices", action="store_true")
    parser.add_argument("--write-h5ad", action="store_true")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    soft_path = args.output_dir / "GSE189926_family.soft.gz"
    download(SOFT_URL, soft_path, minimum_bytes=1_000)
    clinical = parse_soft(soft_path)
    clinical.to_csv(args.output_dir / "GSE189926_clinical_audited.csv", index=False)

    if args.download_matrices:
        for _, sample in clinical.iterrows():
            destination = matrix_path(args.output_dir, sample)
            download(str(sample["matrix_url"]), destination, minimum_bytes=1_000_000)

    missing = [str(matrix_path(args.output_dir, sample)) for _, sample in clinical.iterrows() if not matrix_path(args.output_dir, sample).exists()]
    if missing:
        (args.output_dir / "GSE189926_DOWNLOAD_STATUS.txt").write_text(
            "Clinical metadata are complete. Raw count matrices are not all present.\n"
            "Run with --download-matrices to prepare sample pseudobulk.\n"
            + "\n".join(missing) + "\n",
            encoding="utf-8",
        )
        print(f"Prepared auditable clinical table; {len(missing)} matrices remain to download.")
        return

    sums, adatas = {}, []
    for ordinal, (_, sample) in enumerate(clinical.iterrows(), start=1):
        print(f"[{ordinal}/{len(clinical)}] pseudobulk {sample['sample_id']} ({sample['title']})", flush=True)
        if args.write_h5ad:
            counts, one = read_matrix_as_anndata(matrix_path(args.output_dir, sample), sample)
            sums[str(sample["sample_id"])] = counts
            adatas.append(one)
        else:
            sums[str(sample["sample_id"])] = sum_matrix(matrix_path(args.output_dir, sample))
    expression = pd.DataFrame(sums).fillna(0.0)
    expression.index.name = "gene"
    expression.reset_index().to_csv(args.output_dir / "GSE189926_immune_pseudobulk_counts.csv", index=False)
    h5ad_path = args.output_dir / "GSE189926_CD45_immune_raw_counts.h5ad"
    if args.write_h5ad:
        atlas = ad.concat(adatas, join="outer", merge="same", index_unique=None)
        if not atlas.obs_names.is_unique:
            raise ValueError("GSE189926 sample-qualified cell IDs are not unique.")
        atlas.X = atlas.X.tocsr().astype(np.int32)
        atlas.uns["provenance"] = {
            "geo_accession": ACCESSION,
            "source_page": SOURCE_PAGE,
            "raw_count_matrix_source": RAW_TAR_URL,
            "assay_scope": "CD45-selected immune single-cell suspensions",
            "validation_only": True,
        }
        print(f"[write] {h5ad_path} ({atlas.n_obs} cells, {atlas.n_vars} genes)", flush=True)
        atlas.write_h5ad(h5ad_path, compression="gzip")
    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "dataset_id": ACCESSION,
        "source_page": SOURCE_PAGE,
        "raw_matrix_source": RAW_TAR_URL,
        "n_samples": int(len(clinical)),
        "n_patients": int(clinical["patient_id"].nunique()),
        "raw_count_h5ad": str(h5ad_path) if args.write_h5ad else None,
        "assay_scope": "CD45-selected immune single-cell suspensions",
        "claim_boundary": (
            "Sample-level immune pseudobulk from public raw matrices. It can assess "
            "prespecified immune-state signatures and paired treatment context, but not "
            "CAF/endothelial states or response prediction."
        ),
    }
    (args.output_dir / "GSE189926_PREPARATION_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
