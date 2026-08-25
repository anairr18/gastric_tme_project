#!/usr/bin/env python3
"""Download GSE163558 raw matrices and build a validation-only count H5AD."""

from __future__ import annotations

import argparse
import hashlib
import re
import shutil
import tarfile
import urllib.request
from pathlib import Path

import pandas as pd
import scanpy as sc


GEO_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE163nnn/GSE163558/suppl/GSE163558_RAW.tar"
GEO_DOWNLOAD_FALLBACK = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE163558&format=file&file=GSE163558_RAW.tar"
SAMPLES = [
    ("PT1", "GSM5004180", "Patient_1", "tumor", "primary", "intestinal"),
    ("PT2", "GSM5004181", "Patient_2", "tumor", "primary", "intestinal"),
    ("PT3", "GSM5004182", "Patient_4", "tumor", "primary", "intestinal"),
    ("NT1", "GSM5004183", "Patient_2", "normal", "adjacent", "normal"),
    ("LN1", "GSM5004184", "Patient_4", "metastatic", "lymph_node", "intestinal"),
    ("LN2", "GSM5004185", "Patient_5", "metastatic", "lymph_node", "mixed"),
    ("O1", "GSM5004186", "Patient_3", "metastatic", "ovary", "mixed"),
    ("P1", "GSM5004187", "Patient_6", "metastatic", "peritoneum", "intestinal"),
    ("Li1", "GSM5004188", "Patient_1", "metastatic", "liver", "intestinal"),
    ("Li2", "GSM5004189", "Patient_4", "metastatic", "liver", "intestinal"),
]


def download(url: str, path: Path) -> None:
    if path.exists() and path.stat().st_size > 10_000_000:
        return
    urls = (url, GEO_DOWNLOAD_FALLBACK) if url == GEO_URL else (url,)
    error: Exception | None = None
    for candidate in urls:
        try:
            request = urllib.request.Request(candidate, headers={"User-Agent": "gastric-tme-validation/1.0"})
            print(f"[download] {candidate}", flush=True)
            with urllib.request.urlopen(request) as response, path.open("wb") as handle:
                shutil.copyfileobj(response, handle)
            return
        except Exception as caught:
            error = caught
            path.unlink(missing_ok=True)
    raise RuntimeError(f"Could not download {url}: {type(error).__name__}: {error}")


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1_048_576), b""):
            digest.update(block)
    return digest.hexdigest()


def sample_member(member_name: str, sample: str, gsm: str) -> bool:
    name = Path(member_name).name
    return gsm in name or bool(re.search(rf"(?<![A-Za-z0-9]){re.escape(sample)}(?![A-Za-z0-9])", name))


def role(member_name: str) -> str | None:
    name = Path(member_name).name.lower()
    if "matrix.mtx" in name:
        return "matrix"
    if "barcodes.tsv" in name:
        return "barcodes"
    if "features.tsv" in name or "genes.tsv" in name:
        return "features"
    return None


def suffix(member_name: str) -> str:
    return ".gz" if member_name.lower().endswith(".gz") else ""


def build_h5ad(tar_path: Path, output: Path, work_dir: Path) -> pd.DataFrame:
    work_dir.mkdir(parents=True, exist_ok=True)
    sample_table = pd.DataFrame(
        SAMPLES,
        columns=["sample_id", "gsm_accession", "patient_id", "tissue", "site", "histology"],
    )
    adatas = []
    with tarfile.open(tar_path, "r") as archive:
        members = [member for member in archive.getmembers() if member.isfile()]
        for record in sample_table.itertuples(index=False):
            selected: dict[str, tarfile.TarInfo] = {}
            for member in members:
                member_role = role(member.name)
                if member_role and sample_member(member.name, record.sample_id, record.gsm_accession):
                    if member_role in selected:
                        raise ValueError(f"{record.sample_id}: multiple {member_role} files in GEO archive")
                    selected[member_role] = member
            required = {"matrix", "barcodes", "features"}
            if set(selected) != required:
                available = [member.name for member in members if record.gsm_accession in member.name]
                raise ValueError(f"{record.sample_id}: missing 10x files; GEO members={available[:12]}")
            sample_dir = work_dir / record.sample_id
            sample_dir.mkdir(exist_ok=True)
            destinations = {
                "matrix": sample_dir / f"matrix.mtx{suffix(selected['matrix'].name)}",
                "barcodes": sample_dir / f"barcodes.tsv{suffix(selected['barcodes'].name)}",
                "features": sample_dir / f"features.tsv{suffix(selected['features'].name)}",
            }
            for member_role, member in selected.items():
                with archive.extractfile(member) as source, destinations[member_role].open("wb") as destination:
                    if source is None:
                        raise ValueError(f"Could not read {member.name}")
                    shutil.copyfileobj(source, destination)
            one = sc.read_10x_mtx(sample_dir, var_names="gene_symbols", make_unique=True)
            one.obs_names = pd.Index(record.sample_id + "::" + one.obs_names.astype(str))
            one.obs["sample_id"] = record.sample_id
            one.obs["patient_id"] = record.patient_id
            one.obs["tissue"] = record.tissue
            one.obs["site"] = record.site
            one.obs["histology"] = record.histology
            one.obs["dataset_id"] = "GSE163558"
            adatas.append(one)
    atlas = ad.concat(adatas, join="outer", merge="same", index_unique=None)
    if not atlas.obs_names.is_unique:
        raise ValueError("Sample-qualified GSE163558 cell IDs are not unique")
    atlas.X = atlas.X.astype("float32")
    atlas.uns["provenance"] = {
        "geo_accession": "GSE163558",
        "raw_archive_url": GEO_URL,
        "patient_mapping_source": "Jiang et al., Clin Transl Med 2022, Results section 2.1",
        "validation_only": True,
    }
    atlas.write_h5ad(output, compression="gzip")
    return sample_table.assign(n_cells=atlas.obs["sample_id"].value_counts().reindex(sample_table["sample_id"]).to_numpy())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--keep-extracted", action="store_true")
    args = parser.parse_args()
    root = args.data_root / "external" / "validation" / "GSE163558"
    root.mkdir(parents=True, exist_ok=True)
    archive = root / "GSE163558_RAW.tar"
    output = root / "GSE163558.h5ad"
    metadata_path = root / "GSE163558_sample_metadata.csv"
    work_dir = root / "raw_mtx"
    download(GEO_URL, archive)
    sha256 = sha256sum(archive)
    metadata = build_h5ad(archive, output, work_dir)
    metadata.to_csv(metadata_path, index=False)
    (root / "PROVENANCE.md").write_text(
        f"# GSE163558 Validation Input\n\n- Raw archive: `{GEO_URL}`\n- SHA256: `{sha256}`\n- Patient/sample mapping: Jiang et al., Clin Transl Med 2022, Results 2.1.\n- Role: validation-only; excluded from frozen scVI core.\n",
        encoding="utf-8",
    )
    if not args.keep_extracted:
        shutil.rmtree(work_dir)
    print(f"Wrote {output}")
    print(metadata.to_string(index=False))


if __name__ == "__main__":
    main()
