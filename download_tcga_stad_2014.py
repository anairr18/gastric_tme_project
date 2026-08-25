#!/usr/bin/env python3
"""Prepare open TCGA-STAD 2014 expression and molecular-subtype tables.

This uses the original GDC publication supplements rather than a third-party
portal. The resulting files are inputs to a frozen-signature association, not
to a subtype classifier or survival model.
"""

from __future__ import annotations

import argparse
import json
import shutil
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


MASTER_PATIENT_TABLE = (
    "https://api.gdc.cancer.gov/data/017f7658-4b39-493e-ad16-e00739a56118"
)
RPKM_MATRIX = "https://api.gdc.cancer.gov/data/566f684e-ad4e-41c6-bef0-1915fa6777a7"
SOURCE_PAGE = "https://gdc.cancer.gov/about-data/publications/stad_2014"


def download(url: str, destination: Path, minimum_bytes: int = 1024) -> None:
    if destination.exists() and destination.stat().st_size >= minimum_bytes:
        print(f"[reuse] {destination}", flush=True)
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": "gastric-tme-validation/1.0"})
    print(f"[download] {url}", flush=True)
    with urllib.request.urlopen(request, timeout=300) as response, temporary.open("wb") as handle:
        shutil.copyfileobj(response, handle)
    if temporary.stat().st_size < minimum_bytes:
        temporary.unlink(missing_ok=True)
        raise ValueError(f"Unexpectedly small download: {destination}")
    temporary.replace(destination)


def patient_barcode(sample_id: str) -> str | None:
    fields = str(sample_id).split("-")
    if len(fields) < 4 or fields[0:1] != ["TCGA"]:
        return None
    # The 01 code identifies the primary-tumour aliquots used here.
    if not fields[3].startswith("01"):
        return None
    return "-".join(fields[:3])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("/content/drive/MyDrive/data/external/validation_extensions/TCGA_STAD_SUBTYPES"),
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    master_path = args.output_dir / "TCGA_STAD_2014_MASTER_PATIENT_TABLE.xlsx"
    rpkm_path = args.output_dir / "TCGA_STAD_2014_RPKM.tsv"
    download(MASTER_PATIENT_TABLE, master_path, minimum_bytes=10_000)
    download(RPKM_MATRIX, rpkm_path, minimum_bytes=10_000_000)

    master = pd.read_excel(master_path, sheet_name=0)
    required = {"TCGA barcode", "Molecular Subtype"}
    if missing := required - set(master.columns):
        raise ValueError(f"TCGA master table lacks columns: {sorted(missing)}")
    clinical = master[["TCGA barcode", "Molecular Subtype"]].copy()
    clinical.columns = ["sample_id", "molecular_subtype"]
    clinical["sample_id"] = clinical["sample_id"].astype(str).str.upper()
    clinical["molecular_subtype"] = clinical["molecular_subtype"].astype(str).str.strip()
    clinical = clinical[clinical["molecular_subtype"].ne("") & clinical["molecular_subtype"].ne("nan")]
    clinical = clinical.drop_duplicates("sample_id")

    raw = pd.read_csv(rpkm_path, sep="\t", low_memory=False)
    if "GeneID" not in raw.columns:
        raise ValueError("TCGA RPKM matrix lacks GeneID.")
    selected: dict[str, str] = {}
    for column in raw.columns[1:]:
        patient = patient_barcode(column)
        if patient and patient in set(clinical["sample_id"]) and patient not in selected:
            selected[patient] = str(column)
    if len(selected) < 100:
        raise ValueError(f"Only {len(selected)} TCGA primary tumours matched subtype metadata.")

    expression = raw[["GeneID", *selected.values()]].copy()
    expression.insert(0, "gene", expression.pop("GeneID").astype(str).str.split("|").str[0].str.upper())
    expression = expression[expression["gene"].str.match(r"^[A-Z0-9][A-Z0-9.-]*$", na=False)]
    expression = expression.groupby("gene", as_index=False).mean(numeric_only=True)
    rename = {source: patient for patient, source in selected.items()}
    expression = expression.rename(columns=rename)
    matched_patients = [patient for patient in selected if patient in set(expression.columns)]
    expression = expression[["gene", *matched_patients]]
    clinical = clinical.set_index("sample_id").loc[matched_patients].reset_index()

    expression.to_csv(args.output_dir / "TCGA_STAD_expression.csv", index=False)
    clinical.to_csv(args.output_dir / "TCGA_STAD_clinical.csv", index=False)
    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "source_page": SOURCE_PAGE,
        "master_patient_table": MASTER_PATIENT_TABLE,
        "rpkm_matrix": RPKM_MATRIX,
        "n_primary_tumours_with_subtype": len(matched_patients),
        "n_genes": int(expression.shape[0]),
        "claim_boundary": (
            "Open TCGA-STAD 2014 primary-tumour RPKM data and supplied molecular "
            "subtypes. Frozen signature association only; not a subtype classifier."
        ),
    }
    (args.output_dir / "TCGA_STAD_PREPARATION_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
