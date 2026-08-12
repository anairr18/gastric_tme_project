#!/usr/bin/env python3
"""Reconcile the biological context of every frozen-core cohort.

This script intentionally separates project cell records from source-study sample
and patient counts.  It is designed for the mentor-facing Methods audit, not for
creating new biological conclusions.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


CORE_IDS = {
    "korea_kim2022",
    "kumar2022",
    "sathe2020",
    "zhang2021",
    "diffuse_gc_2021",
    "tcell_exhaustion_2022",
}


def read_csv(path: Path, label: str) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"{label} was not found: {path}")
    table = pd.read_csv(path)
    if table.empty:
        raise ValueError(f"{label} is empty: {path}")
    return table


def locate_column(frame: pd.DataFrame, candidates: list[str]) -> str | None:
    lowered = {str(column).lower(): str(column) for column in frame.columns}
    for candidate in candidates:
        if candidate.lower() in lowered:
            return lowered[candidate.lower()]
    return None


def atlas_project_counts(atlas: Path) -> pd.DataFrame:
    """Read metadata only where possible; never load the expression matrix."""
    if not atlas or not atlas.exists():
        return pd.DataFrame(columns=["dataset_id", "project_atlas_cells", "project_atlas_samples", "project_atlas_patients"])
    try:
        import anndata as ad
    except ImportError as error:
        raise ImportError("Install anndata to audit the atlas metadata.") from error

    data = ad.read_h5ad(atlas, backed="r")
    obs = data.obs.copy()
    data.file.close()
    dataset_col = locate_column(obs, ["dataset_id", "cohort", "canonical_cohort"])
    if dataset_col is None:
        raise ValueError("The atlas lacks dataset_id/cohort metadata; cannot count project cell records.")
    sample_col = locate_column(obs, ["sample_id", "sample", "orig.ident", "library_id", "analysis_unit"])
    patient_col = locate_column(obs, ["patient_id", "patient", "donor_id", "donor"])

    records: list[dict[str, object]] = []
    for dataset_id, subset in obs.groupby(dataset_col, dropna=False, observed=True):
        records.append(
            {
                "dataset_id": str(dataset_id),
                "project_atlas_cells": int(len(subset)),
                "project_atlas_samples": int(subset[sample_col].dropna().astype(str).nunique()) if sample_col else pd.NA,
                "project_atlas_patients": int(subset[patient_col].dropna().astype(str).nunique()) if patient_col else pd.NA,
            }
        )
    return pd.DataFrame.from_records(records)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--registry", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--atlas", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    registry = read_csv(args.registry, "Cohort context registry")
    manifest = read_csv(args.manifest, "Dataset manifest")
    for label, frame in [("registry", registry), ("manifest", manifest)]:
        if "dataset_id" not in frame.columns:
            raise ValueError(f"The {label} needs a dataset_id column.")
        if frame["dataset_id"].duplicated().any():
            duplicates = frame.loc[frame["dataset_id"].duplicated(), "dataset_id"].tolist()
            raise ValueError(f"The {label} has duplicate dataset IDs: {duplicates}")

    registry_ids = set(registry["dataset_id"].astype(str))
    manifest_ids = set(manifest["dataset_id"].astype(str))
    missing_registry = sorted(CORE_IDS - registry_ids)
    missing_manifest = sorted(CORE_IDS - manifest_ids)
    unexpected_registry = sorted(registry_ids - CORE_IDS)
    if missing_registry or missing_manifest or unexpected_registry:
        raise ValueError(
            "Frozen-core registry mismatch. "
            f"Missing registry={missing_registry}; missing manifest={missing_manifest}; "
            f"unexpected registry IDs={unexpected_registry}"
        )

    for required in ["disease_site", "metastatic_site_in_core", "esophageal_or_gej_in_core", "clinical_context", "treatment_context", "claim_boundary"]:
        if required not in registry.columns:
            raise ValueError(f"Cohort context registry is missing {required}.")
    for column in ["metastatic_site_in_core", "esophageal_or_gej_in_core"]:
        values = registry[column].astype(str).str.lower()
        if values.isin(["true", "1", "yes"]).any():
            flagged = registry.loc[values.isin(["true", "1", "yes"]), "dataset_id"].tolist()
            raise ValueError(f"Frozen core was expected to be gastric primary only, but {column} is true for {flagged}.")

    context = registry.merge(manifest, on="dataset_id", how="left", suffixes=("", "_manifest"), validate="one_to_one")
    counts = atlas_project_counts(args.atlas) if args.atlas else pd.DataFrame()
    if not counts.empty:
        context = context.merge(counts, on="dataset_id", how="left", validate="one_to_one")
    else:
        context["project_atlas_cells"] = pd.NA
        context["project_atlas_samples"] = pd.NA
        context["project_atlas_patients"] = pd.NA

    source_cells_col = locate_column(context, ["n_cells_raw", "published_cells", "source_cells"])
    source_patients_col = locate_column(context, ["n_patients", "published_patients"])
    context["source_reported_cells"] = context[source_cells_col] if source_cells_col else pd.NA
    context["source_reported_patients"] = context[source_patients_col] if source_patients_col else pd.NA
    context["cell_count_interpretation"] = (
        "project_atlas_cells are project cell records entering the frozen atlas; "
        "they are not patient counts, sample counts, or necessarily published study totals"
    )
    context["source_project_cell_count_reconciled"] = pd.NA
    both = context["project_atlas_cells"].notna() & context["source_reported_cells"].notna()
    context.loc[both, "source_project_cell_count_reconciled"] = (
        pd.to_numeric(context.loc[both, "project_atlas_cells"], errors="coerce")
        == pd.to_numeric(context.loc[both, "source_reported_cells"], errors="coerce")
    )
    context["reconciliation_action"] = "Use project and source counts with explicit labels; resolve discrepancies against source metadata before manuscript submission."

    columns = [
        "dataset_id", "display_name", "accession_or_doi", "disease_site", "metastatic_site_in_core",
        "esophageal_or_gej_in_core", "clinical_context", "treatment_context", "serial_sampling",
        "primary_comparison", "source_reported_cells", "source_reported_patients", "project_atlas_cells",
        "project_atlas_samples", "project_atlas_patients", "cell_count_interpretation", "source_design_note",
        "claim_boundary", "source_project_cell_count_reconciled", "reconciliation_action",
    ]
    columns = [column for column in columns if column in context.columns]
    slides = context.loc[:, columns].sort_values("dataset_id")
    slides.to_csv(args.output_dir / "COHORT_CONTEXT_FOR_SLIDES.csv", index=False)

    audit = slides.loc[:, ["dataset_id", "disease_site", "metastatic_site_in_core", "esophageal_or_gej_in_core", "treatment_context", "serial_sampling", "primary_comparison", "source_project_cell_count_reconciled", "reconciliation_action"]]
    audit.to_csv(args.output_dir / "COHORT_CONTEXT_AUDIT.csv", index=False)
    reconciliation = slides.loc[:, [column for column in ["dataset_id", "source_reported_cells", "project_atlas_cells", "source_reported_patients", "project_atlas_patients", "source_project_cell_count_reconciled", "reconciliation_action"] if column in slides.columns]]
    reconciliation.to_csv(args.output_dir / "COHORT_CONTEXT_RECONCILIATION.csv", index=False)

    lines = [
        "# Frozen Core Cohort Context",
        "",
        "The frozen scVI core contains six gastric cohorts. It contains no esophageal/GEJ cohort and no metastatic-site cohort.",
        "Project atlas cell counts are cell records, not patients or samples. Source-study counts are retained separately and must not be interchanged.",
        "",
        "## Cohort-specific boundaries",
    ]
    for row in slides.itertuples(index=False):
        lines.append(f"- **{row.dataset_id}**: {row.clinical_context}. {row.treatment_context}. {row.claim_boundary}")
    (args.output_dir / "COHORT_CONTEXT_CLAIM_BOUNDARIES.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    (args.output_dir / "RUN_MANIFEST.json").write_text(json.dumps({"registry": str(args.registry), "manifest": str(args.manifest), "atlas": str(args.atlas) if args.atlas else None, "core_dataset_ids": sorted(CORE_IDS)}, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote cohort context audit: {args.output_dir}")


if __name__ == "__main__":
    main()
