#!/usr/bin/env python3
"""Apply validated reviewed native-state labels and build sample-level composition.

State fractions are calculated within each broad compartment and sample from
the deterministic native-clustering assignment table. This is the appropriate
input to later donor/sample-level condition analyses; cells are never used as
independent observations.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


KEY = ["dataset_id", "broad_label", "native_cluster"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assignments", type=Path, required=True)
    parser.add_argument("--validated-review-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    assignments = pd.read_csv(args.assignments, compression="infer")
    review = pd.read_csv(args.validated_review_dir / "CURATED_NATIVE_CLUSTER_DICTIONARY.csv")
    for table in [assignments, review]:
        table["native_cluster"] = table["native_cluster"].astype(str)
    if review.duplicated(KEY).any():
        raise ValueError("Validated review dictionary has duplicate native-cluster keys.")
    valid = review.loc[review["validation_status"].eq("valid_include")].copy()
    if valid.empty:
        raise ValueError("No valid included native clusters. Complete the review workbook first.")
    state_names = valid["reviewed_state_name"].fillna("").astype(str).str.strip()
    if state_names.eq("").any():
        raise ValueError("A valid included cluster lacks a reviewed state name.")
    merged = assignments.merge(
        valid[KEY + ["reviewed_state_name", "reviewed_marker_rationale", "reviewer"]],
        on=KEY, how="left", validate="many_to_one",
    )
    denominators = (
        assignments.groupby(["dataset_id", "analysis_sample", "broad_label"], observed=True)
        .size().rename("n_cells_broad_compartment").reset_index()
    )
    included = merged.loc[merged["reviewed_state_name"].notna()].copy()
    composition = (
        included.groupby(["dataset_id", "analysis_sample", "broad_label", "reviewed_state_name"], observed=True)
        .size().rename("n_cells_state").reset_index()
        .merge(denominators, on=["dataset_id", "analysis_sample", "broad_label"], how="left", validate="many_to_one")
    )
    composition["state_fraction"] = composition["n_cells_state"] / composition["n_cells_broad_compartment"]
    composition = composition.rename(columns={"analysis_sample": "sample", "broad_label": "compartment", "reviewed_state_name": "reviewed_state_label"})
    composition["analysis_unit"] = "sample; state fraction within broad compartment"
    composition["claim_boundary"] = "Native state composition derives from a stratified raw-cell subset; inferential unit is sample/donor, never cell."
    composition.to_csv(args.output_dir / "CURATED_NATIVE_STATE_SAMPLE_COMPOSITION.csv", index=False)
    mapping_audit = merged.loc[:, ["dataset_id", "analysis_sample", "broad_label", "native_cluster", "reviewed_state_name"]].copy()
    mapping_audit["included_after_manual_review"] = mapping_audit["reviewed_state_name"].notna()
    mapping_audit.to_csv(args.output_dir / "NATIVE_CLUSTER_TO_REVIEWED_STATE_MAPPING_AUDIT.csv", index=False)
    state_inventory = (
        composition.groupby(["compartment", "reviewed_state_label"], observed=True)
        .agg(n_cohorts=("dataset_id", "nunique"), n_samples=("sample", "nunique"), total_state_cells=("n_cells_state", "sum"))
        .reset_index()
    )
    state_inventory["harmonisation_status"] = state_inventory["n_cohorts"].map(
        lambda n: "candidate_cross_cohort" if n >= 3 else "cohort_specific_descriptive"
    )
    state_inventory.to_csv(args.output_dir / "CURATED_NATIVE_STATE_INVENTORY.csv", index=False)
    (args.output_dir / "RUN_MANIFEST.json").write_text(json.dumps({
        "assignments": str(args.assignments),
        "validated_review_dir": str(args.validated_review_dir),
        "n_assignment_rows": int(len(assignments)),
        "n_included_assignment_rows": int(len(included)),
        "n_state_sample_rows": int(len(composition)),
    }, indent=2) + "\n", encoding="utf-8")
    print(f"Reviewed native state composition written to {args.output_dir}")


if __name__ == "__main__":
    main()
