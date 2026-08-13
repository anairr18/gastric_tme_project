#!/usr/bin/env python3
"""Build sample-level composition using the accepted provisional state dictionary."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


KEY = ["dataset_id", "broad_label", "native_cluster"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assignments", type=Path, required=True)
    parser.add_argument("--resolution-ledger", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    assignments = pd.read_csv(args.assignments, compression="infer")
    ledger = pd.read_csv(args.resolution_ledger)
    for frame in (assignments, ledger):
        frame["native_cluster"] = frame["native_cluster"].astype(str)
    if ledger.duplicated(KEY).any():
        raise ValueError("Resolution ledger contains duplicate native-cluster keys.")
    retained = ledger.loc[ledger.analysis_tier.eq("primary_working_state")].copy()
    if retained.empty:
        raise ValueError("No primary working states found in resolution ledger.")
    states = retained[KEY + ["provisional_state_name", "matched_state_markers", "resolution_rationale"]]
    merged = assignments.merge(states, on=KEY, how="left", validate="many_to_one")
    denominators = (
        assignments.groupby(["dataset_id", "analysis_sample", "broad_label"], observed=True)
        .size().rename("n_cells_broad_compartment").reset_index()
    )
    included = merged.loc[merged.provisional_state_name.notna()].copy()
    composition = (
        included.groupby(["dataset_id", "analysis_sample", "broad_label", "provisional_state_name"], observed=True)
        .size().rename("n_cells_state").reset_index()
        .merge(denominators, on=["dataset_id", "analysis_sample", "broad_label"], how="left", validate="many_to_one")
    )
    composition["state_fraction"] = composition.n_cells_state / composition.n_cells_broad_compartment
    composition = composition.rename(columns={
        "analysis_sample": "sample", "broad_label": "compartment", "provisional_state_name": "reviewed_state_label",
    })
    composition["state_test_tier"] = "candidate_confirmatory"
    composition["analysis_unit"] = "sample; within-compartment state fraction"
    composition["claim_boundary"] = (
        "Working computational curation accepted by the project owner; all inferential tests use samples/donors, never cells. "
        "State names remain pending mentor confirmation."
    )
    composition.to_csv(args.output_dir / "PRIMARY_WORKING_STATE_SAMPLE_COMPOSITION.csv", index=False)
    inventory = (
        composition.groupby(["compartment", "reviewed_state_label"], observed=True)
        .agg(n_cohorts=("dataset_id", "nunique"), n_samples=("sample", "nunique"), total_state_cells=("n_cells_state", "sum"))
        .reset_index()
    )
    inventory["crosscohort_candidate"] = inventory.n_cohorts.ge(3)
    inventory.to_csv(args.output_dir / "PRIMARY_WORKING_STATE_INVENTORY.csv", index=False)
    mapping = merged.loc[:, KEY + ["analysis_sample", "provisional_state_name"]].copy()
    mapping["included_primary_working_state"] = mapping.provisional_state_name.notna()
    mapping.to_csv(args.output_dir / "PRIMARY_STATE_ASSIGNMENT_MAPPING_AUDIT.csv", index=False)
    (args.output_dir / "RUN_MANIFEST.json").write_text(json.dumps({
        "assignments": str(args.assignments), "resolution_ledger": str(args.resolution_ledger),
        "n_assignment_rows": int(len(assignments)), "n_primary_assignment_rows": int(len(included)),
        "n_state_sample_rows": int(len(composition)), "status": "working_curation_pending_mentor_review",
    }, indent=2) + "\n", encoding="utf-8")
    print(f"Primary working-state composition written to {args.output_dir}")


if __name__ == "__main__":
    main()
