#!/usr/bin/env python3
"""Validate a completed native-cluster review workbook before harmonisation.

This validates curation fields and produces a state-level inventory. It never
modifies the review sheet or invents labels for unresolved clusters.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


REQUIRED_FOR_INCLUDE = ["reviewed_state_name", "reviewed_marker_rationale", "reviewer"]


def read_review(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in {".xlsx", ".xls"}:
        return pd.read_excel(path, sheet_name="review_all_clusters")
    return pd.read_csv(path)


def nonempty(series: pd.Series) -> pd.Series:
    return series.fillna("").astype(str).str.strip().ne("")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--review-workbook", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--strict", action="store_true", help="Exit non-zero unless all review rows are finalised.")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table = read_review(args.review_workbook)
    required = {"dataset_id", "broad_label", "native_cluster", "review_decision", "curation_recommendation", "n_cells", "n_samples"}
    if missing := required - set(table.columns):
        raise ValueError(f"Review workbook lacks required columns: {sorted(missing)}")
    table["review_decision"] = table["review_decision"].fillna("review").astype(str).str.strip().str.lower()
    legal = {"include", "exclude", "review"}
    invalid = ~table["review_decision"].isin(legal)
    include = table["review_decision"].eq("include")
    missing_include_fields = pd.Series(False, index=table.index)
    for column in REQUIRED_FOR_INCLUDE:
        if column not in table:
            table[column] = ""
        missing_include_fields |= include & ~nonempty(table[column])
    forced_exclusion = table["curation_recommendation"].astype(str).str.startswith("exclude_")
    conflict = include & forced_exclusion
    table["validation_status"] = "valid_exclusion"
    table.loc[table["review_decision"].eq("review"), "validation_status"] = "pending_review"
    table.loc[invalid, "validation_status"] = "invalid_decision"
    table.loc[missing_include_fields, "validation_status"] = "include_missing_required_curation_fields"
    table.loc[conflict, "validation_status"] = "include_conflicts_with_lineage_audit"
    table.loc[include & ~missing_include_fields & ~conflict, "validation_status"] = "valid_include"

    table.to_csv(args.output_dir / "REVIEWED_NATIVE_CLUSTER_VALIDATION.csv", index=False)
    included = table.loc[table["validation_status"].eq("valid_include")].copy()
    included.to_csv(args.output_dir / "CURATED_NATIVE_CLUSTER_DICTIONARY.csv", index=False)
    state_inventory = (
        included.groupby(["broad_label", "reviewed_state_name"], observed=True)
        .agg(
            n_native_clusters=("native_cluster", "size"),
            n_cohorts=("dataset_id", "nunique"),
            cohorts=("dataset_id", lambda values: ";".join(sorted(set(values.astype(str))))),
            total_cells=("n_cells", "sum"),
            median_samples=("n_samples", "median"),
        )
        .reset_index()
    )
    state_inventory["harmonisation_status"] = state_inventory["n_cohorts"].map(
        lambda n: "candidate_cross_cohort" if n >= 3 else "cohort_specific_descriptive"
    )
    state_inventory.to_csv(args.output_dir / "REVIEWED_STATE_HARMONISATION_INVENTORY.csv", index=False)

    gate_summary = pd.DataFrame([
        {"gate": "All decisions finalised", "pass": bool((table["validation_status"] != "pending_review").all()), "detail": f"pending={int(table['validation_status'].eq('pending_review').sum())}"},
        {"gate": "Included clusters documented", "pass": bool((table.loc[include, "validation_status"] == "valid_include").all()), "detail": f"invalid_includes={int((include & ~table['validation_status'].eq('valid_include')).sum())}"},
        {"gate": "Cross-cohort candidate states", "pass": bool((state_inventory["n_cohorts"] >= 3).any()), "detail": f"n_candidates={int((state_inventory['n_cohorts'] >= 3).sum())}"},
    ])
    gate_summary.to_csv(args.output_dir / "REVIEWED_STATE_GATE_SUMMARY.csv", index=False)
    (args.output_dir / "README.md").write_text(
        "# Reviewed-state validation\n\n"
        "Only rows with `validation_status=valid_include` enter `CURATED_NATIVE_CLUSTER_DICTIONARY.csv`. A state becomes a cross-cohort candidate only after the same reviewed state name occurs in at least three cohorts. This inventory still does not establish disease association; it is the input gate for donor/sample-level state-abundance analysis.\n",
        encoding="utf-8",
    )
    (args.output_dir / "RUN_MANIFEST.json").write_text(json.dumps({
        "review_workbook": str(args.review_workbook),
        "n_clusters": int(len(table)),
        "n_valid_includes": int(len(included)),
        "n_pending_review": int(table["validation_status"].eq("pending_review").sum()),
        "n_cross_cohort_candidate_states": int((state_inventory["n_cohorts"] >= 3).sum()),
    }, indent=2) + "\n", encoding="utf-8")
    print(f"Reviewed-state validation written to {args.output_dir}")
    if args.strict and not gate_summary["pass"].all():
        raise RuntimeError("Review workbook has not passed all curation gates.")


if __name__ == "__main__":
    main()
