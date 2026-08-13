#!/usr/bin/env python3
"""Recreate audited native cluster memberships without remaking figures/markers.

The original discovery export contained marker tables and UMAP figures but not
per-cell native cluster membership. This deterministic rerun uses the same raw
inputs, stratified seed, preprocessing, PCA/neighbour graph, and Leiden seed;
it skips UMAP and marker ranking. It fails closed if recreated cluster sizes do
not match the original discovery summary, preventing a review workbook from
being applied to a different partition.
"""

from __future__ import annotations

import argparse
import gc
import json
from pathlib import Path

import pandas as pd

from integration_method_sensitivity_benchmark import frozen_broad_label_lookup
from per_cohort_native_state_discovery import load_native_cohort, native_cluster
from raw_integration_method_audit import CORE_IDS


def canonical_summary(table: pd.DataFrame) -> pd.DataFrame:
    result = table.loc[:, ["dataset_id", "broad_label", "native_cluster", "n_cells", "n_samples"]].copy()
    result["native_cluster"] = result["native_cluster"].astype(str)
    return result.sort_values(["dataset_id", "broad_label", "native_cluster"]).reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--reference-atlas", type=Path, required=True)
    parser.add_argument("--discovery-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--cohorts", nargs="+", choices=CORE_IDS,
        help="Optional subset for interruption-safe cohort checkpoints. A final all-cohort run performs the exact-match audit.",
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    original_path = args.discovery_dir / "NATIVE_PER_COHORT_CLUSTER_SUMMARY.csv"
    if not original_path.exists():
        raise FileNotFoundError(f"Original discovery summary missing: {original_path}")
    original = canonical_summary(pd.read_csv(original_path))
    manifest = pd.read_csv(args.manifest).set_index("dataset_id", drop=False)
    missing = sorted(set(CORE_IDS) - set(manifest.index.astype(str)))
    if missing:
        raise ValueError(f"Manifest is missing core cohorts: {missing}")
    reference = frozen_broad_label_lookup(args.reference_atlas)
    assignment_tables, summary_tables = [], []

    selected_ids = args.cohorts or CORE_IDS
    for ordinal, dataset_id in enumerate(selected_ids):
        print(f"Recreating assignments: {dataset_id}", flush=True)
        data, _ = load_native_cohort(
            dataset_id, manifest.loc[dataset_id], args.data_root, reference, 30000, 17 + ordinal,
        )
        for broad in sorted(data.obs["frozen_broad_label"].astype(str).unique()):
            subset = data[data.obs["frozen_broad_label"].astype(str).eq(broad)].copy()
            if subset.n_obs < 100:
                continue
            clustered = native_cluster(
                subset, resolution=0.6, n_hvgs=3000, compute_umap=False, compute_markers=False,
            )
            summary = (
                clustered.obs.groupby("native_cluster", observed=True)
                .agg(n_cells=("native_cluster", "size"), n_samples=("analysis_sample", "nunique"))
                .reset_index()
                .assign(dataset_id=dataset_id, broad_label=broad)
            )
            assignment = clustered.obs.loc[:, [
                "dataset_id", "frozen_broad_label", "analysis_sample", "source_row_index", "source_cell_id", "native_cluster",
            ]].rename(columns={"frozen_broad_label": "broad_label"})
            assignment_tables.append(assignment)
            summary_tables.append(summary)
            del clustered, subset
            gc.collect()
        del data
        gc.collect()
        if assignment_tables:
            pd.concat(assignment_tables, ignore_index=True).to_csv(
                args.output_dir / "NATIVE_PER_COHORT_CELL_ASSIGNMENTS_PARTIAL.csv.gz", index=False, compression="gzip"
            )

    if not assignment_tables:
        raise RuntimeError("No cohort assignments were created.")

    recreated = canonical_summary(pd.concat(summary_tables, ignore_index=True))
    expected = original.loc[original["dataset_id"].isin(selected_ids)].copy()
    comparison = expected.merge(
        recreated, on=["dataset_id", "broad_label", "native_cluster"], how="outer", suffixes=("_original", "_recreated"), indicator=True,
    )
    comparison["exact_match"] = (
        comparison["_merge"].eq("both")
        & comparison["n_cells_original"].eq(comparison["n_cells_recreated"])
        & comparison["n_samples_original"].eq(comparison["n_samples_recreated"])
    )
    comparison.to_csv(args.output_dir / "NATIVE_CLUSTER_RECREATION_AUDIT.csv", index=False)
    if not comparison["exact_match"].all():
        failures = comparison.loc[~comparison["exact_match"]]
        raise RuntimeError(
            f"Recreated clustering does not exactly match original discovery ({len(failures)} clusters differ). "
            "Do not apply the reviewed workbook; investigate software/input drift."
        )
    assignments = pd.concat(assignment_tables, ignore_index=True)
    assignment_name = (
        "NATIVE_PER_COHORT_CELL_ASSIGNMENTS.csv.gz"
        if args.cohorts is None else f"NATIVE_PER_COHORT_CELL_ASSIGNMENTS_{selected_ids[0]}.csv.gz"
    )
    assignments.to_csv(args.output_dir / assignment_name, index=False, compression="gzip")
    (args.output_dir / "RUN_MANIFEST.json").write_text(json.dumps({
        "original_discovery_dir": str(args.discovery_dir),
        "n_assignment_rows": int(len(assignments)),
        "cluster_recreation_exact_match": True,
        "cohorts": list(selected_ids),
        "method": "same deterministic clustering; UMAP and marker ranking skipped",
    }, indent=2) + "\n", encoding="utf-8")
    print(f"Recreated assignments written to {args.output_dir}")


if __name__ == "__main__":
    main()
