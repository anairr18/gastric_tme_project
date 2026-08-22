#!/usr/bin/env python3
"""Validate and reuse the exact per-cell assignments from native discovery.

Native Leiden labels are implementation-sensitive across Scanpy/igraph/leidenalg
versions. If discovery exported its per-cell assignment table, it is the frozen
source of truth for downstream sample-level inference. This script verifies
that table against the reviewed discovery cluster summary before reusing it.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


CORE_IDS = [
    "korea_kim2022", "kumar2022", "sathe2020", "zhang2021",
    "diffuse_gc_2021", "tcell_exhaustion_2022",
]
KEY = ["dataset_id", "broad_label", "native_cluster"]
REQUIRED_COLUMNS = KEY + ["analysis_sample", "source_row_index", "source_cell_id"]


def cluster_summary(frame: pd.DataFrame) -> pd.DataFrame:
    return (
        frame.groupby(KEY, observed=True)
        .agg(n_cells=("native_cluster", "size"), n_samples=("analysis_sample", "nunique"))
        .reset_index()
        .sort_values(KEY)
        .reset_index(drop=True)
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assignments", type=Path, required=True)
    parser.add_argument("--discovery-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    assignments = pd.read_csv(args.assignments, compression="infer")
    missing = [column for column in REQUIRED_COLUMNS if column not in assignments.columns]
    if missing:
        raise ValueError(f"Frozen assignment table lacks required columns: {missing}")
    assignments = assignments.loc[:, REQUIRED_COLUMNS].copy()
    assignments["dataset_id"] = assignments["dataset_id"].astype(str)
    assignments["broad_label"] = assignments["broad_label"].astype(str)
    assignments["native_cluster"] = assignments["native_cluster"].astype(str)
    assignments["analysis_sample"] = assignments["analysis_sample"].astype(str)

    unexpected = sorted(set(assignments["dataset_id"]) - set(CORE_IDS))
    absent = sorted(set(CORE_IDS) - set(assignments["dataset_id"]))
    if unexpected or absent:
        raise ValueError(f"Frozen assignments have unexpected={unexpected}, absent_core_cohorts={absent}")
    if assignments.duplicated(["dataset_id", "source_row_index"]).any():
        raise ValueError("Frozen assignments have duplicate dataset_id/source_row_index pairs.")

    original = pd.read_csv(args.discovery_dir / "NATIVE_PER_COHORT_CLUSTER_SUMMARY.csv")
    original["dataset_id"] = original["dataset_id"].astype(str)
    original["broad_label"] = original["broad_label"].astype(str)
    original["native_cluster"] = original["native_cluster"].astype(str)
    original = original.loc[:, KEY + ["n_cells", "n_samples"]].sort_values(KEY).reset_index(drop=True)
    recovered = cluster_summary(assignments)
    audit = original.merge(recovered, on=KEY, how="outer", suffixes=("_discovery", "_assignments"), indicator=True)
    audit["exact_match"] = (
        audit["_merge"].eq("both")
        & audit["n_cells_discovery"].eq(audit["n_cells_assignments"])
        & audit["n_samples_discovery"].eq(audit["n_samples_assignments"])
    )
    audit.to_csv(args.output_dir / "FROZEN_NATIVE_ASSIGNMENT_VALIDATION_AUDIT.csv", index=False)
    if not audit["exact_match"].all():
        raise RuntimeError(
            f"Frozen assignment table differs from reviewed discovery summary ({int((~audit['exact_match']).sum())} cluster rows)."
        )

    assignments.to_csv(
        args.output_dir / "NATIVE_PER_COHORT_CELL_ASSIGNMENTS_ALL_COHORTS.csv.gz",
        index=False,
        compression="gzip",
    )
    (args.output_dir / "RUN_MANIFEST.json").write_text(json.dumps({
        "assignment_source": str(args.assignments),
        "original_discovery_dir": str(args.discovery_dir),
        "n_assignment_rows": int(len(assignments)),
        "cluster_partition_exact_match": True,
        "purpose": "reuse frozen native discovery assignments; no reclustering performed",
    }, indent=2) + "\n", encoding="utf-8")
    print(f"Validated frozen assignments written to {args.output_dir}")


if __name__ == "__main__":
    main()
