#!/usr/bin/env python3
"""Combine independently recreated cohort assignments after exact-match checks."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


CORE_IDS = [
    "korea_kim2022", "kumar2022", "sathe2020", "zhang2021",
    "diffuse_gc_2021", "tcell_exhaustion_2022",
]


def summary(frame: pd.DataFrame) -> pd.DataFrame:
    return (
        frame.groupby(["dataset_id", "broad_label", "native_cluster"], observed=True)
        .agg(n_cells=("native_cluster", "size"), n_samples=("analysis_sample", "nunique"))
        .reset_index()
        .assign(native_cluster=lambda table: table.native_cluster.astype(str))
        .sort_values(["dataset_id", "broad_label", "native_cluster"])
        .reset_index(drop=True)
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assignment-dir", type=Path, required=True)
    parser.add_argument("--discovery-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    paths = [args.assignment_dir / f"NATIVE_PER_COHORT_CELL_ASSIGNMENTS_{cohort}.csv.gz" for cohort in CORE_IDS]
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        (args.output_dir / "MISSING_RECREATION_INPUTS.txt").write_text("\n".join(missing) + "\n", encoding="utf-8")
        raise FileNotFoundError("Not all six audited assignment checkpoints exist. See MISSING_RECREATION_INPUTS.txt")
    assignments = pd.concat([pd.read_csv(path, compression="gzip") for path in paths], ignore_index=True)
    assignments["native_cluster"] = assignments["native_cluster"].astype(str)
    original = pd.read_csv(args.discovery_dir / "NATIVE_PER_COHORT_CLUSTER_SUMMARY.csv")
    original["native_cluster"] = original["native_cluster"].astype(str)
    recreated = summary(assignments)
    audit = original.merge(
        recreated, on=["dataset_id", "broad_label", "native_cluster"], how="outer",
        suffixes=("_original", "_recreated"), indicator=True,
    )
    audit["exact_match"] = (
        audit["_merge"].eq("both")
        & audit["n_cells_original"].eq(audit["n_cells_recreated"])
        & audit["n_samples_original"].eq(audit["n_samples_recreated"])
    )
    audit.to_csv(args.output_dir / "NATIVE_CLUSTER_RECREATION_AUDIT_ALL_COHORTS.csv", index=False)
    if not audit.exact_match.all():
        raise RuntimeError("Combined recreation differs from the reviewed native partition; no downstream composition was written.")
    assignments.to_csv(args.output_dir / "NATIVE_PER_COHORT_CELL_ASSIGNMENTS_ALL_COHORTS.csv.gz", index=False, compression="gzip")
    (args.output_dir / "RUN_MANIFEST.json").write_text(json.dumps({
        "n_assignment_rows": int(len(assignments)),
        "cohorts": CORE_IDS,
        "cluster_recreation_exact_match": True,
    }, indent=2) + "\n", encoding="utf-8")
    print(f"Combined exact-checked assignments written to {args.output_dir}")


if __name__ == "__main__":
    main()
