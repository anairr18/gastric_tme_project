"""Measure leave-one-cohort-out marker stability for frozen curated states.

The input profiles are exported from the frozen six-cohort atlas after state
labels have been reviewed. For each cohort/state pair, marker genes are ranked
from the *other* cohorts only, then tested in the held-out cohort against other
states in that same cohort. This is a reproducibility diagnostic, not another
round of clustering or feature selection.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def bh_fdr(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.dropna()
    if valid.empty:
        return result
    order = valid.sort_values().index
    adjusted = valid.loc[order].to_numpy() * len(order) / np.arange(1, len(order) + 1)
    result.loc[order] = np.minimum(np.minimum.accumulate(adjusted[::-1])[::-1], 1.0)
    return result


def state_compartment(state_id: str) -> str:
    return str(state_id).split("::", 1)[0]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort-profiles", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--top-genes", type=int, default=20)
    parser.add_argument("--min-state-cells", type=int, default=30)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    profiles = pd.read_csv(args.cohort_profiles, compression="infer")
    required = {"dataset_id", "state_id", "gene", "mean_log_normalized_expression", "n_state_cells"}
    missing = required - set(profiles.columns)
    if missing:
        raise ValueError(f"Cohort profile table is missing columns: {sorted(missing)}")
    profiles["state_id"] = profiles["state_id"].astype(str)
    profiles["dataset_id"] = profiles["dataset_id"].astype(str)
    profiles["gene"] = profiles["gene"].astype(str).str.upper()
    profiles["expression"] = pd.to_numeric(profiles["mean_log_normalized_expression"], errors="coerce")
    profiles["n_state_cells"] = pd.to_numeric(profiles["n_state_cells"], errors="coerce")
    profiles = profiles.dropna(subset=["expression", "n_state_cells"])

    rows: list[dict[str, object]] = []
    marker_rows: list[pd.DataFrame] = []
    for (cohort, state), heldout in profiles.groupby(["dataset_id", "state_id"], observed=True):
        n_cells = int(heldout["n_state_cells"].iloc[0])
        compartment = state_compartment(state)
        reference = profiles.loc[
            profiles["dataset_id"].ne(cohort) & profiles["state_id"].str.startswith(compartment + "::")
        ].copy()
        same_state = reference.loc[reference["state_id"].eq(state)]
        n_reference_cohorts = same_state["dataset_id"].nunique()
        if n_cells < args.min_state_cells or n_reference_cohorts < 2:
            rows.append(
                {
                    "dataset_id": cohort,
                    "state_id": state,
                    "n_state_cells": n_cells,
                    "n_reference_cohorts": n_reference_cohorts,
                    "status": "insufficient_cells_or_reference_cohorts",
                }
            )
            continue

        state_reference = same_state.groupby("gene", observed=True)["expression"].mean()
        other_reference = reference.loc[~reference["state_id"].eq(state)].groupby("gene", observed=True)["expression"].mean()
        reference_delta = state_reference.sub(other_reference, fill_value=np.nan).dropna().sort_values(ascending=False)
        markers = reference_delta.head(args.top_genes)
        heldout_state = heldout.set_index("gene")["expression"]
        heldout_other = profiles.loc[
            profiles["dataset_id"].eq(cohort)
            & profiles["state_id"].str.startswith(compartment + "::")
            & profiles["state_id"].ne(state)
        ].groupby("gene", observed=True)["expression"].mean()
        evidence = pd.DataFrame({"reference_delta": markers, "heldout_state": heldout_state, "heldout_other": heldout_other}).dropna()
        if len(evidence) < max(5, args.top_genes // 2):
            rows.append(
                {
                    "dataset_id": cohort,
                    "state_id": state,
                    "n_state_cells": n_cells,
                    "n_reference_cohorts": n_reference_cohorts,
                    "n_markers_tested": len(evidence),
                    "status": "insufficient_shared_marker_genes",
                }
            )
            continue
        heldout_delta = evidence["heldout_state"] - evidence["heldout_other"]
        rho, p_value = stats.spearmanr(evidence["reference_delta"], heldout_delta)
        marker_rows.append(
            evidence.assign(
                dataset_id=cohort,
                state_id=state,
                heldout_delta=heldout_delta,
                marker_rank=np.arange(1, len(evidence) + 1),
            ).reset_index(names="gene")
        )
        rows.append(
            {
                "dataset_id": cohort,
                "state_id": state,
                "n_state_cells": n_cells,
                "n_reference_cohorts": n_reference_cohorts,
                "n_markers_tested": len(evidence),
                "spearman_rho": float(rho),
                "spearman_p_value": float(p_value),
                "mean_heldout_marker_delta": float(heldout_delta.mean()),
                "positive_marker_fraction": float((heldout_delta > 0).mean()),
                "status": "tested",
            }
        )

    result_columns = [
        "dataset_id", "state_id", "n_state_cells", "n_reference_cohorts", "n_markers_tested",
        "spearman_rho", "spearman_p_value", "mean_heldout_marker_delta", "positive_marker_fraction",
        "status", "fdr_within_state", "replication_status",
    ]
    results = pd.DataFrame(rows, columns=result_columns)
    if not results.empty:
        tested = results["status"].eq("tested")
        results.loc[tested, "fdr_within_state"] = results.loc[tested].groupby("state_id", group_keys=False)[
            "spearman_p_value"
        ].apply(bh_fdr)
        results["replication_status"] = np.where(
            tested
            & results["spearman_rho"].ge(0.30)
            & results["mean_heldout_marker_delta"].gt(0)
            & results["positive_marker_fraction"].ge(0.60),
            "marker_consistent",
            np.where(tested, "marker_inconsistent", "not_tested"),
        )
    results.to_csv(args.output_dir / "STATE_MARKER_REPLICATION_BY_COHORT.csv", index=False)
    marker_table = pd.concat(marker_rows, ignore_index=True) if marker_rows else pd.DataFrame()
    marker_table.to_csv(args.output_dir / "STATE_MARKER_REPLICATION_EVIDENCE.csv", index=False)
    summary = (
        results.loc[results.get("status", pd.Series(dtype=str)).eq("tested")]
        .groupby("state_id", observed=True)
        .agg(
            n_tested_cohorts=("dataset_id", "nunique"),
            n_marker_consistent=("replication_status", lambda value: int((value == "marker_consistent").sum())),
            median_spearman_rho=("spearman_rho", "median"),
            median_positive_marker_fraction=("positive_marker_fraction", "median"),
        )
        .reset_index()
    )
    if not summary.empty:
        summary["state_marker_gate"] = np.where(
            (summary["n_tested_cohorts"] >= 3)
            & (summary["n_marker_consistent"] >= 3)
            & (summary["median_spearman_rho"] >= 0.30),
            "pass",
            "review_or_descriptive",
        )
    summary = summary.reindex(columns=[
        "state_id", "n_tested_cohorts", "n_marker_consistent", "median_spearman_rho",
        "median_positive_marker_fraction", "state_marker_gate",
    ])
    summary.to_csv(args.output_dir / "STATE_MARKER_REPLICATION_SUMMARY.csv", index=False)
    print(f"State-marker replication results: {args.output_dir}")
    if not summary.empty:
        print(summary.sort_values(["state_marker_gate", "median_spearman_rho"], ascending=[True, False]).to_string(index=False))


if __name__ == "__main__":
    main()
