#!/usr/bin/env python3
"""Resolve every native cluster to an auditable analysis disposition.

This intentionally separates a *cluster disposition* from a *biological state
label*.  All clusters receive a disposition, while only marker-supported,
lineage-audited clusters enter the primary working state set.  This prevents
the common but misleading practice of assigning a precise biological identity
to technical, mixed-lineage, or underpowered Leiden clusters.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dictionary", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table = pd.read_csv(args.dictionary).copy()
    required = {
        "dataset_id", "broad_label", "native_cluster", "n_cells", "n_samples",
        "curation_recommendation", "provisional_decision", "provisional_state_name",
        "state_marker_score", "runner_up_state_marker_score", "matched_state_markers",
        "top_30_native_markers",
    }
    if missing := required - set(table.columns):
        raise ValueError(f"Dictionary missing required columns: {sorted(missing)}")

    table["analysis_tier"] = "excluded"
    table["cluster_disposition"] = "unresolved"
    table["state_label_for_primary_analysis"] = ""
    table["candidate_state_for_descriptive_inventory"] = ""
    table["resolution_rationale"] = ""

    primary = table.provisional_decision.eq("include_provisional")
    table.loc[primary, "analysis_tier"] = "primary_working_state"
    table.loc[primary, "cluster_disposition"] = "retain_provisional_state"
    table.loc[primary, "state_label_for_primary_analysis"] = table.loc[primary, "provisional_state_name"]
    table.loc[primary, "candidate_state_for_descriptive_inventory"] = table.loc[primary, "provisional_state_name"]
    table.loc[primary, "resolution_rationale"] = (
        "Passed lineage audit, support threshold, and transparent raw-marker panel. "
        "Pending mentor confirmation before final manuscript claims."
    )

    cross = table.curation_recommendation.eq("exclude_suspected_cross_lineage_or_misassignment")
    table.loc[cross, "cluster_disposition"] = "exclude_cross_lineage_or_misassignment"
    table.loc[cross, "resolution_rationale"] = (
        "Excluded by lineage audit because discordant marker evidence makes a fine-state label unreliable."
    )

    low = table.curation_recommendation.eq("exclude_insufficient_cell_or_sample_support")
    table.loc[low, "cluster_disposition"] = "exclude_insufficient_cell_or_sample_support"
    table.loc[low, "resolution_rationale"] = (
        "Excluded by prespecified cell/sample support rule; retain only in the supplemental cluster log."
    )

    manual = table.curation_recommendation.eq("manual_review_required")
    clear_candidate = (
        manual
        & table.provisional_state_name.ne("unresolved")
        & table.state_marker_score.ge(2)
        & table.state_marker_score.gt(table.runner_up_state_marker_score)
    )
    table.loc[clear_candidate, "analysis_tier"] = "secondary_candidate"
    table.loc[clear_candidate, "cluster_disposition"] = "candidate_requires_mentor_confirmation"
    table.loc[clear_candidate, "candidate_state_for_descriptive_inventory"] = table.loc[clear_candidate, "provisional_state_name"]
    table.loc[clear_candidate, "resolution_rationale"] = (
        "Marker panel suggests a state, but the initial lineage audit required manual review. "
        "Excluded from primary inference until mentor confirmation."
    )

    competing = manual & ~clear_candidate & table.state_marker_score.ge(2)
    table.loc[competing, "cluster_disposition"] = "unresolved_competing_state_programs"
    table.loc[competing, "resolution_rationale"] = (
        "More than one state program is similarly supported by native markers; no state label is assigned."
    )

    weak = manual & ~clear_candidate & ~competing
    table.loc[weak, "cluster_disposition"] = "unresolved_inadequate_state_marker_evidence"
    table.loc[weak, "resolution_rationale"] = (
        "No sufficiently specific native marker program supports a reproducible fine-state label."
    )

    eligible_ambiguous = (
        table.curation_recommendation.eq("eligible_after_manual_state_annotation")
        & ~primary
    )
    table.loc[eligible_ambiguous, "cluster_disposition"] = "unresolved_eligible_marker_ambiguity"
    table.loc[eligible_ambiguous, "resolution_rationale"] = (
        "Broad lineage passed the audit but no decisive fine-state program was detected; excluded conservatively."
    )

    table["review_status"] = "computationally_resolved_pending_mentor_review"
    table["claim_boundary"] = (
        "A complete disposition is not equivalent to expert annotation. Primary working states may enter planned "
        "donor/sample-level analyses but remain provisional until mentor review; secondary candidates and exclusions "
        "cannot support confirmatory biological claims."
    )
    table = table.sort_values(["dataset_id", "broad_label", "native_cluster"])
    table.to_csv(args.output_dir / "ALL_505_NATIVE_CLUSTER_RESOLUTION_LEDGER.csv", index=False)

    summary = (
        table.groupby(["analysis_tier", "cluster_disposition"], observed=True)
        .agg(n_clusters=("native_cluster", "size"), n_cells=("n_cells", "sum"), n_cohorts=("dataset_id", "nunique"))
        .reset_index()
        .sort_values(["analysis_tier", "cluster_disposition"])
    )
    summary.to_csv(args.output_dir / "ALL_505_RESOLUTION_SUMMARY.csv", index=False)
    candidates = table.loc[table.analysis_tier.eq("secondary_candidate")].copy()
    candidates.to_csv(args.output_dir / "SECONDARY_CANDIDATES_REQUIRING_MENTOR_CONFIRMATION.csv", index=False)

    manifest = {
        "source_dictionary": str(args.dictionary.resolve()),
        "n_clusters": int(len(table)),
        "primary_working_clusters": int(primary.sum()),
        "secondary_candidates": int(clear_candidate.sum()),
        "excluded_or_unresolved": int((~primary & ~clear_candidate).sum()),
        "status": "all_clusters_dispositioned_pending_mentor_review",
    }
    (args.output_dir / "ALL_505_RESOLUTION_MANIFEST.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    (args.output_dir / "README.md").write_text(
        "# Full native-cluster disposition ledger\n\n"
        "All 505 native clusters have a reproducible disposition. `primary_working_state` clusters are the only clusters "
        "that may enter preplanned donor/sample-level state analyses. `secondary_candidate` clusters are marker-supported "
        "but await mentor confirmation. All other rows are explicitly excluded or unresolved rather than assigned a false "
        "fine-state label.\n",
        encoding="utf-8",
    )
    print(summary.to_string(index=False))
    print(f"Wrote complete resolution ledger: {args.output_dir}")


if __name__ == "__main__":
    main()
