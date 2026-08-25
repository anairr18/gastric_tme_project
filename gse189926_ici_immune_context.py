#!/usr/bin/env python3
"""Evaluate frozen immune-state signatures in the public GSE189926 ICI cohort.

The cohort is CD45 selected. Consequently, this is an immune-context analysis
only, with sample/patient as the unit of inference. It does not assess CAF or
endothelial signatures and never trains a response predictor.
"""

from __future__ import annotations

import argparse
import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from frozen_bulk_clinical_validation import signature_panels


IMMUNE_PREFIXES = ("Myeloid::", "T/NK::")
RESPONDER = {"CR", "PR"}
NONRESPONDER = {"PD"}


def fdr(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame.copy()
    result["fdr_bh"] = np.nan
    tested = result["status"].eq("tested") & result["p_value"].notna()
    if tested.any():
        result.loc[tested, "fdr_bh"] = multipletests(result.loc[tested, "p_value"], method="fdr_bh")[1]
    return result


def response_group(value: object) -> str | None:
    value = str(value).strip().upper()
    if value in RESPONDER:
        return "responder"
    if value in NONRESPONDER:
        return "nonresponder"
    return None


def group_test(data: pd.DataFrame, state_id: str, value_column: str) -> dict[str, object]:
    frame = data[["response_group", value_column]].dropna()
    positive = frame.loc[frame["response_group"].eq("responder"), value_column].to_numpy(float)
    negative = frame.loc[frame["response_group"].eq("nonresponder"), value_column].to_numpy(float)
    row: dict[str, object] = {
        "state_id": state_id,
        "n_responders": int(len(positive)),
        "n_nonresponders": int(len(negative)),
    }
    if len(positive) < 3 or len(negative) < 3:
        row["status"] = "insufficient_response_groups"
        return row
    test = stats.mannwhitneyu(positive, negative, alternative="two-sided")
    row.update({
        "status": "tested",
        "mean_responder_score": float(positive.mean()),
        "mean_nonresponder_score": float(negative.mean()),
        "difference_responder_minus_nonresponder": float(positive.mean() - negative.mean()),
        "p_value": float(test.pvalue),
    })
    return row


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expression", required=True, type=Path)
    parser.add_argument("--clinical", required=True, type=Path)
    parser.add_argument("--reference-profiles", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    expression = pd.read_csv(args.expression).set_index("gene")
    expression.index = expression.index.astype(str).str.upper()
    expression = expression.groupby(level=0).sum()
    clinical = pd.read_csv(args.clinical)
    required = {"sample_id", "patient_id", "timepoint", "response", "assay_scope"}
    if missing := required - set(clinical.columns):
        raise ValueError(f"GSE189926 clinical input lacks: {sorted(missing)}")
    clinical["sample_id"] = clinical["sample_id"].astype(str)
    shared = [sample for sample in clinical["sample_id"] if sample in expression.columns]
    if len(shared) < 10:
        raise ValueError(f"Only {len(shared)} count pseudobulks matched clinical metadata.")
    clinical = clinical.set_index("sample_id").loc[shared].reset_index()
    expression = expression.loc[:, shared]

    reference = pd.read_csv(args.reference_profiles, compression="infer")
    all_panels, panel_table = signature_panels(reference, top_n=12)
    panels = {state: genes for state, genes in all_panels.items() if state.startswith(IMMUNE_PREFIXES)}
    if not panels:
        raise ValueError("Frozen reference profiles did not yield any immune state panels.")
    ranks = expression.rank(axis=0, pct=True, method="average")
    scores, coverage = {}, []
    for state_id, markers in panels.items():
        available = [gene for gene in markers if gene in ranks.index]
        coverage.append({
            "state_id": state_id,
            "n_requested": len(markers),
            "n_available": len(available),
            "coverage_fraction": len(available) / max(1, len(markers)),
            "available_genes": ";".join(available),
        })
        if len(available) >= max(3, int(np.ceil(0.6 * len(markers)))):
            scores[state_id] = ranks.loc[available].mean(axis=0)
    if not scores:
        raise RuntimeError("No frozen immune-state panel met the 60% gene-coverage gate.")

    sample_scores = pd.DataFrame(scores, index=shared).reset_index(names="sample_id")
    sample_scores = sample_scores.merge(clinical, on="sample_id", how="inner", validate="one_to_one")
    sample_scores["response_group"] = sample_scores["response"].map(response_group)
    sample_scores.to_csv(args.output_dir / "GSE189926_FROZEN_IMMUNE_STATE_SAMPLE_SCORES.csv", index=False)
    panel_table[panel_table["state_id"].isin(panels)].to_csv(args.output_dir / "GSE189926_FROZEN_IMMUNE_STATE_PANELS.csv", index=False)
    pd.DataFrame(coverage).to_csv(args.output_dir / "GSE189926_FROZEN_IMMUNE_STATE_COVERAGE.csv", index=False)

    state_columns = list(scores)
    patient_time = (
        sample_scores.groupby(["patient_id", "timepoint", "response", "response_group", "assay_scope"], observed=True)[state_columns]
        .mean()
        .reset_index()
    )
    patient_time.to_csv(args.output_dir / "GSE189926_FROZEN_IMMUNE_STATE_PATIENT_TIMEPOINT_SCORES.csv", index=False)

    baseline = patient_time[patient_time["timepoint"].eq("baseline")].copy()
    baseline_rows = [group_test(baseline, state_id, state_id) for state_id in state_columns]
    baseline_results = fdr(pd.DataFrame(baseline_rows))
    baseline_results["analysis"] = "baseline responder versus PD nonresponder"
    baseline_results.to_csv(args.output_dir / "GSE189926_BASELINE_RESPONSE_ASSOCIATIONS.csv", index=False)

    pivot = patient_time.pivot_table(index=["patient_id", "response", "response_group"], columns="timepoint", values=state_columns, aggfunc="mean")
    paired_rows, interaction_rows = [], []
    for state_id in state_columns:
        if state_id not in pivot.columns.get_level_values(0):
            continue
        wide = pivot[state_id].copy()
        if not {"baseline", "on_treatment"}.issubset(wide.columns):
            paired_rows.append({"state_id": state_id, "status": "no_auditable_paired_timepoints"})
            interaction_rows.append({"state_id": state_id, "status": "no_auditable_paired_timepoints"})
            continue
        paired = wide.dropna(subset=["baseline", "on_treatment"]).reset_index()
        paired["delta_on_treatment_minus_baseline"] = paired["on_treatment"] - paired["baseline"]
        if len(paired) < 4:
            paired_rows.append({"state_id": state_id, "status": "insufficient_pairs", "n_pairs": len(paired)})
        else:
            test = stats.wilcoxon(paired["delta_on_treatment_minus_baseline"], alternative="two-sided", method="auto")
            paired_rows.append({
                "state_id": state_id,
                "status": "tested",
                "n_pairs": int(len(paired)),
                "mean_delta_on_treatment_minus_baseline": float(paired["delta_on_treatment_minus_baseline"].mean()),
                "p_value": float(test.pvalue),
            })
        interaction_rows.append(group_test(paired, state_id, "delta_on_treatment_minus_baseline"))
    paired_results = fdr(pd.DataFrame(paired_rows))
    paired_results["analysis"] = "within-patient paired baseline to on-treatment change"
    paired_results.to_csv(args.output_dir / "GSE189926_PAIRED_TREATMENT_CHANGES.csv", index=False)
    interaction_results = fdr(pd.DataFrame(interaction_rows))
    interaction_results["analysis"] = "response-group difference in within-patient change"
    interaction_results.to_csv(args.output_dir / "GSE189926_RESPONSE_BY_TREATMENT_CHANGE_ASSOCIATIONS.csv", index=False)

    paired_patient_count = int(
        patient_time.groupby("patient_id", observed=True)["timepoint"]
        .agg(lambda values: {str(value) for value in values})
        .map(lambda values: {"baseline", "on_treatment"}.issubset(values))
        .sum()
    )

    manifest = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "n_samples": int(sample_scores["sample_id"].nunique()),
        "n_patients": int(patient_time["patient_id"].nunique()),
        "n_baseline_samples": int(baseline["patient_id"].nunique()),
        "n_auditable_pairs": paired_patient_count,
        "assay_scope": "CD45-selected immune cells",
        "claim_boundary": (
            "Frozen immune-state signature associations in a small, CD45-selected, "
            "single-arm chemoimmunotherapy cohort. Neither response prediction nor "
            "CAF/endothelial inference is supported."
        ),
    }
    (args.output_dir / "GSE189926_ICI_CONTEXT_MANIFEST.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
