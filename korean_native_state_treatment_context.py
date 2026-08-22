#!/usr/bin/env python3
"""Describe patient-paired Korean native-state changes during treatment.

This is a longitudinal context analysis. It is not a response classifier and
does not substitute a related Korean state for an absent cross-cohort state.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


TARGETS = [
    "Endothelial::PLVAP_CA4_capillary_endothelial",
    "Endothelial::activated_endothelial",
    "Fibroblast::CFD_DCN_resident_fibroblast",
    "Fibroblast::PDPN_CTSK_inflammatory_CAF_candidate",
    "Myeloid::C1QC_APOC1_macrophage",
    "Myeloid::CD1C_antigen_presenting_myeloid",
    "T/NK::activated_cytotoxic_T",
    "T/NK::cytotoxic_NK",
]


def find_column(frame: pd.DataFrame, candidates: list[str]) -> str | None:
    lookup = {str(column).lower(): str(column) for column in frame.columns}
    return next((lookup[name.lower()] for name in candidates if name.lower() in lookup), None)


def patient_and_timepoint(sample: str) -> tuple[str | None, str | None]:
    text = str(sample)
    patient = re.match(r"^(E\d+)", text)
    if text.endswith("_B") or "_BL" in text:
        return (patient.group(1) if patient else None), "Baseline"
    if "_F1" in text:
        return (patient.group(1) if patient else None), "FU1"
    if "_F2" in text:
        return (patient.group(1) if patient else None), "FU2"
    return (patient.group(1) if patient else None), None


def clinical_table(path: Path) -> pd.DataFrame:
    source = pd.read_csv(path)
    patient_col = find_column(source, ["patient", "patient_id"])
    response_col = find_column(source, ["best_overall_response"])
    binary_col = find_column(source, ["response_binary"])
    if patient_col is None or response_col is None:
        raise ValueError("Clinical table requires patient/patient_id and best_overall_response.")
    result = source[[patient_col, response_col]].copy().rename(columns={patient_col: "patient_id", response_col: "best_overall_response"})
    result["patient_id"] = result.patient_id.astype(str)
    response = result.best_overall_response.astype(str).str.upper().str.strip()
    result["response_binary"] = np.where(response.isin(["CR", "PR"]), 1, np.where(response.isin(["SD", "PD"]), 0, np.nan))
    if binary_col is not None:
        declared = pd.to_numeric(source[binary_col], errors="raise").to_numpy()
        auditable = pd.Series(result["response_binary"]).notna().to_numpy()
        if not np.array_equal(declared[auditable], result.loc[auditable, "response_binary"].astype(int).to_numpy()):
            raise ValueError("Clinical response_binary disagrees with CR/PR versus SD/PD coding.")
    result = result.dropna(subset=["response_binary"]).drop_duplicates("patient_id")
    result["response_binary"] = result.response_binary.astype(int)
    if result.response_binary.nunique() != 2:
        raise ValueError("Audited Korean clinical data need both responder and non-responder labels.")
    return result


def permutation_difference(delta: np.ndarray, response: np.ndarray, *, seed: int, n_permutations: int) -> tuple[float, float]:
    observed = float(delta[response == 1].mean() - delta[response == 0].mean())
    rng = np.random.default_rng(seed)
    extreme = 0
    for _ in range(n_permutations):
        shuffled = rng.permutation(response)
        value = float(delta[shuffled == 1].mean() - delta[shuffled == 0].mean())
        extreme += abs(value) >= abs(observed)
    return observed, float((extreme + 1) / (n_permutations + 1))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--composition", required=True, type=Path)
    parser.add_argument("--clinical", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--permutations", type=int, default=10000)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    composition = pd.read_csv(args.composition)
    required = {"dataset_id", "sample", "compartment", "reviewed_state_label", "state_fraction", "zero_filled_absent_state"}
    if missing := required - set(composition.columns):
        raise ValueError(f"Composition is missing required columns: {sorted(missing)}")
    korea = composition.loc[composition.dataset_id.astype(str).eq("korea_kim2022")].copy()
    korea["state_id"] = korea.compartment.astype(str) + "::" + korea.reviewed_state_label.astype(str)
    patient_time = korea.sample.map(patient_and_timepoint)
    korea[["patient_id", "timepoint"]] = pd.DataFrame(patient_time.tolist(), index=korea.index)
    korea = korea.dropna(subset=["patient_id", "timepoint"])
    clinical = clinical_table(args.clinical)
    availability = pd.DataFrame({"state_id": TARGETS})
    availability["available_in_korean_native_dictionary"] = availability.state_id.isin(korea.state_id.unique())
    availability["interpretation"] = np.where(
        availability.available_in_korean_native_dictionary,
        "Eligible for longitudinal descriptive analysis if paired samples exist.",
        "Not tested: no equivalent Korean native state was curated; no substitution is made.",
    )
    availability.to_csv(args.output_dir / "KOREAN_TARGET_STATE_AVAILABILITY.csv", index=False)

    rows = []
    for state_id, state in korea.groupby("state_id", observed=True):
        values = state.groupby(["patient_id", "timepoint"], observed=True).state_fraction.mean().unstack("timepoint")
        for followup, role in [("FU1", "primary"), ("FU2", "exploratory")]:
            if "Baseline" not in values or followup not in values:
                continue
            paired = values[["Baseline", followup]].dropna().merge(clinical, left_index=True, right_on="patient_id", how="inner")
            if paired.empty:
                continue
            delta = (paired[followup] - paired["Baseline"]).to_numpy(float)
            response = paired.response_binary.to_numpy(int)
            row = {
                "state_id": state_id, "followup": followup, "analysis_role": role,
                "n_pairs": int(len(paired)), "n_responders": int(response.sum()),
                "n_nonresponders": int((response == 0).sum()), "mean_within_patient_change": float(delta.mean()),
                "responder_mean_change": float(delta[response == 1].mean()) if (response == 1).any() else np.nan,
                "nonresponder_mean_change": float(delta[response == 0].mean()) if (response == 0).any() else np.nan,
            }
            if response.sum() >= 5 and (response == 0).sum() >= 5:
                effect, p_value = permutation_difference(delta, response, seed=711 + len(rows), n_permutations=args.permutations)
                row.update({"response_group_difference": effect, "response_group_permutation_p_value": p_value, "test_status": "tested"})
            else:
                row.update({"response_group_difference": np.nan, "response_group_permutation_p_value": np.nan, "test_status": "insufficient_response_groups"})
            rows.append(row)
    results = pd.DataFrame(rows)
    if results.empty:
        raise RuntimeError("No Korean baseline-follow-up state pairs were available.")
    results["response_group_fdr_bh"] = np.nan
    for followup in results.followup.unique():
        index = results.followup.eq(followup) & results.response_group_permutation_p_value.notna()
        if index.any():
            results.loc[index, "response_group_fdr_bh"] = multipletests(results.loc[index, "response_group_permutation_p_value"], method="fdr_bh")[1]
    results["claim_boundary"] = "Patient-paired longitudinal context; not a response-prediction model."
    results.to_csv(args.output_dir / "KOREAN_NATIVE_STATE_TREATMENT_CONTEXT.csv", index=False)
    primary = results.loc[results.followup.eq("FU1")].sort_values("response_group_fdr_bh", na_position="last")
    fig, ax = plt.subplots(figsize=(8.5, max(3.5, 0.38 * len(primary))))
    ax.scatter(primary["mean_within_patient_change"], np.arange(len(primary)), color="#1f4e79")
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_yticks(np.arange(len(primary)), primary.state_id)
    ax.set_xlabel("Mean baseline-to-FU1 within-patient change")
    ax.set_title("Korean serial-treatment native-state context")
    fig.tight_layout()
    fig.savefig(args.output_dir / "KOREAN_NATIVE_STATE_TREATMENT_CONTEXT.png", dpi=300, bbox_inches="tight")
    fig.savefig(args.output_dir / "KOREAN_NATIVE_STATE_TREATMENT_CONTEXT.pdf", bbox_inches="tight")
    plt.close(fig)
    print(results.to_string(index=False))


if __name__ == "__main__":
    main()
