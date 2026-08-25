#!/usr/bin/env python3
"""Audit and describe available Korean pathology linkage for frozen native states."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

from korean_native_state_treatment_context import patient_and_timepoint


FIELDS = {
    "lauren_subtype": ("lauren", "lauren_subtype", "histology"),
    "stage": ("stage", "pathologic_stage", "ajcc_stage"),
    "tumour_regression_grade": ("tumor_regression_grade", "tumour_regression_grade", "trg"),
    "metastatic_site": ("metastatic_site", "metastasis_site", "site"),
    "treatment_regimen": ("treatment_regimen", "regimen", "therapy"),
}


def find_column(frame: pd.DataFrame, aliases: tuple[str, ...]) -> str | None:
    lookup = {str(column).lower(): str(column) for column in frame.columns}
    return next((lookup[item.lower()] for item in aliases if item.lower() in lookup), None)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--composition", required=True, type=Path)
    parser.add_argument("--clinical", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    composition = pd.read_csv(args.composition)
    korea = composition.loc[composition.dataset_id.astype(str).eq("korea_kim2022")].copy()
    parsed = korea["sample"].map(patient_and_timepoint)
    korea[["patient_id", "timepoint"]] = pd.DataFrame(parsed.tolist(), index=korea.index)
    korea = korea.loc[korea.timepoint.eq("Baseline")].copy()
    korea["state_id"] = korea.compartment.astype(str) + "::" + korea.reviewed_state_label.astype(str)
    values = korea.pivot_table(index="patient_id", columns="state_id", values="state_fraction", aggfunc="mean", fill_value=0)

    clinical = pd.read_csv(args.clinical)
    patient = find_column(clinical, ("patient", "patient_id"))
    if patient is None:
        raise ValueError("Clinical table has no patient/patient_id column.")
    clinical = clinical.copy()
    clinical["patient_id"] = clinical[patient].astype(str)
    clinical = clinical.drop_duplicates("patient_id").set_index("patient_id")
    data = values.join(clinical, how="inner")
    audit, tests = [], []
    for label, aliases in FIELDS.items():
        column = find_column(data, aliases)
        if column is None:
            audit.append({"pathology_field": label, "status": "not_available", "source_column": ""})
            continue
        groups = data[column].dropna().astype(str)
        counts = groups.value_counts()
        eligible = counts[counts >= 3]
        audit.append({"pathology_field": label, "status": "descriptive_only" if len(eligible) < 2 else "exploratory_testable", "source_column": column, "n_patients": int(groups.index.nunique()), "groups": ";".join(f"{key}:{value}" for key, value in counts.items())})
        if len(eligible) < 2:
            continue
        for state_id in values.columns:
            samples = [data.loc[groups.index[groups.eq(group)], state_id].to_numpy(float) for group in eligible.index]
            if any(len(sample) < 3 for sample in samples):
                continue
            p_value = stats.kruskal(*samples).pvalue
            tests.append({"pathology_field": label, "source_column": column, "state_id": state_id, "n_groups": len(samples), "p_value": p_value})
    audit_table, test_table = pd.DataFrame(audit), pd.DataFrame(tests)
    if not test_table.empty:
        test_table["fdr_bh"] = multipletests(test_table.p_value, method="fdr_bh")[1]
    test_table["claim_boundary"] = "Exploratory baseline association with available Korean clinical/pathology metadata; no causal or predictive claim."
    audit_table.to_csv(args.output_dir / "KOREAN_PATHOLOGY_METADATA_AUDIT.csv", index=False)
    test_table.to_csv(args.output_dir / "KOREAN_BASELINE_STATE_PATHOLOGY_ASSOCIATIONS.csv", index=False)
    print(audit_table.to_string(index=False))


if __name__ == "__main__":
    main()
