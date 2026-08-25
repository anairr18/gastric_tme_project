#!/usr/bin/env python3
"""Test frozen native-state signatures in an external clinical bulk cohort.

The script never trains a predictor.  It only scores marker panels derived from
the already frozen native-state reference profiles and tests prespecified
clinical metadata supplied with the external cohort.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests


TARGETS = (
    "Endothelial::PLVAP_CA4_capillary_endothelial",
    "Fibroblast::CFD_DCN_resident_fibroblast",
    "Fibroblast::PDPN_CTSK_inflammatory_CAF_candidate",
    "Myeloid::IL1B_FCN1_inflammatory_myeloid",
    "Myeloid::C1QC_APOC1_macrophage",
    "T/NK::GZMK_effector_memory_T",
)
EXCLUDED_PREFIXES = ("RPL", "RPS", "MT-", "MALAT1", "HBA", "HBB")


def find_column(frame: pd.DataFrame, candidates: tuple[str, ...]) -> str | None:
    lookup = {str(column).lower(): str(column) for column in frame.columns}
    return next((lookup[item.lower()] for item in candidates if item.lower() in lookup), None)


def read_expression(path: Path) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t" if path.suffix.lower() == ".tsv" else ",")
    if table.shape[1] < 3:
        raise ValueError("Expression table must be gene-by-sample with at least two samples.")
    gene_column = find_column(table, ("gene", "gene_symbol", "symbol", "hugo_symbol")) or str(table.columns[0])
    expression = table.set_index(gene_column).apply(pd.to_numeric, errors="coerce")
    expression.index = expression.index.astype(str).str.upper()
    expression = expression.groupby(level=0).mean().dropna(axis=1, how="all")
    if expression.empty:
        raise ValueError("Expression table contained no numeric sample columns after parsing.")
    return expression


def signature_panels(reference: pd.DataFrame, top_n: int) -> tuple[dict[str, list[str]], pd.DataFrame]:
    required = {"state_id", "gene", "mean_log_normalized_expression"}
    if missing := required - set(reference.columns):
        raise ValueError(f"Reference profile is missing columns: {sorted(missing)}")
    profile = reference.copy()
    profile["gene"] = profile["gene"].astype(str).str.upper()
    wide = profile.pivot_table(index="state_id", columns="gene", values="mean_log_normalized_expression", aggfunc="mean").fillna(0.0)
    panels, rows = {}, []
    for state_id in TARGETS:
        if state_id not in wide.index:
            rows.append({"state_id": state_id, "status": "missing_from_reference"})
            continue
        compartment = state_id.split("::", 1)[0]
        peers = [state for state in wide.index if state.startswith(compartment + "::") and state != state_id]
        if not peers:
            rows.append({"state_id": state_id, "status": "no_within_compartment_peer_states"})
            continue
        contrast = wide.loc[state_id] - wide.loc[peers].mean(axis=0)
        genes = [
            gene for gene in contrast.sort_values(ascending=False).index
            if contrast[gene] > 0 and not gene.startswith(EXCLUDED_PREFIXES)
        ][:top_n]
        panels[state_id] = genes
        rows.append({"state_id": state_id, "status": "derived", "n_markers": len(genes), "markers": ";".join(genes)})
    return panels, pd.DataFrame(rows)


def response_binary(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    if numeric.dropna().isin([0, 1]).all() and numeric.notna().any():
        return numeric
    text = values.astype(str).str.upper().str.strip()
    return pd.Series(np.where(text.isin(["CR", "PR", "RESPONDER", "RESPONSE", "1"]), 1,
                              np.where(text.isin(["SD", "PD", "NONRESPONDER", "NON-RESPONDER", "0"]), 0, np.nan)), index=values.index)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expression", required=True, type=Path)
    parser.add_argument("--clinical", required=True, type=Path)
    parser.add_argument("--reference-profiles", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--mode", choices=("ici_response", "molecular_subtype"), required=True)
    parser.add_argument("--top-markers", type=int, default=12)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    expression = read_expression(args.expression)
    clinical = pd.read_csv(args.clinical)
    sample_column = find_column(clinical, ("sample_id", "sample", "patient_id", "patient", "barcode"))
    if sample_column is None:
        raise ValueError("Clinical table needs an exact sample_id/sample/patient_id column matching expression columns.")
    clinical = clinical.copy()
    clinical["sample_id"] = clinical[sample_column].astype(str)
    clinical = clinical.drop_duplicates("sample_id")
    shared = [sample for sample in expression.columns.astype(str) if sample in set(clinical.sample_id)]
    if len(shared) < 10:
        raise ValueError(f"Only {len(shared)} expression samples matched clinical metadata; require at least 10 exact IDs.")
    expression = expression.loc[:, shared]
    clinical = clinical.set_index("sample_id").loc[shared].reset_index()

    reference = pd.read_csv(args.reference_profiles, compression="infer")
    panels, panel_table = signature_panels(reference, args.top_markers)
    ranks = expression.rank(axis=0, pct=True, method="average")
    scores, coverage = {}, []
    for state_id, markers in panels.items():
        present = [gene for gene in markers if gene in ranks.index]
        coverage.append({"state_id": state_id, "n_requested": len(markers), "n_available": len(present), "coverage_fraction": len(present) / max(len(markers), 1), "available_genes": ";".join(present)})
        if len(present) >= max(3, int(np.ceil(0.6 * len(markers)))):
            scores[state_id] = ranks.loc[present].mean(axis=0)
    score_table = pd.DataFrame(scores, index=shared).reset_index(names="sample_id").merge(clinical, on="sample_id", how="left", validate="one_to_one")
    if score_table.shape[1] <= clinical.shape[1]:
        raise RuntimeError("No frozen state panel passed the 60% coverage requirement.")

    rows = []
    if args.mode == "ici_response":
        response_column = find_column(clinical, ("best_overall_response", "response", "recist", "response_binary"))
        if response_column is None:
            raise ValueError("ICI response mode requires best_overall_response/response/recist/response_binary.")
        response = response_binary(score_table[response_column])
        for state_id in scores:
            values = pd.to_numeric(score_table[state_id], errors="coerce")
            eligible = values.notna() & response.notna()
            positive, negative = values[eligible & response.eq(1)], values[eligible & response.eq(0)]
            if len(positive) < 5 or len(negative) < 5:
                rows.append({"state_id": state_id, "status": "insufficient_response_groups", "n_responders": len(positive), "n_nonresponders": len(negative)})
                continue
            test = stats.mannwhitneyu(positive, negative, alternative="two-sided")
            rows.append({"state_id": state_id, "status": "tested", "n_responders": len(positive), "n_nonresponders": len(negative), "mean_responder_score": positive.mean(), "mean_nonresponder_score": negative.mean(), "difference": positive.mean() - negative.mean(), "p_value": test.pvalue})
        boundary = "Frozen signature association with observed ICI response; no model training, AUROC, or predictive claim."
    else:
        subtype_column = find_column(clinical, ("molecular_subtype", "tcga_subtype", "subtype"))
        if subtype_column is None:
            raise ValueError("Molecular-subtype mode requires molecular_subtype/tcga_subtype/subtype.")
        for state_id in scores:
            values = pd.to_numeric(score_table[state_id], errors="coerce")
            groups = [group[state_id].dropna().to_numpy(float) for _, group in score_table.assign(**{state_id: values}).groupby(subtype_column, observed=True) if len(group[state_id].dropna()) >= 5]
            if len(groups) < 3:
                rows.append({"state_id": state_id, "status": "insufficient_subtype_groups", "n_eligible_subtypes": len(groups)})
                continue
            test = stats.kruskal(*groups)
            rows.append({"state_id": state_id, "status": "tested", "n_eligible_subtypes": len(groups), "p_value": test.pvalue})
        boundary = "Frozen signature association with supplied molecular subtype labels; association only, not subtype prediction."

    results = pd.DataFrame(rows)
    results["fdr_bh"] = np.nan
    tested = results.status.eq("tested") & results.p_value.notna()
    if tested.any():
        results.loc[tested, "fdr_bh"] = multipletests(results.loc[tested, "p_value"], method="fdr_bh")[1]
    results["claim_boundary"] = boundary
    panel_table.to_csv(args.output_dir / "FROZEN_STATE_MARKER_PANELS.csv", index=False)
    pd.DataFrame(coverage).to_csv(args.output_dir / "FROZEN_STATE_PANEL_COVERAGE.csv", index=False)
    score_table.to_csv(args.output_dir / "FROZEN_STATE_CLINICAL_SCORES.csv", index=False)
    results.to_csv(args.output_dir / "FROZEN_STATE_CLINICAL_ASSOCIATIONS.csv", index=False)
    print(results.to_string(index=False))


if __name__ == "__main__":
    main()
