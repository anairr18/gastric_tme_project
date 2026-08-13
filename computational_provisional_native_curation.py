#!/usr/bin/env python3
"""Create an auditable provisional state dictionary from native cluster markers.

This is deliberately conservative.  It accepts only clusters that passed the
existing lineage audit, labels them using transparent marker panels, and leaves
clusters unresolved when the evidence is weak or mixed.  The output is a
working dictionary for planned analyses, not a replacement for domain-expert
review or a source of final biological claims.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

import pandas as pd


# Each label has a small, interpretable marker program.  These are used only on
# within-cohort raw-expression marker ranks produced by the native discovery.
STATE_PANELS: dict[str, dict[str, set[str]]] = {
    "B/Plasma": {
        "antigen_presenting_B": {"MS4A1", "CD79A", "CD79B", "CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD37", "CD22", "FCER2"},
        "activated_B": {"MS4A1", "CD79A", "CD74", "CD69", "NR4A1", "NR4A2", "JUNB", "FOS", "CXCR4", "GPR18"},
        "IgA_plasma": {"MZB1", "XBP1", "DERL3", "JCHAIN", "IGHA1", "IGHA2", "SSR4", "SSR3", "PRDM1", "SDC1", "TNFRSF17"},
        "IgG_plasma": {"MZB1", "XBP1", "DERL3", "IGHG1", "IGHG2", "IGHG3", "IGHG4", "JCHAIN", "SSR4", "PRDM1"},
    },
    "Endothelial": {
        "PLVAP_CA4_capillary_endothelial": {"PLVAP", "CA4", "RGCC", "GNG11", "KDR", "EMCN", "FABP5", "RBP7", "AQP1", "CD36", "SLC14A1"},
        "ACKR1_venular_endothelial": {"ACKR1", "SELP", "MADCAM1", "VCAM1", "ICAM1", "CLDN5", "MMRN1", "CYP1B1", "IL1R1"},
        "lymphatic_endothelial": {"CCL21", "PROX1", "LYVE1", "PDPN", "FLT4", "CCL14", "MRC1"},
        "activated_endothelial": {"ICAM1", "VCAM1", "SELE", "SELP", "JUN", "FOS", "DUSP1", "GADD45B"},
    },
    "Epithelial": {
        "gastric_mucous_epithelial": {"MUC5AC", "TFF1", "TFF2", "GKN1", "GKN2", "PGC", "SPINK1", "CA2", "PSCA", "MSMB"},
        "intestinal_metaplastic_epithelial": {"TFF3", "REG4", "FABP1", "FABP2", "KRT20", "MUC2", "MUC17", "ANPEP", "PHGR1", "OLFM4", "LGALS4", "DMBT1", "PIGR", "ALDOB"},
        "glandular_secretory_epithelial": {"MUC6", "AQP5", "REG1A", "REG3A", "LTF", "BPIFB1", "ZG16B", "PRR4", "WFDC2", "LCN2"},
        "epithelial_stress_or_tumour_candidate": {"KRT7", "KRT17", "KRT19", "KRT8", "EPCAM", "CLDN3", "CLDN4", "S100P", "MSLN", "GDF15", "CEACAM6", "TNFRSF12A"},
    },
    "Fibroblast": {
        "CFD_DCN_resident_fibroblast": {"CFD", "DCN", "LUM", "C7", "C1R", "C1S", "FBLN1", "FBLN2", "PTGDS", "ABCA8", "ADH1B", "COL15A1", "DPT"},
        "PDPN_CTSK_inflammatory_CAF_candidate": {"PDPN", "CTSK", "CXCL12", "CXCL14", "CCL2", "CCL11", "IL6", "WNT5A", "F3", "PTGS1", "COL3A1", "POSTN"},
        "POSTN_CTHRC1_matrix_CAF_candidate": {"POSTN", "CTHRC1", "COL1A1", "COL1A2", "COL3A1", "THBS2", "FAP", "BGN", "SPARC", "COL8A1", "FN1", "INHBA"},
    },
    "Myeloid": {
        "IL1B_FCN1_inflammatory_myeloid": {"IL1B", "IL1RN", "FCN1", "S100A8", "S100A9", "S100A12", "CXCL8", "CXCL2", "CXCL3", "C15orf48", "SOD2", "PLAUR", "G0S2"},
        "C1QC_APOC1_macrophage": {"C1QA", "C1QB", "C1QC", "APOE", "APOC1", "LIPA", "LILRB1", "MSR1", "SEPP1", "LGMN", "CTSD"},
        "SPP1_TREM2_macrophage_candidate": {"SPP1", "TREM2", "GPNMB", "APOC1", "APOE", "HAMP", "MARCO", "LGALS3", "FN1", "ACP5", "FCGR2B"},
        "CD1C_antigen_presenting_myeloid": {"CD1C", "CD1E", "CLEC10A", "FCER1A", "HLA-DRA", "HLA-DPA1", "HLA-DPB1", "CD74", "FSCN1", "LAMP3"},
        "mast_cell": {"TPSAB1", "TPSB2", "CPA3", "KIT", "HDC", "HPGDS", "MS4A2"},
        "pDC_candidate": {"GZMB", "GZMB", "LILRA4", "GZMB", "TCF4", "IL3RA", "GZMB", "PLAC8", "GZMB", "SCT"},
    },
    "Stromal": {
        "contractile_perivascular_stromal": {"RGS5", "MCAM", "CSPG4", "ACTA2", "TAGLN", "MYH11", "MYL9", "CNN1", "RERGL", "PLN", "PDGFRB", "NOTCH3"},
    },
    "T/NK": {
        "GZMK_effector_memory_T": {"GZMK", "GZMA", "CCL5", "NKG7", "CXCR3", "TRAC", "CD3D", "CD3E", "DUSP2", "CMC1"},
        "cytotoxic_NK": {"NKG7", "GNLY", "GZMB", "GZMH", "PRF1", "FGFBP2", "KLRD1", "FCGR3A", "CTSW", "XCL1", "XCL2", "SPON2"},
        "activated_CD4_Treg_candidate": {"FOXP3", "IL2RA", "CTLA4", "TIGIT", "TNFRSF4", "TNFRSF18", "ICOS", "CCR7", "MAGEH1"},
        "CD4_memory_T": {"IL7R", "LTB", "CCR7", "KLRB1", "MAL", "LST1", "S100A4", "IL32", "CD3D", "CD3E"},
        "activated_cytotoxic_T": {"CD8A", "CCL5", "GZMB", "HOPX", "KLRC1", "NKG7", "PRF1", "IFNG", "CD69"},
    },
}

MIN_CELLS = 50
MIN_SAMPLES = 3


def key(row: pd.Series) -> tuple[str, str, int]:
    return (str(row.dataset_id), str(row.broad_label), int(row.native_cluster))


def marker_sets(markers: pd.DataFrame) -> dict[tuple[str, str, int], list[str]]:
    required = {"dataset_id", "broad_label", "native_cluster", "names"}
    missing = required - set(markers.columns)
    if missing:
        raise ValueError(f"Marker table missing columns: {sorted(missing)}")
    rank_col = "rank" if "rank" in markers.columns else None
    out: dict[tuple[str, str, int], list[str]] = {}
    for group_key, group in markers.groupby(["dataset_id", "broad_label", "native_cluster"], observed=True):
        if rank_col:
            group = group.sort_values(rank_col)
        names = (
            group["names"].astype(str).str.upper().str.strip().head(30).tolist()
        )
        out[(str(group_key[0]), str(group_key[1]), int(group_key[2]))] = names
    return out


def label_cluster(broad_label: str, markers: list[str]) -> tuple[str, list[str], int, int]:
    panels = STATE_PANELS.get(broad_label, {})
    scores: list[tuple[int, str, list[str]]] = []
    marker_set = set(markers)
    for state, panel in panels.items():
        hits = sorted(marker_set & panel)
        scores.append((len(hits), state, hits))
    scores.sort(key=lambda value: (-value[0], value[1]))
    if not scores or scores[0][0] < 2:
        return "unresolved", [], scores[0][0] if scores else 0, scores[1][0] if len(scores) > 1 else 0
    best_score, best_state, best_hits = scores[0]
    runner_up = scores[1][0] if len(scores) > 1 else 0
    # A tied program often reflects a genuine transitional phenotype or an
    # overbroad panel; retain rather than invent a false fine-state label.
    if best_score == runner_up and best_score < 4:
        return "unresolved", best_hits, best_score, runner_up
    return best_state, best_hits, best_score, runner_up


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--markers", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    audit = pd.read_csv(args.audit)
    by_cluster = marker_sets(pd.read_csv(args.markers))
    rows = []
    for _, row in audit.iterrows():
        cluster_markers = by_cluster.get(key(row), [])
        state, hits, score, runner_up = label_cluster(str(row.broad_label), cluster_markers)
        eligible = row.curation_recommendation == "eligible_after_manual_state_annotation"
        support = int(row.n_cells) >= MIN_CELLS and int(row.n_samples) >= MIN_SAMPLES
        if not eligible:
            decision = "exclude"
            reason = str(row.curation_recommendation)
        elif not support:
            decision = "exclude"
            reason = "below_provisional_minimum_support"
        elif state == "unresolved":
            decision = "exclude"
            reason = "insufficient_or_ambiguous_state_marker_support"
        else:
            decision = "include_provisional"
            reason = "transparent_marker_panel_support"
        rows.append({
            **row.to_dict(),
            "top_30_native_markers": ";".join(cluster_markers),
            "provisional_state_name": state,
            "matched_state_markers": ";".join(hits),
            "state_marker_score": score,
            "runner_up_state_marker_score": runner_up,
            "provisional_decision": decision,
            "provisional_rationale": reason,
            "reviewer": "Aadi Nair; computational proposal generated by Codex",
            "review_status": "provisional_pending_mentor_review",
            "claim_boundary": (
                "Computational working label based on within-cohort raw-expression marker ranks. "
                "Not expert-curated, not a final state taxonomy, and not evidence of disease association."
            ),
        })
    dictionary = pd.DataFrame(rows)
    dictionary.to_csv(args.output_dir / "COMPUTATIONAL_PROVISIONAL_STATE_DICTIONARY.csv", index=False)

    included = dictionary.loc[dictionary.provisional_decision.eq("include_provisional")].copy()
    inventory = (
        included.groupby(["broad_label", "provisional_state_name"], observed=True)
        .agg(
            n_clusters=("native_cluster", "size"),
            n_cohorts=("dataset_id", "nunique"),
            cohorts=("dataset_id", lambda values: ";".join(sorted(set(values)))),
            n_cells=("n_cells", "sum"),
            minimum_samples=("n_samples", "min"),
        )
        .reset_index()
        .sort_values(["broad_label", "provisional_state_name"])
    )
    inventory["eligible_for_crosscohort_harmonisation"] = inventory.n_cohorts.ge(3)
    inventory["claim_boundary"] = (
        "Candidate recurrence only. Formal tumour-normal or progression testing requires "
        "mentor approval and donor/sample-level inference."
    )
    inventory.to_csv(args.output_dir / "PROVISIONAL_STATE_RECURRENCE_INVENTORY.csv", index=False)
    excluded = dictionary.loc[~dictionary.provisional_decision.eq("include_provisional")].copy()
    excluded.to_csv(args.output_dir / "PROVISIONAL_CLUSTER_EXCLUSION_LOG.csv", index=False)

    summary = pd.DataFrame([
        {"metric": "native_clusters_total", "value": len(dictionary)},
        {"metric": "included_provisional_clusters", "value": len(included)},
        {"metric": "excluded_or_unresolved_clusters", "value": len(excluded)},
        {"metric": "recurrent_candidate_states_ge_3_cohorts", "value": int(inventory.eligible_for_crosscohort_harmonisation.sum())},
    ])
    summary.to_csv(args.output_dir / "COMPUTATIONAL_CURATION_SUMMARY.csv", index=False)
    manifest = {
        "input_audit": str(args.audit.resolve()),
        "input_markers": str(args.markers.resolve()),
        "minimum_cells": MIN_CELLS,
        "minimum_samples": MIN_SAMPLES,
        "state_panels": {compartment: {state: sorted(panel) for state, panel in states.items()} for compartment, states in STATE_PANELS.items()},
        "status": "provisional_pending_mentor_review",
    }
    (args.output_dir / "COMPUTATIONAL_CURATION_MANIFEST.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    (args.output_dir / "README.md").write_text(
        "# Computational provisional curation\n\n"
        "This folder records a conservative working state dictionary derived from the completed native per-cohort marker tables. "
        "Clusters flagged by the lineage audit, below support thresholds, or without decisive marker-panel evidence are excluded. "
        "All included labels are provisional and require mentor review before statistical testing or manuscript claims.\n",
        encoding="utf-8",
    )
    print(summary.to_string(index=False))
    print(f"Wrote provisional curation: {args.output_dir}")


if __name__ == "__main__":
    main()
