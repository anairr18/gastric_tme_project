#!/usr/bin/env python3
"""Create a reviewer-facing curation workbook for native per-cohort clusters.

Candidate labels are generated from transparent marker panels only to focus
review. They are never accepted automatically. A cluster may enter a later
cross-cohort harmonisation only after a reviewer sets ``review_decision`` to
``include`` and records a reviewed state name and rationale.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


STATE_PANELS: dict[str, dict[str, set[str]]] = {
    "B/Plasma": {
        "HLAII_B_cell": {"MS4A1", "CD79A", "CD79B", "CD74", "HLA-DRA", "CD37", "CD22", "CD19"},
        "secretory_plasma": {"MZB1", "XBP1", "DERL3", "JCHAIN", "SDC1", "TNFRSF17", "SSR4"},
        "cycling_B_cell": {"MKI67", "TOP2A", "STMN1", "HMGB2", "TYMS", "TUBB"},
    },
    "Endothelial": {
        "capillary_endothelial": {"CA4", "PLVAP", "RGCC", "GNG11", "KDR", "BTNL9", "FABP4"},
        "lymphatic_endothelial": {"CCL21", "PROX1", "LYVE1", "PDPN", "FLT4", "CCL14"},
        "arterial_endothelial": {"EFNB2", "GJA4", "SOX17", "SEMA3G", "HEY1"},
        "activated_endothelial": {"ICAM1", "VCAM1", "SELE", "SELP", "ACKR1", "HLA-DRA", "CXCL10"},
    },
    "Epithelial": {
        "gastric_mucous_epithelial": {"MUC5AC", "TFF1", "TFF2", "GKN1", "GKN2", "PGC", "LIPF"},
        "intestinal_metaplastic_epithelial": {"TFF3", "KRT20", "FABP1", "ANPEP", "MUC2", "LGALS4", "REG4", "PHGR1"},
        "parietal_epithelial": {"ATP4A", "ATP4B", "GIF", "KCNE2", "CKB"},
        "tumour_epithelial_candidate": {"EPCAM", "KRT7", "KRT19", "MMP7", "TM4SF1", "CEACAM5", "KRT17", "CLDN7", "ERBB2"},
        "cycling_epithelial": {"MKI67", "TOP2A", "STMN1", "TYMS", "TUBA1B"},
    },
    "Fibroblast": {
        "matrix_CAF": {"COL1A1", "COL1A2", "COL3A1", "POSTN", "THBS2", "COL6A3", "FAP", "CTHRC1"},
        "inflammatory_CAF": {"CXCL12", "CXCL14", "CCL11", "IL6", "CCL2", "PDPN", "CTSK"},
        "resident_fibroblast": {"DCN", "LUM", "CFD", "C7", "C1R", "C1S", "PTGDS", "SFRP2"},
        "perivascular_stromal": {"RGS5", "MCAM", "CSPG4", "ACTA2", "TAGLN", "PDGFRB"},
    },
    "Myeloid": {
        "C1QC_macrophage": {"C1QA", "C1QB", "C1QC", "APOE", "APOC1", "LGMN", "CTSD"},
        "SPP1_macrophage": {"SPP1", "GPNMB", "TREM2", "APOC1", "LGALS3", "HMOX1", "HAMP"},
        "inflammatory_monocyte": {"S100A8", "S100A9", "FCN1", "IL1B", "CXCL8", "S100A12", "VCAN"},
        "antigen_presenting_myeloid": {"HLA-DRA", "HLA-DPA1", "HLA-DPB1", "CD74", "CLEC10A", "CD1C"},
        "mast_cell": {"TPSAB1", "TPSB2", "CPA3", "HDC", "KIT", "MS4A2"},
        "pDC": {"GZMB", "GZMB", "LILRA4", "GZMB", "TCF4", "IL3RA", "GZMB"},
    },
    "Stromal": {
        "pericyte_smooth_muscle": {"RGS5", "MCAM", "CSPG4", "ACTA2", "TAGLN", "MYH11", "CNN1", "DES"},
        "matrix_stromal": {"COL4A1", "COL4A2", "COL1A1", "COL1A2", "SPARC", "THY1"},
    },
    "T/NK": {
        "CD4_memory_T": {"IL7R", "LTB", "CCR7", "MALAT1", "LST1", "KLRB1", "MAL"},
        "Treg_candidate": {"FOXP3", "IL2RA", "CTLA4", "TIGIT", "TNFRSF4", "TNFRSF18", "IKZF2"},
        "GZMK_effector_memory_T": {"GZMK", "GZMA", "CCL5", "NKG7", "CXCR3", "TRAC"},
        "cytotoxic_NK_T": {"NKG7", "GNLY", "GZMB", "GZMH", "PRF1", "FGFBP2", "KLRD1", "XCL1"},
        "activated_T": {"CD69", "NR4A1", "NR4A2", "NR4A3", "TNFRSF9", "HLA-DRA", "FOS"},
        "cycling_T_NK": {"MKI67", "TOP2A", "STMN1", "TYMS", "TUBB", "HMGB2"},
    },
}


def marker_candidates(group: pd.DataFrame, broad: str) -> tuple[str, int, str]:
    ranked = group.sort_values(["pvals_adj", "scores"], ascending=[True, False])
    genes = ranked["names"].astype(str).head(30).tolist()
    programs = STATE_PANELS.get(broad, {})
    hits = {name: [gene for gene in genes if gene in panel] for name, panel in programs.items()}
    if not hits:
        return "unresolved", 0, ""
    winner, win_hits = max(hits.items(), key=lambda item: len(item[1]))
    if not win_hits:
        return "unresolved", 0, ""
    return winner, len(win_hits), ";".join(win_hits)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--discovery-dir", type=Path, required=True)
    parser.add_argument("--audit-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    markers = pd.read_csv(args.discovery_dir / "NATIVE_PER_COHORT_CLUSTER_MARKERS.csv")
    audit = pd.read_csv(args.audit_dir / "NATIVE_CLUSTER_LINEAGE_AUDIT.csv")
    key = ["dataset_id", "broad_label", "native_cluster"]
    markers["native_cluster"] = markers["native_cluster"].astype(str)
    audit["native_cluster"] = audit["native_cluster"].astype(str)

    candidates = []
    for values, group in markers.groupby(key, observed=True):
        dataset_id, broad, cluster = values
        proposal, n_hits, hits = marker_candidates(group, broad)
        ranked = group.sort_values(["pvals_adj", "scores"], ascending=[True, False])
        candidates.append({
            "dataset_id": dataset_id,
            "broad_label": broad,
            "native_cluster": str(cluster),
            "marker_panel_proposal": proposal,
            "proposal_marker_hits": hits,
            "n_proposal_marker_hits": n_hits,
            "top_15_markers": ";".join(ranked["names"].astype(str).head(15)),
            "top_15_logfoldchanges": ";".join(
                pd.to_numeric(ranked["logfoldchanges"], errors="coerce").head(15).round(2).astype(str)
            ),
        })
    review = audit.merge(pd.DataFrame(candidates), on=key, how="left", validate="one_to_one")
    forced_exclude = review["curation_recommendation"].str.startswith("exclude_")
    review["recommended_initial_decision"] = np.where(forced_exclude, "exclude", "review")
    review["review_decision"] = review["recommended_initial_decision"]
    review["reviewed_state_name"] = ""
    review["reviewed_marker_rationale"] = ""
    review["reviewer"] = ""
    review["review_date"] = ""
    review["eligible_for_cross_cohort_harmonisation"] = False
    review = review.sort_values(
        ["recommended_initial_decision", "dataset_id", "broad_label", "n_samples", "n_cells"],
        ascending=[True, True, True, False, False],
    )
    review_csv = args.output_dir / "NATIVE_CLUSTER_REVIEW_WORKSHEET.csv"
    review.to_csv(review_csv, index=False)

    summary = (
        review.groupby(["dataset_id", "broad_label", "recommended_initial_decision"], observed=True)
        .agg(n_clusters=("native_cluster", "size"), n_cells=("n_cells", "sum"), median_samples=("n_samples", "median"))
        .reset_index()
    )
    summary.to_csv(args.output_dir / "NATIVE_CLUSTER_REVIEW_SUMMARY.csv", index=False)

    workbook_path = args.output_dir / "NATIVE_CLUSTER_REVIEW_WORKBOOK.xlsx"
    with pd.ExcelWriter(workbook_path, engine="xlsxwriter") as writer:
        review.to_excel(writer, sheet_name="review_all_clusters", index=False)
        summary.to_excel(writer, sheet_name="review_summary", index=False)
        for broad, subset in review.groupby("broad_label", observed=True):
            subset.to_excel(writer, sheet_name=broad.replace("/", "_")[:31], index=False)
        notes = pd.DataFrame({"Review instructions": [
            "This workbook does not create accepted biological states automatically.",
            "Start with rows marked review and inspect their native UMAP plus full Wilcoxon marker table.",
            "Keep review_decision=exclude for cross-lineage, low-sample, erythroid, technical, or unresolved clusters.",
            "For an included cluster, fill reviewed_state_name, reviewed_marker_rationale, reviewer, and review_date.",
            "A later script must validate the signed worksheet before cross-cohort harmonisation.",
        ]})
        notes.to_excel(writer, sheet_name="instructions", index=False)
        workbook = writer.book
        header = workbook.add_format({"bold": True, "bg_color": "#D9EAF7", "text_wrap": True})
        for worksheet in writer.sheets.values():
            worksheet.freeze_panes(1, 0)
            worksheet.set_row(0, 28, header)
            worksheet.set_default_row(18)
            worksheet.set_column(0, 2, 19)
            worksheet.set_column(3, 5, 12)
            worksheet.set_column(6, 12, 24)
            worksheet.set_column(13, 22, 22)

    (args.output_dir / "README.md").write_text(
        "# Native cluster review package\n\n"
        "`NATIVE_CLUSTER_REVIEW_WORKBOOK.xlsx` is the human curation record. The marker-panel proposal is a transparent triage suggestion, not an annotation. Do not edit the original marker or lineage-audit tables; record decisions only in the review worksheet.\n\n"
        "A cluster becomes eligible for cross-cohort harmonisation only after it is marked `include`, has a reviewed state name and rationale, and passes the later worksheet validator.\n",
        encoding="utf-8",
    )
    (args.output_dir / "RUN_MANIFEST.json").write_text(json.dumps({
        "discovery_dir": str(args.discovery_dir),
        "lineage_audit_dir": str(args.audit_dir),
        "n_clusters": int(len(review)),
        "n_forced_exclusion": int(forced_exclude.sum()),
        "purpose": "manual curation before cross-cohort state harmonisation",
    }, indent=2) + "\n", encoding="utf-8")
    print(f"Reviewer package written to {args.output_dir}")


if __name__ == "__main__":
    main()
