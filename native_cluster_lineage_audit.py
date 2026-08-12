#!/usr/bin/env python3
"""Audit provisional native clusters before cross-cohort state harmonisation.

This consumes the native per-cohort discovery marker tables rather than
re-using integrated coordinates. It flags clusters whose raw marker evidence
is inconsistent with their frozen broad-compartment correspondence, so they
cannot quietly become biological states in later abundance analyses.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


PANELS: dict[str, set[str]] = {
    "B/Plasma": {
        "MS4A1", "CD79A", "CD79B", "CD74", "HLA-DRA", "CD37", "CD22",
        "CD83", "CD19", "CD40", "MZB1", "JCHAIN", "XBP1", "SDC1",
        "DERL3", "IGKC", "IGHM", "IGHA1", "IGHG1",
    },
    "Endothelial": {
        "VWF", "KDR", "EMCN", "PECAM1", "CLDN5", "RAMP2", "PLVAP", "CA4",
        "ESAM", "ENG", "RGCC", "ACKR1", "KLF2",
    },
    "Epithelial": {
        "EPCAM", "KRT8", "KRT18", "KRT19", "KRT20", "KRT7", "MUC1", "MUC5AC",
        "MUC6", "TFF1", "TFF2", "TFF3", "CLDN18", "KRT17", "KRT23",
    },
    "Fibroblast": {
        "COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2", "DCN", "LUM",
        "PDGFRA", "COL5A1", "COL5A2", "C7", "CFD", "COL15A1", "FAP", "POSTN",
    },
    "Myeloid": {
        "LYZ", "LST1", "TYROBP", "FCER1G", "CTSS", "AIF1", "LILRB1", "LILRB2",
        "C1QA", "C1QB", "C1QC", "APOE", "S100A8", "S100A9", "FCN1", "LGALS3",
        "SIGLEC10", "PLAC8", "BCL2A1",
    },
    "Stromal": {
        "RGS5", "CSPG4", "MCAM", "ACTA2", "TAGLN", "DES", "PDGFRB", "NOTCH3",
        "MYH11", "CNN1", "COL4A1", "COL4A2",
    },
    "T/NK": {
        "CD3D", "CD3E", "TRAC", "TRBC1", "TRBC2", "CD247", "LCK", "NKG7",
        "KLRD1", "KLRB1", "GNLY", "GZMB", "GZMH", "PRF1", "CCL5", "IL7R",
        "MALAT1", "TRAT1", "CD2",
    },
}

TECHNICAL = {
    "MKI67", "TOP2A", "HMGB2", "STMN1", "TUBA1B", "TUBB", "H2AFZ", "PCNA",
    "HBB", "HBA1", "HBA2", "ALAS2", "SLC25A37", "HSP90AA1", "HSPA1A", "HSPA1B",
    "RPL", "RPS", "MALAT1", "MT-",
}


def technical_gene(gene: str) -> bool:
    return any(gene == token or gene.startswith(token) for token in TECHNICAL)


def gene_hits(genes: list[str], panel: set[str]) -> list[str]:
    return [gene for gene in genes if gene in panel]


def audit_cluster(group: pd.DataFrame, n_cells: int, n_samples: int) -> dict:
    dataset_id = str(group["dataset_id"].iloc[0])
    broad = str(group["broad_label"].iloc[0])
    cluster = str(group["native_cluster"].iloc[0])
    ranked = group.sort_values(["pvals_adj", "scores"], ascending=[True, False])
    top_genes = ranked["names"].astype(str).head(25).tolist()
    supports = {label: gene_hits(top_genes, panel) for label, panel in PANELS.items()}
    expected_hits = supports.get(broad, [])
    other = {label: hits for label, hits in supports.items() if label != broad}
    discordant_label, discordant_hits = max(other.items(), key=lambda item: len(item[1]))
    nontechnical = [gene for gene in top_genes if not technical_gene(gene)]

    if len(expected_hits) >= 2 and len(discordant_hits) <= 1 and n_cells >= 50 and n_samples >= 3:
        recommendation = "eligible_after_manual_state_annotation"
    elif len(expected_hits) == 0 and len(discordant_hits) >= 2:
        recommendation = "exclude_suspected_cross_lineage_or_misassignment"
    elif n_cells < 30 or n_samples < 3:
        recommendation = "exclude_insufficient_cell_or_sample_support"
    elif len(nontechnical) < 5:
        recommendation = "exclude_technical_or_low_information_cluster"
    else:
        recommendation = "manual_review_required"

    return {
        "dataset_id": dataset_id,
        "broad_label": broad,
        "native_cluster": cluster,
        "n_cells": int(n_cells),
        "n_samples": int(n_samples),
        "expected_panel_hits": ";".join(expected_hits),
        "n_expected_panel_hits": len(expected_hits),
        "discordant_panel": discordant_label,
        "discordant_panel_hits": ";".join(discordant_hits),
        "n_discordant_panel_hits": len(discordant_hits),
        "top_nontechnical_markers": ";".join(nontechnical[:12]),
        "curation_recommendation": recommendation,
        "manual_state_name": "",
        "manual_marker_rationale": "",
        "reviewer": "",
        "review_status": "pending_review",
        "claim_boundary": (
            "A retained cluster is a within-cohort transcriptional state candidate. "
            "It is not a cross-cohort state or disease-associated result until manual review, "
            "cross-cohort harmonisation, and donor/sample-level inference are complete."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--discovery-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    markers_path = args.discovery_dir / "NATIVE_PER_COHORT_CLUSTER_MARKERS.csv"
    summary_path = args.discovery_dir / "NATIVE_PER_COHORT_CLUSTER_SUMMARY.csv"
    if not markers_path.exists() or not summary_path.exists():
        raise FileNotFoundError("Discovery directory lacks final marker and cluster-summary tables.")
    markers = pd.read_csv(markers_path)
    summary = pd.read_csv(summary_path)
    required_markers = {"dataset_id", "broad_label", "native_cluster", "names", "scores", "pvals_adj"}
    if missing := required_markers - set(markers.columns):
        raise ValueError(f"Marker table lacks: {sorted(missing)}")

    summary["native_cluster"] = summary["native_cluster"].astype(str)
    markers["native_cluster"] = markers["native_cluster"].astype(str)
    sizes = summary.set_index(["dataset_id", "broad_label", "native_cluster"])
    rows = []
    for key, group in markers.groupby(["dataset_id", "broad_label", "native_cluster"], observed=True):
        if key not in sizes.index:
            raise ValueError(f"Missing cluster size summary for {key}")
        values = sizes.loc[key]
        rows.append(audit_cluster(group, int(values["n_cells"]), int(values["n_samples"])))
    audit = pd.DataFrame(rows).sort_values(
        ["curation_recommendation", "dataset_id", "broad_label", "n_cells"],
        ascending=[True, True, True, False],
    )
    audit.to_csv(args.output_dir / "NATIVE_CLUSTER_LINEAGE_AUDIT.csv", index=False)
    audit.loc[audit["curation_recommendation"].ne("eligible_after_manual_state_annotation")].to_csv(
        args.output_dir / "NATIVE_CLUSTER_REVIEW_QUEUE.csv", index=False
    )
    qc = (
        audit.groupby(["dataset_id", "broad_label", "curation_recommendation"], observed=True)
        .agg(n_clusters=("native_cluster", "size"), n_cells=("n_cells", "sum"))
        .reset_index()
    )
    qc.to_csv(args.output_dir / "NATIVE_CLUSTER_LINEAGE_AUDIT_SUMMARY.csv", index=False)

    overview = (
        audit.pivot_table(index="dataset_id", columns="curation_recommendation", values="native_cluster", aggfunc="count", fill_value=0)
        .reindex(columns=sorted(audit["curation_recommendation"].unique()), fill_value=0)
    )
    fig, ax = plt.subplots(figsize=(10, 5))
    overview.plot(kind="bar", stacked=True, ax=ax, colormap="tab20")
    ax.set_ylabel("Provisional native clusters")
    ax.set_xlabel("Cohort")
    ax.set_title("Native-cluster lineage audit: review status")
    ax.legend(title="Audit recommendation", bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(args.output_dir / "NATIVE_CLUSTER_LINEAGE_AUDIT_OVERVIEW.png", dpi=300)
    fig.savefig(args.output_dir / "NATIVE_CLUSTER_LINEAGE_AUDIT_OVERVIEW.pdf")
    plt.close(fig)

    (args.output_dir / "README_AND_REVIEW_INSTRUCTIONS.md").write_text(
        "# Native Cluster Lineage Audit\n\n"
        "This is a marker-panel triage, not an automatic annotation system. Review every retained and manual-review cluster using the complete Wilcoxon marker table and UMAP. Clusters flagged as cross-lineage must not enter state abundance, interaction, or clinical analyses unless a reviewer documents a correction with evidence.\n\n"
        "For example, a cluster assigned to B/Plasma but dominated by myeloid markers such as SIGLEC10, PLAC8, and BCL2A1 is a candidate cross-lineage/misassignment and is excluded by default.\n",
        encoding="utf-8",
    )
    (args.output_dir / "RUN_MANIFEST.json").write_text(json.dumps({
        "source_discovery_dir": str(args.discovery_dir),
        "n_provisional_clusters": int(len(audit)),
        "n_cohorts": int(audit["dataset_id"].nunique()),
        "purpose": "marker-panel triage before manual state annotation and cross-cohort harmonisation",
    }, indent=2) + "\n", encoding="utf-8")
    print(f"Native cluster lineage audit written to {args.output_dir}")


if __name__ == "__main__":
    main()
