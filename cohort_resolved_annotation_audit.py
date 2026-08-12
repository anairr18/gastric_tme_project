#!/usr/bin/env python3
"""Audit broad and curated-cell-state annotations separately in every cohort.

The same frozen taxonomy is applied to all cohorts.  This is deliberately not a
new pooled differential-expression analysis: it reports where labels, markers,
and biological comparisons are actually supported before hypotheses are chosen.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import sparse


CORE_IDS = ["korea_kim2022", "kumar2022", "sathe2020", "zhang2021", "diffuse_gc_2021", "tcell_exhaustion_2022"]
BROAD_COLUMNS = ["cell_type_coarse", "broad_cell_type", "cell_type", "cell_type_coarse_fine"]
DATASET_COLUMNS = ["dataset_id", "cohort", "canonical_cohort"]
SAMPLE_COLUMNS = ["sample_id", "sample", "orig.ident", "library_id", "analysis_unit"]
MARKERS = {
    "T/NK": ["CD3D", "CD3E", "TRBC2", "NKG7"],
    "Myeloid": ["LYZ", "FCER1G", "TYROBP", "C1QC"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM"],
    "Endothelial": ["KDR", "EMCN", "VWF", "CLDN5"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "MUC1"],
    "B/Plasma": ["CD74", "MS4A1", "CD79A", "MZB1"],
}


def col(frame: pd.DataFrame, candidates: list[str]) -> str | None:
    lookup = {str(name).lower(): str(name) for name in frame.columns}
    return next((lookup[name.lower()] for name in candidates if name.lower() in lookup), None)


def canonical_broad(series: pd.Series) -> pd.Series:
    # Pandas Categoricals cannot receive a new value with fillna until converted.
    value = series.astype("string").fillna("Unassigned").astype(str)
    lowered = value.str.lower()
    out = value.copy()
    out.loc[lowered.str.contains("t.?cell|nk|lymph", regex=True)] = "T/NK"
    out.loc[lowered.str.contains("myeloid|macro|mono|mast|dendritic", regex=True)] = "Myeloid"
    out.loc[lowered.str.contains("fibro|caf")] = "Fibroblast"
    out.loc[lowered.str.contains("endo|vascular")] = "Endothelial"
    out.loc[lowered.str.contains("epithelial|tumou?r|malignant")] = "Epithelial"
    out.loc[lowered.str.contains("b.?cell|plasma")] = "B/Plasma"
    out.loc[~out.isin(MARKERS)] = "Other/Unassigned"
    return out


def load_state_labels(
    path: Path | None,
    obs: pd.DataFrame,
    dataset_col: str,
    sample_col: str | None,
) -> tuple[pd.Series | None, str]:
    if path is None or not path.exists():
        return None, "State-level audit blocked: frozen cell-level state labels were not supplied."
    labels = pd.read_csv(path, compression="infer")
    id_col = col(labels, ["cell_uid", "cell_id", "obs_name", "cell", "barcode"])
    state_col = col(labels, ["state_id", "curated_state", "state", "reviewed_state"])
    if id_col is None or state_col is None:
        return None, "State-level audit blocked: state-label table needs cell_id/obs_name and state_id columns."
    if labels[id_col].duplicated().any():
        return None, "State-level audit blocked: state-label table has duplicate cell IDs."
    lookup = labels.set_index(id_col)[state_col].astype(str)
    direct = pd.Series(obs.index.astype(str), index=obs.index).map(lookup)
    candidates = [("bare atlas cell IDs", direct)]
    if sample_col:
        qualified = (
            obs[dataset_col].astype(str) + "::" + obs[sample_col].astype(str) + "::" + obs.index.astype(str)
        )
        candidates.append(("dataset::sample::cell IDs", pd.Series(qualified.to_numpy(), index=obs.index).map(lookup)))
    scheme, aligned = max(candidates, key=lambda item: int(item[1].notna().sum()))
    matched = int(aligned.notna().sum())
    if matched == 0:
        return None, f"State-level audit blocked: no atlas cells match frozen state labels using {scheme}."
    if matched != len(obs):
        return aligned, (
            f"State-level audit has partial coverage: {matched}/{len(obs)} atlas cells match frozen state labels using {scheme}. "
            "Unmatched cells are excluded from the reviewed state-test taxonomy and are reported explicitly."
        )
    return aligned, f"State-level audit passed: all atlas cells matched unique frozen state labels using {scheme}."


def marker_support(data, obs: pd.DataFrame, dataset_col: str, broad: pd.Series) -> pd.DataFrame:
    var_lookup = {str(gene).upper(): index for index, gene in enumerate(data.var_names.astype(str))}
    matrix = data.X
    rows: list[dict[str, object]] = []
    for dataset_id in CORE_IDS:
        cohort_indices = np.flatnonzero(obs[dataset_col].astype(str).to_numpy() == dataset_id)
        for state, genes in MARKERS.items():
            indices = [var_lookup[gene] for gene in genes if gene in var_lookup]
            selected = cohort_indices[broad.iloc[cohort_indices].to_numpy() == state]
            if not len(selected) or not indices:
                rows.append({"dataset_id": dataset_id, "broad_cell_type": state, "n_cells": int(len(selected)), "n_marker_genes_available": len(indices), "mean_marker_expression": np.nan, "fraction_cells_expressing_any_marker": np.nan, "marker_support_status": "not_assessable"})
                continue
            values = matrix[selected][:, indices]
            mean = float(values.mean())
            if sparse.issparse(values):
                fraction = float(np.asarray((values.sum(axis=1) > 0)).mean())
            else:
                fraction = float((np.asarray(values).sum(axis=1) > 0).mean())
            rows.append({"dataset_id": dataset_id, "broad_cell_type": state, "n_cells": int(len(selected)), "n_marker_genes_available": len(indices), "mean_marker_expression": mean, "fraction_cells_expressing_any_marker": fraction, "marker_support_status": "assessed"})
    return pd.DataFrame(rows)


def figure_composition(frame: pd.DataFrame, destination: Path) -> None:
    pivot = frame.pivot(index="dataset_id", columns="broad_cell_type", values="cell_fraction").fillna(0).reindex(CORE_IDS)
    pivot.plot(kind="bar", stacked=True, figsize=(12, 6), colormap="tab20")
    plt.ylabel("Cell fraction within cohort")
    plt.xlabel("")
    plt.legend(title="Frozen broad taxonomy", bbox_to_anchor=(1.02, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(destination.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.savefig(destination.with_suffix(".pdf"), bbox_inches="tight")
    plt.close()


def figure_state_coverage(frame: pd.DataFrame, destination: Path) -> None:
    if frame.empty:
        return
    pivot = frame.pivot(index="state_id", columns="dataset_id", values="n_cells").fillna(0).reindex(columns=CORE_IDS, fill_value=0)
    pivot = pivot.loc[pivot.sum(axis=1).sort_values(ascending=False).head(40).index]
    plt.figure(figsize=(11, max(4, 0.28 * len(pivot))))
    sns.heatmap(np.log10(1 + pivot), cmap="viridis", cbar_kws={"label": "log10(1 + cells)"})
    plt.xlabel("Cohort")
    plt.ylabel("Frozen curated state")
    plt.tight_layout()
    plt.savefig(destination.with_suffix(".png"), dpi=300, bbox_inches="tight")
    plt.savefig(destination.with_suffix(".pdf"), bbox_inches="tight")
    plt.close()


def write_questions(context: pd.DataFrame, state_status: str, destination: Path) -> None:
    rows = context.set_index("dataset_id")
    questions = [
        "# Data-Grounded Research Questions",
        "",
        "These are analysis questions justified by the frozen cohort designs. They are not conclusions.",
        "",
        "1. **Cross-cohort tumour-normal state abundance:** estimate a donor/sample-level effect only in cohorts whose audited metadata explicitly map tumour and normal tissue; pool only comparable effects.",
        "2. **Korean longitudinal treatment context:** describe within-patient baseline-to-follow-up state changes, separated by documented response labels after source metadata audit. This is not an external response-prediction analysis.",
        "3. **Zhang disease continuum:** test pre-specified state abundance and programme trends across NAG, CAG, IM, and EGC within this cohort only; label it exploratory until independently replicated.",
        "4. **Diffuse gastric depth:** test normal versus superficial/deep tumour state variation only after reconciling Jeong donor and tissue labels.",
        "5. **Interaction hypotheses:** only test state pairs that show marker support and adequate sample coverage across cohorts. Ligand-receptor coexpression remains hypothesis-generating.",
        "",
        "## Boundaries",
        f"- {state_status}",
        "- The frozen core has no esophageal/GEJ or metastatic-site cohort, so it cannot answer cross-site or primary-versus-metastatic questions.",
        "- Only the Korean cohort is serial chemoimmunotherapy; it cannot establish general treatment-response prediction.",
    ]
    destination.write_text("\n".join(questions) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, required=True)
    parser.add_argument("--cohort-context", type=Path, required=True)
    parser.add_argument("--state-labels", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    import anndata as ad
    data = ad.read_h5ad(args.atlas, backed="r")
    try:
        obs = data.obs.copy()
        dataset_col = col(obs, DATASET_COLUMNS)
        broad_col = col(obs, BROAD_COLUMNS)
        sample_col = col(obs, SAMPLE_COLUMNS)
        if dataset_col is None or broad_col is None:
            raise ValueError("Atlas needs dataset_id/cohort and a broad cell-type annotation column.")
        obs["dataset_id"] = obs[dataset_col].astype(str)
        unexpected = sorted(set(obs["dataset_id"]) - set(CORE_IDS))
        if unexpected:
            raise ValueError(f"Atlas includes unexpected cohorts for this frozen-core audit: {unexpected}")
        obs["broad_cell_type"] = canonical_broad(obs[broad_col])
        obs["analysis_unit"] = obs[sample_col].astype(str) if sample_col else "MISSING_SAMPLE"
        broad_counts = obs.groupby(["dataset_id", "broad_cell_type"], observed=True).size().rename("n_cells").reset_index()
        broad_counts["cohort_cells"] = broad_counts.groupby("dataset_id")["n_cells"].transform("sum")
        broad_counts["cell_fraction"] = broad_counts["n_cells"] / broad_counts["cohort_cells"]
        broad_counts.to_csv(args.output_dir / "PER_COHORT_BROAD_ANNOTATION_COUNTS.csv", index=False)
        sample_counts = obs.groupby(["dataset_id", "analysis_unit", "broad_cell_type"], observed=True).size().rename("n_cells").reset_index()
        sample_counts.to_csv(args.output_dir / "PER_COHORT_SAMPLE_BROAD_ANNOTATION_COUNTS.csv", index=False)
        annotation_qc = obs.groupby("dataset_id", observed=True).agg(n_cells=("dataset_id", "size"), n_samples=("analysis_unit", "nunique"), n_broad_types=("broad_cell_type", "nunique")).reset_index()
        annotation_qc["broad_annotation_column"] = broad_col
        annotation_qc["sample_metadata_available"] = sample_col is not None
        annotation_qc.to_csv(args.output_dir / "PER_COHORT_ANNOTATION_QC.csv", index=False)
        support = marker_support(data, obs, dataset_col, obs["broad_cell_type"])
        support.to_csv(args.output_dir / "PER_COHORT_BROAD_MARKER_SUPPORT.csv", index=False)

        states, state_status = load_state_labels(args.state_labels, obs, dataset_col, sample_col)
        state_counts = pd.DataFrame(columns=["dataset_id", "analysis_unit", "state_id", "n_cells"])
        if states is not None:
            obs["state_id"] = states.to_numpy()
            labeled = obs.loc[obs["state_id"].notna()].copy()
            state_counts = labeled.groupby(["dataset_id", "analysis_unit", "state_id"], observed=True).size().rename("n_cells").reset_index()
            state_counts.to_csv(args.output_dir / "PER_COHORT_STATE_ANNOTATION_COUNTS.csv", index=False)
            state_coverage_audit = obs.groupby("dataset_id", observed=True).agg(total_atlas_cells=("dataset_id", "size"), state_labeled_cells=("state_id", "count")).reset_index()
            state_coverage_audit["state_label_coverage_fraction"] = state_coverage_audit["state_labeled_cells"] / state_coverage_audit["total_atlas_cells"]
            state_coverage_audit.to_csv(args.output_dir / "PER_COHORT_STATE_LABEL_COVERAGE_AUDIT.csv", index=False)
            coverage = state_counts.groupby(["dataset_id", "state_id"], observed=True).agg(n_cells=("n_cells", "sum"), n_samples=("analysis_unit", "nunique")).reset_index()
            coverage["n_labeled_cohorts"] = coverage.groupby("state_id")["dataset_id"].transform("nunique")
            coverage.to_csv(args.output_dir / "STATE_PRESENCE_AND_COVERAGE.csv", index=False)
            figure_state_coverage(coverage, args.output_dir / "STATE_PRESENCE_HEATMAP")
        else:
            pd.DataFrame({"state_audit_status": [state_status]}).to_csv(args.output_dir / "PER_COHORT_STATE_ANNOTATION_STATUS.csv", index=False)
        figure_composition(broad_counts, args.output_dir / "PER_COHORT_BROAD_COMPOSITION")
        context = pd.read_csv(args.cohort_context)
        write_questions(context, state_status, args.output_dir / "DATA_DRIVEN_RESEARCH_QUESTIONS.md")
        (args.output_dir / "RUN_MANIFEST.json").write_text(json.dumps({"atlas": str(args.atlas), "broad_annotation_column": broad_col, "sample_column": sample_col, "state_labels": str(args.state_labels) if args.state_labels else None, "state_audit_status": state_status}, indent=2) + "\n", encoding="utf-8")
    finally:
        data.file.close()
    print(f"Per-cohort annotation audit written to {args.output_dir}")


if __name__ == "__main__":
    main()
