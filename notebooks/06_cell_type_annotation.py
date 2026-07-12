"""
06_cell_type_annotation.py
Cell type annotation of the integrated meta-dataset.

Two-pass strategy:
  Pass 1 — CellTypist automated annotation:
      Uses pretrained models trained on human immune + epithelial atlases.
      Models used: "Immune_All_High.pkl" (immune) + "Human_Colorectal_Cancer.pkl" (GI epithelial)
      Produces majority-vote labels per Leiden cluster.

  Pass 2 — Gastric TME marker refinement:
      Known marker gene sets for gastric cancer cell types are scored
      (sc.tl.score_genes) and used to verify / correct CellTypist assignments.
      Coarse cell type hierarchy is written to obs["cell_type_coarse"],
      fine-grained labels to obs["cell_type_fine"].

Outputs:
  data/processed/integrated/gastric_meta_annotated.h5ad
  outputs/annotation/umap_cell_type_coarse.png
  outputs/annotation/umap_cell_type_fine.png
  outputs/annotation/dotplot_markers_coarse.png
  outputs/annotation/cluster_annotation_table.csv

Usage:
    python 06_cell_type_annotation.py
    python 06_cell_type_annotation.py --input data/processed/integrated/gastric_meta_integrated.h5ad
    python 06_cell_type_annotation.py --mode markers-only   # skip CellTypist, use markers only
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", line_buffering=True)
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", line_buffering=True)

import os
import argparse
import warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import scanpy as sc
import celltypist
from celltypist import models

warnings.filterwarnings("ignore")

BASE     = os.path.expanduser("~/gastric_tme_project")
INT_DIR  = os.path.join(BASE, "data/processed/integrated")
ANN_DIR  = os.path.join(BASE, "outputs/annotation")
MANIFEST = os.path.join(BASE, "data/external/dataset_manifest.csv")

os.makedirs(ANN_DIR, exist_ok=True)
sc.settings.figdir = ANN_DIR


# ─── Gastric TME marker gene sets ─────────────────────────────────────────────
# Evidence-based markers from Kumar 2022, Sathe 2020, Zhang 2021
GASTRIC_MARKERS = {
    # Epithelial compartment
    "Epithelial / Tumor":       ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "MUC5AC", "TFF1", "CLDN4"],
    "Cancer-associated Fibro":  ["FAP", "ACTA2", "MMP11", "PDPN", "POSTN", "COL1A1", "COL1A2", "THY1"],

    # Myeloid compartment
    "Macrophage":               ["CD68", "CD163", "MRC1", "MARCO", "MSR1", "APOE", "C1QA", "C1QB"],
    "M1 Macrophage":            ["CD80", "CD86", "IL1B", "TNF", "CXCL9", "CXCL10", "FCGR1A"],
    "M2 Macrophage":            ["MRC1", "CD163", "IL10", "TGM2", "FOLR2", "LYVE1", "SELENOP"],
    "Dendritic Cell":           ["ITGAX", "CD1C", "FCER1A", "CLEC9A", "IRF8", "BATF3", "XCR1"],
    "Monocyte":                 ["CD14", "FCGR3A", "LYZ", "S100A8", "S100A9", "CST3", "VCAN"],
    "Mast Cell":                ["TPSAB1", "TPSB2", "CPA3", "MS4A2", "FCER1A", "KIT", "HDC"],

    # T cell compartment
    "CD4+ T cell":              ["CD3D", "CD3E", "CD4", "IL7R", "LDHB", "TCF7", "CCR7"],
    "CD8+ T cell":              ["CD3D", "CD3E", "CD8A", "CD8B", "GZMA", "GZMB", "PRF1"],
    "Treg":                     ["CD4", "FOXP3", "CTLA4", "IL2RA", "IKZF2", "TNFRSF9"],
    "Exhausted CD8+ T":         ["CD8A", "PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "ENTPD1"],
    "NK Cell":                  ["NCAM1", "NKG7", "GNLY", "KLRD1", "KLRB1", "FCGR3A", "XCL1"],
    "Cycling T / NK":           ["MKI67", "TOP2A", "CDK1", "PCNA", "STMN1", "TUBB"],

    # B cell compartment
    "B Cell":                   ["CD19", "MS4A1", "CD79A", "CD79B", "IGHM", "PAX5"],
    "Plasma Cell":               ["MZB1", "IGHG1", "IGHG2", "SDC1", "XBP1", "JCHAIN", "PRDM1"],

    # Stromal / vascular
    "Endothelial":              ["PECAM1", "VWF", "CDH5", "CLDN5", "PLVAP", "ERG", "ENG"],
    "Smooth Muscle / Pericyte": ["RGS5", "MCAM", "PDGFRB", "CSPG4", "ACTA2", "MYH11", "CNN1"],
}

COARSE_MAP = {
    "Epithelial / Tumor":       "Epithelial",
    "Cancer-associated Fibro":  "Fibroblast",
    "Macrophage":               "Myeloid",
    "M1 Macrophage":            "Myeloid",
    "M2 Macrophage":            "Myeloid",
    "Dendritic Cell":           "Myeloid",
    "Monocyte":                 "Myeloid",
    "Mast Cell":                "Myeloid",
    "CD4+ T cell":              "T/NK",
    "CD8+ T cell":              "T/NK",
    "Treg":                     "T/NK",
    "Exhausted CD8+ T":         "T/NK",
    "NK Cell":                  "T/NK",
    "Cycling T / NK":           "T/NK",
    "B Cell":                   "B/Plasma",
    "Plasma Cell":              "B/Plasma",
    "Endothelial":              "Endothelial",
    "Smooth Muscle / Pericyte": "Stromal",
}


# ─── helpers ──────────────────────────────────────────────────────────────────

def run_celltypist(adata):
    """Run CellTypist with two complementary models and return predictions."""
    print("\n--- CellTypist annotation ---")

    # CellTypist needs log-normalized counts in .X (10k normalization, log1p)
    # The integrated object has raw counts in .X after concatenation for scVI.
    # We need to normalize a copy for CellTypist.
    adata_ct = adata.copy()
    if adata_ct.raw is not None:
        adata_ct.X = adata_ct.raw.to_adata().X.copy()
    sc.pp.normalize_total(adata_ct, target_sum=1e4)
    sc.pp.log1p(adata_ct)

    results = {}

    for model_name in ["Immune_All_High.pkl", "Human_Colorectal_Cancer.pkl"]:
        print(f"  loading model: {model_name} ...")
        try:
            model = models.Model.load(model=model_name)
        except Exception:
            print(f"  downloading {model_name} ...")
            models.download_models(model=model_name)
            model = models.Model.load(model=model_name)

        print(f"  predicting ...")
        pred = celltypist.annotate(
            adata_ct,
            model=model,
            majority_voting=True,
            over_clustering="leiden",
        )
        key = model_name.replace(".pkl", "")
        adata.obs[f"celltypist_{key}"]              = pred.predicted_labels["predicted_labels"].values
        adata.obs[f"celltypist_{key}_majority"]     = pred.predicted_labels["majority_voting"].values
        adata.obs[f"celltypist_{key}_conf"]         = pred.predicted_labels["conf_score"].values
        results[model_name] = pred
        print(f"  done — {adata.obs[f'celltypist_{key}_majority'].nunique()} unique labels")

    return adata, results


def score_marker_genes(adata):
    """Score all marker gene sets; assign best-scoring label per cell."""
    print("\n--- Marker gene scoring ---")

    # score_genes requires log-normalized X; integrated object has raw counts.
    # Normalize in-place for scoring and restore from layers["counts"] after.
    x_check = adata.X
    max_val = float(x_check.data.max()) if sp.issparse(x_check) and x_check.nnz > 0 else float(x_check.max())
    needs_norm = max_val > 50
    if needs_norm:
        print("  normalizing X for scoring (raw counts in .X); will restore after ...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    score_cols = []
    for cell_type, markers in GASTRIC_MARKERS.items():
        valid_markers = [g for g in markers if g in adata.var_names]
        if len(valid_markers) < 2:
            print(f"  SKIP {cell_type}: only {len(valid_markers)} markers in gene set")
            continue
        key = f"score_{cell_type.replace(' ', '_').replace('/', '_')}"
        sc.tl.score_genes(adata, gene_list=valid_markers, score_name=key)
        score_cols.append((cell_type, key))
        print(f"  {cell_type}: scored {len(valid_markers)}/{len(markers)} markers")

    if score_cols:
        score_matrix = np.column_stack([adata.obs[k].values for _, k in score_cols])
        best_idx     = np.argmax(score_matrix, axis=1)
        cell_labels  = [score_cols[i][0] for i in best_idx]
        adata.obs["cell_type_marker"]       = cell_labels
        adata.obs["cell_type_coarse_marker"] = [COARSE_MAP[l] for l in cell_labels]
        print(f"\n  Cell type distribution (marker-based, coarse):")
        print(adata.obs["cell_type_coarse_marker"].value_counts().to_string())

    if needs_norm and "counts" in adata.layers:
        adata.X = adata.layers["counts"]
        print("  restored raw counts to .X")

    return adata, score_cols


def assign_final_labels(adata):
    """
    Combine CellTypist majority voting with marker scores per Leiden cluster.
    Priority: if CellTypist confidence is high (>0.6) → use CellTypist;
              else fall back to marker-based label.
    Final labels written to obs["cell_type_fine"] and obs["cell_type_coarse"].
    """
    print("\n--- Assigning final labels ---")

    # Try to use Immune_All_High majority voting as primary
    ct_col = None
    for col in ["celltypist_Immune_All_High_majority",
                "celltypist_Human_Colorectal_Cancer_majority"]:
        if col in adata.obs.columns:
            ct_col = col
            break

    if ct_col and "cell_type_marker" in adata.obs.columns:
        conf_col = ct_col.replace("majority", "conf")
        conf = adata.obs.get(conf_col, pd.Series(0.0, index=adata.obs.index))
        fine = np.where(conf > 0.6,
                        adata.obs[ct_col].values,
                        adata.obs["cell_type_marker"].values)
        adata.obs["cell_type_fine"] = fine
    elif ct_col:
        adata.obs["cell_type_fine"] = adata.obs[ct_col]
    elif "cell_type_marker" in adata.obs.columns:
        adata.obs["cell_type_fine"] = adata.obs["cell_type_marker"]
    else:
        adata.obs["cell_type_fine"] = "Unknown"

    # Coarse: map fine label through COARSE_MAP (best effort, unknown otherwise)
    def to_coarse(label):
        if label in COARSE_MAP:
            return COARSE_MAP[label]
        # partial match
        for k, v in COARSE_MAP.items():
            if k.lower() in label.lower() or label.lower() in k.lower():
                return v
        return "Other"

    adata.obs["cell_type_coarse"] = adata.obs["cell_type_fine"].map(to_coarse)

    print("\nCoarse cell type distribution:")
    print(adata.obs["cell_type_coarse"].value_counts().to_string())

    return adata


def save_annotation_plots(adata, score_cols):
    print(f"\nSaving annotation plots to {ANN_DIR} ...")

    sc.pl.umap(adata, color=["cell_type_coarse", "cell_type_fine"],
               ncols=1, save="_cell_type_coarse_fine.png", show=False)
    sc.pl.umap(adata, color=["leiden", "dataset_id"],
               ncols=2, save="_leiden_dataset.png", show=False)

    # Dotplot of top markers per coarse type
    coarse_types = adata.obs["cell_type_coarse"].unique().tolist()
    marker_genes_for_dot = []
    seen = set()
    for cell_type, markers in GASTRIC_MARKERS.items():
        for g in markers[:3]:
            if g in adata.var_names and g not in seen:
                marker_genes_for_dot.append(g)
                seen.add(g)
    if marker_genes_for_dot:
        sc.pl.dotplot(adata, var_names=marker_genes_for_dot,
                      groupby="cell_type_coarse",
                      save="_markers_coarse.png", show=False)

    # Cluster annotation table
    cluster_table = (adata.obs.groupby("leiden")
                     .agg(
                         n_cells=("cell_type_coarse", "count"),
                         cell_type_coarse=("cell_type_coarse", lambda x: x.value_counts().index[0]),
                         cell_type_fine=("cell_type_fine", lambda x: x.value_counts().index[0]),
                     )
                     .reset_index())
    cluster_table.to_csv(os.path.join(ANN_DIR, "cluster_annotation_table.csv"), index=False)
    print(f"  Cluster annotation table saved.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=os.path.join(INT_DIR, "gastric_meta_integrated.h5ad"))
    parser.add_argument("--output", default=os.path.join(INT_DIR, "gastric_meta_annotated.h5ad"))
    parser.add_argument("--mode", choices=["full", "markers-only", "celltypist-only"],
                        default="full")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    if os.path.exists(args.output) and not args.force:
        print(f"Annotated h5ad already exists: {args.output}")
        print("Use --force to redo.")
        return

    if not os.path.exists(args.input):
        print(f"Input not found: {args.input}")
        print("Run 05_integration.py first.")
        sys.exit(1)

    print(f"Loading integrated object from {args.input} ...")
    adata = sc.read_h5ad(args.input)
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    score_cols = []

    # Pass 1: CellTypist
    if args.mode in ("full", "celltypist-only"):
        try:
            adata, _ = run_celltypist(adata)
        except Exception as e:
            print(f"CellTypist failed: {e}")
            print("Falling back to marker-only annotation.")
            args.mode = "markers-only"

    # Pass 2: marker scores
    if args.mode in ("full", "markers-only"):
        adata, score_cols = score_marker_genes(adata)

    # Combine
    adata = assign_final_labels(adata)

    # Plots
    save_annotation_plots(adata, score_cols)

    # Save
    print(f"\nSaving annotated object to {args.output} ...")
    adata.write_h5ad(args.output)
    size_gb = os.path.getsize(args.output) / 1e9
    print(f"Saved ({size_gb:.1f} GB)")

    print(f"\nNext step: run 07_meta_analysis.py")


if __name__ == "__main__":
    main()
