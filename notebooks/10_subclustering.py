"""
10_subclustering.py
Subcluster CD8/T/NK cells and Myeloid cells from the Korean annotated cohort.
Identify exhaustion subpopulations and macrophage polarization states.

Outputs in outputs/subclustering/:
  umap_tcell_subclusters.png
  umap_myeloid_subclusters.png
  dotplot_tcell_markers.png
  dotplot_myeloid_markers.png
  tcell_subcluster_by_progression.png
  myeloid_subcluster_by_progression.png
  korea_kim2022_subclustered.h5ad  (adds tcell_subtype, myeloid_subtype to obs)
"""
import sys, os
import argparse
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
import scipy.sparse as sp

warnings.filterwarnings("ignore")

BASE    = os.path.expanduser("~/gastric_tme_project")
PER_DIR = os.path.join(BASE, "data/processed/per_dataset")
OUT_DIR = os.path.join(BASE, "outputs/subclustering")
os.makedirs(OUT_DIR, exist_ok=True)

# Marker genes for T cell subtype annotation
TCELL_MARKERS = {
    "Naive_T":              ["CCR7", "SELL", "TCF7", "LEF1", "KLF2"],
    "Effector_CD8":         ["GZMB", "GZMK", "IFNG", "PRF1", "CX3CR1", "KLRG1"],
    "Progenitor_Exhausted": ["TCF7", "PDCD1", "CXCR5", "TOX", "SLAMF6", "CCR7"],
    "Terminal_Exhausted":   ["PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "TOX2",
                             "ENTPD1", "CTLA4", "NR4A1"],
    "CD4_Helper":           ["CD4", "IL7R", "CXCR5", "ICOS", "BCL6"],
    "Treg":                 ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18"],
    "NK":                   ["GNLY", "NKG7", "KLRD1", "NCAM1", "FCGR3A"],
    "Cycling_T":            ["MKI67", "TOP2A", "PCNA", "TYMS"],
}

# Marker genes for macrophage subtype annotation
MYELOID_MARKERS = {
    "SPP1_Macro":        ["SPP1", "APOC1", "APOE", "C1QB", "TREM2"],
    "FOLR2_Macro":       ["FOLR2", "LYVE1", "SELENOP", "MRC1", "CD163"],
    "Inflammatory_M1":   ["IL1B", "TNF", "CXCL9", "CXCL10", "NOS2", "CD80"],
    "Immunosuppressive": ["IL10", "TGM2", "CCL18", "ARG1", "PPARG", "CD200R1"],
    "Monocyte":          ["S100A8", "S100A9", "LYZ", "FCN1", "CD14", "VCAN"],
    "cDC1":              ["CLEC9A", "XCR1", "CADM1", "IDO1"],
    "cDC2":              ["CD1C", "FCER1A", "CLEC10A", "SIRPA"],
    "pDC":               ["LILRA4", "CLEC4C", "IL3RA", "JCHAIN"],
}


def subcluster(adata_sub, name, n_pcs=20, resolution=0.5):
    """Re-PCA, re-neighbors, re-Leiden on a subset."""
    print(f"  Subclustering {name}: {adata_sub.n_obs:,} cells ...")
    sc.pp.highly_variable_genes(adata_sub, n_top_genes=3000, flavor="seurat", subset=False)
    sc.tl.pca(adata_sub, n_comps=n_pcs, use_highly_variable=True)
    sc.pp.neighbors(adata_sub, n_pcs=n_pcs, n_neighbors=20)
    sc.tl.umap(adata_sub)
    sc.tl.leiden(adata_sub, resolution=resolution, key_added=f"{name}_leiden")
    return adata_sub


def score_subtypes(adata_sub, marker_dict, key_prefix):
    """Score each subtype gene program and assign best-scoring label."""
    var_names = set(adata_sub.var_names)
    score_cols = []
    for subtype, markers in marker_dict.items():
        valid = [g for g in markers if g in var_names]
        if len(valid) < 2:
            continue
        key = f"{key_prefix}_{subtype}"
        sc.tl.score_genes(adata_sub, gene_list=valid, score_name=key)
        score_cols.append((subtype, key))
        print(f"    {subtype}: scored {len(valid)}/{len(markers)} markers")

    if score_cols:
        score_matrix = np.column_stack([adata_sub.obs[k].values for _, k in score_cols])
        best_idx = np.argmax(score_matrix, axis=1)
        adata_sub.obs[f"{key_prefix}_subtype"] = [score_cols[i][0] for i in best_idx]
    return adata_sub, score_cols


def plot_subclusters(adata_sub, name, subtype_col, score_cols, out_dir):
    """UMAP + dotplot for subclusters."""
    sc.settings.figdir = out_dir

    # UMAP colored by subtype
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    sc.pl.umap(adata_sub, color=subtype_col, ax=axes[0], show=False,
               title=f"{name} subtypes", legend_loc="on data", legend_fontsize=7)
    leiden_key = f"{name.lower().replace(' ', '_')}_leiden"
    if leiden_key in adata_sub.obs.columns:
        sc.pl.umap(adata_sub, color=leiden_key, ax=axes[1], show=False,
                   title=f"{name} Leiden clusters")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"umap_{name.lower().replace(' ', '_')}_subtypes.png"),
                dpi=150, bbox_inches="tight")
    plt.close()
    print(f"    saved: umap_{name.lower().replace(' ', '_')}_subtypes.png")

    # Dotplot
    if score_cols:
        marker_genes = []
        for subtype, markers in (TCELL_MARKERS if "T" in name else MYELOID_MARKERS).items():
            marker_genes.extend(markers[:3])
        marker_genes = [g for g in list(dict.fromkeys(marker_genes))
                        if g in adata_sub.var_names][:30]
        if marker_genes and subtype_col in adata_sub.obs.columns:
            try:
                sc.pl.dotplot(adata_sub, var_names=marker_genes,
                              groupby=subtype_col, show=False,
                              save=f"_{name.lower().replace(' ', '_')}_markers.png")
                print(f"    saved: dotplot_{name.lower().replace(' ', '_')}_markers.png")
            except Exception as e:
                print(f"    dotplot failed: {e}")


def composition_by_progression(adata_sub, subtype_col, name, out_dir):
    """Bar chart: subtype fraction in slow vs fast progressors."""
    if "progression_category" not in adata_sub.obs.columns:
        return
    obs = adata_sub.obs[[subtype_col, "progression_category"]].dropna()
    comp = (obs.groupby(["progression_category", subtype_col]).size()
               .unstack(subtype_col, fill_value=0))
    comp = comp.div(comp.sum(axis=1), axis=0)

    fig, ax = plt.subplots(figsize=(max(8, len(comp.columns) * 0.8), 5))
    comp.T.plot(kind="bar", ax=ax, width=0.7)
    ax.set_title(f"{name} subtype composition: Slow vs Fast progressors")
    ax.set_ylabel("Fraction of cells")
    ax.set_xlabel("Subtype")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, f"{name.lower().replace(' ', '_')}_by_progression.png"),
                dpi=150, bbox_inches="tight")
    plt.close()
    print(f"    saved: {name.lower().replace(' ', '_')}_by_progression.png")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=os.path.join(
        PER_DIR, "korea_kim2022_annotated_scored.h5ad"))
    parser.add_argument("--output", default=os.path.join(
        PER_DIR, "korea_kim2022_subclustered.h5ad"))
    parser.add_argument("--tcell-res",   type=float, default=0.6)
    parser.add_argument("--myeloid-res", type=float, default=0.5)
    args = parser.parse_args()

    print(f"Loading {args.input} ...")
    adata = sc.read_h5ad(args.input)
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # Use existing HVG matrix (.X has 3000 HVGs). All major marker genes
    # (PDCD1, LAG3, HAVCR2, TOX, SPP1 etc.) are in the HVG set.
    print(f"  Using {adata.n_vars:,} HVG matrix for subclustering")

    ct_col = "cell_type_coarse"
    if ct_col not in adata.obs.columns:
        adata.obs[ct_col] = "cluster_" + adata.obs["leiden"].astype(str)

    adata.obs["tcell_subtype"]   = "Other"
    adata.obs["myeloid_subtype"] = "Other"

    # ── T/NK subclustering ────────────────────────────────────────────────────
    print("\n[T/NK subclustering]")
    tcell_mask = adata.obs[ct_col].isin(["T/NK", "T cell", "CD8+ T cell", "CD4+ T cell", "NK"])
    if tcell_mask.sum() > 500:
        adata_t = adata[tcell_mask].copy()
        adata_t = subcluster(adata_t, "tcell", resolution=args.tcell_res)
        adata_t, t_scores = score_subtypes(adata_t, TCELL_MARKERS, "tc")
        plot_subclusters(adata_t, "T cell", "tc_subtype", t_scores, OUT_DIR)
        composition_by_progression(adata_t, "tc_subtype", "Tcell", OUT_DIR)

        # Write back subtype labels to full object
        adata.obs.loc[tcell_mask, "tcell_subtype"] = adata_t.obs["tc_subtype"].values

        # Print distribution
        print(f"\n  T cell subtype distribution:")
        print(adata_t.obs["tc_subtype"].value_counts().to_string())

        # Save exhaustion trajectory data
        exhaustion_genes = ["PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "ENTPD1"]
        valid_ex = [g for g in exhaustion_genes if g in adata_t.var_names]
        if valid_ex:
            sc.tl.score_genes(adata_t, valid_ex, score_name="exhaustion_score_detailed")
        if "progression_category" in adata_t.obs.columns and "tc_subtype" in adata_t.obs.columns:
            ex_by_subtype = adata_t.obs.groupby("tc_subtype")["exhaustion_score_detailed"].mean()
            print(f"\n  Mean exhaustion score by T cell subtype:")
            print(ex_by_subtype.sort_values(ascending=False).to_string())
            ex_by_subtype.to_csv(os.path.join(OUT_DIR, "exhaustion_by_subtype.csv"))
    else:
        print(f"  Only {tcell_mask.sum()} T/NK cells — skipping")

    # ── Myeloid subclustering ─────────────────────────────────────────────────
    print("\n[Myeloid subclustering]")
    myeloid_mask = adata.obs[ct_col].isin(["Myeloid", "Macrophage", "Monocyte", "DC"])
    if myeloid_mask.sum() > 200:
        adata_m = adata[myeloid_mask].copy()
        adata_m = subcluster(adata_m, "myeloid", resolution=args.myeloid_res)
        adata_m, m_scores = score_subtypes(adata_m, MYELOID_MARKERS, "my")
        plot_subclusters(adata_m, "Myeloid", "my_subtype", m_scores, OUT_DIR)
        composition_by_progression(adata_m, "my_subtype", "Myeloid", OUT_DIR)

        adata.obs.loc[myeloid_mask, "myeloid_subtype"] = adata_m.obs["my_subtype"].values

        print(f"\n  Myeloid subtype distribution:")
        print(adata_m.obs["my_subtype"].value_counts().to_string())

        # SPP1 vs FOLR2 macro ratio by progression
        if "progression_category" in adata_m.obs.columns:
            spp1_mask  = adata_m.obs["my_subtype"] == "SPP1_Macro"
            folr2_mask = adata_m.obs["my_subtype"] == "FOLR2_Macro"
            by_prog = adata_m.obs.groupby("progression_category")["my_subtype"].value_counts(normalize=True)
            by_prog.to_csv(os.path.join(OUT_DIR, "myeloid_subtype_by_progression.csv"))
            print("\n  Myeloid subtypes by progression:")
            print(by_prog.unstack().to_string())
    else:
        print(f"  Only {myeloid_mask.sum()} myeloid cells — skipping")

    # ── Save ──────────────────────────────────────────────────────────────────
    print(f"\nSaving subclustered object to {args.output} ...")
    adata.write_h5ad(args.output)
    print(f"  Saved.")

    print(f"\n{'='*60}")
    print(f"Subclustering complete. Outputs in: {OUT_DIR}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
