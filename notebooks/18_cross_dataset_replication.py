"""
18_cross_dataset_replication.py
Validate SPP1+ macrophage and exhaustion findings across all 5 scRNA-seq datasets.

Strategy:
  1. Load integrated AnnData (gastric_meta_annotated.h5ad)
  2. Per dataset: compute SPP1/APOE/APOC1/TREM2 co-expression in Myeloid cells
     → proxy for SPP1_Macro abundance
  3. Per dataset: compute PDCD1/LAG3/HAVCR2/TOX co-expression in T cells
     → proxy for exhaustion level
  4. Correlate SPP1_Macro% with CAF marker expression
  5. Plot cross-dataset comparison bar charts

Also checks whether subclustered Korean cohort SPP1_Macro proportions
translate to marker gene expression in other datasets.

Outputs in outputs/cross_dataset/
"""
import os, warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc

warnings.filterwarnings("ignore")
sc.settings.verbosity = 0

BASE    = os.path.expanduser("~/gastric_tme_project")
INT_DIR = os.path.join(BASE, "data/processed/integrated")
OUT_DIR = os.path.join(BASE, "outputs/cross_dataset")
os.makedirs(OUT_DIR, exist_ok=True)

# Marker genes for each signature
SPP1_MACRO_MARKERS  = ["SPP1", "APOC1", "APOE", "C1QB", "TREM2", "GPNMB", "FN1"]
EXHAUSTION_MARKERS  = ["PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX"]
EFFECTOR_MARKERS    = ["GZMB", "PRF1", "IFNG", "GNLY"]
CAF_MARKERS         = ["ACTA2", "FAP", "PDGFRB", "POSTN", "MMP2"]
TREG_MARKERS        = ["FOXP3", "IL2RA", "IKZF2"]

DATASET_ORDER = ["korea_kim2022", "kumar2022", "zhang2021",
                 "diffuse_gc_2021", "tcell_exhaustion_2022"]
DATASET_LABELS = {
    "korea_kim2022":       "Korea Kim 2022\n(n=429k cells)",
    "kumar2022":           "Kumar 2022\n(n=159k cells)",
    "zhang2021":           "Zhang 2021\n(n=44k cells)",
    "diffuse_gc_2021":     "Diffuse GC 2021\n(n=23k cells)",
    "tcell_exhaustion_2022":"T-exh 2022\n(n=111k cells)",
}

MYELOID_KEYWORDS = ["myeloid", "macrophage", "monocyte", "dendritic", "dc", "macro"]
TCELL_KEYWORDS   = ["t_nk", "t/nk", "t cell", "nk", "cd8", "cd4", "treg", "exhaust"]
FIBRO_KEYWORDS   = ["fibroblast", "stromal", "caf", "mes"]


def score_genes(X, genes_in_data, markers):
    """Mean expression of marker genes in each cell (sparse-safe)."""
    present = [g for g in markers if g in genes_in_data]
    if not present:
        return np.zeros(X.shape[0])
    idx = [genes_in_data.index(g) for g in present]
    if sp.issparse(X):
        sub = np.asarray(X[:, idx].todense())
    else:
        sub = X[:, idx]
    return sub.mean(axis=1)


def fraction_detectable(scores):
    """Fraction of cells with detectable (>0) marker co-expression.

    Avoids two failure modes seen during development:
      - a per-dataset percentile cutoff is tautological (always ~25% of
        that dataset by construction, regardless of true abundance)
      - a z-score cutoff calibrated on one cohort and applied to another
        breaks down on zero-inflated/skewed distributions (e.g. zhang2021
        myeloid SPP1 scores are 76% exact zeros with a long right tail,
        so the reference z-threshold sat below the zero floor and
        "detected" 100% of cells)
    "Percent expressing" (score > 0) is calibration-free and standard
    in scRNA-seq marker analysis, so it's used here instead.
    """
    return (scores > 0).mean()


def main():
    # Load integrated data
    adata_path = os.path.join(INT_DIR, "gastric_meta_annotated.h5ad")
    if not os.path.exists(adata_path):
        print(f"ERROR: {adata_path} not found. Run annotation step first.")
        return

    print("Loading integrated AnnData ...")
    adata = sc.read_h5ad(adata_path)
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # Identify cell type column
    ct_col = None
    for col in ["cell_type", "cell_type_coarse", "cell_type_fine", "leiden", "cluster_label"]:
        if col in adata.obs.columns:
            ct_col = col
            break
    if ct_col is None:
        print("  ERROR: No cell type column found")
        return
    print(f"  Cell type column: {ct_col}")

    # Dataset column
    ds_col = None
    for col in ["dataset", "batch", "sample", "source"]:
        if col in adata.obs.columns:
            ds_col = col
            break
    if ds_col is None:
        print("  ERROR: No dataset column found")
        return
    print(f"  Dataset column: {ds_col}")
    print(f"  Datasets: {adata.obs[ds_col].value_counts().to_dict()}")

    # Gene names
    genes = adata.var_names.str.upper().tolist()
    print(f"  Genes available: {len(genes):,}")

    # Results storage
    results = []
    ct_vals = adata.obs[ct_col].str.lower().values
    ds_vals = adata.obs[ds_col].values

    for ds_id in DATASET_ORDER:
        ds_mask = ds_vals == ds_id
        n_cells = ds_mask.sum()
        if n_cells < 100:
            print(f"  {ds_id}: {n_cells} cells — skipping")
            continue
        print(f"\n  {ds_id}: {n_cells:,} cells")

        X_ds = adata.X[ds_mask]
        ct_ds = ct_vals[ds_mask]

        # Identify myeloid cells
        myeloid_mask = np.array([any(kw in ct for kw in MYELOID_KEYWORDS) for ct in ct_ds])
        n_myeloid = myeloid_mask.sum()
        print(f"    Myeloid: {n_myeloid:,} cells")

        if n_myeloid > 50:
            X_m = X_ds[myeloid_mask]
            spp1_scores = score_genes(X_m, genes, SPP1_MACRO_MARKERS)
            spp1_frac = fraction_detectable(spp1_scores)

            # Mean SPP1 expression
            spp1_idx = genes.index("SPP1") if "SPP1" in genes else None
            if spp1_idx is not None:
                if sp.issparse(X_m):
                    spp1_expr = np.asarray(X_m[:, spp1_idx].todense()).flatten().mean()
                else:
                    spp1_expr = X_m[:, spp1_idx].mean()
            else:
                spp1_expr = 0.0
        else:
            spp1_frac = np.nan
            spp1_expr = np.nan

        # Identify T/NK cells
        tcell_mask = np.array([any(kw in ct for kw in TCELL_KEYWORDS) for ct in ct_ds])
        n_tcell = tcell_mask.sum()
        print(f"    T/NK: {n_tcell:,} cells")

        if n_tcell > 50:
            X_t = X_ds[tcell_mask]
            exhaust_scores  = score_genes(X_t, genes, EXHAUSTION_MARKERS)
            effector_scores = score_genes(X_t, genes, EFFECTOR_MARKERS)
            exhaust_frac  = fraction_detectable(exhaust_scores)
            ex_ratio = exhaust_scores.mean() / (effector_scores.mean() + 1e-8)
        else:
            exhaust_frac = np.nan
            ex_ratio = np.nan

        # Fibroblast CAF signature
        fibro_mask = np.array([any(kw in ct for kw in FIBRO_KEYWORDS) for ct in ct_ds])
        n_fibro = fibro_mask.sum()
        if n_fibro > 20:
            X_f = X_ds[fibro_mask]
            caf_scores = score_genes(X_f, genes, CAF_MARKERS)
            caf_mean = caf_scores.mean()
        else:
            caf_mean = np.nan

        results.append({
            "dataset": ds_id,
            "label": DATASET_LABELS.get(ds_id, ds_id),
            "n_cells": int(n_cells),
            "n_myeloid": int(n_myeloid) if n_myeloid > 0 else 0,
            "n_tcell": int(n_tcell) if n_tcell > 0 else 0,
            "spp1_macro_frac": spp1_frac,
            "spp1_mean_expr": spp1_expr,
            "exhaustion_frac": exhaust_frac,
            "exhaustion_ratio": ex_ratio,
            "caf_mean": caf_mean,
            "frac_myeloid": n_myeloid / n_cells,
        })

    df = pd.DataFrame(results)
    df.to_csv(os.path.join(OUT_DIR, "cross_dataset_summary.csv"), index=False)
    print(f"\n{df[['dataset', 'n_cells', 'spp1_macro_frac', 'exhaustion_frac']].to_string(index=False)}")

    # Chi-square test: does SPP1+ macrophage / exhausted-T-cell abundance
    # (using the fixed, Korean-cohort-calibrated threshold) differ across cohorts?
    from scipy.stats import chi2_contingency
    for label, frac_col, n_col in [("SPP1+ macrophage", "spp1_macro_frac", "n_myeloid"),
                                    ("Exhausted T/NK", "exhaustion_frac", "n_tcell")]:
        valid = df.dropna(subset=[frac_col, n_col])
        valid = valid[valid[n_col] > 0]
        if len(valid) >= 2:
            high = (valid[frac_col] * valid[n_col]).round().astype(int)
            low  = valid[n_col].astype(int) - high
            table = np.vstack([high, low])
            chi2, p, dof, _ = chi2_contingency(table)
            print(f"  Chi-square ({label} fraction across cohorts): chi2={chi2:.2f}, dof={dof}, p={p:.2e}")

    # ---- Plots ----
    datasets_present = df["dataset"].tolist()
    labels = df["label"].tolist()
    n = len(df)
    x = np.arange(n)
    w = 0.6

    # SPP1+ macrophage fraction
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    ax = axes[0]
    vals = df["spp1_macro_frac"].fillna(0).values
    bars = ax.bar(x, vals, width=w, color="firebrick", alpha=0.8)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=7, rotation=15, ha="right")
    ax.set_ylabel("Fraction of Myeloid cells with detectable\nSPP1+ macrophage marker co-expression")
    ax.set_title("SPP1+ Macrophage Signature\nacross datasets")
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, val + 0.002, f"{val:.3f}",
                ha="center", va="bottom", fontsize=7)

    ax = axes[1]
    vals = df["exhaustion_frac"].fillna(0).values
    bars = ax.bar(x, vals, width=w, color="steelblue", alpha=0.8)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=7, rotation=15, ha="right")
    ax.set_ylabel("Fraction of T/NK cells with detectable\nexhaustion marker co-expression")
    ax.set_title("T Cell Exhaustion Signature\nacross datasets")
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, val + 0.002, f"{val:.3f}",
                ha="center", va="bottom", fontsize=7)

    ax = axes[2]
    vals = df["frac_myeloid"].fillna(0).values
    bars = ax.bar(x, vals, width=w, color="forestgreen", alpha=0.8)
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=7, rotation=15, ha="right")
    ax.set_ylabel("Myeloid cell fraction")
    ax.set_title("Overall Myeloid Proportion\nacross datasets")
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, val + 0.002, f"{val:.3f}",
                ha="center", va="bottom", fontsize=7)

    plt.suptitle("Cross-dataset replication: SPP1+ macrophage and exhaustion signatures\n"
                 "(5 gastric cancer scRNA-seq cohorts, n=766k cells total)",
                 fontsize=11, fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "cross_dataset_comparison.png"), dpi=150, bbox_inches="tight")
    plt.savefig(os.path.join(OUT_DIR, "cross_dataset_comparison.pdf"), bbox_inches="tight")
    plt.close()
    print("  saved: cross_dataset_comparison.png")

    # Correlation: SPP1_Macro vs myeloid fraction
    fig, ax = plt.subplots(figsize=(6, 5))
    valid = df.dropna(subset=["spp1_macro_frac", "frac_myeloid"])
    if len(valid) >= 3:
        ax.scatter(valid["frac_myeloid"], valid["spp1_macro_frac"],
                   s=100, c="firebrick", alpha=0.8, zorder=5)
        for _, row in valid.iterrows():
            ax.annotate(row["dataset"].replace("_", "\n"),
                        (row["frac_myeloid"], row["spp1_macro_frac"]),
                        fontsize=7, xytext=(5, 5), textcoords="offset points")
        from numpy.polynomial import polynomial as P
        coeffs = np.polyfit(valid["frac_myeloid"], valid["spp1_macro_frac"], 1)
        xfit = np.linspace(valid["frac_myeloid"].min(), valid["frac_myeloid"].max(), 100)
        ax.plot(xfit, np.polyval(coeffs, xfit), "k--", alpha=0.5)
        from scipy.stats import pearsonr
        r, p = pearsonr(valid["frac_myeloid"], valid["spp1_macro_frac"])
        ax.set_title(f"Myeloid fraction vs SPP1+ Macro proportion\nr={r:.2f}, p={p:.3f}")
    ax.set_xlabel("Myeloid cell fraction"); ax.set_ylabel("SPP1+ Macro fraction")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "spp1_myeloid_correlation.png"), dpi=150)
    plt.close()
    print("  saved: spp1_myeloid_correlation.png")

    print(f"\n{'='*60}")
    print(f"Cross-dataset replication complete. Outputs in: {OUT_DIR}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
