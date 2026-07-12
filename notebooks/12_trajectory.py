"""
12_trajectory.py
Pseudotime / trajectory analysis on T/NK cells using PAGA + diffusion pseudotime.
Tracks the naive → effector → exhausted arc and asks where responders fall.

Outputs in outputs/trajectory/:
  paga_graph_tcells.png
  umap_pseudotime.png
  pseudotime_by_progression.png
  pseudotime_by_timepoint.png
  pseudotime_vs_exhaustion.png
  paga_connectivity.csv
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
from scipy import stats

warnings.filterwarnings("ignore")

BASE    = os.path.expanduser("~/gastric_tme_project")
PER_DIR = os.path.join(BASE, "data/processed/per_dataset")
OUT_DIR = os.path.join(BASE, "outputs/trajectory")
os.makedirs(OUT_DIR, exist_ok=True)

NAIVE_MARKERS       = ["CCR7", "SELL", "TCF7", "LEF1"]
EXHAUSTION_MARKERS  = ["PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "ENTPD1", "CTLA4"]
EFFECTOR_MARKERS    = ["GZMB", "PRF1", "IFNG", "CX3CR1"]


def prepare_tcells(adata):
    """Subset to T/NK cells and re-process."""
    ct_col = "cell_type_coarse" if "cell_type_coarse" in adata.obs.columns else None
    if ct_col:
        tcell_labels = ["T/NK", "T cell", "CD8+ T cell", "CD4+ T cell",
                        "NK", "Terminal_Exhausted", "Progenitor_Exhausted",
                        "Effector_CD8", "Naive_T", "Treg"]
        mask = adata.obs[ct_col].isin(tcell_labels)
        # Also use subtype if available
        if "tcell_subtype" in adata.obs.columns:
            mask |= adata.obs["tcell_subtype"] != "Other"
    else:
        # Fall back to marker-based rough selection
        genes = ["CD3D", "CD3E", "CD3G", "TRAC"]
        valid = [g for g in genes if g in adata.var_names]
        if valid:
            import scipy.sparse as sp
            X = adata[:, valid].X
            X = X.toarray() if sp.issparse(X) else X
            mask = (X > 0).any(axis=1)
        else:
            mask = pd.Series(True, index=adata.obs_names)

    adata_t = adata[mask].copy()
    print(f"  T/NK subset: {adata_t.n_obs:,} cells")

    if adata_t.n_obs < 200:
        return None

    # Re-process
    sc.pp.highly_variable_genes(adata_t, n_top_genes=2000, flavor="seurat", subset=False)
    sc.tl.pca(adata_t, n_comps=20, use_highly_variable=True)
    sc.pp.neighbors(adata_t, n_pcs=20, n_neighbors=15)
    sc.tl.umap(adata_t)

    return adata_t


def score_programs(adata_t):
    """Score naive, effector, exhaustion programs."""
    var = set(adata_t.var_names)
    for name, genes in [("naive",      NAIVE_MARKERS),
                         ("effector",   EFFECTOR_MARKERS),
                         ("exhaustion", EXHAUSTION_MARKERS)]:
        valid = [g for g in genes if g in var]
        if len(valid) >= 2:
            sc.tl.score_genes(adata_t, valid, score_name=f"{name}_score")
            print(f"  {name}: scored {len(valid)}/{len(genes)} genes")
    return adata_t


def run_paga(adata_t, cluster_col):
    """Run PAGA on T cell subset."""
    if cluster_col not in adata_t.obs.columns:
        sc.tl.leiden(adata_t, resolution=0.4, key_added=cluster_col)
    sc.tl.paga(adata_t, groups=cluster_col)
    return adata_t


def run_dpt(adata_t, cluster_col):
    """Run diffusion pseudotime starting from the most naive-like cluster."""
    if "naive_score" not in adata_t.obs.columns:
        return adata_t

    # Find root cell = cell with highest naive score (use positional index, not obs_name)
    root_int = int(np.argmax(adata_t.obs["naive_score"].values))
    adata_t.uns["iroot"] = root_int

    sc.tl.diffmap(adata_t)
    sc.tl.dpt(adata_t)
    print("  Diffusion pseudotime computed.")
    return adata_t


def plot_trajectory(adata_t, cluster_col, out_dir):
    """PAGA graph + UMAP colored by pseudotime + progression."""
    sc.settings.figdir = out_dir

    # PAGA graph
    try:
        fig, ax = plt.subplots(figsize=(7, 6))
        sc.pl.paga(adata_t, ax=ax, show=False, node_size_scale=2,
                   title="PAGA — T cell trajectory")
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, "paga_graph_tcells.png"), dpi=150)
        plt.close()
        print("  saved: paga_graph_tcells.png")
    except Exception as e:
        print(f"  PAGA plot failed: {e}")

    # Save connectivity
    if "paga" in adata_t.uns:
        conn = pd.DataFrame(adata_t.uns["paga"]["connectivities"].toarray(),
                            index=adata_t.obs[cluster_col].cat.categories,
                            columns=adata_t.obs[cluster_col].cat.categories)
        conn.to_csv(os.path.join(out_dir, "paga_connectivity.csv"))

    color_keys = ["dpt_pseudotime", "naive_score", "exhaustion_score",
                  "effector_score", cluster_col]
    color_keys = [c for c in color_keys if c in adata_t.obs.columns]

    if color_keys:
        n = len(color_keys)
        fig, axes = plt.subplots(1, n, figsize=(5 * n, 5))
        if n == 1:
            axes = [axes]
        for ax, key in zip(axes, color_keys):
            sc.pl.umap(adata_t, color=key, ax=ax, show=False,
                       title=key.replace("_", " "), frameon=False)
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, "umap_pseudotime.png"), dpi=150, bbox_inches="tight")
        plt.close()
        print("  saved: umap_pseudotime.png")


def pseudotime_by_group(adata_t, out_dir):
    """Violin plots of pseudotime by progression category and timepoint."""
    if "dpt_pseudotime" not in adata_t.obs.columns:
        return

    for group_col, fname in [("progression_category", "pseudotime_by_progression.png"),
                               ("timepoint_label",      "pseudotime_by_timepoint.png")]:
        if group_col not in adata_t.obs.columns:
            continue
        obs = adata_t.obs[[group_col, "dpt_pseudotime", "exhaustion_score"]].dropna()
        groups = obs[group_col].unique()
        if len(groups) < 2:
            continue

        fig, axes = plt.subplots(1, 2, figsize=(10, 5))
        for ax, col in zip(axes, ["dpt_pseudotime", "exhaustion_score"]):
            if col not in obs.columns:
                continue
            order = sorted(groups)
            data = [obs[obs[group_col] == g][col].values for g in order]
            ax.violinplot(data, positions=range(len(order)), showmedians=True)
            ax.set_xticks(range(len(order)))
            ax.set_xticklabels(order, rotation=30, ha="right")
            ax.set_ylabel(col.replace("_", " "))
            ax.set_title(f"{col} by {group_col}")

            # Mann-Whitney p-values for pairs
            if len(order) == 2:
                _, p = stats.mannwhitneyu(data[0], data[1], alternative="two-sided")
                ax.set_xlabel(f"p={p:.4f}")

        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, fname), dpi=150)
        plt.close()
        print(f"  saved: {fname}")

    # Scatter: pseudotime vs exhaustion score
    if "exhaustion_score" in adata_t.obs.columns:
        obs2 = adata_t.obs[["dpt_pseudotime", "exhaustion_score"]].dropna()
        rho, p = stats.spearmanr(obs2["dpt_pseudotime"], obs2["exhaustion_score"])
        fig, ax = plt.subplots(figsize=(6, 5))
        ax.hexbin(obs2["dpt_pseudotime"], obs2["exhaustion_score"],
                  gridsize=50, cmap="YlOrRd", mincnt=1)
        ax.set_xlabel("Diffusion pseudotime")
        ax.set_ylabel("Exhaustion score")
        ax.set_title(f"Pseudotime vs Exhaustion\nSpearman ρ={rho:.3f}, p={p:.2e}")
        plt.colorbar(ax.collections[0], ax=ax, label="Cell count")
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, "pseudotime_vs_exhaustion.png"), dpi=150)
        plt.close()
        print("  saved: pseudotime_vs_exhaustion.png")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=os.path.join(
        PER_DIR, "korea_kim2022_annotated_scored.h5ad"))
    parser.add_argument("--cluster-col", default="tc_subtype",
                        help="Column to use for PAGA groups (falls back to new Leiden)")
    args = parser.parse_args()

    print(f"Loading {args.input} ...")
    adata = sc.read_h5ad(args.input)
    print(f"  {adata.n_obs:,} cells")

    # Use the existing HVG matrix (.X has 3000 HVGs) — marker genes for T cell
    # exhaustion (PDCD1, LAG3, TOX etc.) are included in the HVG set

    print("\n[Preparing T/NK subset]")
    adata_t = prepare_tcells(adata)
    if adata_t is None:
        print("Too few T/NK cells — exiting.")
        return

    print("\n[Scoring gene programs]")
    adata_t = score_programs(adata_t)

    cluster_col = args.cluster_col
    if cluster_col not in adata_t.obs.columns:
        cluster_col = "trajectory_leiden"

    print("\n[PAGA]")
    adata_t = run_paga(adata_t, cluster_col)

    print("\n[Diffusion pseudotime]")
    adata_t = run_dpt(adata_t, cluster_col)

    print("\n[Plotting]")
    plot_trajectory(adata_t, cluster_col, OUT_DIR)
    pseudotime_by_group(adata_t, OUT_DIR)

    # Save T cell object with pseudotime
    out_path = os.path.join(PER_DIR, "korea_tcell_trajectory.h5ad")
    adata_t.write_h5ad(out_path)
    print(f"\nSaved T cell trajectory object: {out_path}")

    print(f"\n{'='*60}")
    print(f"Trajectory analysis complete. Outputs in: {OUT_DIR}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
