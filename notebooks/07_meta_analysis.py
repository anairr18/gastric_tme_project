"""
07_meta_analysis.py
Meta-analysis of the integrated, annotated gastric TME atlas.

Analyses:
  A. TME composition
     - Cell type fractions per patient per dataset
     - Composition differences: Slow vs Fast progressors (Korean cohort)
     - Composition differences: immunotherapy response (PR vs PD vs SD)

  B. Pseudobulk differential expression (pydeseq2)
     - Per cell type: Slow Progressor vs Fast Progressor
     - Per cell type: baseline vs FU1 vs FU2 (longitudinal)

  C. CD8+ T cell exhaustion trajectory
     - Score exhaustion gene program per cell
     - Compare exhaustion across response groups and timepoints

  D. Macrophage polarization
     - M1/M2 score per macrophage
     - Compare polarization in responders vs non-responders

  E. Biomarker summary
     - Top features predicting immunotherapy response (across cohorts with treatment data)
     - Survival-associated cell type abundances (PFS/OS correlation)

Outputs in outputs/meta_analysis/:
  tme_composition_heatmap.png
  composition_slow_vs_fast.png
  cd8_exhaustion_scores.png
  macrophage_polarization.png
  pseudobulk_de_<cell_type>.csv (one per major cell type)
  biomarker_summary.png

Usage:
    python 07_meta_analysis.py
    python 07_meta_analysis.py --analyses composition de exhaustion macrophage
"""
import sys, os
import argparse
import warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.stats as stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc

warnings.filterwarnings("ignore")

BASE      = os.path.expanduser("~/gastric_tme_project")
INT_DIR   = os.path.join(BASE, "data/processed/integrated")
MA_DIR    = os.path.join(BASE, "outputs/meta_analysis")
MANIFEST  = os.path.join(BASE, "data/external/dataset_manifest.csv")

os.makedirs(MA_DIR, exist_ok=True)
sc.settings.figdir = MA_DIR

# Exhaustion gene program (TOX-driven, from published gastric / pan-cancer refs)
EXHAUSTION_GENES = ["PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "TOX2",
                     "ENTPD1", "CD160", "CTLA4", "BTLA", "CD244"]

M1_GENES = ["CD80", "CD86", "IL1B", "TNF", "CXCL9", "CXCL10", "NOS2", "IL6",
             "FCGR1A", "FCGR3A", "HLA-DRA", "HLA-DRB1"]
M2_GENES = ["MRC1", "CD163", "IL10", "TGM2", "FOLR2", "LYVE1", "SELENOP",
             "CCL18", "CCL22", "CD200R1", "PPARG", "ARG1"]


# ─── A: TME Composition ───────────────────────────────────────────────────────

def compute_composition(adata, cell_type_col="cell_type_coarse",
                         patient_col="patient", dataset_col="dataset_id"):
    """
    Return a DataFrame: rows = patients, columns = cell types, values = fraction.
    """
    # Deduplicate groupby keys (patient_col may equal dataset_col when no patient metadata)
    group_cols = list(dict.fromkeys([patient_col, dataset_col]))
    obs = adata.obs[group_cols + [cell_type_col]].copy()
    obs = obs.dropna(subset=[patient_col, cell_type_col])

    comp = (obs.groupby(group_cols + [cell_type_col])
               .size()
               .unstack(cell_type_col, fill_value=0))
    comp = comp.div(comp.sum(axis=1), axis=0)
    comp.columns.name = "cell_type"
    return comp.reset_index()


def plot_composition_heatmap(comp, outdir, row_col="patient"):
    meta_cols = {"patient", "dataset_id", row_col}
    ct_cols = [c for c in comp.columns if c not in meta_cols]
    row_col = row_col if row_col in comp.columns else comp.columns[0]
    mat = comp.set_index(row_col)[ct_cols]
    plt.figure(figsize=(max(10, len(ct_cols) * 0.8), max(8, len(mat) * 0.2)))
    sns.heatmap(mat, cmap="YlOrRd", vmin=0, vmax=1, yticklabels=True,
                linewidths=0.2, linecolor="white")
    plt.title("TME Composition per Sample")
    plt.xlabel("Cell Type")
    plt.ylabel(row_col)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "tme_composition_heatmap.png"), dpi=150)
    plt.close()
    print(f"  saved: tme_composition_heatmap.png")


def plot_composition_by_response(adata, comp, outdir, patient_col="patient"):
    if "progression_category" not in adata.obs.columns:
        print("  no progression_category column — skipping response comparison")
        return
    if patient_col not in adata.obs.columns:
        print(f"  no {patient_col} column — skipping response comparison")
        return

    prog = (adata.obs[[patient_col, "progression_category"]]
            .drop_duplicates(patient_col)
            .set_index(patient_col))
    row_col = patient_col if patient_col in comp.columns else comp.columns[0]
    comp2 = comp.set_index(row_col).join(prog, how="left").reset_index()
    comp2 = comp2.dropna(subset=["progression_category"])

    meta_cols = {patient_col, "dataset_id", "progression_category", row_col}
    ct_cols = [c for c in comp2.columns if c not in meta_cols]

    fig, axes = plt.subplots(2, max(1, (len(ct_cols) + 1) // 2),
                              figsize=(16, 8), sharey=False)
    axes = axes.flatten()
    for i, ct in enumerate(ct_cols):
        if i >= len(axes):
            break
        grp = comp2.groupby("progression_category")[ct]
        slow = grp.get_group("Slow") if "Slow" in grp.groups else pd.Series(dtype=float)
        fast = grp.get_group("Fast") if "Fast" in grp.groups else pd.Series(dtype=float)
        ax = axes[i]
        ax.boxplot([slow.values, fast.values], labels=["Slow", "Fast"])
        ax.set_title(ct, fontsize=8)
        ax.set_ylabel("Fraction")
        if len(slow) > 2 and len(fast) > 2:
            _, pval = stats.mannwhitneyu(slow, fast, alternative="two-sided")
            ax.set_xlabel(f"p={pval:.3f}", fontsize=7)

    for j in range(i + 1, len(axes)):
        axes[j].set_visible(False)

    plt.suptitle("TME Composition: Slow vs Fast Progressors", y=1.01)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "composition_slow_vs_fast.png"), dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  saved: composition_slow_vs_fast.png")


# ─── B: Pseudobulk DE ────────────────────────────────────────────────────────

def pseudobulk_de(adata, cell_type_col="cell_type_coarse",
                   patient_col="patient", condition_col="progression_category",
                   target_cell_types=None, outdir=MA_DIR):
    """
    Pseudobulk differential expression using pydeseq2.
    For each cell type, aggregate counts per patient, then run DESeq2.
    """
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError:
        print("  pydeseq2 not available — skipping DE analysis")
        return {}

    if condition_col not in adata.obs.columns:
        print(f"  {condition_col} not in obs — skipping DE")
        return {}

    de_outdir = os.path.join(outdir, "pseudobulk_de")
    os.makedirs(de_outdir, exist_ok=True)

    if target_cell_types is None:
        target_cell_types = adata.obs[cell_type_col].value_counts().head(8).index.tolist()

    results = {}
    for ct in target_cell_types:
        print(f"  DE for cell type: {ct} ...")
        ct_mask = adata.obs[cell_type_col] == ct
        adata_ct = adata[ct_mask].copy()

        if adata_ct.n_obs < 50:
            print(f"    SKIP: only {adata_ct.n_obs} cells")
            continue

        # Pseudobulk: sum raw counts per patient
        obs_ct = adata_ct.obs[[patient_col, condition_col]].copy()
        obs_ct = obs_ct.dropna()
        valid_patients = obs_ct[patient_col].unique()

        if len(valid_patients) < 4:
            print(f"    SKIP: only {len(valid_patients)} patients with data")
            continue

        # Get raw counts — stay sparse to avoid OOM on large cell type subsets
        if adata_ct.raw is not None:
            X_sparse = adata_ct.raw.to_adata().X
        else:
            X_sparse = adata_ct.X
        if not sp.issparse(X_sparse):
            X_sparse = sp.csr_matrix(X_sparse)
        else:
            X_sparse = X_sparse.tocsr()

        # Sum per patient without materializing the full dense matrix
        pb_counts = {}
        pb_meta   = {}
        pat_col_vals = adata_ct.obs[patient_col].values
        for pat in valid_patients:
            mask = pat_col_vals == pat
            row_sum = np.asarray(X_sparse[mask].sum(axis=0)).flatten().astype(int)
            pb_counts[pat] = row_sum
            cond = obs_ct.loc[obs_ct[patient_col] == pat, condition_col].iloc[0]
            pb_meta[pat] = {condition_col: cond}

        count_df = pd.DataFrame(pb_counts, index=adata_ct.var_names).T
        meta_df  = pd.DataFrame(pb_meta).T
        meta_df[condition_col] = meta_df[condition_col].astype(str)

        # Filter low-count genes
        count_df = count_df.loc[:, count_df.sum(axis=0) > 10]
        if count_df.shape[1] < 100:
            print(f"    SKIP: too few genes after filtering ({count_df.shape[1]})")
            continue

        try:
            dds = DeseqDataSet(
                counts=count_df,
                metadata=meta_df,
                design_factors=condition_col,
            )
            dds.deseq2()
            # pydeseq2 >=0.4 requires explicit contrast
            conds = meta_df[condition_col].unique().tolist()
            contrast = (condition_col, conds[0], conds[1]) if len(conds) >= 2 else None
            ds = DeseqStats(dds, contrast=contrast) if contrast else DeseqStats(dds)
            ds.summary()
            de_df = ds.results_df.sort_values("padj")
            safe_ct = ct.replace("/", "_").replace(" ", "_")
            out_path = os.path.join(de_outdir, f"de_{safe_ct}_slow_vs_fast.csv")
            de_df.to_csv(out_path)
            results[ct] = de_df
            n_sig = (de_df["padj"] < 0.05).sum()
            print(f"    {n_sig} significant genes (padj < 0.05)")
        except Exception as e:
            print(f"    DE failed: {e}")

    return results


# ─── C: CD8+ T cell exhaustion ───────────────────────────────────────────────

def _ensure_lognorm(adata):
    """Normalize X in-place if it contains raw counts. Returns True if normalized."""
    x = adata.X
    max_val = float(x.data.max()) if sp.issparse(x) and x.nnz > 0 else float(x.max())
    if max_val > 50:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        return True
    return False


def _restore_counts(adata, was_normalized):
    """Restore raw counts to X from layers["counts"] if normalization was applied."""
    if was_normalized and "counts" in adata.layers:
        adata.X = adata.layers["counts"]


def score_exhaustion(adata):
    valid_genes = [g for g in EXHAUSTION_GENES if g in adata.var_names]
    if len(valid_genes) < 3:
        print(f"  only {len(valid_genes)} exhaustion genes in dataset — skipping")
        return adata
    was_norm = _ensure_lognorm(adata)
    sc.tl.score_genes(adata, gene_list=valid_genes, score_name="exhaustion_score")
    _restore_counts(adata, was_norm)
    print(f"  exhaustion scored using {len(valid_genes)}/{len(EXHAUSTION_GENES)} genes")
    return adata


def plot_exhaustion(adata, outdir):
    if "exhaustion_score" not in adata.obs.columns:
        return

    cd8_mask = adata.obs.get("cell_type_coarse", pd.Series("", index=adata.obs.index)).str.contains("T/NK")
    adata_cd8 = adata[cd8_mask.values]

    if adata_cd8.n_obs < 10:
        print("  too few T/NK cells for exhaustion plot")
        return

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    # By progression category
    if "progression_category" in adata_cd8.obs.columns:
        prog_data = adata_cd8.obs.groupby("progression_category")["exhaustion_score"]
        groups = [prog_data.get_group(g) for g in prog_data.groups]
        labels = list(prog_data.groups.keys())
        axes[0].violinplot(groups, positions=range(len(groups)), showmedians=True)
        axes[0].set_xticks(range(len(groups)))
        axes[0].set_xticklabels(labels)
        axes[0].set_title("CD8 Exhaustion by Progression")
        axes[0].set_ylabel("Exhaustion Score")

    # By timepoint
    if "timepoint_label" in adata_cd8.obs.columns:
        tp_order = ["Baseline", "FU1", "FU2"]
        tp_data  = adata_cd8.obs.groupby("timepoint_label")["exhaustion_score"]
        groups   = [tp_data.get_group(t) for t in tp_order if t in tp_data.groups]
        labels   = [t for t in tp_order if t in tp_data.groups]
        axes[1].violinplot(groups, positions=range(len(groups)), showmedians=True)
        axes[1].set_xticks(range(len(groups)))
        axes[1].set_xticklabels(labels)
        axes[1].set_title("CD8 Exhaustion by Timepoint")
        axes[1].set_ylabel("Exhaustion Score")

    plt.suptitle("T Cell Exhaustion Dynamics")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "cd8_exhaustion_scores.png"), dpi=150)
    plt.close()
    print(f"  saved: cd8_exhaustion_scores.png")


# ─── D: Macrophage polarization ──────────────────────────────────────────────

def score_macrophage_polarization(adata):
    valid_m1 = [g for g in M1_GENES if g in adata.var_names]
    valid_m2 = [g for g in M2_GENES if g in adata.var_names]
    was_norm = _ensure_lognorm(adata)
    if len(valid_m1) >= 3:
        sc.tl.score_genes(adata, gene_list=valid_m1, score_name="M1_score")
    if len(valid_m2) >= 3:
        sc.tl.score_genes(adata, gene_list=valid_m2, score_name="M2_score")
    _restore_counts(adata, was_norm)
    if "M1_score" in adata.obs and "M2_score" in adata.obs:
        adata.obs["M1_M2_ratio"] = adata.obs["M1_score"] - adata.obs["M2_score"]
    print(f"  M1 scored: {len(valid_m1)} genes; M2 scored: {len(valid_m2)} genes")
    return adata


def plot_macrophage_polarization(adata, outdir):
    mac_mask = adata.obs.get("cell_type_coarse", pd.Series("", index=adata.obs.index)) == "Myeloid"
    adata_mac = adata[mac_mask.values]
    if adata_mac.n_obs < 10 or "M1_M2_ratio" not in adata_mac.obs.columns:
        return

    fig, ax = plt.subplots(figsize=(8, 5))
    if "progression_category" in adata_mac.obs.columns:
        groups = adata_mac.obs.groupby("progression_category")["M1_M2_ratio"]
        for i, (name, g) in enumerate(groups):
            ax.violinplot([g.values], positions=[i], showmedians=True)
        ax.set_xticks(range(len(groups)))
        ax.set_xticklabels([n for n, _ in groups])
        ax.set_xlabel("Progression Category")
    ax.set_ylabel("M1 - M2 Score")
    ax.set_title("Macrophage Polarization (M1 vs M2)")
    ax.axhline(0, color="gray", linestyle="--", linewidth=0.8)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "macrophage_polarization.png"), dpi=150)
    plt.close()
    print(f"  saved: macrophage_polarization.png")


# ─── E: Biomarker summary ─────────────────────────────────────────────────────

def survival_correlation(adata, outdir):
    """Correlate cell type abundances with PFS and OS (spearman)."""
    required = ["patient", "PFS_days", "OS_days", "cell_type_coarse"]
    if not all(c in adata.obs.columns for c in required):
        print("  survival correlation: missing columns — skipping")
        return

    comp = compute_composition(adata)
    comp = comp.set_index("patient")
    ct_cols = [c for c in comp.columns if c not in ["dataset_id"]]

    surv = (adata.obs[["patient", "PFS_days", "OS_days"]]
            .drop_duplicates("patient").set_index("patient"))
    merged = comp.join(surv, how="inner").dropna()

    rows = []
    for ct in ct_cols:
        for endpoint in ["PFS_days", "OS_days"]:
            if endpoint not in merged.columns:
                continue
            valid = merged[[ct, endpoint]].dropna()
            if len(valid) < 5:
                continue
            rho, p = stats.spearmanr(valid[ct], valid[endpoint])
            rows.append({"cell_type": ct, "endpoint": endpoint,
                         "spearman_rho": rho, "p_value": p})

    if not rows:
        return

    corr_df = pd.DataFrame(rows).sort_values("p_value")
    corr_df.to_csv(os.path.join(outdir, "survival_correlation.csv"), index=False)

    # Dot plot
    corr_pivot = corr_df.pivot(index="cell_type", columns="endpoint", values="spearman_rho")
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(corr_pivot, cmap="RdBu_r", center=0, annot=True, fmt=".2f",
                linewidths=0.5, ax=ax)
    ax.set_title("Cell Type Abundance vs Survival (Spearman ρ)")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, "survival_correlation.png"), dpi=150)
    plt.close()
    print(f"  saved: survival_correlation.png")


# ─── main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=os.path.join(INT_DIR, "gastric_meta_annotated.h5ad"),
                        help="annotated integrated h5ad")
    parser.add_argument("--analyses", nargs="+",
                        choices=["composition", "de", "exhaustion", "macrophage", "survival"],
                        default=["composition", "exhaustion", "macrophage", "survival"],
                        help="which analyses to run")
    parser.add_argument("--include-de", action="store_true",
                        help="include pseudobulk DE (slow, can take 30+ min)")
    args = parser.parse_args()

    # Fallback to unannotated integrated object if annotated doesn't exist
    if not os.path.exists(args.input):
        fallback = os.path.join(INT_DIR, "gastric_meta_integrated.h5ad")
        if os.path.exists(fallback):
            print(f"Annotated h5ad not found; falling back to {fallback}")
            args.input = fallback
        else:
            # Final fallback: run on single Korean cohort processed object
            fallback2 = os.path.join(BASE, "data/processed/per_dataset/korea_kim2022_processed.h5ad")
            if os.path.exists(fallback2):
                print(f"Using single-dataset fallback: {fallback2}")
                args.input = fallback2
            else:
                fallback3 = os.path.join(BASE, "data/processed/gastric_processed.h5ad")
                if os.path.exists(fallback3):
                    print(f"Using original processed object: {fallback3}")
                    args.input = fallback3
                else:
                    print(f"No processed h5ad found. Run earlier pipeline steps first.")
                    sys.exit(1)

    print(f"Loading {args.input} ...")
    adata = sc.read_h5ad(args.input)
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print(f"  obs columns: {adata.obs.columns.tolist()}")

    # If only HVGs are in .X but .raw has the full gene set, restore it for scoring
    if adata.raw is not None and adata.raw.n_vars > adata.n_vars:
        print(f"  restoring full gene set from .raw ({adata.raw.n_vars:,} genes)")
        obs_backup = adata.obs.copy()
        adata = adata.raw.to_adata()
        for col in obs_backup.columns:
            if col not in adata.obs.columns:
                adata.obs[col] = obs_backup[col].values
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        print(f"  gene set restored and log-normalised: {adata.n_vars:,} genes")

    # Flatten any obs columns that got stored as 2D arrays during h5ad round-trip
    for col in adata.obs.columns:
        if adata.obs[col].values.ndim > 1:
            adata.obs[col] = adata.obs[col].values.ravel()

    # Map cell_type_coarse if annotation ran; otherwise use leiden as proxy
    if "cell_type_coarse" not in adata.obs.columns:
        print("  cell_type_coarse not found — using leiden clusters as proxy cell types")
        adata.obs["cell_type_coarse"] = "cluster_" + adata.obs["leiden"].astype(str)

    # ── score genes (exhaustion + macrophage) first so plots can use them ────
    if "exhaustion" in args.analyses or "de" in args.analyses:
        print("\n[Exhaustion scoring]")
        adata = score_exhaustion(adata)

    if "macrophage" in args.analyses:
        print("\n[Macrophage polarization scoring]")
        adata = score_macrophage_polarization(adata)

    # ── A: Composition ───────────────────────────────────────────────────────
    if "composition" in args.analyses:
        print("\n[TME Composition]")
        patient_col = "patient" if "patient" in adata.obs.columns else "dataset_id"
        comp = compute_composition(adata, patient_col=patient_col)
        plot_composition_heatmap(comp, MA_DIR, row_col=patient_col)
        plot_composition_by_response(adata, comp, MA_DIR, patient_col=patient_col)

    # ── B: DE ────────────────────────────────────────────────────────────────
    if args.include_de or "de" in args.analyses:
        print("\n[Pseudobulk DE]")
        pseudobulk_de(adata)

    # ── C: Exhaustion ────────────────────────────────────────────────────────
    if "exhaustion" in args.analyses:
        print("\n[CD8 Exhaustion plots]")
        plot_exhaustion(adata, MA_DIR)

    # ── D: Macrophage ────────────────────────────────────────────────────────
    if "macrophage" in args.analyses:
        print("\n[Macrophage polarization plots]")
        plot_macrophage_polarization(adata, MA_DIR)

    # ── E: Survival ──────────────────────────────────────────────────────────
    if "survival" in args.analyses:
        print("\n[Survival correlation]")
        survival_correlation(adata, MA_DIR)

    # Save object with scores added
    scored_out = args.input.replace(".h5ad", "_scored.h5ad")
    print(f"\nSaving scored object to {scored_out} ...")
    adata.write_h5ad(scored_out)

    print(f"\n{'='*60}")
    print("Meta-analysis complete.")
    print(f"Outputs in: {MA_DIR}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
