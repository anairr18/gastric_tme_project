import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", line_buffering=True)
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", line_buffering=True)

import os
import warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import scanpy as sc
import scrublet as scr

warnings.filterwarnings("ignore")

# paths
BASE      = os.path.expanduser("~/gastric_tme_project")
RAW_H5AD  = os.path.join(BASE, "data/raw/ge_korea_raw_data_count_matricies_raw_combined.h5ad")
CKPT_H5AD = os.path.join(BASE, "data/processed/checkpoint_post_scrublet.h5ad")
OUT_H5AD  = os.path.join(BASE, "data/processed/gastric_processed.h5ad")
QC_DIR    = os.path.join(BASE, "outputs/qc_plots")
TABLE_S1  = os.path.join(BASE, "data/raw/Table S1.xlsx")
TABLE_S2  = os.path.join(BASE, "data/raw/Table S2.xlsx")

os.makedirs(QC_DIR, exist_ok=True)
os.makedirs(os.path.join(BASE, "data/processed"), exist_ok=True)
sc.settings.figdir = QC_DIR
sc.settings.verbosity = 2

# skip steps 1-5 if checkpoint exists
if os.path.exists(CKPT_H5AD):
    print("loading post-scrublet checkpoint ...")
    adata = sc.read_h5ad(CKPT_H5AD)
    print(f"loaded: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")
    print(f"obs columns: {adata.obs.columns.tolist()}")

else:
    # step 1: load raw data
    print("\nstep 1: loading raw data")
    adata = sc.read_h5ad(RAW_H5AD)
    print(f"shape: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")
    print(f"obs columns: {adata.obs.columns.tolist()}")
    print(f"var columns: {adata.var.columns.tolist()}")
    print(f"uns keys: {list(adata.uns.keys())}")
    print(f"\nfirst 5 obs rows:\n{adata.obs.head()}")
    print(f"\nfirst 5 var rows:\n{adata.var.head()}")

    # step 2: QC metrics
    print("\nstep 2: QC metrics")
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    n_mt = adata.var["mt"].sum()
    print(f"mitochondrial genes: {n_mt}")

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )
    print(f"n_genes_by_counts -- median: {adata.obs['n_genes_by_counts'].median():.0f}, "
          f"mean: {adata.obs['n_genes_by_counts'].mean():.0f}")
    print(f"total_counts      -- median: {adata.obs['total_counts'].median():.0f}, "
          f"mean: {adata.obs['total_counts'].mean():.0f}")
    print(f"pct_counts_mt     -- median: {adata.obs['pct_counts_mt'].median():.2f}%, "
          f"mean: {adata.obs['pct_counts_mt'].mean():.2f}%")

    # step 3: QC plots
    print("\nstep 3: saving QC plots")
    sc.pl.violin(
        adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4, save="_qc_violin.png", show=False,
    )
    sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt",
                  save="_mt_scatter.png", show=False)
    sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts",
                  save="_genes_scatter.png", show=False)

    # step 4: filter cells & genes
    print("\nstep 4: filtering")
    n_before = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_cells(adata, max_genes=6000)
    adata = adata[adata.obs["pct_counts_mt"] < 20].copy()
    n_after_cells = adata.n_obs
    print(f"cells: {n_before:,} -> {n_after_cells:,} (removed {n_before - n_after_cells:,})")

    g_before = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=3)
    g_after = adata.n_vars
    print(f"genes: {g_before:,} -> {g_after:,} (removed {g_before - g_after:,})")

    # step 5: scrublet doublet detection (per sample)
    print("\nstep 5: scrublet (per sample)")
    doublet_scores = np.zeros(adata.n_obs, dtype=float)
    predicted_dbl  = np.zeros(adata.n_obs, dtype=bool)

    samples = adata.obs["sample"].unique()
    print(f"running on {len(samples)} samples ...")

    for samp in samples:
        mask = (adata.obs["sample"] == samp).values
        n_cells = mask.sum()
        if n_cells < 50:
            print(f"  [{samp}]  {n_cells} cells -- skipping (< 50)")
            continue
        X_samp = adata.X[mask]
        if sp.issparse(X_samp):
            X_samp = X_samp.toarray()
        try:
            scrub = scr.Scrublet(X_samp, random_state=42)
            scores, pred = scrub.scrub_doublets(verbose=False)
            doublet_scores[mask] = scores
            predicted_dbl[mask]  = pred
            n_pred = pred.sum()
            print(f"  [{samp}]  {n_cells} cells, {n_pred} predicted doublets "
                  f"({100*n_pred/n_cells:.1f}%)")
        except Exception as e:
            print(f"  [{samp}]  scrublet failed: {e}")

    adata.obs["doublet_score"]     = doublet_scores
    adata.obs["predicted_doublet"] = predicted_dbl
    print(f"total doublets: {predicted_dbl.sum():,} / {adata.n_obs:,} ({100*predicted_dbl.mean():.2f}%)")

    print(f"\nsaving checkpoint to {CKPT_H5AD} ...")
    adata.write_h5ad(CKPT_H5AD)
    print("checkpoint saved.")

# step 6: normalize, log1p, HVGs
print("\nstep 6: normalize -> log1p -> HVGs")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
print("normalized to 10k counts/cell, log1p-transformed.")

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=3000,
    batch_key="sample",
    flavor="seurat_v3",
    subset=False,
)
n_hvg = adata.var["highly_variable"].sum()
print(f"HVGs: {n_hvg:,}")

# step 7: scale, PCA, neighbors, UMAP, leiden
print("\nstep 7: scale -> PCA -> neighbors -> UMAP -> leiden")

adata.raw = adata  # stash log-normalized counts before HVG subset

print(f"subsetting to {adata.var.highly_variable.sum():,} HVGs ...")
adata = adata[:, adata.var.highly_variable].copy()
print(f"shape: {adata.shape[0]:,} x {adata.shape[1]:,}")

# densify as float32 (~4.8 GB) to avoid float64 peak (~9.6 GB) from sc.pp.scale
print("densifying as float32 ...")
if sp.issparse(adata.X):
    adata.X = adata.X.toarray().astype(np.float32)
else:
    adata.X = adata.X.astype(np.float32)
print(f"matrix in memory: {adata.X.nbytes / 1e9:.1f} GB")

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50, svd_solver="randomized", random_state=42)
sc.pl.pca_variance_ratio(adata, n_pcs=50, save="_pca_variance.png", show=False)

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, flavor="igraph", n_iterations=2)
n_clusters = adata.obs["leiden"].nunique()
print(f"leiden: {n_clusters} clusters at resolution 0.5")

# step 8: save
print("\nstep 8: saving processed AnnData")
adata.write_h5ad(OUT_H5AD)
print(f"saved: {OUT_H5AD}")

# step 9: clinical metadata
print("\nstep 9: clinical metadata")

s1 = pd.read_excel(TABLE_S1, sheet_name=0)
print(f"Table S1 shape: {s1.shape}")

s1 = s1.rename(columns={
    "Subject No":                                                                    "patient",
    "Best overall response":                                                         "best_overall_response",
    "Response caegory (0 = Slow Progressor, 1 = fast progressor)":                  "_prog_code",
    "PFS (days)":                                                                    "PFS_days",
    "OS (days)":                                                                     "OS_days",
    "HER2 ":                                                                         "HER2_status",
    "EBV":                                                                           "EBV_status",
    "MSI/MSS":                                                                       "MSI_status",
    "TCGA Subtype (0 = Indeterminate, 1 = GS, 2 = CIN, 3 = MSI-High, 4 = EBV":    "TCGA_subtype_code",
    "PDL1 22C3 score (Baseline)":                                                   "PDL1_baseline_CPS",
})

s1["progression_category"] = s1["_prog_code"].map({0: "Slow", 1: "Fast"})
tcga_map = {0: "Indeterminate", 1: "GS", 2: "CIN", 3: "MSI-High", 4: "EBV"}
s1["TCGA_subtype"] = s1["TCGA_subtype_code"].map(tcga_map)

META_COLS_S1 = [
    "patient", "best_overall_response", "progression_category",
    "PFS_days", "OS_days", "HER2_status", "EBV_status",
    "MSI_status", "TCGA_subtype", "PDL1_baseline_CPS",
]
meta_s1 = s1[META_COLS_S1].drop_duplicates("patient")
print(f"patients in Table S1: {meta_s1['patient'].nunique()}")

s2 = pd.read_excel(TABLE_S2, sheet_name=0)
print(f"Table S2 shape: {s2.shape}")
s2 = s2.rename(columns={"Patient": "patient", "Progression": "progression_s2"})
meta_s2 = s2[["patient", "progression_s2"]].drop_duplicates("patient")

obs = adata.obs.copy()
for col in obs.select_dtypes("category").columns:
    obs[col] = obs[col].astype(str).replace("nan", np.nan)

tp_map = {"B": "Baseline", "F1": "FU1", "F2": "FU2"}
obs["timepoint_label"] = obs["timepoint"].map(tp_map).fillna(obs["timepoint"])
print(f"timepoints: {obs['timepoint'].value_counts().to_dict()}")

obs = obs.merge(meta_s1, on="patient", how="left")
obs = obs.merge(meta_s2, on="patient", how="left")
obs.index = adata.obs.index
adata.obs = obs
print("merge done.")

# step 10: summary
print("\nstep 10: summary")
print(adata)
print(f"\nshape: {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")
print(f"leiden clusters: {adata.obs['leiden'].nunique()}")
print(f"patients: {adata.obs['patient'].nunique()}")
print(f"\nall .obs columns: {adata.obs.columns.tolist()}")

META_CHECK = [
    "progression_category", "best_overall_response",
    "MSI_status", "EBV_status", "HER2_status",
    "TCGA_subtype", "PDL1_baseline_CPS",
    "PFS_days", "OS_days", "timepoint_label",
]
print("\nmetadata fill rates:")
for col in META_CHECK:
    if col in adata.obs.columns:
        n_filled = adata.obs[col].notna().sum()
        pct = 100 * n_filled / adata.n_obs
        print(f"  {col:<28}  {n_filled:>7,} / {adata.n_obs:,} ({pct:.1f}%)")
    else:
        print(f"  {col:<28}  *** MISSING ***")

adata.write_h5ad(OUT_H5AD)
print(f"\nfinal object saved to: {OUT_H5AD}")
print("done.")
