"""
04_standardized_qc.py
Standardized QC and preprocessing pipeline — runs on any dataset in the manifest.
Produces a consistent processed h5ad in data/processed/per_dataset/<dataset_id>_processed.h5ad.

Usage:
    python 04_standardized_qc.py                       # process all downloaded datasets
    python 04_standardized_qc.py korea_kim2022         # process one dataset
    python 04_standardized_qc.py --list                # show processing status

QC parameters (same thresholds as 01_preprocessing.py for the Korean cohort):
    min_genes  : 200     cells with fewer genes are removed (empty droplets)
    max_genes  : 6000    cells with more genes are removed (likely doublets)
    max_mt_pct : 20      cells with >20% MT reads are removed (dying cells)
    min_cells  : 3       genes expressed in fewer cells are removed
    n_hvgs     : 3000    highly variable genes for PCA
    n_pcs      : 50      PCA components
    n_neighbors: 15      kNN graph
    leiden_res : 0.5     Leiden clustering resolution
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
import scrublet as scr

warnings.filterwarnings("ignore")

BASE      = os.path.expanduser("~/gastric_tme_project")
EXT_DIR   = os.path.join(BASE, "data/external")
PER_DS    = os.path.join(BASE, "data/processed/per_dataset")
QC_BASE   = os.path.join(BASE, "outputs/qc_plots")
MANIFEST  = os.path.join(EXT_DIR, "dataset_manifest.csv")

os.makedirs(PER_DS, exist_ok=True)

# ─── QC parameters ────────────────────────────────────────────────────────────
QC_PARAMS = dict(
    min_genes   = 200,
    max_genes   = 6000,
    max_mt_pct  = 20.0,
    min_cells   = 3,
    n_hvgs      = 3000,
    n_pcs       = 50,
    n_neighbors = 15,
    leiden_res  = 0.5,
)


def get_raw_h5ad(row):
    """Return path to raw h5ad for this dataset (local or downloaded)."""
    if pd.notna(row.get("local_h5ad")) and row["local_h5ad"]:
        return row["local_h5ad"]
    return os.path.join(EXT_DIR, row["dataset_id"], f"{row['dataset_id']}_raw.h5ad")


def get_processed_h5ad(dataset_id):
    return os.path.join(PER_DS, f"{dataset_id}_processed.h5ad")


def get_checkpoint(dataset_id):
    return os.path.join(PER_DS, f"{dataset_id}_ckpt_scrublet.h5ad")


def infer_sample_key(adata):
    """Pick the best obs column to use as sample/batch key."""
    candidates = ["sample", "Sample", "library_id", "batch", "donor", "patient",
                  "orig.ident", "orig_ident", "sample_id", "SampleID", "sample_file"]
    for c in candidates:
        if c in adata.obs.columns and adata.obs[c].nunique() > 1:
            return c
    # fallback: use geo_accession or dataset_id if present
    for c in ["geo_accession", "dataset_id"]:
        if c in adata.obs.columns:
            return c
    return None


def run_scrublet(adata, sample_key):
    """Per-sample doublet detection using Scrublet."""
    doublet_scores = np.zeros(adata.n_obs, dtype=float)
    predicted_dbl  = np.zeros(adata.n_obs, dtype=bool)

    if sample_key and sample_key in adata.obs.columns:
        samples = adata.obs[sample_key].unique()
        print(f"  running Scrublet on {len(samples)} samples (key='{sample_key}') ...")
    else:
        samples = ["__all__"]
        print(f"  running Scrublet on full dataset (no sample key found) ...")

    for samp in samples:
        if samp == "__all__":
            mask = np.ones(adata.n_obs, dtype=bool)
        else:
            mask = (adata.obs[sample_key] == samp).values

        n_cells = mask.sum()
        if n_cells < 50:
            print(f"    [{samp}] {n_cells} cells — skipping (< 50)")
            continue

        X_samp = adata.X[mask]
        if sp.issparse(X_samp):
            X_samp = X_samp.toarray()

        try:
            scrub = scr.Scrublet(X_samp, random_state=42)
            scores, pred = scrub.scrub_doublets(verbose=False)
            doublet_scores[mask] = scores
            predicted_dbl[mask]  = pred
            print(f"    [{samp}] {n_cells:,} cells, {pred.sum()} doublets ({100*pred.mean():.1f}%)")
        except Exception as e:
            print(f"    [{samp}] Scrublet failed: {e}")

    adata.obs["doublet_score"]     = doublet_scores
    adata.obs["predicted_doublet"] = predicted_dbl
    total_dbl = predicted_dbl.sum()
    print(f"  total doublets: {total_dbl:,} / {adata.n_obs:,} ({100*predicted_dbl.mean():.2f}%)")
    return adata


def process_dataset(row, params=QC_PARAMS, force=False):
    dataset_id  = row["dataset_id"]
    raw_h5ad    = get_raw_h5ad(row)
    out_h5ad    = get_processed_h5ad(dataset_id)
    ckpt_h5ad   = get_checkpoint(dataset_id)
    qc_dir      = os.path.join(QC_BASE, dataset_id)
    os.makedirs(qc_dir, exist_ok=True)
    sc.settings.figdir = qc_dir

    print(f"\n{'='*70}")
    print(f"Processing: {dataset_id}")
    print(f"{'='*70}")

    if not os.path.exists(raw_h5ad):
        print(f"  SKIP: raw h5ad not found at {raw_h5ad}")
        print(f"        Run 03_dataset_download.py first for this dataset.")
        return False

    if os.path.exists(out_h5ad) and not force:
        size_gb = os.path.getsize(out_h5ad) / 1e9
        print(f"  Already processed ({size_gb:.1f} GB). Use --force to rerun.")
        return True

    # ── load ────────────────────────────────────────────────────────────────
    if os.path.exists(ckpt_h5ad) and not force:
        print(f"  Loading scrublet checkpoint ...")
        adata = sc.read_h5ad(ckpt_h5ad)
        print(f"  loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    else:
        print(f"  Loading raw data ...")
        adata = sc.read_h5ad(raw_h5ad)
        print(f"  shape: {adata.n_obs:,} cells x {adata.n_vars:,} genes")

        # tag source
        adata.obs["dataset_id"] = dataset_id
        if "geo_accession" not in adata.obs.columns:
            acc = row.get("geo_accession", "")
            adata.obs["geo_accession"] = acc if pd.notna(acc) else ""

        # ensure gene names are strings
        adata.var_names = adata.var_names.astype(str)
        adata.var_names_make_unique()

        # ── QC metrics ──────────────────────────────────────────────────────
        print(f"\n  QC metrics ...")
        adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
        sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None,
                                   log1p=False, inplace=True)
        print(f"    genes/cell  median={adata.obs['n_genes_by_counts'].median():.0f}  "
              f"mean={adata.obs['n_genes_by_counts'].mean():.0f}")
        print(f"    UMI/cell    median={adata.obs['total_counts'].median():.0f}  "
              f"mean={adata.obs['total_counts'].mean():.0f}")
        print(f"    MT%         median={adata.obs['pct_counts_mt'].median():.2f}  "
              f"mean={adata.obs['pct_counts_mt'].mean():.2f}")

        # ── QC plots ────────────────────────────────────────────────────────
        print(f"  Saving QC plots to {qc_dir} ...")
        sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
                     jitter=0.4, save=f"_{dataset_id}_violin.png", show=False)
        sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt",
                      save=f"_{dataset_id}_mt_scatter.png", show=False)
        sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts",
                      save=f"_{dataset_id}_gene_scatter.png", show=False)

        # ── cell / gene filtering ────────────────────────────────────────────
        print(f"\n  Filtering ...")
        n0 = adata.n_obs
        sc.pp.filter_cells(adata, min_genes=params["min_genes"])
        sc.pp.filter_cells(adata, max_genes=params["max_genes"])
        adata = adata[adata.obs["pct_counts_mt"] < params["max_mt_pct"]].copy()
        n1 = adata.n_obs
        print(f"    cells: {n0:,} → {n1:,} (removed {n0-n1:,}, {100*(n0-n1)/n0:.1f}%)")

        g0 = adata.n_vars
        sc.pp.filter_genes(adata, min_cells=params["min_cells"])
        g1 = adata.n_vars
        print(f"    genes: {g0:,} → {g1:,} (removed {g0-g1:,})")

        if adata.n_obs == 0 or adata.n_vars == 0:
            print(f"  SKIP: 0 cells or 0 genes after filtering — dataset not suited for QC thresholds")
            return False

        # ── doublet detection ────────────────────────────────────────────────
        print(f"\n  Doublet detection (Scrublet) ...")
        sample_key = infer_sample_key(adata)
        adata = run_scrublet(adata, sample_key)

        print(f"  Saving scrublet checkpoint ...")
        adata.write_h5ad(ckpt_h5ad)

    # ── normalization + HVGs ──────────────────────────────────────────────────
    print(f"\n  Normalize → log1p → HVGs ...")
    # stash raw counts before normalization so 05_integration.py can pass them to scVI
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sample_key = infer_sample_key(adata)
    # seurat_v3 with batch_key can fail on tiny batches (LOESS svd error);
    # fall back to seurat flavor (log-normalized) if that happens.
    try:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=params["n_hvgs"],
            batch_key=sample_key,
            flavor="seurat_v3",
            subset=False,
        )
    except Exception as e:
        print(f"  seurat_v3 HVG failed ({e}); retrying with seurat flavor ...")
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=params["n_hvgs"],
            batch_key=sample_key,
            flavor="seurat",
            subset=False,
        )
    n_hvg = adata.var["highly_variable"].sum()
    print(f"  HVGs: {n_hvg:,}")

    # ── dimensionality reduction ──────────────────────────────────────────────
    print(f"\n  Scale → PCA → kNN → UMAP → Leiden ...")
    adata.raw = adata  # stash log-normalized before HVG subset
    adata_hvg = adata[:, adata.var.highly_variable].copy()
    print(f"  HVG subset: {adata_hvg.n_obs:,} x {adata_hvg.n_vars:,}")

    if sp.issparse(adata_hvg.X):
        adata_hvg.X = adata_hvg.X.toarray().astype(np.float32)
    else:
        adata_hvg.X = adata_hvg.X.astype(np.float32)

    sc.pp.scale(adata_hvg, max_value=10)
    sc.tl.pca(adata_hvg, n_comps=params["n_pcs"], svd_solver="randomized", random_state=42)
    sc.pl.pca_variance_ratio(adata_hvg, n_pcs=params["n_pcs"],
                             save=f"_{dataset_id}_pca_variance.png", show=False)
    sc.pp.neighbors(adata_hvg, n_neighbors=params["n_neighbors"],
                    n_pcs=params["n_pcs"])
    sc.tl.umap(adata_hvg)
    sc.tl.leiden(adata_hvg, resolution=params["leiden_res"],
                 flavor="igraph", n_iterations=2)

    n_clusters = adata_hvg.obs["leiden"].nunique()
    print(f"  Leiden clusters: {n_clusters} @ resolution {params['leiden_res']}")

    # copy embeddings back to full adata
    adata.obsm["X_pca"]  = adata_hvg.obsm["X_pca"]
    adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"]
    adata.obs["leiden"]  = adata_hvg.obs["leiden"]
    adata.obsp           = adata_hvg.obsp
    adata.uns            = {**adata_hvg.uns, **{k: v for k, v in adata.uns.items()
                                                 if k not in adata_hvg.uns}}

    # ── save ──────────────────────────────────────────────────────────────────
    print(f"\n  Saving processed h5ad ...")
    adata.write_h5ad(out_h5ad)
    size_gb = os.path.getsize(out_h5ad) / 1e9
    print(f"  Saved: {out_h5ad} ({size_gb:.1f} GB)")

    # ── summary ───────────────────────────────────────────────────────────────
    print(f"\n  Summary for {dataset_id}:")
    print(f"    cells: {adata.n_obs:,}  genes: {adata.n_vars:,}  clusters: {n_clusters}")
    if "predicted_doublet" in adata.obs.columns:
        dbl_n = adata.obs["predicted_doublet"].sum()
        print(f"    doublets flagged (not removed): {dbl_n:,} ({100*dbl_n/adata.n_obs:.2f}%)")

    return True


# ─── main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("dataset_id", nargs="?", help="specific dataset_id to process")
    parser.add_argument("--list",  action="store_true", help="show processing status")
    parser.add_argument("--force", action="store_true", help="reprocess even if output exists")
    args = parser.parse_args()

    if not os.path.exists(MANIFEST):
        print("Manifest not found. Run 02_dataset_manifest.py first.")
        sys.exit(1)

    manifest = pd.read_csv(MANIFEST)

    if args.list:
        print("\n--- Processing status ---")
        for _, row in manifest.iterrows():
            did = row["dataset_id"]
            raw = get_raw_h5ad(row)
            out = get_processed_h5ad(did)
            raw_ok = os.path.exists(raw)
            out_ok = os.path.exists(out)
            if out_ok:
                size_gb = os.path.getsize(out) / 1e9
                print(f"  {did:<25}  processed ({size_gb:.1f} GB)")
            elif raw_ok:
                print(f"  {did:<25}  raw available, not yet processed")
            else:
                print(f"  {did:<25}  raw NOT found — download first")
        return

    if args.dataset_id:
        rows = manifest[manifest["dataset_id"] == args.dataset_id]
        if rows.empty:
            print(f"dataset_id '{args.dataset_id}' not in manifest.")
            sys.exit(1)
    else:
        rows = manifest

    results = {}
    for _, row in rows.iterrows():
        ok = process_dataset(row, force=args.force)
        results[row["dataset_id"]] = "OK" if ok else "SKIPPED/FAILED"

    print("\n\n=== Processing summary ===")
    for did, status in results.items():
        print(f"  {did:<25}  {status}")


if __name__ == "__main__":
    main()
