"""
05_integration.py
Multi-dataset integration using scVI (scverse/scvi-tools).

Inputs : data/processed/per_dataset/<dataset_id>_processed.h5ad (one per dataset)
Output : data/processed/integrated/gastric_meta_integrated.h5ad
         data/processed/integrated/scvi_model/          (saved scVI model)
         outputs/integration/                           (UMAP plots, QC)

Strategy:
  1. Load each processed per-dataset h5ad.
  2. Find the shared gene universe (intersection of highly variable genes
     across all datasets, then top N by overall variance).
  3. Concatenate into a single AnnData, tagging batch = dataset_id.
  4. Train scVI (VAE) on raw counts (use adata.raw or counts layer).
  5. Extract the scVI latent embedding (obsm["X_scVI"]).
  6. Compute neighbors in latent space, UMAP, Leiden clustering.
  7. Save the integrated object.

Usage:
    python 05_integration.py                 # run full integration
    python 05_integration.py --skip-training # load saved model and re-embed
    python 05_integration.py --n-hvgs 4000   # adjust shared HVG count
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", line_buffering=True)
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", line_buffering=True)

import gc
import os
import argparse
import warnings
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import scanpy as sc
import scvi
import h5py

warnings.filterwarnings("ignore")

BASE      = os.path.expanduser("~/gastric_tme_project")
PER_DS    = os.path.join(BASE, "data/processed/per_dataset")
INT_DIR   = os.path.join(BASE, "data/processed/integrated")
MODEL_DIR = os.path.join(INT_DIR, "scvi_model")
OUT_DIR   = os.path.join(BASE, "outputs/integration")
MANIFEST  = os.path.join(BASE, "data/external/dataset_manifest.csv")

os.makedirs(INT_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)
sc.settings.figdir = OUT_DIR


MIN_GENES_FOR_INTEGRATION = 5000  # datasets with fewer genes are excluded from main integration


def _read_h5ad_var_names(path):
    """Read var_names and raw var_names from an h5ad file without loading X into RAM."""
    with h5py.File(path, 'r') as f:
        # var names: stored at var/_index in newer anndata; some old files use different key
        var_names = None
        if 'var' in f:
            vg = f['var']
            if isinstance(vg, h5py.Group):
                if '_index' in vg:
                    var_names = vg['_index'].asstr()[:]
                else:
                    # fall back to first key
                    k = list(vg.keys())[0]
                    var_names = vg[k].asstr()[:] if vg[k].dtype.kind in ('S', 'O', 'U') else list(vg.keys())
            elif isinstance(vg, h5py.Dataset):
                var_names = vg.asstr()[:]
        var_names = list(var_names) if var_names is not None else []

        raw_var_names = None
        if 'raw' in f and 'var' in f['raw']:
            rvg = f['raw']['var']
            if isinstance(rvg, h5py.Group) and '_index' in rvg:
                raw_var_names = list(rvg['_index'].asstr()[:])

        hvg_flags = None
        dispersions_norm = None
        if 'var' in f and isinstance(f['var'], h5py.Group):
            vg = f['var']
            if 'highly_variable' in vg:
                hvg_flags = vg['highly_variable'][:].astype(bool).tolist()
            if 'dispersions_norm' in vg:
                dispersions_norm = vg['dispersions_norm'][:].tolist()

    return var_names, raw_var_names, hvg_flags, dispersions_norm


def load_gene_metadata(manifest):
    """
    Read only var_names and HVG metadata from each processed h5ad — no X loading.
    Returns a dict: dataset_id → {path, var_names, raw_var_names, hvg_flags, dispersions_norm}.
    Datasets excluded due to too few genes are silently dropped.
    """
    print("Reading gene metadata (no X loading) ...")
    meta = {}
    excluded = []
    for _, row in manifest.iterrows():
        did  = row["dataset_id"]
        path = os.path.join(PER_DS, f"{did}_processed.h5ad")
        if not os.path.exists(path):
            print(f"  SKIP {did}: processed h5ad not found")
            continue
        size_gb = os.path.getsize(path) / 1e9
        try:
            var_names, raw_var_names, hvg_flags, dispersions_norm = _read_h5ad_var_names(path)
            n_effective = len(raw_var_names) if raw_var_names and len(raw_var_names) > len(var_names) else len(var_names)
            if n_effective < MIN_GENES_FOR_INTEGRATION:
                print(f"  SKIP {did} ({size_gb:.1f} GB): only {n_effective:,} genes (< {MIN_GENES_FOR_INTEGRATION:,})")
                excluded.append(did)
                continue
            meta[did] = dict(
                path=path,
                size_gb=size_gb,
                var_names=var_names,
                raw_var_names=raw_var_names,
                hvg_flags=hvg_flags,
                dispersions_norm=dispersions_norm,
            )
            n_total = len(raw_var_names) if raw_var_names and len(raw_var_names) > len(var_names) else len(var_names)
            n_hvg   = sum(hvg_flags) if hvg_flags else len(var_names)
            print(f"  {did} ({size_gb:.1f} GB): {n_hvg:,} HVGs | {n_total:,} total genes")
        except Exception as e:
            print(f"  FAILED to read {did}: {e}")

    if excluded:
        print(f"\n  NOTE: excluded from integration (too few genes): {excluded}")
    return meta


def find_shared_hvgs(meta, n_hvgs=4000):
    """
    Build shared HVG set from gene metadata (no loaded AnnData required).
    Union of per-dataset HVGs, intersected with genes present in ALL datasets,
    then ranked by mean dispersions_norm.
    """
    print(f"\nSelecting {n_hvgs} shared HVGs across {len(meta)} datasets ...")
    hvg_sets = []
    all_gene_sets = []
    for did, m in meta.items():
        effective = m['raw_var_names'] if (m['raw_var_names'] and len(m['raw_var_names']) > len(m['var_names'])) else m['var_names']
        all_gene_set = set(effective)
        all_gene_sets.append(all_gene_set)
        ds_hvgs = set(g for g, h in zip(m['var_names'], m['hvg_flags']) if h) if m['hvg_flags'] else all_gene_set
        hvg_sets.append(ds_hvgs)
        print(f"  {did}: {len(ds_hvgs):,} HVGs  ({len(all_gene_set):,} total genes)")

    union_genes  = set.union(*hvg_sets)
    common_genes = set.intersection(*all_gene_sets)
    candidates   = list(union_genes & common_genes)
    print(f"  union HVGs: {len(union_genes):,}")
    print(f"  genes in all datasets: {len(common_genes):,}")
    print(f"  candidates (union ∩ all-datasets): {len(candidates):,}")

    if len(candidates) <= n_hvgs:
        print(f"  using all {len(candidates):,} candidates (< requested {n_hvgs})")
        return candidates

    disp_dict = {}
    for did, m in meta.items():
        if m['dispersions_norm']:
            s = pd.Series(dict(zip(m['var_names'], m['dispersions_norm']))).reindex(candidates).fillna(0)
            disp_dict[did] = s.values
    if disp_dict:
        mean_disp = np.mean(list(disp_dict.values()), axis=0)
        ranked    = [candidates[i] for i in np.argsort(mean_disp)[::-1]]
        shared    = ranked[:n_hvgs]
    else:
        shared = candidates[:n_hvgs]

    print(f"  final shared HVG set: {len(shared):,}")
    return shared


def _lognorm_to_approx_counts(X):
    if sp.issparse(X):
        out = X.copy().astype(np.float64)
        out.data = np.round(np.expm1(out.data)).astype(np.float32)
        return out.astype(np.float32)
    return np.round(np.expm1(X.astype(np.float64))).astype(np.float32)


def _load_and_prepare_one(did, m, shared_genes):
    """
    Load a single h5ad, subset to shared_genes, ensure raw counts in .X.
    Deletes the full AnnData from RAM before returning the prepared subset.
    """
    print(f"  loading {did} ({m['size_gb']:.1f} GB) ...")
    a = sc.read_h5ad(m['path'])
    a.obs["dataset_id"] = did

    # Subset to shared genes — prefer .raw if it has more genes than .var_names
    genes_in_var = [g for g in shared_genes if g in a.var_names]
    if len(genes_in_var) == len(shared_genes):
        a_sub = a[:, shared_genes].copy()
    elif a.raw is not None:
        raw_full = a.raw.to_adata()
        available = [g for g in shared_genes if g in raw_full.var_names]
        if len(available) < len(shared_genes):
            print(f"    WARNING: {did}: only {len(available)}/{len(shared_genes)} shared genes found")
        a_sub = raw_full[:, available].copy()
        a_sub.obs = a.obs.copy()
        del raw_full
    else:
        a_sub = a[:, genes_in_var].copy()
        if len(genes_in_var) < len(shared_genes):
            print(f"    WARNING: {did}: only {len(genes_in_var)}/{len(shared_genes)} shared genes in .var_names")

    del a  # free the full h5ad from RAM immediately
    gc.collect()

    # Restore raw counts for scVI (prefers layers["counts"])
    if "counts" in a_sub.layers:
        a_sub.X = a_sub.layers["counts"].astype(np.float32)
        print(f"    {did}: raw counts from layers['counts']")
    elif a_sub.raw is not None:
        raw_full2 = a_sub.raw.to_adata()
        avail2    = [g for g in shared_genes if g in raw_full2.var_names]
        approx    = _lognorm_to_approx_counts(raw_full2[:, avail2].X)
        a_sub.X   = approx
        a_sub.layers["counts"] = approx.copy()
        del raw_full2
        print(f"    {did}: back-transformed log1p→counts from .raw")
    else:
        approx = _lognorm_to_approx_counts(a_sub.X)
        a_sub.X = approx
        a_sub.layers["counts"] = approx.copy()
        print(f"    {did}: back-transformed log1p→counts from .X")

    a_sub.obs["dataset_id"] = did
    gc.collect()
    return a_sub


def prepare_concat(meta, shared_genes):
    """
    Load each dataset one at a time, prepare it (subset + counts), return dict of AnnDatas.
    Peak RAM ≈ largest single dataset (Korea ~10 GB) + accumulated prepared subsets.
    """
    prepared = {}
    for did, m in meta.items():
        a_sub = _load_and_prepare_one(did, m, shared_genes)
        prepared[did] = a_sub
    return prepared


def concatenate_datasets(prepared):
    """Concatenate all prepared AnnDatas into one."""
    dids   = list(prepared.keys())
    adatas = list(prepared.values())

    print(f"\nConcatenating {len(adatas)} datasets ...")
    combined = sc.concat(
        adatas,
        label="dataset_id",
        keys=dids,
        merge="same",
        uns_merge="same",
    )
    # obs_names may collide across datasets — make unique
    combined.obs_names_make_unique()
    combined.obs["batch"] = combined.obs["dataset_id"].astype(str)

    print(f"Combined: {combined.n_obs:,} cells x {combined.n_vars:,} genes")
    print(f"Datasets: {combined.obs['dataset_id'].value_counts().to_dict()}")
    return combined


def _prevent_sleep():
    """Prevent Windows from sleeping during training."""
    try:
        import ctypes
        ES_CONTINUOUS        = 0x80000000
        ES_SYSTEM_REQUIRED   = 0x00000001
        ES_AWAYMODE_REQUIRED = 0x00000040
        ctypes.windll.kernel32.SetThreadExecutionState(
            ES_CONTINUOUS | ES_SYSTEM_REQUIRED | ES_AWAYMODE_REQUIRED
        )
        print("  sleep prevention: active (Windows SetThreadExecutionState)")
    except Exception as e:
        print(f"  sleep prevention: unavailable ({e})")


def _allow_sleep():
    try:
        import ctypes
        ctypes.windll.kernel32.SetThreadExecutionState(0x80000000)  # ES_CONTINUOUS only
    except Exception:
        pass


class _PeriodicSave:
    """Lightning-compatible callback that saves scVI model weights every N epochs."""
    def __init__(self, model_ref, save_dir, every_n=10):
        self._model  = model_ref
        self._dir    = save_dir
        self._every  = every_n

    def __getattr__(self, name):
        return lambda *args, **kwargs: None

    def on_train_epoch_end(self, trainer, pl_module):
        import torch, pickle
        epoch = trainer.current_epoch + 1
        if epoch % self._every == 0:
            ckpt_dir = self._dir + f"_ckpt_ep{epoch:03d}"
            os.makedirs(ckpt_dir, exist_ok=True)
            # Save weights directly — model.save() may require _is_trained=True
            torch.save(
                {"model_state_dict": self._model.module.state_dict(), "epoch": epoch},
                os.path.join(ckpt_dir, "model.pt"),
            )
            # Save var_names so checkpoint is self-contained
            var_names = self._model.adata.var_names.tolist()
            with open(os.path.join(ckpt_dir, "var_names.txt"), "w") as f:
                f.write("\n".join(var_names))
            print(f"\n  [checkpoint] epoch {epoch} weights → {ckpt_dir}")


def train_scvi(combined, n_latent=30, n_layers=2, n_epochs=None, model_dir=MODEL_DIR):
    """
    Set up and train scVI model on the concatenated object.
    scVI expects raw (non-normalized) counts in adata.X.
    """
    print(f"\nSetting up scVI model ...")
    print(f"  latent dims: {n_latent}   layers: {n_layers}")

    _prevent_sleep()

    scvi.model.SCVI.setup_anndata(
        combined,
        layer=None,          # use .X
        batch_key="batch",
    )

    model = scvi.model.SCVI(
        combined,
        n_latent=n_latent,
        n_layers=n_layers,
        gene_likelihood="nb",   # negative binomial — appropriate for scRNA-seq counts
        dispersion="gene-batch",
    )

    if n_epochs is None:
        n_epochs = max(400, int(np.ceil(1e5 / combined.n_obs) * 100))
        n_epochs = min(n_epochs, 400)
    print(f"  training epochs: {n_epochs}")

    ckpt_cb = _PeriodicSave(model, model_dir, every_n=10)

    model.train(
        max_epochs=n_epochs,
        early_stopping=True,
        early_stopping_patience=30,
        plan_kwargs={"lr": 1e-3},
        callbacks=[ckpt_cb],
    )

    _allow_sleep()

    print(f"  saving model to {model_dir} ...")
    os.makedirs(model_dir, exist_ok=True)
    model.save(model_dir, overwrite=True)
    print(f"  model saved.")
    return model


def embed_and_cluster(combined, model, n_pcs=50, n_neighbors=20, leiden_res=0.5):
    """Extract scVI latent embedding, compute UMAP and Leiden on it."""
    print(f"\nExtracting scVI latent representation ...")
    combined.obsm["X_scVI"] = model.get_latent_representation()
    print(f"  X_scVI shape: {combined.obsm['X_scVI'].shape}")

    print(f"  computing kNN graph in scVI latent space ...")
    sc.pp.neighbors(combined, use_rep="X_scVI", n_neighbors=n_neighbors)
    sc.tl.umap(combined)
    sc.tl.leiden(combined, resolution=leiden_res, flavor="igraph", n_iterations=2)
    n_cl = combined.obs["leiden"].nunique()
    print(f"  Leiden clusters (res={leiden_res}): {n_cl}")
    return combined


def save_plots(combined):
    print(f"\nSaving integration plots to {OUT_DIR} ...")

    color_keys = ["dataset_id", "leiden"]
    # add more coloring if columns are present
    for col in ["progression_category", "best_overall_response", "timepoint_label", "MSI_status"]:
        if col in combined.obs.columns:
            color_keys.append(col)

    sc.pl.umap(combined, color=color_keys, ncols=2,
               save="_integration_overview.png", show=False)

    # per-dataset fraction plot
    fig_data = (combined.obs.groupby(["leiden", "dataset_id"])
                .size().unstack(fill_value=0))
    fig_data_pct = fig_data.div(fig_data.sum(axis=1), axis=0)
    ax = fig_data_pct.plot(kind="bar", stacked=True, figsize=(14, 5))
    ax.set_title("Dataset composition per Leiden cluster")
    ax.set_xlabel("Leiden cluster")
    ax.set_ylabel("Fraction of cells")
    ax.legend(bbox_to_anchor=(1.01, 1), loc="upper left", fontsize=8)
    import matplotlib.pyplot as plt
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "dataset_composition_per_cluster.png"), dpi=150)
    plt.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-training", action="store_true",
                        help="load saved scVI model instead of retraining")
    parser.add_argument("--n-hvgs",   type=int, default=4000, help="shared HVG count")
    parser.add_argument("--n-latent", type=int, default=30,   help="scVI latent dims")
    parser.add_argument("--leiden-res", type=float, default=0.5, help="Leiden resolution")
    parser.add_argument("--force",    action="store_true",
                        help="redo integration even if output exists")
    parser.add_argument("--exclude", nargs="+", default=[],
                        help="dataset_ids to exclude from integration")
    parser.add_argument("--use-checkpoint", default=None, metavar="CKPT_DIR",
                        help="load weights from a partial checkpoint dir and skip training")
    args = parser.parse_args()

    OUT_H5AD = os.path.join(INT_DIR, "gastric_meta_integrated.h5ad")

    if os.path.exists(OUT_H5AD) and not args.force and not args.skip_training:
        print(f"Integrated object already exists: {OUT_H5AD}")
        print("Use --force to redo, or --skip-training to re-embed with saved model.")
        return

    if not os.path.exists(MANIFEST):
        print("Manifest not found. Run 02_dataset_manifest.py first.")
        sys.exit(1)

    manifest = pd.read_csv(MANIFEST)
    if args.exclude:
        manifest = manifest[~manifest["dataset_id"].isin(args.exclude)]
        print(f"Excluding from integration: {args.exclude}")

    # ── pass 1: read gene metadata only (no X loading) ────────────────────────
    meta = load_gene_metadata(manifest)

    if len(meta) < 2:
        print(f"\nOnly {len(meta)} dataset(s) available for integration.")
        print("At least 2 datasets are needed. Run 04_standardized_qc.py for more datasets.")
        return

    # ── shared HVGs (from metadata, no AnnData needed) ────────────────────────
    shared_genes = find_shared_hvgs(meta, n_hvgs=args.n_hvgs)

    # ── pass 2: load each dataset one at a time, prepare, then concatenate ────
    print("\nPreparing datasets one at a time (subset to shared HVGs + restore counts) ...")
    prepared = prepare_concat(meta, shared_genes)
    combined = concatenate_datasets(prepared)

    # ── train or load scVI model ───────────────────────────────────────────────
    if args.use_checkpoint and os.path.exists(args.use_checkpoint):
        import torch
        ckpt_path = os.path.join(args.use_checkpoint, "model.pt")
        print(f"\nLoading weights from checkpoint: {ckpt_path}")
        scvi.model.SCVI.setup_anndata(combined, batch_key="batch")
        model = scvi.model.SCVI(
            combined, n_latent=args.n_latent, n_layers=2,
            gene_likelihood="nb", dispersion="gene-batch",
        )
        state = torch.load(ckpt_path, map_location="cpu", weights_only=False)
        sd = state.get("model_state_dict", state)
        # strict=False: checkpoint may contain pyro_param_store or other non-VAE keys
        model.module.load_state_dict(sd, strict=False)
        model.is_trained_ = True
        epoch = state.get("epoch", "?")
        print(f"  loaded weights from epoch {epoch}")
    elif args.skip_training and os.path.exists(MODEL_DIR):
        print(f"\nLoading saved scVI model from {MODEL_DIR} ...")
        scvi.model.SCVI.setup_anndata(combined, batch_key="batch")
        model = scvi.model.SCVI.load(MODEL_DIR, adata=combined)
    else:
        model = train_scvi(combined, n_latent=args.n_latent)

    # ── embed + cluster ───────────────────────────────────────────────────────
    combined = embed_and_cluster(combined, model, leiden_res=args.leiden_res)

    # ── plots ─────────────────────────────────────────────────────────────────
    save_plots(combined)

    # ── save ──────────────────────────────────────────────────────────────────
    print(f"\nSaving integrated AnnData to {OUT_H5AD} ...")
    combined.write_h5ad(OUT_H5AD)
    size_gb = os.path.getsize(OUT_H5AD) / 1e9
    print(f"Saved ({size_gb:.1f} GB)")

    print(f"\n{'='*60}")
    print("Integration complete.")
    print(f"  Cells: {combined.n_obs:,}")
    print(f"  Genes: {combined.n_vars:,}")
    print(f"  Leiden clusters: {combined.obs['leiden'].nunique()}")
    print(f"  Datasets: {combined.obs['dataset_id'].nunique()}")
    print(f"\nNext step: run 06_cell_type_annotation.py")


if __name__ == "__main__":
    main()
