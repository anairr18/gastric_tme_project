"""
11_cell_communication.py
Ligand-receptor cell-cell communication analysis using LIANA.
Focuses on Myeloid → T/NK signalling differences between slow and fast progressors.

Outputs in outputs/cell_communication/:
  liana_dotplot_all.png
  liana_dotplot_myeloid_tcell.png
  differential_interactions.csv
  interaction_heatmap_slow.png
  interaction_heatmap_fast.png
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
import liana
from liana.method import singlecellsignalr, connectome, cellphonedb, natmi, logfc
from liana.method import rank_aggregate

warnings.filterwarnings("ignore")

BASE    = os.path.expanduser("~/gastric_tme_project")
PER_DIR = os.path.join(BASE, "data/processed/per_dataset")
OUT_DIR = os.path.join(BASE, "outputs/cell_communication")
os.makedirs(OUT_DIR, exist_ok=True)

KEY_INTERACTIONS = [
    # Myeloid → T cell immunosuppression
    ("Myeloid", "T/NK", "PDCD1LG2", "PDCD1"),    # PDL2-PD1
    ("Myeloid", "T/NK", "CD274", "PDCD1"),         # PDL1-PD1
    ("Myeloid", "T/NK", "LGALS9", "HAVCR2"),       # Galectin9-TIM3
    ("Myeloid", "T/NK", "IL10", "IL10RA"),          # IL10 immunosuppression
    ("Myeloid", "T/NK", "MIF", "CD74"),             # MIF immunosuppression
    # T cell → Myeloid recruitment
    ("T/NK", "Myeloid", "IFNG", "IFNGR1"),          # IFNg signaling
    ("T/NK", "Myeloid", "TNF", "TNFRSF1A"),         # TNF
]


def run_liana(adata, ct_col, groupby=None, min_cells=10, max_cells_per_type=500):
    """Run LIANA rank aggregate on an AnnData object."""
    # Filter to cell types with enough cells
    ct_counts = adata.obs[ct_col].value_counts()
    valid_cts  = ct_counts[ct_counts >= min_cells].index.tolist()
    adata_filt = adata[adata.obs[ct_col].isin(valid_cts)].copy()

    if adata_filt.n_obs < 100:
        print(f"  Too few cells after filtering: {adata_filt.n_obs}")
        return None

    # Downsample to avoid OOM — LIANA permutation test allocates O(cells * pairs)
    if max_cells_per_type and adata_filt.n_obs > max_cells_per_type * len(valid_cts):
        rng = np.random.default_rng(42)
        keep_mask = np.zeros(adata_filt.n_obs, dtype=bool)
        ct_arr = adata_filt.obs[ct_col].values
        for ct in valid_cts:
            pos = np.where(ct_arr == ct)[0]
            sel = rng.choice(pos, min(len(pos), max_cells_per_type), replace=False)
            keep_mask[sel] = True
        adata_filt = adata_filt[keep_mask].copy()

    print(f"  Running LIANA on {adata_filt.n_obs:,} cells, {len(valid_cts)} cell types ...")

    try:
        rank_aggregate(adata_filt, groupby=ct_col, use_raw=False, verbose=False,
                       n_perms=100, seed=42)
        return adata_filt.uns.get("liana_res", None)
    except Exception as e:
        print(f"  LIANA rank_aggregate failed: {e}")
        # Try with fewer perms
        try:
            cellphonedb(adata_filt, groupby=ct_col, use_raw=False, verbose=False,
                        n_perms=50, seed=42)
            return adata_filt.uns.get("liana_res", None)
        except Exception as e2:
            print(f"  cellphonedb also failed: {e2}")
            return None


def dotplot_interactions(res_df, title, outpath, top_n=30,
                         source_col="source", target_col="target",
                         ligand_col="ligand_complex", receptor_col="receptor_complex",
                         score_col="magnitude_rank"):
    """Dot plot of top ligand-receptor interactions."""
    if res_df is None or len(res_df) == 0:
        return

    # Filter to Myeloid ↔ T/NK interactions
    mask = (
        (res_df[source_col].str.contains("Myeloid|Macro|Mono", case=False) &
         res_df[target_col].str.contains("T|NK", case=False)) |
        (res_df[source_col].str.contains("T|NK", case=False) &
         res_df[target_col].str.contains("Myeloid|Macro|Mono", case=False))
    )
    df = res_df[mask].copy() if mask.sum() > 0 else res_df.copy()
    df = df.sort_values(score_col).head(top_n)

    if len(df) == 0:
        return

    df["interaction"] = df[ligand_col] + " → " + df[receptor_col]
    df["pair"] = df[source_col] + " → " + df[target_col]

    pivot = df.pivot_table(index="interaction", columns="pair",
                           values=score_col, aggfunc="min")

    fig, ax = plt.subplots(figsize=(max(8, len(pivot.columns) * 2),
                                    max(6, len(pivot) * 0.35)))
    sns.heatmap(pivot, cmap="RdYlBu_r", ax=ax, linewidths=0.3,
                cbar_kws={"label": score_col})
    ax.set_title(title, fontsize=11)
    plt.xticks(rotation=45, ha="right", fontsize=8)
    plt.yticks(fontsize=7)
    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()


def differential_interactions(res_slow, res_fast,
                               source_col="source", target_col="target",
                               ligand_col="ligand_complex", receptor_col="receptor_complex",
                               score_col="magnitude_rank"):
    """Find interactions enriched in slow vs fast progressors."""
    if res_slow is None or res_fast is None:
        return None

    def make_key(df):
        return (df[source_col] + "|" + df[target_col] + "|" +
                df[ligand_col] + "|" + df[receptor_col])

    slow = res_slow.copy(); slow["key"] = make_key(slow)
    fast = res_fast.copy(); fast["key"] = make_key(fast)

    merged = slow[["key", score_col]].rename(columns={score_col: "score_slow"}).merge(
             fast[["key", score_col]].rename(columns={score_col: "score_fast"}),
             on="key", how="outer").fillna(1.0)

    merged["delta"] = merged["score_fast"] - merged["score_slow"]
    merged[["source", "target", "ligand", "receptor"]] = (
        merged["key"].str.split("|", expand=True)
    )
    merged = merged.sort_values("delta")
    return merged


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=os.path.join(
        PER_DIR, "korea_kim2022_annotated_scored.h5ad"))
    parser.add_argument("--min-cells", type=int, default=20)
    args = parser.parse_args()

    print(f"Loading {args.input} ...")
    adata = sc.read_h5ad(args.input)
    print(f"  {adata.n_obs:,} cells")

    # LIANA only needs the L-R gene set (~2000 genes) — all within the 3000 HVGs
    # Do NOT restore raw here: expanding to 32k genes causes OOM in permutation tests

    ct_col = "cell_type_coarse"
    if ct_col not in adata.obs.columns:
        adata.obs[ct_col] = "cluster_" + adata.obs["leiden"].astype(str)

    # ── Full cohort LIANA ──────────────────────────────────────────────────────
    print("\n[LIANA — full cohort]")
    res_all = run_liana(adata, ct_col, min_cells=args.min_cells)
    if res_all is not None:
        res_all.to_csv(os.path.join(OUT_DIR, "liana_all.csv"), index=False)
        dotplot_interactions(res_all, "Cell-cell interactions — all samples",
                             os.path.join(OUT_DIR, "liana_dotplot_myeloid_tcell.png"))
        print("  saved: liana_dotplot_myeloid_tcell.png, liana_all.csv")

    # ── Stratified by progression ──────────────────────────────────────────────
    if "progression_category" in adata.obs.columns:
        for prog in ["Slow", "Fast"]:
            mask = adata.obs["progression_category"] == prog
            if mask.sum() < 200:
                print(f"  Too few {prog} cells: {mask.sum()}")
                continue
            print(f"\n[LIANA — {prog} progressors ({mask.sum():,} cells)]")
            adata_prog = adata[mask].copy()
            res = run_liana(adata_prog, ct_col, min_cells=args.min_cells)
            if res is not None:
                res.to_csv(os.path.join(OUT_DIR, f"liana_{prog.lower()}.csv"), index=False)
                dotplot_interactions(
                    res, f"Cell-cell interactions — {prog} progressors",
                    os.path.join(OUT_DIR, f"interaction_heatmap_{prog.lower()}.png"))
                print(f"  saved: interaction_heatmap_{prog.lower()}.png")

        # Differential
        slow_path = os.path.join(OUT_DIR, "liana_slow.csv")
        fast_path = os.path.join(OUT_DIR, "liana_fast.csv")
        if os.path.exists(slow_path) and os.path.exists(fast_path):
            res_slow = pd.read_csv(slow_path)
            res_fast = pd.read_csv(fast_path)
            diff = differential_interactions(res_slow, res_fast)
            if diff is not None:
                diff.to_csv(os.path.join(OUT_DIR, "differential_interactions.csv"), index=False)
                print("\n  Top interactions enriched in Fast progressors (immunosuppression):")
                print(diff.tail(10)[["source", "target", "ligand", "receptor", "delta"]]
                      .to_string(index=False))
                print("\n  Top interactions enriched in Slow progressors (anti-tumor):")
                print(diff.head(10)[["source", "target", "ligand", "receptor", "delta"]]
                      .to_string(index=False))
    else:
        print("  No progression_category column — skipping stratified analysis")

    print(f"\n{'='*60}")
    print(f"Cell communication complete. Outputs in: {OUT_DIR}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
