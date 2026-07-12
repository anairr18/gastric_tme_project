"""
16_pathway_enrichment.py
Pathway enrichment analysis on pseudobulk DE results (Slow vs Fast progressors).
Also runs enrichment on our marker gene sets.

Steps:
  1. Load all pseudobulk DE CSVs from outputs/meta_analysis/pseudobulk_de/
  2. Enrichr on significant genes per cell type
  3. Preranked GSEA on full T cell and Myeloid gene lists
  4. Heatmap of top DE genes across patients
  5. Volcano plots per cell type
  6. SPP1 pathway deepdive: extract SPP1+ macrophage L-R interactions from LIANA

Outputs in outputs/pathway_enrichment/
"""
import os, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import gseapy as gp

warnings.filterwarnings("ignore")

BASE    = os.path.expanduser("~/gastric_tme_project")
DE_DIR  = os.path.join(BASE, "outputs/meta_analysis/pseudobulk_de")
OUT_DIR = os.path.join(BASE, "outputs/pathway_enrichment")
os.makedirs(OUT_DIR, exist_ok=True)

GENE_SETS_ENRICHR = [
    "MSigDB_Hallmark_2020",
    "KEGG_2021_Human",
    "GO_Biological_Process_2021",
    "Reactome_2022",
    "WikiPathways_2019_Human",
]

CELL_TYPES_OF_INTEREST = ["Myeloid", "T_NK", "Fibroblast", "Epithelial"]


def load_de_results():
    """Load all pseudobulk DE CSVs."""
    de = {}
    for fname in os.listdir(DE_DIR):
        if not fname.endswith(".csv"):
            continue
        ct = fname.replace("de_", "").replace("_slow_vs_fast.csv", "")
        df = pd.read_csv(os.path.join(DE_DIR, fname), index_col=0)
        de[ct] = df
        n_sig = (df["padj"] < 0.05).sum() if "padj" in df.columns else 0
        print(f"  {ct}: {len(df):,} genes, {n_sig} sig (padj<0.05)")
    return de


def volcano_plot(de_df, title, out_path, padj_thresh=0.05, lfc_thresh=0.5):
    """Volcano plot: -log10(p) vs log2FC."""
    df = de_df.copy().dropna(subset=["pvalue", "log2FoldChange"])
    df["-log10p"] = -np.log10(df["pvalue"].clip(lower=1e-300))

    fig, ax = plt.subplots(figsize=(7, 6))

    # Not significant
    mask_ns = (df["padj"] >= padj_thresh) | (df["padj"].isna())
    ax.scatter(df[mask_ns]["log2FoldChange"], df[mask_ns]["-log10p"],
               s=3, alpha=0.3, color="gray")

    # Significant
    mask_sig = (df["padj"] < padj_thresh) & (~df["padj"].isna())
    mask_up   = mask_sig & (df["log2FoldChange"] >  lfc_thresh)
    mask_down = mask_sig & (df["log2FoldChange"] < -lfc_thresh)
    ax.scatter(df[mask_up]["log2FoldChange"],   df[mask_up]["-log10p"],
               s=8, alpha=0.8, color="firebrick",  label=f"Up in Fast ({mask_up.sum()})")
    ax.scatter(df[mask_down]["log2FoldChange"], df[mask_down]["-log10p"],
               s=8, alpha=0.8, color="steelblue", label=f"Up in Slow ({mask_down.sum()})")

    # Label top genes
    top = df[mask_sig].nlargest(10, "-log10p")
    for _, row in top.iterrows():
        ax.annotate(row.name, (row["log2FoldChange"], row["-log10p"]),
                    fontsize=6, xytext=(3, 3), textcoords="offset points")

    ax.axvline(lfc_thresh,  color="black", linestyle="--", linewidth=0.8)
    ax.axvline(-lfc_thresh, color="black", linestyle="--", linewidth=0.8)
    ax.axhline(-np.log10(padj_thresh), color="black", linestyle=":", linewidth=0.8)
    ax.set_xlabel("log2 Fold Change (Fast / Slow)")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title(title)
    ax.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


def enrichr_analysis(gene_list, label, out_dir, gene_sets=None):
    """Run Enrichr on a gene list."""
    if gene_sets is None:
        gene_sets = GENE_SETS_ENRICHR
    if len(gene_list) < 3:
        print(f"  SKIP enrichr for {label}: only {len(gene_list)} genes")
        return None

    try:
        enr = gp.enrichr(
            gene_list=gene_list,
            gene_sets=gene_sets[:3],   # top 3 databases
            organism="human",
            outdir=None,
            verbose=False,
        )
        res = enr.res2d
        if res is not None and len(res) > 0:
            top = res[res["Adjusted P-value"] < 0.05].head(20)
            if len(top) > 0:
                top.to_csv(os.path.join(out_dir, f"enrichr_{label}.csv"), index=False)
                # Dot plot of top 15 terms
                plot_enrichr(top, label, os.path.join(out_dir, f"enrichr_{label}.png"))
            return res
    except Exception as e:
        print(f"  Enrichr {label} failed: {e}")
    return None


def plot_enrichr(df, title, outpath, top_n=15):
    df = df.head(top_n).copy()
    df["-log10padj"] = -np.log10(df["Adjusted P-value"].clip(lower=1e-50))
    df["Term"] = df["Term"].str[:60]
    df = df.sort_values("-log10padj")

    fig, ax = plt.subplots(figsize=(8, max(4, len(df) * 0.35)))
    colors = plt.cm.RdYlGn(df["-log10padj"] / df["-log10padj"].max())
    bars = ax.barh(df["Term"], df["-log10padj"], color=colors)
    ax.set_xlabel("-log10(Adjusted P-value)")
    ax.set_title(f"Enrichr: {title}", fontsize=10)
    for bar, val in zip(bars, df["-log10padj"]):
        ax.text(val + 0.05, bar.get_y() + bar.get_height()/2,
                f"{val:.2f}", va="center", fontsize=7)
    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()


def preranked_gsea(de_df, label, out_dir):
    """Preranked GSEA using -log10(p) * sign(LFC) as rank metric."""
    if "pvalue" not in de_df.columns or "log2FoldChange" not in de_df.columns:
        return
    df = de_df.dropna(subset=["pvalue", "log2FoldChange"]).copy()
    df["rank_metric"] = -np.log10(df["pvalue"].clip(lower=1e-300)) * np.sign(df["log2FoldChange"])
    rnk = df["rank_metric"].sort_values(ascending=False)
    rnk.index = df.index

    try:
        pre_res = gp.prerank(
            rnk=rnk,
            gene_sets=["MSigDB_Hallmark_2020"],
            min_size=5, max_size=500,
            permutation_num=100,
            outdir=None, verbose=False, seed=42,
        )
        result = pre_res.res2d
        if result is not None and len(result) > 0:
            sig = result[result["FDR q-val"] < 0.25].sort_values("NES", ascending=False)
            sig.to_csv(os.path.join(out_dir, f"gsea_preranked_{label}.csv"), index=False)
            print(f"  {label} GSEA: {len(sig)} significant pathways (FDR<0.25)")
            if len(sig) > 0:
                print(sig[["Term", "NES", "FDR q-val"]].head(10).to_string(index=False))
    except Exception as e:
        print(f"  Preranked GSEA {label} failed: {e}")


def heatmap_top_genes(de_results, out_dir, n_per_ct=10):
    """Heatmap of top DE genes across cell types."""
    sig_genes_per_ct = {}
    for ct, df in de_results.items():
        if "padj" not in df.columns:
            continue
        sig = df[df["padj"] < 0.05].copy()
        if len(sig) == 0:
            sig = df.nsmallest(n_per_ct, "pvalue")
        sig_genes_per_ct[ct] = sig.nsmallest(n_per_ct, "padj") if "padj" in sig.columns else sig.head(n_per_ct)

    if not sig_genes_per_ct:
        return

    # Build LFC matrix: genes x cell_types
    all_genes = list(set(g for df in sig_genes_per_ct.values() for g in df.index))
    lfc_matrix = pd.DataFrame(index=all_genes, columns=list(de_results.keys()), dtype=float)

    for ct, df in de_results.items():
        if "log2FoldChange" in df.columns:
            for gene in all_genes:
                if gene in df.index:
                    lfc_matrix.loc[gene, ct] = df.loc[gene, "log2FoldChange"]

    lfc_matrix = lfc_matrix.dropna(how="all").fillna(0)
    if len(lfc_matrix) == 0:
        return

    fig, ax = plt.subplots(figsize=(max(8, lfc_matrix.shape[1] * 1.5),
                                    max(6, lfc_matrix.shape[0] * 0.3)))
    sns.heatmap(lfc_matrix.astype(float), cmap="RdBu_r", center=0,
                ax=ax, linewidths=0.3, yticklabels=True,
                cbar_kws={"label": "log2 FC (Fast / Slow)"})
    ax.set_title("Top DE genes: Fast vs Slow progressors\n(Positive = higher in Fast)", fontsize=11)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "de_heatmap_top_genes.png"), dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  saved: de_heatmap_top_genes.png ({lfc_matrix.shape[0]} genes)")


def spp1_liana_deepdive(out_dir):
    """Extract SPP1+ macrophage L-R interactions from LIANA differential results."""
    diff_path = os.path.join(BASE, "outputs/cell_communication/differential_interactions.csv")
    if not os.path.exists(diff_path):
        print("  No differential_interactions.csv found")
        return

    diff = pd.read_csv(diff_path)
    print(f"\n[SPP1 L-R deepdive] ({len(diff)} interactions)")

    # Focus on Myeloid source interactions
    myeloid_src = diff[diff["source"].str.contains("Myeloid|Macro|Mono", case=False, na=False)]
    if len(myeloid_src) == 0:
        myeloid_src = diff.head(50)

    # Top enriched in Slow (anti-tumor) vs Fast (immunosuppressive)
    top_slow = diff.nsmallest(10, "delta")   # negative delta = enriched in Slow
    top_fast = diff.nlargest(10, "delta")    # positive delta = enriched in Fast

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for ax, df_sub, title, color in [
        (axes[0], top_slow, "Enriched in Slow (anti-tumor)", "steelblue"),
        (axes[1], top_fast, "Enriched in Fast (immunosuppressive)", "firebrick"),
    ]:
        if "ligand" in df_sub.columns and "receptor" in df_sub.columns:
            labels = (df_sub["source"].str[:8] + ": " +
                      df_sub["ligand"].fillna("?") + " -> " +
                      df_sub["receptor"].fillna("?"))
        else:
            labels = df_sub.index.astype(str)
        ax.barh(range(len(df_sub)), df_sub["delta"].abs(), color=color)
        ax.set_yticks(range(len(df_sub)))
        ax.set_yticklabels(labels.values, fontsize=7)
        ax.set_xlabel("|delta score|")
        ax.set_title(title, fontsize=10)

    plt.suptitle("Differential L-R interactions: Slow vs Fast progressors", fontsize=11)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "liana_spp1_deepdive.png"), dpi=150, bbox_inches="tight")
    plt.close()
    print("  saved: liana_spp1_deepdive.png")


def main():
    print("[Loading DE results]")
    de_results = load_de_results()

    print("\n[Volcano plots]")
    for ct, df in de_results.items():
        if "pvalue" in df.columns and "log2FoldChange" in df.columns:
            volcano_plot(df, f"{ct}: Fast vs Slow progressors",
                         os.path.join(OUT_DIR, f"volcano_{ct}.png"))
            print(f"  saved: volcano_{ct}.png")

    print("\n[Enrichr analysis]")
    for ct, df in de_results.items():
        if "padj" not in df.columns:
            continue
        sig = df[df["padj"] < 0.05].index.tolist()
        if len(sig) < 3:
            # Use top 50 by p-value as background
            sig = df.nsmallest(50, "pvalue").index.tolist()
        print(f"  {ct}: {len(sig)} genes -> Enrichr")
        enrichr_analysis(sig, ct, OUT_DIR)

    print("\n[Preranked GSEA — Myeloid + T/NK]")
    for ct in ["Myeloid", "T_NK", "Fibroblast"]:
        if ct in de_results:
            preranked_gsea(de_results[ct], ct, OUT_DIR)

    print("\n[DE heatmap]")
    heatmap_top_genes(de_results, OUT_DIR)

    print("\n[SPP1 L-R deepdive]")
    spp1_liana_deepdive(OUT_DIR)

    print(f"\n{'='*60}")
    print(f"Pathway enrichment complete. Outputs in: {OUT_DIR}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
