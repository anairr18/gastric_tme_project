"""
17_figure_assembly.py
Assemble publication-quality composite figures from analysis outputs.
Each figure corresponds to a main paper figure for Gut submission.

Figure layout:
  Fig 1 — Cohort overview: composition + KM top hits
  Fig 2 — SPP1+ macrophage: subtype UMAP + bar + KM
  Fig 3 — T cell trajectory: PAGA + pseudotime + exhaustion
  Fig 4 — Cell communication: heatmap + SPP1 deepdive
  Fig 5 — TCGA validation: Cox forest + KM curves
  Fig 6 — TME predictor: ROC + C-index + coefficients
"""
import os, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.image import imread
import seaborn as sns

warnings.filterwarnings("ignore")

BASE = os.path.expanduser("~/gastric_tme_project")

# Output dirs
OUT_DIR   = os.path.join(BASE, "outputs/figures")
SURV_DIR  = os.path.join(BASE, "outputs/survival")
META_DIR  = os.path.join(BASE, "outputs/meta_analysis")
SUB_DIR   = os.path.join(BASE, "outputs/subclustering")
TRAJ_DIR  = os.path.join(BASE, "outputs/trajectory")
COMM_DIR  = os.path.join(BASE, "outputs/cell_communication")
TCGA_DIR  = os.path.join(BASE, "outputs/tcga_validation")
PRED_DIR  = os.path.join(BASE, "outputs/tme_predictor")
PATH_DIR  = os.path.join(BASE, "outputs/pathway_enrichment")
GEO_DIR   = os.path.join(BASE, "outputs/geo_validation")
os.makedirs(OUT_DIR, exist_ok=True)

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8,
    "figure.dpi": 300,
})


def load_img(path):
    """Load image or return None."""
    if os.path.exists(path):
        return imread(path)
    return None


def show_img(ax, path, title=None):
    """Show image in axis with no ticks."""
    img = load_img(path)
    if img is not None:
        ax.imshow(img)
    else:
        ax.text(0.5, 0.5, f"Missing:\n{os.path.basename(path)}", ha="center", va="center",
                transform=ax.transAxes, fontsize=7, color="red")
        ax.set_facecolor("#fff8f8")
    ax.axis("off")
    if title:
        ax.set_title(title, fontsize=9, fontweight="bold")


def label_panel(ax, letter, x=-0.05, y=1.08):
    """Add panel label (A, B, C, ...) to axis."""
    ax.text(x, y, letter, transform=ax.transAxes,
            fontsize=14, fontweight="bold", va="top", ha="right")


# ============================================================
# FIGURE 1: Cohort overview
# ============================================================
def figure1():
    fig = plt.figure(figsize=(18, 12))
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.3)

    # A: TME composition heatmap (Slow vs Fast)
    ax_a = fig.add_subplot(gs[0, 0])
    show_img(ax_a, os.path.join(META_DIR, "tme_composition_heatmap.png"),
             "A. TME composition by progression")
    label_panel(ax_a, "A")

    # B: Composition stacked bar (Slow vs Fast)
    ax_b = fig.add_subplot(gs[0, 1])
    show_img(ax_b, os.path.join(META_DIR, "composition_slow_vs_fast.png"),
             "B. Cell type proportions")
    label_panel(ax_b, "B")

    # C: Top KM curve (M2 macrophage or exhaustion score)
    ax_c = fig.add_subplot(gs[0, 2])
    km_path = os.path.join(SURV_DIR, "km_PFS_score_M2_Macrophage.png")
    if not os.path.exists(km_path):
        km_path = os.path.join(SURV_DIR, "km_PFS_score_Exhausted_CD8+_T.png")
    show_img(ax_c, km_path, "C. M2 Macrophage score — PFS")
    label_panel(ax_c, "C")

    # D: PDL1 vs TME scatter
    ax_d = fig.add_subplot(gs[1, 0])
    show_img(ax_d, os.path.join(SURV_DIR, "pdl1_vs_tme.png"),
             "D. PDL1 vs TME predictors")
    label_panel(ax_d, "D")

    # E: Exhaustion score KM
    ax_e = fig.add_subplot(gs[1, 1])
    show_img(ax_e, os.path.join(SURV_DIR, "km_PFS_exhaustion_score.png"),
             "E. Exhaustion score — PFS")
    label_panel(ax_e, "E")

    # F: Macrophage polarization
    ax_f = fig.add_subplot(gs[1, 2])
    show_img(ax_f, os.path.join(META_DIR, "macrophage_polarization.png"),
             "F. Macrophage M1/M2 balance")
    label_panel(ax_f, "F")

    fig.suptitle("Figure 1: Gastric cancer TME landscape in immunotherapy response\n"
                 "Korean cohort (n=33 patients)", fontsize=12, fontweight="bold", y=0.98)
    plt.savefig(os.path.join(OUT_DIR, "Fig1_cohort_overview.pdf"), bbox_inches="tight")
    plt.savefig(os.path.join(OUT_DIR, "Fig1_cohort_overview.png"), bbox_inches="tight", dpi=300)
    plt.close()
    print("  saved: Fig1_cohort_overview")


# ============================================================
# FIGURE 2: SPP1+ macrophage deep dive
# ============================================================
def figure2():
    fig = plt.figure(figsize=(18, 12))
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.3)

    # A: Myeloid subtype UMAP
    ax_a = fig.add_subplot(gs[0, 0:2])
    show_img(ax_a, os.path.join(SUB_DIR, "umap_myeloid_subtypes.png"),
             "A. Myeloid subtypes (UMAP)")
    label_panel(ax_a, "A")

    # B: SPP1 dotplot
    ax_b = fig.add_subplot(gs[0, 2])
    show_img(ax_b, os.path.join(SUB_DIR, "dotplot__myeloid_markers.png"),
             "B. Myeloid marker expression")
    label_panel(ax_b, "B")

    # C: SPP1 Macro by progression
    ax_c = fig.add_subplot(gs[1, 0])
    show_img(ax_c, os.path.join(SUB_DIR, "myeloid_by_progression.png"),
             "C. SPP1_Macro abundance by progression")
    label_panel(ax_c, "C")

    # D: KM curve for M2 macrophage
    ax_d = fig.add_subplot(gs[1, 1])
    show_img(ax_d, os.path.join(SURV_DIR, "km_PFS_score_M2_Macrophage.png"),
             "D. M2 Macrophage score — PFS")
    label_panel(ax_d, "D")

    # E: Myeloid DE volcano
    ax_e = fig.add_subplot(gs[1, 2])
    show_img(ax_e, os.path.join(PATH_DIR, "volcano_Myeloid.png"),
             "E. Myeloid: Slow vs Fast (DEG)")
    label_panel(ax_e, "E")

    # Add key finding text
    fig.text(0.5, 0.01,
             "SPP1+ macrophages: 16.9% in Slow vs 7.3% in Fast progressors (p<0.05)\n"
             "Key immunosuppressive niche enriched in poor responders",
             ha="center", va="bottom", fontsize=9, style="italic",
             bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))

    fig.suptitle("Figure 2: SPP1+ macrophages define immunosuppressive niche in poor responders",
                 fontsize=12, fontweight="bold", y=0.98)
    plt.savefig(os.path.join(OUT_DIR, "Fig2_SPP1_macrophage.pdf"), bbox_inches="tight")
    plt.savefig(os.path.join(OUT_DIR, "Fig2_SPP1_macrophage.png"), bbox_inches="tight", dpi=300)
    plt.close()
    print("  saved: Fig2_SPP1_macrophage")


# ============================================================
# FIGURE 3: T cell exhaustion trajectory
# ============================================================
def figure3():
    fig = plt.figure(figsize=(18, 12))
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.3)

    ax_a = fig.add_subplot(gs[0, 0])
    show_img(ax_a, os.path.join(TRAJ_DIR, "paga_graph_tcells.png"),
             "A. T cell PAGA graph")
    label_panel(ax_a, "A")

    ax_b = fig.add_subplot(gs[0, 1])
    show_img(ax_b, os.path.join(TRAJ_DIR, "umap_pseudotime.png"),
             "B. Pseudotime (diffusion)")
    label_panel(ax_b, "B")

    ax_c = fig.add_subplot(gs[0, 2])
    show_img(ax_c, os.path.join(TRAJ_DIR, "pseudotime_vs_exhaustion.png"),
             "C. Pseudotime vs exhaustion")
    label_panel(ax_c, "C")

    ax_d = fig.add_subplot(gs[1, 0])
    show_img(ax_d, os.path.join(TRAJ_DIR, "pseudotime_by_progression.png"),
             "D. Pseudotime by progression")
    label_panel(ax_d, "D")

    ax_e = fig.add_subplot(gs[1, 1])
    show_img(ax_e, os.path.join(SUB_DIR, "umap_t_cell_subtypes.png"),
             "E. T cell subtypes (UMAP)")
    label_panel(ax_e, "E")

    ax_f = fig.add_subplot(gs[1, 2])
    show_img(ax_f, os.path.join(META_DIR, "cd8_exhaustion_scores.png"),
             "F. CD8 exhaustion by progression")
    label_panel(ax_f, "F")

    fig.suptitle("Figure 3: Accelerated T cell exhaustion trajectory in fast progressors",
                 fontsize=12, fontweight="bold", y=0.98)
    plt.savefig(os.path.join(OUT_DIR, "Fig3_trajectory.pdf"), bbox_inches="tight")
    plt.savefig(os.path.join(OUT_DIR, "Fig3_trajectory.png"), bbox_inches="tight", dpi=300)
    plt.close()
    print("  saved: Fig3_trajectory")


# ============================================================
# FIGURE 4: Cell-cell communication
# ============================================================
def figure4():
    fig = plt.figure(figsize=(18, 10))
    gs = gridspec.GridSpec(1, 3, figure=fig, hspace=0.3, wspace=0.35)

    ax_a = fig.add_subplot(gs[0, 0])
    show_img(ax_a, os.path.join(COMM_DIR, "interaction_heatmap_slow.png"),
             "A. Interactions enriched in Slow")
    label_panel(ax_a, "A")

    ax_b = fig.add_subplot(gs[0, 1])
    show_img(ax_b, os.path.join(COMM_DIR, "interaction_heatmap_fast.png"),
             "B. Interactions enriched in Fast")
    label_panel(ax_b, "B")

    ax_c = fig.add_subplot(gs[0, 2])
    show_img(ax_c, os.path.join(PATH_DIR, "liana_spp1_deepdive.png"),
             "C. SPP1+ L-R interactions differential")
    label_panel(ax_c, "C")

    fig.suptitle("Figure 4: Differential cell-cell communication in Slow vs Fast progressors",
                 fontsize=12, fontweight="bold", y=1.02)
    plt.savefig(os.path.join(OUT_DIR, "Fig4_communication.pdf"), bbox_inches="tight")
    plt.savefig(os.path.join(OUT_DIR, "Fig4_communication.png"), bbox_inches="tight", dpi=300)
    plt.close()
    print("  saved: Fig4_communication")


# ============================================================
# FIGURE 5: TCGA validation
# ============================================================
def figure5():
    fig = plt.figure(figsize=(18, 12))
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.35, wspace=0.3)

    # A: Cox forest plot (ssGSEA)
    ax_a = fig.add_subplot(gs[0, 0:2])
    show_img(ax_a, os.path.join(TCGA_DIR, "tcga_cox_forest.png"),
             "A. TCGA-STAD: Multivariate Cox regression")
    label_panel(ax_a, "A")

    # B: CAF KM (most significant: HR=2.31, p=0.033)
    ax_b = fig.add_subplot(gs[0, 2])
    show_img(ax_b, os.path.join(TCGA_DIR, "km_tcga_CAF.png"),
             "B. CAF signature — OS (TCGA-STAD)")
    label_panel(ax_b, "B")

    # C: NNLS deconv forest plot
    ax_c = fig.add_subplot(gs[1, 0:2])
    show_img(ax_c, os.path.join(TCGA_DIR, "tcga_deconv_cox_forest.png"),
             "C. TCGA-STAD: NNLS deconvolution Cox")
    label_panel(ax_c, "C")

    # D: FOLR2 Macro KM
    ax_d = fig.add_subplot(gs[1, 2])
    show_img(ax_d, os.path.join(TCGA_DIR, "km_tcga_FOLR2_Macro.png"),
             "D. FOLR2_Macro signature — OS (TCGA-STAD)")
    label_panel(ax_d, "D")

    fig.text(0.5, 0.01,
             "TCGA-STAD validation (n=397): CAF signature HR=2.31 (p=0.033)\n"
             "FOLR2_Macro HR=2.17 (p=0.059) | Multivariate Cox C-index=0.577",
             ha="center", va="bottom", fontsize=9, style="italic",
             bbox=dict(boxstyle="round", facecolor="lightcyan", alpha=0.8))

    fig.suptitle("Figure 5: TCGA-STAD (n=397) validates TME prognostic signatures",
                 fontsize=12, fontweight="bold", y=0.98)
    plt.savefig(os.path.join(OUT_DIR, "Fig5_TCGA_validation.pdf"), bbox_inches="tight")
    plt.savefig(os.path.join(OUT_DIR, "Fig5_TCGA_validation.png"), bbox_inches="tight", dpi=300)
    plt.close()
    print("  saved: Fig5_TCGA_validation")


# ============================================================
# FIGURE 6: TME predictor
# ============================================================
def figure6():
    fig = plt.figure(figsize=(18, 10))
    gs = gridspec.GridSpec(1, 3, figure=fig, hspace=0.3, wspace=0.35)

    ax_a = fig.add_subplot(gs[0, 0])
    show_img(ax_a, os.path.join(PRED_DIR, "roc_response.png"),
             "A. ROC: TME features vs response")
    label_panel(ax_a, "A")

    ax_b = fig.add_subplot(gs[0, 1])
    show_img(ax_b, os.path.join(PRED_DIR, "cindex_bar.png"),
             "B. C-index: composite vs PDL1 CPS")
    label_panel(ax_b, "B")

    ax_c = fig.add_subplot(gs[0, 2])
    # Show LASSO coefficients as bar chart
    coef_path = os.path.join(PRED_DIR, "predictor_coefficients.csv")
    if os.path.exists(coef_path):
        coefs = pd.read_csv(coef_path, index_col=0)
        coefs_nonzero = coefs[coefs["coef"] != 0].sort_values("coef")
        if len(coefs_nonzero) > 0:
            colors = ["firebrick" if x > 0 else "steelblue" for x in coefs_nonzero["coef"]]
            coefs_nonzero["coef"].plot(kind="barh", ax=ax_c, color=colors)
            ax_c.axvline(0, color="black", linewidth=0.8)
            ax_c.set_xlabel("LASSO coefficient")
            ax_c.set_title("C. LASSO predictor coefficients")
            ax_c.set_ylabel("")
        else:
            show_img(ax_c, None, "C. LASSO coefficients (no non-zero)")
    else:
        show_img(ax_c, None, "C. LASSO coefficients (missing)")
    label_panel(ax_c, "C")

    fig.text(0.5, 0.01,
             "LASSO predictor: M2_Macro (-0.22), M1_M2_ratio (+0.09), Exhaustion (+0.20)\n"
             "Composite TME score outperforms PDL1 CPS for PFS prediction",
             ha="center", va="bottom", fontsize=9, style="italic",
             bbox=dict(boxstyle="round", facecolor="lightyellow", alpha=0.8))

    fig.suptitle("Figure 6: Composite TME predictor outperforms PDL1 CPS",
                 fontsize=12, fontweight="bold", y=1.02)
    plt.savefig(os.path.join(OUT_DIR, "Fig6_predictor.pdf"), bbox_inches="tight")
    plt.savefig(os.path.join(OUT_DIR, "Fig6_predictor.png"), bbox_inches="tight", dpi=300)
    plt.close()
    print("  saved: Fig6_predictor")


# ============================================================
# SUPPLEMENTARY FIGURES
# ============================================================
def supp_figure_km():
    """All 24 KM curves as supplementary panel."""
    km_files = sorted([f for f in os.listdir(SURV_DIR) if f.startswith("km_PFS_") and f.endswith(".png")])
    if not km_files:
        print("  No KM files found")
        return

    n = len(km_files)
    ncols = 6
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 3.5 * nrows))
    axes = axes.flatten()

    for i, fname in enumerate(km_files):
        show_img(axes[i], os.path.join(SURV_DIR, fname),
                 fname.replace("km_PFS_", "").replace(".png", "").replace("_", " "))

    for j in range(i + 1, len(axes)):
        axes[j].axis("off")

    plt.suptitle("Supplementary: KM curves for all TME features (Korean cohort, PFS)",
                 fontsize=12, fontweight="bold", y=1.01)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "SuppFig_KM_all.pdf"), bbox_inches="tight")
    plt.savefig(os.path.join(OUT_DIR, "SuppFig_KM_all.png"), bbox_inches="tight", dpi=200)
    plt.close()
    print(f"  saved: SuppFig_KM_all ({n} curves)")


def supp_figure_de():
    """Volcano plots + DE heatmap + pathway enrichment."""
    fig = plt.figure(figsize=(20, 16))
    gs = gridspec.GridSpec(3, 4, figure=fig, hspace=0.4, wspace=0.3)

    ct_order = ["Myeloid", "T_NK", "Fibroblast", "Epithelial", "Endothelial", "Stromal", "B_Plasma"]
    for i, ct in enumerate(ct_order):
        row, col = divmod(i, 4)
        ax = fig.add_subplot(gs[row, col])
        show_img(ax, os.path.join(PATH_DIR, f"volcano_{ct}.png"), f"{ct}")

    # Heatmap in last slot
    ax_hm = fig.add_subplot(gs[1, 3])
    show_img(ax_hm, os.path.join(PATH_DIR, "de_heatmap_top_genes.png"), "DE heatmap")

    plt.suptitle("Supplementary: Pseudobulk DE (Slow vs Fast) and pathway analysis",
                 fontsize=12, fontweight="bold")
    plt.savefig(os.path.join(OUT_DIR, "SuppFig_DE_pathways.pdf"), bbox_inches="tight")
    plt.savefig(os.path.join(OUT_DIR, "SuppFig_DE_pathways.png"), bbox_inches="tight", dpi=200)
    plt.close()
    print("  saved: SuppFig_DE_pathways")


def supp_figure_gsea():
    """GSEA preranked results for Myeloid and Fibroblast."""
    for ct in ["Myeloid", "Fibroblast"]:
        csv = os.path.join(PATH_DIR, f"gsea_preranked_{ct}.csv")
        if not os.path.exists(csv):
            continue
        df = pd.read_csv(csv)
        if "FDR q-val" not in df.columns:
            continue
        df = df[df["FDR q-val"] < 0.25].sort_values("NES", ascending=False).head(20)
        if len(df) == 0:
            continue
        df["Term"] = df["Term"].str.replace("MSigDB_Hallmark_2020__", "").str[:50]

        fig, ax = plt.subplots(figsize=(10, max(4, len(df) * 0.4)))
        colors = ["firebrick" if n > 0 else "steelblue" for n in df["NES"]]
        ax.barh(range(len(df)), df["NES"], color=colors)
        ax.set_yticks(range(len(df)))
        ax.set_yticklabels(df["Term"].values, fontsize=7)
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_xlabel("NES (Normalized Enrichment Score)")
        ax.set_title(f"{ct}: MSigDB Hallmark pathway enrichment (Slow vs Fast progressors)\n"
                     f"(positive = enriched in Fast, negative = enriched in Slow)", fontsize=9)
        for i, (_, row) in enumerate(df.iterrows()):
            ax.text(row["NES"] + 0.03, i, f"q={row['FDR q-val']:.3f}", va="center", fontsize=6)
        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, f"SuppFig_GSEA_{ct}.pdf"), bbox_inches="tight")
        plt.savefig(os.path.join(OUT_DIR, f"SuppFig_GSEA_{ct}.png"), bbox_inches="tight", dpi=200)
        plt.close()
        print(f"  saved: SuppFig_GSEA_{ct} ({len(df)} pathways)")


def supp_figure_geo():
    """GEO validation KM curves (if available)."""
    geo_files = [f for f in os.listdir(GEO_DIR) if f.startswith("km_") and f.endswith(".png")]
    if not geo_files:
        print("  No GEO KM curves available yet")
        return

    n = len(geo_files)
    ncols = min(n, 4)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 4.5 * nrows))
    if n == 1:
        axes = [axes]
    else:
        axes = axes.flatten() if nrows > 1 else axes

    for i, fname in enumerate(geo_files):
        show_img(axes[i], os.path.join(GEO_DIR, fname),
                 fname.replace("km_", "").replace(".png", "").replace("_", " ")[:40])

    for j in range(i + 1, len(axes)):
        axes[j].axis("off")

    plt.suptitle("Supplementary: GEO cohort validation (ssGSEA + KM)",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "SuppFig_GEO_validation.pdf"), bbox_inches="tight")
    plt.savefig(os.path.join(OUT_DIR, "SuppFig_GEO_validation.png"), bbox_inches="tight", dpi=200)
    plt.close()
    print(f"  saved: SuppFig_GEO_validation ({n} plots)")


def main():
    print("[Figure 1: Cohort overview]")
    figure1()

    print("[Figure 2: SPP1+ macrophage]")
    figure2()

    print("[Figure 3: T cell trajectory]")
    figure3()

    print("[Figure 4: Cell communication]")
    figure4()

    print("[Figure 5: TCGA validation]")
    figure5()

    print("[Figure 6: TME predictor]")
    figure6()

    print("[Supplementary: All KM curves]")
    supp_figure_km()

    print("[Supplementary: DE + Volcano + Pathways]")
    supp_figure_de()

    print("[Supplementary: GSEA results]")
    supp_figure_gsea()

    print("[Supplementary: GEO validation]")
    supp_figure_geo()

    print(f"\n{'='*60}")
    print(f"Figure assembly complete. Outputs in: {OUT_DIR}")
    print(f"{'='*60}")
    print("\nFigures generated:")
    for f in sorted(os.listdir(OUT_DIR)):
        size_kb = os.path.getsize(os.path.join(OUT_DIR, f)) // 1024
        print(f"  {f} ({size_kb} KB)")


if __name__ == "__main__":
    main()
