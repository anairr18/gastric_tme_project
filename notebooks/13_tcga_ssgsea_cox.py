"""
13_tcga_ssgsea_cox.py
ssGSEA scoring of TCGA-STAD using our cell type marker gene sets,
followed by multivariate Cox proportional-hazards regression.

Gene sets:
  SPP1_Macro, FOLR2_Macro, Exhausted_CD8, Effector_CD8,
  Treg, Inflammatory_M1, CAF, Mast

Outputs in outputs/tcga_validation/:
  tcga_ssgsea_scores.csv
  tcga_cox_forest.png
  tcga_cox_results.csv
  km_tcga_SPP1_Macro.png
  km_tcga_Exhausted_CD8.png
  tcga_multivariate_cox.csv
"""
import os, sys, json, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import gseapy as gp
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
import requests, mygene

warnings.filterwarnings("ignore")

BASE    = os.path.expanduser("~/gastric_tme_project")
RAW_DIR = os.path.join(BASE, "data/external/tcga_stad")
OUT_DIR = os.path.join(BASE, "outputs/tcga_validation")
os.makedirs(OUT_DIR, exist_ok=True)

GENE_SETS = {
    "SPP1_Macro":       ["SPP1", "APOC1", "APOE", "C1QB", "TREM2", "GPNMB", "FN1"],
    "FOLR2_Macro":      ["FOLR2", "LYVE1", "SELENOP", "MRC1", "CD163", "STAB1"],
    "Exhausted_CD8":    ["PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "ENTPD1", "CTLA4", "NR4A1"],
    "Effector_CD8":     ["GZMB", "PRF1", "IFNG", "CX3CR1", "KLRG1", "GNLY"],
    "Treg":             ["FOXP3", "IL2RA", "IKZF2", "TNFRSF18", "LAYN"],
    "Inflammatory_M1":  ["IL1B", "TNF", "CXCL9", "CXCL10", "IDO1", "CXCL11"],
    "Immunosuppressive":["IL10", "TGM2", "CCL18", "ARG1", "PPARG"],
    "CAF":              ["ACTA2", "FAP", "PDGFRB", "THY1", "POSTN", "MMP2"],
}


def load_and_convert_tcga_counts():
    """Load cached raw counts and convert ENSEMBL -> symbol, caching result."""
    sym_cache = os.path.join(RAW_DIR, "tcga_stad_counts_symbols.parquet")
    if os.path.exists(sym_cache):
        print("Loading cached symbol-mapped TCGA counts ...")
        return pd.read_parquet(sym_cache)

    raw_cache = os.path.join(RAW_DIR, "tcga_stad_counts.csv")
    if not os.path.exists(raw_cache):
        print("ERROR: tcga_stad_counts.csv not found. Run 09_tcga_validation.py first.")
        return None

    print("Loading raw TCGA counts ...")
    counts = pd.read_csv(raw_cache, index_col=0)
    # Strip ENSEMBL version
    counts.index = counts.index.str.split(".").str[0]

    print(f"  Converting {len(counts):,} ENSEMBL IDs to gene symbols via mygene ...")
    mg = mygene.MyGeneInfo()
    result = mg.querymany(counts.index.tolist(), scopes="ensembl.gene",
                          fields="symbol", species="human",
                          as_dataframe=True, verbose=False)
    if "symbol" in result.columns:
        mapping = result["symbol"].dropna().to_dict()
    else:
        mapping = {}
    counts.index = [mapping.get(g, g) for g in counts.index]
    counts = counts.groupby(counts.index).mean()
    print(f"  {len(counts):,} genes after symbol mapping")
    counts.to_parquet(sym_cache)
    return counts


def run_ssgsea(counts_df):
    """Run ssGSEA and return samples x gene_sets score matrix."""
    print("Running ssGSEA ...")
    # Uppercase gene names to match gene set symbols
    counts_up = counts_df.copy()
    counts_up.index = counts_up.index.str.upper()
    counts_up = counts_up.groupby(counts_up.index).mean()

    # Gene sets already uppercase
    result = gp.ssgsea(
        data=counts_up,
        gene_sets=GENE_SETS,
        min_size=2,
        outdir=None,
        sample_norm_method="rank",
        no_plot=True,
        verbose=False,
    )
    res = result.res2d
    # gseapy >= 1.0 returns long-format; older versions return wide (gene_sets x samples)
    if "Name" in res.columns and "Term" in res.columns:
        score_col = "NES" if "NES" in res.columns else "ES"
        scores = res.pivot(index="Name", columns="Term", values=score_col)
        scores.columns.name = None
    else:
        scores = res.T
    scores.index.name = "sample_id"
    print(f"  ssGSEA scores: {scores.shape[0]} samples x {scores.shape[1]} gene sets")
    return scores


def load_tcga_clinical():
    """Load TCGA clinical and map to barcode."""
    clin_path = os.path.join(RAW_DIR, "tcga_stad_clinical.csv")
    if not os.path.exists(clin_path):
        print("Clinical cache not found.")
        return None
    return pd.read_csv(clin_path)


def load_tcga_manifest():
    """Return UUID -> TCGA barcode mapping from cached manifest."""
    mf_path = os.path.join(RAW_DIR, "tcga_stad_file_manifest.json")
    if not os.path.exists(mf_path):
        return {}
    with open(mf_path) as f:
        hits = json.load(f)
    mapping = {}
    for h in hits:
        fid = h.get("file_id", "")
        cases = h.get("cases", [{}])
        sub = cases[0].get("submitter_id", "") if cases else ""
        if fid and sub:
            mapping[fid] = sub
    return mapping


def merge_scores_clinical(scores_df, clin_df, uuid_map):
    """Align ssGSEA scores with clinical OS data."""
    # Map sample IDs (file UUIDs or already barcodes)
    def to_barcode(sid):
        barcode = uuid_map.get(sid, sid)
        return barcode[:12] if barcode else sid

    scores_df = scores_df.copy()
    scores_df.index = [to_barcode(i) for i in scores_df.index]
    scores_df["sample_key"] = scores_df.index

    clin_df = clin_df.copy()
    clin_df["sample_key"] = clin_df["submitter_id"].str[:12]

    merged = scores_df.merge(
        clin_df[["sample_key", "OS_days", "os_event"]].dropna(subset=["OS_days"]),
        on="sample_key", how="inner"
    )
    print(f"  Merged: {len(merged)} samples with OS data")
    return merged


def km_plot(df, score_col, out_dir, min_events=5):
    """KM curves stratified by median ssGSEA score."""
    valid = df[[score_col, "OS_days", "os_event"]].dropna()
    if valid["os_event"].sum() < min_events:
        print(f"  SKIP KM {score_col}: too few events ({valid['os_event'].sum()})")
        return None

    med = valid[score_col].median()
    low  = valid[valid[score_col] <= med]
    high = valid[valid[score_col] >  med]

    res = logrank_test(low["OS_days"], high["OS_days"],
                       event_observed_A=low["os_event"],
                       event_observed_B=high["os_event"])
    p = res.p_value

    fig, ax = plt.subplots(figsize=(6, 5))
    kmf = KaplanMeierFitter()
    kmf.fit(low["OS_days"],  event_observed=low["os_event"],  label=f"Low {score_col} (n={len(low)})")
    kmf.plot_survival_function(ax=ax, ci_show=True, color="steelblue")
    kmf.fit(high["OS_days"], event_observed=high["os_event"], label=f"High {score_col} (n={len(high)})")
    kmf.plot_survival_function(ax=ax, ci_show=True, color="firebrick")

    sig = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "ns"))
    ax.set_title(f"TCGA-STAD: {score_col}\nlog-rank p={p:.4f} {sig}", fontsize=10)
    ax.set_xlabel("Days"); ax.set_ylabel("Overall Survival probability")
    plt.tight_layout()
    fname = f"km_tcga_{score_col}.png"
    plt.savefig(os.path.join(out_dir, fname), dpi=150)
    plt.close()
    print(f"  saved: {fname}  (p={p:.4f} {sig})")
    return p


def cox_analysis(merged):
    """Univariate + multivariate Cox on ssGSEA scores."""
    gs_cols = [c for c in GENE_SETS.keys() if c in merged.columns]
    rows = []

    # Univariate
    for col in gs_cols:
        df = merged[["OS_days", "os_event", col]].dropna()
        if df["os_event"].sum() < 5:
            continue
        # Normalize to 0-1
        rng = df[col].max() - df[col].min()
        if rng > 0:
            df[col] = (df[col] - df[col].min()) / rng
        cph = CoxPHFitter(penalizer=0.1)
        try:
            cph.fit(df, duration_col="OS_days", event_col="os_event")
            row = cph.summary.loc[col].to_dict()
            row["feature"] = col
            rows.append(row)
        except Exception as e:
            print(f"  Cox {col} failed: {e}")

    uni_df = pd.DataFrame(rows)
    if len(uni_df) > 0:
        uni_df = uni_df.sort_values("exp(coef)")
        uni_df.to_csv(os.path.join(OUT_DIR, "tcga_cox_univariate.csv"), index=False)

    # Multivariate: top features
    if len(gs_cols) >= 2:
        df_multi = merged[["OS_days", "os_event"] + gs_cols].dropna()
        for col in gs_cols:
            rng = df_multi[col].max() - df_multi[col].min()
            if rng > 0:
                df_multi[col] = (df_multi[col] - df_multi[col].min()) / rng
        cph_multi = CoxPHFitter(penalizer=0.1)
        try:
            cph_multi.fit(df_multi, duration_col="OS_days", event_col="os_event")
            cph_multi.summary.to_csv(os.path.join(OUT_DIR, "tcga_multivariate_cox.csv"))
            print(f"  Multivariate Cox: {len(gs_cols)} features, C-index={cph_multi.concordance_index_:.3f}")
            forest_plot(cph_multi.summary, "TCGA-STAD Multivariate Cox (OS)",
                        os.path.join(OUT_DIR, "tcga_cox_forest.png"))
        except Exception as e:
            print(f"  Multivariate Cox failed: {e}")

    return uni_df


def forest_plot(summary_df, title, outpath):
    df = summary_df.copy().sort_values("exp(coef)")
    fig, ax = plt.subplots(figsize=(8, max(4, len(df) * 0.45)))
    y = range(len(df))
    ax.errorbar(df["exp(coef)"], y,
                xerr=[df["exp(coef)"] - df["exp(coef) lower 95%"],
                      df["exp(coef) upper 95%"] - df["exp(coef)"]],
                fmt="o", color="navy", ecolor="gray", capsize=4, markersize=8)
    ax.axvline(1, color="red", linestyle="--", linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels(df.index, fontsize=9)
    ax.set_xlabel("Hazard Ratio (95% CI)", fontsize=11)
    ax.set_title(title, fontsize=11)
    for i, (idx, row) in enumerate(df.iterrows()):
        p_str = f"p={row['p']:.3f}" if row['p'] >= 0.001 else "p<0.001"
        star = " ***" if row['p'] < 0.001 else (" **" if row['p'] < 0.01 else (" *" if row['p'] < 0.05 else ""))
        ax.text(df["exp(coef) upper 95%"].max() * 1.05, i,
                f"{p_str}{star}", va="center", fontsize=8)
    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  saved: {os.path.basename(outpath)}")


def main():
    # Load + convert counts
    counts_df = load_and_convert_tcga_counts()
    if counts_df is None:
        return

    # ssGSEA
    scores_df = run_ssgsea(counts_df)
    scores_df.to_csv(os.path.join(OUT_DIR, "tcga_ssgsea_scores.csv"))
    print("  saved: tcga_ssgsea_scores.csv")

    # Merge with clinical
    clin_df  = load_tcga_clinical()
    uuid_map = load_tcga_manifest()
    if clin_df is None:
        print("No clinical data — exiting.")
        return
    merged = merge_scores_clinical(scores_df, clin_df, uuid_map)
    merged.to_csv(os.path.join(OUT_DIR, "tcga_ssgsea_clinical.csv"), index=False)

    # KM curves
    print("\n[KM curves]")
    for gs in GENE_SETS:
        km_plot(merged, gs, OUT_DIR)

    # Cox
    print("\n[Cox regression]")
    uni_df = cox_analysis(merged)
    if len(uni_df) > 0:
        print("\nUnivariate Cox summary:")
        print(uni_df[["feature", "exp(coef)", "exp(coef) lower 95%",
                       "exp(coef) upper 95%", "p"]].to_string(index=False))

    print(f"\n{'='*60}")
    print(f"TCGA ssGSEA Cox complete. Outputs in: {OUT_DIR}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
