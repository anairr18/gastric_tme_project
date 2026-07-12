"""
08_survival_analysis.py
Kaplan-Meier survival analysis, Cox regression, and PDL1 comparison
on the Korean immunotherapy cohort (korea_kim2022_annotated_scored.h5ad).

Outputs in outputs/survival/:
  km_pfs_<cell_type>.png   - KM curves for PFS stratified by cell type abundance
  km_os_<cell_type>.png    - KM curves for OS
  cox_forest_pfs.png       - Forest plot of Cox HR (PFS)
  cox_forest_os.png        - Forest plot of Cox HR (OS)
  pdl1_vs_tme.png          - C-index comparison: PDL1 CPS vs TME model
  survival_summary.csv     - All HR, CI, p-values
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
from scipy import stats
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test
from lifelines.utils import concordance_index

warnings.filterwarnings("ignore")

BASE    = os.path.expanduser("~/gastric_tme_project")
PER_DIR = os.path.join(BASE, "data/processed/per_dataset")
OUT_DIR = os.path.join(BASE, "outputs/survival")
os.makedirs(OUT_DIR, exist_ok=True)


def build_patient_table(adata):
    """Aggregate per-cell data to per-patient level."""
    obs = adata.obs.copy()

    # Clinical columns expected
    clin_cols = ["patient", "PFS_days", "OS_days", "progression_s2",
                 "best_overall_response", "progression_category",
                 "HER2_status", "EBV_status", "MSI_status",
                 "TCGA_subtype", "PDL1_baseline_CPS"]
    clin_cols = [c for c in clin_cols if c in obs.columns]

    # Cell type fractions per patient
    ct_col = "cell_type_coarse"
    ct_fracs = (obs.groupby(["patient", ct_col]).size()
                   .unstack(ct_col, fill_value=0))
    ct_fracs = ct_fracs.div(ct_fracs.sum(axis=1), axis=0)
    ct_fracs.columns = [f"frac_{c.replace(' ', '_').replace('/', '_')}"
                        for c in ct_fracs.columns]

    # Mean exhaustion score per patient (all cells)
    score_cols = [c for c in obs.columns if c.startswith("score_") or
                  c in ("exhaustion_score", "M1_score", "M2_score", "M1_M2_ratio")]
    mean_scores = obs.groupby("patient")[score_cols].mean() if score_cols else pd.DataFrame()

    # Clinical table (one row per patient)
    clin = (obs[clin_cols].drop_duplicates("patient").set_index("patient")
            if "patient" in clin_cols else pd.DataFrame())

    pt = ct_fracs.join(mean_scores, how="left").join(clin, how="left")
    pt = pt.reset_index().rename(columns={"index": "patient"})

    # Event indicators (1=event occurred, 0=censored)
    # progression_s2 / progression_category store "Fast"/"Slow" string labels
    for prog_col in ["progression_s2", "progression_category"]:
        if prog_col in pt.columns and "pfs_event" not in pt.columns:
            raw = pt[prog_col].astype(str).str.strip()
            # "Fast" progression = event, "Slow" = censored
            # Also handle numeric/bool encodings
            pt["pfs_event"] = raw.map(
                {"Fast": 1, "Slow": 0, "True": 1, "False": 0, "1": 1, "0": 0, "1.0": 1, "0.0": 0}
            ).fillna(pd.to_numeric(raw, errors="coerce").fillna(0)).astype(int)

    if "best_overall_response" in pt.columns:
        # PD = progressive disease = event for OS proxy
        pt["os_event"] = pt["best_overall_response"].isin(["PD"]).astype(int)
        if pt["os_event"].sum() == 0:
            # All responders — use progression as OS proxy
            pt["os_event"] = pt.get("pfs_event", pd.Series(0, index=pt.index))

    return pt


def km_plot(pt, duration_col, event_col, strat_col, label, outpath, cut=None):
    """KM curve with log-rank test, stratified by median of strat_col."""
    valid = pt[[duration_col, event_col, strat_col]].dropna()
    if len(valid) < 10 or valid[event_col].sum() < 3:
        print(f"  SKIP {label}: too few events ({valid[event_col].sum()})")
        return None

    cut = cut or valid[strat_col].median()
    low  = valid[valid[strat_col] <= cut]
    high = valid[valid[strat_col] >  cut]
    if len(low) < 5 or len(high) < 5:
        return None

    result = logrank_test(low[duration_col], high[duration_col],
                          event_observed_A=low[event_col],
                          event_observed_B=high[event_col])
    p = result.p_value

    fig, ax = plt.subplots(figsize=(7, 5))
    kmf = KaplanMeierFitter()
    kmf.fit(low[duration_col],  event_observed=low[event_col],  label=f"Low  (n={len(low)})")
    kmf.plot_survival_function(ax=ax, ci_show=True, color="steelblue")
    kmf.fit(high[duration_col], event_observed=high[event_col], label=f"High (n={len(high)})")
    kmf.plot_survival_function(ax=ax, ci_show=True, color="firebrick")

    ax.set_title(f"{label}\nlog-rank p={p:.4f}", fontsize=11)
    ax.set_xlabel("Days"); ax.set_ylabel("Survival probability")
    ax.text(0.02, 0.05, f"p = {p:.4f}", transform=ax.transAxes, fontsize=10,
            color="black" if p > 0.05 else "red")
    plt.tight_layout()
    plt.savefig(outpath, dpi=150); plt.close()
    return p


def cox_analysis(pt, duration_col, event_col, feature_cols, covariate_cols, label):
    """Multivariate Cox PH regression. Returns summary DataFrame."""
    cols = [duration_col, event_col] + feature_cols + covariate_cols
    cols = [c for c in cols if c in pt.columns]
    df = pt[cols].dropna()
    if len(df) < 20 or df[event_col].sum() < 5:
        print(f"  SKIP Cox {label}: too few events")
        return None

    # Encode categoricals
    cat_cols = [c for c in covariate_cols if c in df.columns and df[c].dtype == object]
    for c in cat_cols:
        df[c] = pd.Categorical(df[c]).codes.astype(float)

    # Normalize continuous features to 0-1 for interpretable HRs
    for fc in feature_cols:
        if fc in df.columns:
            rng = df[fc].max() - df[fc].min()
            if rng > 0:
                df[fc] = (df[fc] - df[fc].min()) / rng

    cph = CoxPHFitter(penalizer=0.1)
    try:
        cph.fit(df, duration_col=duration_col, event_col=event_col)
    except Exception as e:
        print(f"  Cox fit failed: {e}")
        return None

    return cph.summary


def forest_plot(summary_df, title, outpath):
    """Forest plot from lifelines Cox summary."""
    df = summary_df.copy()
    df = df.sort_values("exp(coef)")
    fig, ax = plt.subplots(figsize=(8, max(4, len(df) * 0.4)))
    y = range(len(df))
    ax.errorbar(df["exp(coef)"], y,
                xerr=[df["exp(coef)"] - df["exp(coef) lower 95%"],
                      df["exp(coef) upper 95%"] - df["exp(coef)"]],
                fmt="o", color="navy", ecolor="gray", capsize=4)
    ax.axvline(1, color="red", linestyle="--", linewidth=1)
    ax.set_yticks(y)
    ax.set_yticklabels(df.index, fontsize=8)
    ax.set_xlabel("Hazard Ratio (95% CI)")
    ax.set_title(title)

    # Annotate p-values
    for i, (idx, row) in enumerate(df.iterrows()):
        p_str = f"p={row['p']:.3f}" if row['p'] >= 0.001 else "p<0.001"
        ax.text(df["exp(coef) upper 95%"].max() * 1.05, i, p_str, va="center", fontsize=7)

    plt.tight_layout()
    plt.savefig(outpath, dpi=150, bbox_inches="tight")
    plt.close()


def pdl1_comparison(pt, duration_col, event_col, feature_cols):
    """Compare c-index of TME model vs PDL1 CPS alone."""
    results = {}
    base = pt[[duration_col, event_col, "PDL1_baseline_CPS"]].dropna()
    if len(base) >= 20 and base[event_col].sum() >= 5:
        ci_pdl1 = concordance_index(base[duration_col], -base["PDL1_baseline_CPS"], base[event_col])
        results["PDL1_CPS"] = ci_pdl1

    tme_cols = [c for c in feature_cols if c in pt.columns]
    combined = pt[[duration_col, event_col] + tme_cols + ["PDL1_baseline_CPS"]].dropna()
    if len(combined) >= 20 and combined[event_col].sum() >= 5:
        from sklearn.linear_model import LogisticRegression
        X = combined[tme_cols].values
        y = combined[event_col].values
        # Simple linear combination score
        score = X.mean(axis=1)
        ci_tme = concordance_index(combined[duration_col], -score, combined[event_col])
        results["TME_model"] = ci_tme
        # Combined
        X_comb = np.column_stack([X, combined["PDL1_baseline_CPS"].values])
        score_comb = X_comb.mean(axis=1)
        ci_comb = concordance_index(combined[duration_col], -score_comb, combined[event_col])
        results["TME+PDL1"] = ci_comb

    if not results:
        return

    fig, ax = plt.subplots(figsize=(6, 4))
    bars = ax.bar(list(results.keys()), list(results.values()),
                  color=["steelblue", "firebrick", "forestgreen"][:len(results)])
    ax.axhline(0.5, color="gray", linestyle="--", label="Random (0.5)")
    ax.set_ylim(0, 1)
    ax.set_ylabel("C-index (concordance)")
    ax.set_title("Predictive performance: TME vs PDL1 CPS")
    for bar, val in zip(bars, results.values()):
        ax.text(bar.get_x() + bar.get_width() / 2, val + 0.01, f"{val:.3f}",
                ha="center", va="bottom", fontsize=10)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "pdl1_vs_tme.png"), dpi=150)
    plt.close()
    print(f"  saved: pdl1_vs_tme.png  {results}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", default=os.path.join(
        PER_DIR, "korea_kim2022_annotated_scored.h5ad"))
    args = parser.parse_args()

    print(f"Loading {args.input} ...")
    adata = sc.read_h5ad(args.input)
    print(f"  {adata.n_obs:,} cells | obs: {adata.obs.columns.tolist()}")

    # Restore full gene set if needed
    if adata.raw is not None and adata.raw.n_vars > adata.n_vars:
        obs_bak = adata.obs.copy()
        adata = adata.raw.to_adata()
        for col in obs_bak.columns:
            if col not in adata.obs.columns:
                adata.obs[col] = obs_bak[col].values
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    if "cell_type_coarse" not in adata.obs.columns:
        adata.obs["cell_type_coarse"] = "cluster_" + adata.obs["leiden"].astype(str)

    if "patient" not in adata.obs.columns:
        print("No patient column — cannot run survival analysis"); return

    print("Building per-patient table ...")
    pt = build_patient_table(adata)
    print(f"  {len(pt)} patients")
    pt.to_csv(os.path.join(OUT_DIR, "patient_features.csv"), index=False)

    frac_cols = [c for c in pt.columns if c.startswith("frac_")]
    score_cols = [c for c in pt.columns if "exhaustion" in c or "M1_M2" in c or "score_" in c]
    feature_cols = frac_cols + score_cols
    cov_cols = ["HER2_status", "EBV_status", "MSI_status", "TCGA_subtype", "PDL1_baseline_CPS"]
    cov_cols = [c for c in cov_cols if c in pt.columns]

    rows = []

    # ── KM curves ─────────────────────────────────────────────────────────────
    print("\n[Kaplan-Meier curves]")
    for fc in feature_cols:
        for endpoint, dur, evt in [("PFS", "PFS_days", "pfs_event"),
                                    ("OS",  "OS_days",  "os_event")]:
            if dur not in pt.columns or evt not in pt.columns:
                continue
            label = f"{fc} | {endpoint}"
            outpath = os.path.join(OUT_DIR, f"km_{endpoint}_{fc}.png")
            p = km_plot(pt, dur, evt, fc, label, outpath)
            if p is not None:
                rows.append({"feature": fc, "endpoint": endpoint, "km_logrank_p": p})
                sig = " ***" if p < 0.001 else (" *" if p < 0.05 else "")
                print(f"  {fc} | {endpoint}: p={p:.4f}{sig} -> saved")

    # ── Cox regression ────────────────────────────────────────────────────────
    print("\n[Cox regression — PFS]")
    if "PFS_days" in pt.columns and "pfs_event" in pt.columns:
        cox_pfs = cox_analysis(pt, "PFS_days", "pfs_event", feature_cols, cov_cols, "PFS")
        if cox_pfs is not None:
            forest_plot(cox_pfs, "Cox PH — PFS", os.path.join(OUT_DIR, "cox_forest_pfs.png"))
            cox_pfs.to_csv(os.path.join(OUT_DIR, "cox_pfs.csv"))
            print("  saved: cox_forest_pfs.png, cox_pfs.csv")

    print("\n[Cox regression — OS]")
    if "OS_days" in pt.columns and "os_event" in pt.columns:
        cox_os = cox_analysis(pt, "OS_days", "os_event", feature_cols, cov_cols, "OS")
        if cox_os is not None:
            forest_plot(cox_os, "Cox PH — OS", os.path.join(OUT_DIR, "cox_forest_os.png"))
            cox_os.to_csv(os.path.join(OUT_DIR, "cox_os.csv"))
            print("  saved: cox_forest_os.png, cox_os.csv")

    # ── PDL1 comparison ───────────────────────────────────────────────────────
    print("\n[PDL1 vs TME model comparison]")
    if "PDL1_baseline_CPS" in pt.columns:
        pdl1_comparison(pt, "PFS_days", "pfs_event", frac_cols)
    else:
        print("  no PDL1_baseline_CPS column — skipping")

    # ── Summary table ─────────────────────────────────────────────────────────
    if rows:
        summary = pd.DataFrame(rows).sort_values("km_logrank_p")
        summary.to_csv(os.path.join(OUT_DIR, "survival_summary.csv"), index=False)
        print(f"\nTop hits:")
        print(summary.head(10).to_string(index=False))

    print(f"\n{'='*60}")
    print(f"Survival analysis complete. Outputs in: {OUT_DIR}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
