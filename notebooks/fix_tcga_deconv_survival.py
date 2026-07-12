"""
fix_tcga_deconv_survival.py
Re-merge NNLS deconvolution fractions with corrected clinical (os_event now has real events)
and regenerate KM survival plots.
"""
import os, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

warnings.filterwarnings("ignore")

BASE    = os.path.expanduser("~/gastric_tme_project")
RAW_DIR = os.path.join(BASE, "data/external/tcga_stad")
OUT_DIR = os.path.join(BASE, "outputs/tcga_validation")

# Load deconv fractions
frac_path = os.path.join(OUT_DIR, "tcga_deconv_fractions.csv")
frac_df = pd.read_csv(frac_path, index_col=0)
print(f"Deconv fractions: {frac_df.shape}")

# Load corrected clinical
clin_df = pd.read_csv(os.path.join(RAW_DIR, "tcga_stad_clinical.csv"))
print(f"Clinical: {len(clin_df)} cases, {clin_df['os_event'].sum()} deaths")

# Need UUID -> barcode mapping
import json
mf_path = os.path.join(RAW_DIR, "tcga_stad_file_manifest.json")
uuid_map = {}
if os.path.exists(mf_path):
    with open(mf_path) as f:
        hits = json.load(f)
    for h in hits:
        fid = h.get("file_id", "")
        cases = h.get("cases", [{}])
        sub = cases[0].get("submitter_id", "") if cases else ""
        if fid and sub:
            uuid_map[fid] = sub
    print(f"UUID map: {len(uuid_map)} entries")

# Convert frac_df index: UUID -> barcode
frac_df.index = [uuid_map.get(i, i) for i in frac_df.index]
frac_df["sample_key"] = frac_df.index.str[:12]

clin_df["sample_key"] = clin_df["submitter_id"].str[:12]
clin_valid = clin_df.dropna(subset=["OS_days"])[["sample_key", "OS_days", "os_event"]]

merged = frac_df.merge(clin_valid, on="sample_key", how="inner")
print(f"Merged: {len(merged)} samples, {merged['os_event'].sum():.0f} events")
merged.to_csv(os.path.join(OUT_DIR, "tcga_deconv_clinical.csv"), index=False)
print("  saved: tcga_deconv_clinical.csv")

# KM plots for each cell type fraction
ct_cols = [c for c in frac_df.columns if c != "sample_key"]
cox_rows = []
for ct in ct_cols:
    df = merged[["OS_days", "os_event", ct]].dropna()
    if df["os_event"].sum() < 5:
        continue
    med = df[ct].median()
    low  = df[df[ct] <= med]
    high = df[df[ct] >  med]
    res  = logrank_test(low["OS_days"], high["OS_days"],
                        event_observed_A=low["os_event"],
                        event_observed_B=high["os_event"])
    p = res.p_value
    sig = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "ns"))

    fig, ax = plt.subplots(figsize=(5, 4))
    kmf = KaplanMeierFitter()
    kmf.fit(low["OS_days"],  event_observed=low["os_event"],  label=f"Low {ct} (n={len(low)})")
    kmf.plot_survival_function(ax=ax, ci_show=False, color="steelblue")
    kmf.fit(high["OS_days"], event_observed=high["os_event"], label=f"High {ct} (n={len(high)})")
    kmf.plot_survival_function(ax=ax, ci_show=False, color="firebrick")
    ax.set_title(f"TCGA-STAD: {ct} fraction -> OS\np={p:.4f} {sig}", fontsize=9)
    ax.set_xlabel("Days"); ax.set_ylabel("OS probability")
    plt.tight_layout()
    safe_ct = ct.replace(' ', '_').replace('/', '_')
    fname = f"km_tcga_deconv_{safe_ct}.png"
    plt.savefig(os.path.join(OUT_DIR, fname), dpi=150)
    plt.close()
    print(f"  {ct}: p={p:.4f} {sig}")

    # Cox
    rng = df[ct].max() - df[ct].min()
    if rng > 0:
        df[ct] = (df[ct] - df[ct].min()) / rng
    try:
        cph = CoxPHFitter(penalizer=0.1)
        cph.fit(df, duration_col="OS_days", event_col="os_event")
        row = cph.summary.loc[ct].to_dict()
        row["feature"] = ct
        cox_rows.append(row)
    except Exception:
        pass

if cox_rows:
    cox_df = pd.DataFrame(cox_rows).sort_values("exp(coef)")
    cox_df.to_csv(os.path.join(OUT_DIR, "tcga_deconv_cox.csv"), index=False)
    print("\nCox summary (NNLS deconv):")
    print(cox_df[["feature", "exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]].to_string(index=False))

    # Updated forest plot
    fig, ax = plt.subplots(figsize=(8, max(4, len(cox_df) * 0.5)))
    y = range(len(cox_df))
    ax.errorbar(cox_df["exp(coef)"], y,
                xerr=[cox_df["exp(coef)"] - cox_df["exp(coef) lower 95%"],
                      cox_df["exp(coef) upper 95%"] - cox_df["exp(coef)"]],
                fmt="o", color="navy", ecolor="gray", capsize=4, markersize=8)
    ax.axvline(1, color="red", linestyle="--")
    ax.set_yticks(y); ax.set_yticklabels(cox_df["feature"], fontsize=9)
    ax.set_xlabel("Hazard Ratio (95% CI)")
    ax.set_title("TCGA-STAD: Cell type fraction Cox OS (NNLS deconvolution)")
    for i, (_, row) in enumerate(cox_df.iterrows()):
        p_str = "p<0.001" if row["p"] < 0.001 else f"p={row['p']:.3f}"
        ax.text(cox_df["exp(coef) upper 95%"].max() * 1.05, i, p_str, va="center", fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "tcga_deconv_cox_forest.png"), dpi=150, bbox_inches="tight")
    plt.close()
    print("  saved: tcga_deconv_cox_forest.png")

print("\nDone.")
