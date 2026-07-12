"""
15_tme_predictor.py
Build a composite TME-based response predictor from the Korean cohort
and compare against PDL1 CPS.

Steps:
  1. Build features: cell type fractions + exhaustion score from patient_features.csv
  2. ROC analysis for response (PR/CR vs SD/PD) in Korean cohort
  3. C-index comparison: SPP1_score, exhaustion_score, composite, PDL1_CPS
  4. Apply composite signature to TCGA (using ssGSEA scores from script 13)
  5. Clinical subgroup analysis: MSI/HER2/EBV enrichment in SPP1-high tumors

Outputs in outputs/tme_predictor/
  roc_response.png
  cindex_bar.png
  subgroup_heatmap.png
  tcga_score_km.png
  predictor_summary.csv
"""
import os, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.utils import concordance_index

try:
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score, roc_curve
    SKLEARN = True
except ImportError:
    SKLEARN = False

warnings.filterwarnings("ignore")

BASE    = os.path.expanduser("~/gastric_tme_project")
OUT_DIR = os.path.join(BASE, "outputs/tme_predictor")
os.makedirs(OUT_DIR, exist_ok=True)


def load_korean_features():
    path = os.path.join(BASE, "outputs/survival/patient_features.csv")
    if not os.path.exists(path):
        print("patient_features.csv not found — run 08_survival_analysis.py first")
        return None
    pt = pd.read_csv(path)
    print(f"Korean cohort: {len(pt)} patients")
    print(f"  Columns: {pt.columns.tolist()}")
    return pt


def build_response_label(pt):
    """Binary response: PR/CR=1, SD/PD=0."""
    if "best_overall_response" not in pt.columns:
        return pt
    pt = pt.copy()
    pt["response_binary"] = pt["best_overall_response"].isin(["CR", "PR"]).astype(int)
    print(f"  Response: {pt['response_binary'].sum()} responders / {(~pt['response_binary'].astype(bool)).sum()} non-responders")
    return pt


def roc_analysis(pt, feature_cols, out_dir):
    """ROC curves for each feature predicting response."""
    if "response_binary" not in pt.columns:
        return
    y = pt["response_binary"]
    if y.sum() < 3 or (1 - y).sum() < 3:
        print("  Too few of one class for ROC analysis")
        return

    fig, ax = plt.subplots(figsize=(7, 6))
    aucs = {}
    colors = plt.cm.tab10(np.linspace(0, 1, len(feature_cols)))

    for col, color in zip(feature_cols, colors):
        if col not in pt.columns:
            continue
        valid = pt[[col, "response_binary"]].dropna()
        if len(valid) < 10:
            continue
        try:
            fpr, tpr, _ = roc_curve(valid["response_binary"], valid[col])
            auc = roc_auc_score(valid["response_binary"], valid[col])
            aucs[col] = auc
            ax.plot(fpr, tpr, color=color, linewidth=1.5,
                    label=f"{col} (AUC={auc:.2f})")
        except Exception:
            continue

    ax.plot([0, 1], [0, 1], "k--", linewidth=0.8, label="Random")
    ax.set_xlabel("1 - Specificity (FPR)", fontsize=11)
    ax.set_ylabel("Sensitivity (TPR)", fontsize=11)
    ax.set_title("ROC: TME features predicting immunotherapy response\n(Korean cohort)", fontsize=11)
    ax.legend(fontsize=7, loc="lower right")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "roc_response.png"), dpi=150)
    plt.close()
    print(f"  saved: roc_response.png")
    pd.Series(aucs).sort_values(ascending=False).to_csv(
        os.path.join(out_dir, "roc_aucs.csv"), header=["AUC"])
    return aucs


def build_composite_score(pt, feature_cols, response_col="response_binary"):
    """LASSO-based composite predictor or simple mean of top features."""
    valid = pt[feature_cols + [response_col]].dropna()
    X = valid[feature_cols].values
    y = valid[response_col].values

    if not SKLEARN or len(valid) < 10 or len(feature_cols) < 2:
        # Simple mean
        pt["composite_score"] = pt[feature_cols].mean(axis=1)
        return pt, feature_cols

    scaler = StandardScaler()
    X_sc = scaler.fit_transform(X)
    lr = LogisticRegression(penalty="l1", solver="liblinear", C=0.5, max_iter=500)
    try:
        lr.fit(X_sc, y)
        coefs = pd.Series(lr.coef_[0], index=feature_cols)
        print(f"\n  LASSO coefficients:\n{coefs[coefs != 0].sort_values()}")
        # Score all patients
        all_valid = pt[feature_cols].dropna()
        X_all = scaler.transform(all_valid.values)
        pt.loc[all_valid.index, "composite_score"] = lr.predict_proba(X_all)[:, 1]
        top_features = coefs[coefs != 0].index.tolist()
        coefs.to_csv(os.path.join(OUT_DIR, "predictor_coefficients.csv"), header=["coef"])
    except Exception as e:
        print(f"  LASSO failed: {e}, using mean score")
        pt["composite_score"] = pt[feature_cols].mean(axis=1)
        top_features = feature_cols

    return pt, top_features


def cindex_comparison(pt, duration_col, event_col, feature_cols, out_dir):
    """Bar chart comparing C-index of each feature + PDL1 CPS."""
    results = {}
    for col in feature_cols + ["PDL1_baseline_CPS", "composite_score"]:
        if col not in pt.columns:
            continue
        valid = pt[[duration_col, event_col, col]].dropna()
        if valid[event_col].sum() < 3:
            continue
        try:
            ci = concordance_index(valid[duration_col], -valid[col], valid[event_col])
            results[col] = ci
        except Exception:
            pass

    if not results:
        print("  No valid C-index comparisons")
        return

    fig, ax = plt.subplots(figsize=(max(6, len(results) * 1.2), 5))
    cols_sorted = sorted(results, key=results.get, reverse=True)
    vals = [results[c] for c in cols_sorted]
    colors_bar = ["gold" if c == "PDL1_baseline_CPS" else
                  ("forestgreen" if c == "composite_score" else "steelblue")
                  for c in cols_sorted]
    bars = ax.bar(range(len(cols_sorted)), vals, color=colors_bar)
    ax.axhline(0.5, color="gray", linestyle="--", label="Random (0.5)")
    ax.set_xticks(range(len(cols_sorted)))
    ax.set_xticklabels([c.replace("_", "\n") for c in cols_sorted], fontsize=8)
    ax.set_ylabel("C-index (concordance)"); ax.set_ylim(0, 1)
    ax.set_title("Predictive performance: TME features vs PDL1 CPS\n(Korean cohort, PFS)")
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width()/2, val + 0.01, f"{val:.3f}",
                ha="center", va="bottom", fontsize=9)
    ax.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "cindex_bar.png"), dpi=150)
    plt.close()
    print(f"  saved: cindex_bar.png")
    pd.Series(results, name="c_index").to_csv(os.path.join(out_dir, "cindex_results.csv"))
    return results


def subgroup_analysis(pt, score_col, out_dir):
    """Enrichment of molecular subtypes in score-high vs score-low patients."""
    subgroup_cols = [c for c in ["MSI_status", "HER2_status", "EBV_status", "TCGA_subtype",
                                 "progression_category", "best_overall_response"]
                     if c in pt.columns]
    if not subgroup_cols or score_col not in pt.columns:
        return

    med = pt[score_col].median()
    pt_sg = pt.copy()
    pt_sg["score_group"] = np.where(pt_sg[score_col] > med, "High", "Low")

    n_cols = len(subgroup_cols)
    fig, axes = plt.subplots(1, n_cols, figsize=(3.5 * n_cols, 4))
    if n_cols == 1:
        axes = [axes]

    for ax, col in zip(axes, subgroup_cols):
        ct = pd.crosstab(pt_sg["score_group"], pt_sg[col])
        ct_norm = ct.div(ct.sum(axis=1), axis=0)
        ct_norm.plot(kind="bar", ax=ax, stacked=True, colormap="tab10", legend=True)
        ax.set_title(f"{col}", fontsize=9)
        ax.set_xlabel(""); ax.set_ylabel("Fraction")
        ax.legend(fontsize=6, bbox_to_anchor=(1, 1))
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=0)

    plt.suptitle(f"{score_col}: High vs Low — molecular subgroup enrichment", fontsize=10, y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "subgroup_heatmap.png"), dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  saved: subgroup_heatmap.png")


def validate_in_tcga(composite_feature_cols, out_dir):
    """Apply TME signature to TCGA ssGSEA scores if available."""
    tcga_path = os.path.join(BASE, "outputs/tcga_validation/tcga_ssgsea_clinical.csv")
    if not os.path.exists(tcga_path):
        print("  TCGA ssGSEA scores not found — skipping TCGA validation")
        return

    tcga = pd.read_csv(tcga_path)
    avail = [c for c in composite_feature_cols if c in tcga.columns]
    if not avail:
        return

    tcga["composite_score_tcga"] = tcga[avail].mean(axis=1)
    valid = tcga[["OS_days", "os_event", "composite_score_tcga"]].dropna()
    if valid["os_event"].sum() < 5:
        return

    med = valid["composite_score_tcga"].median()
    low  = valid[valid["composite_score_tcga"] <= med]
    high = valid[valid["composite_score_tcga"] >  med]
    res  = logrank_test(low["OS_days"], high["OS_days"],
                        event_observed_A=low["os_event"],
                        event_observed_B=high["os_event"])
    p = res.p_value
    sig = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "ns"))

    fig, ax = plt.subplots(figsize=(6, 5))
    kmf = KaplanMeierFitter()
    kmf.fit(low["OS_days"],  event_observed=low["os_event"],  label=f"Low score (n={len(low)})")
    kmf.plot_survival_function(ax=ax, ci_show=True, color="steelblue")
    kmf.fit(high["OS_days"], event_observed=high["os_event"], label=f"High score (n={len(high)})")
    kmf.plot_survival_function(ax=ax, ci_show=True, color="firebrick")
    ax.set_title(f"TCGA-STAD: Composite TME score -> OS\np={p:.4f} {sig}", fontsize=10)
    ax.set_xlabel("Days"); ax.set_ylabel("Overall Survival probability")
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, "tcga_composite_score_km.png"), dpi=150)
    plt.close()
    print(f"  TCGA composite score KM: p={p:.4f} {sig}")
    print(f"  saved: tcga_composite_score_km.png")


def main():
    pt = load_korean_features()
    if pt is None:
        return

    pt = build_response_label(pt)

    # Feature columns for prediction
    frac_cols  = [c for c in pt.columns if c.startswith("frac_")]
    score_cols = [c for c in pt.columns if any(k in c for k in
                  ["exhaustion", "M1_M2", "score_Exhausted", "score_M1", "score_M2",
                   "score_Macrophage", "score_CD8", "score_CD4", "score_NK"])]
    feature_cols = frac_cols + score_cols
    print(f"\n  Using {len(feature_cols)} features: {feature_cols[:6]} ...")

    # ROC analysis
    print("\n[ROC Analysis]")
    aucs = roc_analysis(pt, feature_cols, OUT_DIR)

    # Build composite predictor
    print("\n[Composite Predictor]")
    pt, top_feats = build_composite_score(pt, feature_cols)

    # C-index comparison
    print("\n[C-index Comparison]")
    for dur, evt, label in [("PFS_days", "pfs_event", "PFS"), ("OS_days", "os_event", "OS")]:
        if dur in pt.columns and evt in pt.columns:
            print(f"\n  {label}:")
            cindex_comparison(pt, dur, evt, feature_cols + ["composite_score"], OUT_DIR)

    # Subgroup analysis
    print("\n[Subgroup Analysis]")
    subgroup_analysis(pt, "composite_score", OUT_DIR)

    # Summary table
    summary_rows = []
    if aucs:
        for k, v in aucs.items():
            summary_rows.append({"feature": k, "metric": "AUC", "value": v})
    pd.DataFrame(summary_rows).to_csv(os.path.join(OUT_DIR, "predictor_summary.csv"), index=False)

    # TCGA validation
    print("\n[TCGA Validation]")
    validate_in_tcga(["SPP1_Macro", "Exhausted_CD8", "Effector_CD8", "Treg"], OUT_DIR)

    print(f"\n{'='*60}")
    print(f"TME predictor complete. Outputs in: {OUT_DIR}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
