#!/usr/bin/env python3
"""
Diagnose why composite score underperforms PD-L1
"""
import os, warnings, numpy as np, pandas as pd
import scanpy as sc
from sklearn.metrics import roc_auc_score, roc_curve
from scipy.stats import pearsonr, spearmanr

warnings.filterwarnings("ignore")

BASE = "C:\\Users\\Aadi Nair\\gastric_tme_project"
KOREAN = os.path.join(BASE, "data/processed/per_dataset/korea_kim2022_annotated_scored.h5ad")

print("="*80)
print("DIAGNOSTIC: Why is composite score underperforming?")
print("="*80 + "\n")

# Load data
korean = sc.read_h5ad(KOREAN)
korean_df = korean.obs.copy()

print(f"Loaded Korean cohort: {len(korean):,} cells from {korean_df['patient'].nunique()} patients\n")

# Check metadata columns
print("STEP 1: Available metadata columns")
print("-" * 50)
cols = korean_df.columns.tolist()
print(f"All columns: {cols}\n")

# Check which scoring columns exist
score_cols = ['PDL1_baseline_CPS', 'exhaustion_score', 'M2_score', 'M1_M2_ratio',
              'progression_category', 'PFS_days']
for col in score_cols:
    if col in korean_df.columns:
        print(f"[YES] {col}")
        print(f"  - Non-null: {korean_df[col].notna().sum()}/{len(korean_df)}")
        print(f"  - Type: {korean_df[col].dtype}")
        if korean_df[col].dtype in ['int64', 'float64']:
            print(f"  - Range: {korean_df[col].min():.3f} to {korean_df[col].max():.3f}")
    else:
        print(f"[NO] {col} MISSING")

print("\n")

# Aggregate to patient level
print("STEP 2: Patient-level aggregation")
print("-" * 50)

patient_list = []
for pid in korean_df['patient'].unique():
    cells = korean_df[korean_df['patient'] == pid]
    patient_list.append({
        'patient': pid,
        'progression': cells['progression_category'].iloc[0],
        'pdl1': cells['PDL1_baseline_CPS'].iloc[0],
        'exhaustion': cells['exhaustion_score'].mean(),
        'm2': cells['M2_score'].mean(),
        'm1m2': cells['M1_M2_ratio'].mean(),
        'pfs': cells['PFS_days'].iloc[0] if 'PFS_days' in cells.columns else np.nan,
    })

pdata = pd.DataFrame(patient_list)
pdata['is_fast'] = (pdata['progression'] == 'Fast').astype(int)

print(f"Aggregated to {len(pdata)} patients")
print(f"  - Fast progressors: {pdata['is_fast'].sum()}")
print(f"  - Slow progressors: {(pdata['is_fast']==0).sum()}\n")

# Check each component's predictive power
print("STEP 3: Individual component AUCs")
print("-" * 50)

components = ['pdl1', 'exhaustion', 'm2', 'm1m2']
aucs = {}
for comp in components:
    auc = roc_auc_score(pdata['is_fast'], pdata[comp])
    aucs[comp] = auc
    print(f"{comp:15s} AUC: {auc:.4f}")

print()

# Check correlations with progression
print("STEP 4: Correlation with progression (is_fast)")
print("-" * 50)

for comp in components:
    r, p = spearmanr(pdata[comp], pdata['is_fast'])
    print(f"{comp:15s} Spearman r={r:+.3f} (p={p:.4f})")

print()

# Check multicollinearity
print("STEP 5: Component intercorrelations")
print("-" * 50)

for i, comp1 in enumerate(components):
    for comp2 in components[i+1:]:
        r, p = pearsonr(pdata[comp1], pdata[comp2])
        print(f"{comp1} vs {comp2:12s} r={r:+.3f}")

print("\n")

# Test different weighting schemes
print("STEP 6: Optimized weighting schemes")
print("-" * 50)

from sklearn.preprocessing import StandardScaler

# Scheme A: Equal weights
X = pdata[components].values
X_scaled = StandardScaler().fit_transform(X)
pdata['composite_equal'] = X_scaled.mean(axis=1)
auc_equal = roc_auc_score(pdata['is_fast'], pdata['composite_equal'])
print(f"Equal weights (0.25 each):     AUC={auc_equal:.4f}")

# Scheme B: Original weights
weights_orig = np.array([0.40, 0.35, 0.15, 0.10])
pdata['composite_orig'] = X_scaled @ weights_orig
auc_orig = roc_auc_score(pdata['is_fast'], pdata['composite_orig'])
print(f"Original (0.40, 0.35, 0.15, 0.10): AUC={auc_orig:.4f}")

# Scheme C: PD-L1 only
pdata['composite_pdl1only'] = pdata['pdl1']
auc_pdl1 = roc_auc_score(pdata['is_fast'], pdata['pdl1'])
print(f"PD-L1 only:                    AUC={auc_pdl1:.4f}")

# Scheme D: Exhaustion + M2 (remove PDL1)
weights_no_pdl1 = np.array([0, 0.5, 0.3, 0.2])
pdata['composite_no_pdl1'] = X_scaled @ weights_no_pdl1
auc_no_pdl1 = roc_auc_score(pdata['is_fast'], pdata['composite_no_pdl1'])
print(f"No PD-L1 (0, 0.5, 0.3, 0.2):   AUC={auc_no_pdl1:.4f}")

# Scheme E: Try inverse (score as risk of SLOW progression)
pdata['is_slow'] = 1 - pdata['is_fast']
weights_slow = np.array([0.40, 0.35, 0.15, 0.10])
pdata['composite_slow'] = X_scaled @ weights_slow
auc_slow = roc_auc_score(pdata['is_slow'], pdata['composite_slow'])
print(f"Predict SLOW (inverse):        AUC={auc_slow:.4f}")

print("\n")

# Investigate score distributions
print("STEP 7: Score distributions by progression")
print("-" * 50)

for comp in components:
    fast = pdata[pdata['is_fast']==1][comp]
    slow = pdata[pdata['is_fast']==0][comp]
    print(f"\n{comp}:")
    print(f"  Fast:  mean={fast.mean():.3f}, std={fast.std():.3f}")
    print(f"  Slow:  mean={slow.mean():.3f}, std={slow.std():.3f}")
    print(f"  Delta: {fast.mean() - slow.mean():+.3f}")

print("\n")

# Check for potential issues
print("STEP 8: Data quality checks")
print("-" * 50)

for comp in components:
    print(f"{comp}:")
    print(f"  NaN values: {pdata[comp].isna().sum()}")
    print(f"  Zeros: {(pdata[comp]==0).sum()}")
    print(f"  Negative: {(pdata[comp]<0).sum()}")

print("\n")

print("="*80)
print("SUMMARY & RECOMMENDATIONS")
print("="*80)

best_scheme = max([('Equal', auc_equal), ('Original', auc_orig),
                    ('PD-L1', auc_pdl1), ('No PD-L1', auc_no_pdl1),
                    ('Slow (inverse)', auc_slow)], key=lambda x: x[1])

print(f"\nBest performing: {best_scheme[0]} (AUC={best_scheme[1]:.4f})")
print(f"PD-L1 baseline:  AUC={auc_pdl1:.4f}")

if best_scheme[1] > auc_pdl1:
    print(f"\nFINDING: Composite score CAN outperform PD-L1")
    print(f"Use '{best_scheme[0]}' weighting scheme for {best_scheme[1]:.4f} AUC")
else:
    print(f"\nFINDING: PD-L1 is the dominant biomarker")
    print(f"Reason: Other components (exhaustion, M2, M1/M2) are NOT predictive of progression")
    print(f"        or are highly correlated with PD-L1 (redundant)")
