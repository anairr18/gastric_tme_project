"""
09_tcga_validation.py
Download TCGA STAD bulk RNA-seq, deconvolve cell type fractions using our
scRNA-seq reference (NNLS), then replicate Korean cohort survival findings.

Steps:
  1. Build pseudo-bulk reference from korea_kim2022_annotated_scored.h5ad
  2. Download TCGA-STAD HTSeq counts + clinical data via GDC API
  3. NNLS deconvolution of each TCGA sample
  4. Kaplan-Meier + Cox survival in TCGA using deconvolved fractions

Outputs in outputs/tcga_validation/:
  tcga_deconv_fractions.csv
  km_tcga_pfs_<cell_type>.png
  cox_tcga_os.png / .csv
  validation_comparison.png
"""
import sys, os, json, time, io
import argparse
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
import requests
from scipy.optimize import nnls
from scipy import stats
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

warnings.filterwarnings("ignore")

BASE    = os.path.expanduser("~/gastric_tme_project")
PER_DIR = os.path.join(BASE, "data/processed/per_dataset")
OUT_DIR = os.path.join(BASE, "outputs/tcga_validation")
RAW_DIR = os.path.join(BASE, "data/external/tcga_stad")
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(RAW_DIR, exist_ok=True)

GDC_FILES  = "https://api.gdc.cancer.gov/files"
GDC_DATA   = "https://api.gdc.cancer.gov/data"
GDC_CASES  = "https://api.gdc.cancer.gov/cases"


def build_reference(adata, ct_col="cell_type_coarse", n_genes=5000):
    """Build genes × cell_types pseudo-bulk reference (mean log-normalized expression)."""
    print("Building scRNA-seq pseudo-bulk reference ...")

    # Use full gene set if raw available
    if adata.raw is not None and adata.raw.n_vars > adata.n_vars:
        obs_bak = adata.obs.copy()
        adata = adata.raw.to_adata()
        for col in obs_bak.columns:
            if col not in adata.obs.columns:
                adata.obs[col] = obs_bak[col].values
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

    cell_types = adata.obs[ct_col].unique()
    ref = {}
    for ct in cell_types:
        mask = adata.obs[ct_col] == ct
        if mask.sum() < 10:
            continue
        sub = adata[mask]
        import scipy.sparse as sp
        if sp.issparse(sub.X):
            # avoid materializing the dense matrix; sparse mean is memory-efficient
            ref[ct] = np.asarray(sub.X.mean(axis=0)).flatten()
        else:
            ref[ct] = sub.X.mean(axis=0)

    ref_df = pd.DataFrame(ref, index=adata.var_names)

    # Select top variable genes across cell types
    ref_std = ref_df.std(axis=1)
    top_genes = ref_std.nlargest(n_genes).index
    ref_df = ref_df.loc[top_genes]

    print(f"  Reference: {ref_df.shape[0]} genes × {ref_df.shape[1]} cell types")
    ref_df.to_csv(os.path.join(OUT_DIR, "scrnaseq_reference.csv"))
    return ref_df


def fetch_tcga_stad_files():
    """Query GDC API for TCGA-STAD HTSeq-Counts files."""
    cache = os.path.join(RAW_DIR, "tcga_stad_file_manifest.json")
    if os.path.exists(cache):
        print("  Loading cached file manifest ...")
        with open(cache) as f:
            return json.load(f)

    print("Querying GDC API for TCGA-STAD HTSeq-Counts ...")
    filters = {
        "op": "and",
        "content": [
            {"op": "=", "content": {"field": "cases.project.project_id", "value": "TCGA-STAD"}},
            {"op": "=", "content": {"field": "data_type", "value": "Gene Expression Quantification"}},
            {"op": "=", "content": {"field": "analysis.workflow_type", "value": "STAR - Counts"}},
        ]
    }
    params = {
        "filters": json.dumps(filters),
        "fields": "file_id,file_name,cases.case_id,cases.submitter_id",
        "format": "JSON",
        "size": "600"
    }
    r = requests.get(GDC_FILES, params=params, timeout=60)
    r.raise_for_status()
    hits = r.json()["data"]["hits"]
    print(f"  Found {len(hits)} files")
    with open(cache, "w") as f:
        json.dump(hits, f)
    return hits


def download_tcga_counts(hits, max_files=400):
    """Batch-download TCGA counts files via GDC data endpoint."""
    counts_cache = os.path.join(RAW_DIR, "tcga_stad_counts.csv")
    if os.path.exists(counts_cache):
        print("  Loading cached TCGA counts ...")
        return pd.read_csv(counts_cache, index_col=0)

    print(f"Downloading {min(len(hits), max_files)} TCGA-STAD count files ...")
    file_ids = [h["file_id"] for h in hits[:max_files]]

    # Download in batches of 50
    all_counts = {}
    batch_size = 50
    for i in range(0, len(file_ids), batch_size):
        batch = file_ids[i:i + batch_size]
        payload = {"ids": batch}
        r = requests.post(GDC_DATA, json=payload, timeout=300)
        if r.status_code != 200:
            print(f"  Batch {i//batch_size} failed: {r.status_code}")
            continue

        # Parse tar.gz response
        import tarfile
        try:
            with tarfile.open(fileobj=io.BytesIO(r.content)) as tar:
                for member in tar.getmembers():
                    if member.name.endswith(".tsv") or member.name.endswith(".counts"):
                        f = tar.extractfile(member)
                        if f is None:
                            continue
                        try:
                            df = pd.read_csv(f, sep="\t", comment="#",
                                             header=None, names=["gene_id", "count"],
                                             usecols=[0, 3])
                            # Filter PAR_Y and keep ENSG
                            df = df[df["gene_id"].str.startswith("ENSG")]
                            sample_id = member.name.split("/")[0]
                            all_counts[sample_id] = df.set_index("gene_id")["count"]
                        except Exception:
                            continue
        except Exception as e:
            print(f"  Parse error batch {i//batch_size}: {e}")
            continue

        print(f"  Downloaded batch {i//batch_size+1}/{(len(file_ids)-1)//batch_size+1} "
              f"({len(all_counts)} samples so far)")
        time.sleep(0.5)

    if not all_counts:
        print("  No counts downloaded.")
        return None

    counts_df = pd.DataFrame(all_counts).fillna(0)
    counts_df.to_csv(counts_cache)
    print(f"  TCGA counts: {counts_df.shape[0]} genes × {counts_df.shape[1]} samples")
    return counts_df


def fetch_tcga_clinical():
    """Fetch TCGA-STAD clinical data via GDC cases endpoint."""
    cache = os.path.join(RAW_DIR, "tcga_stad_clinical.csv")
    if os.path.exists(cache):
        return pd.read_csv(cache)

    print("Fetching TCGA-STAD clinical data ...")
    filters = {
        "op": "=",
        "content": {"field": "project.project_id", "value": "TCGA-STAD"}
    }
    params = {
        "filters": json.dumps(filters),
        "fields": ("submitter_id,case_id,"
                   "diagnoses.days_to_death,diagnoses.vital_status,"
                   "diagnoses.days_to_last_follow_up"),
        "format": "JSON",
        "size": "600"
    }
    r = requests.get(GDC_CASES, params=params, timeout=60)
    r.raise_for_status()
    hits = r.json()["data"]["hits"]

    rows = []
    for h in hits:
        diag = h.get("diagnoses", [{}])[0]
        rows.append({
            "case_id": h.get("case_id"),
            "submitter_id": h.get("submitter_id"),
            "vital_status": diag.get("vital_status"),
            "days_to_death": diag.get("days_to_death"),
            "days_to_last_follow_up": diag.get("days_to_last_follow_up"),
        })

    df = pd.DataFrame(rows)
    df["OS_days"] = np.where(df["vital_status"] == "Dead",
                             df["days_to_death"],
                             df["days_to_last_follow_up"])
    df["os_event"] = (df["vital_status"] == "Dead").astype(int)
    df.to_csv(cache, index=False)
    print(f"  Clinical: {len(df)} cases")
    return df


def convert_ensembl_to_symbol(counts_df):
    """Map ENSEMBL IDs to gene symbols using mygene."""
    # Always strip ENSEMBL version numbers first
    stripped = counts_df.index.str.split(".").str[0]
    counts_df.index = stripped

    try:
        import mygene
        mg = mygene.MyGeneInfo()
        ids = counts_df.index.tolist()
        print(f"  Querying mygene for {len(ids):,} ENSEMBL IDs ...")
        result = mg.querymany(ids, scopes="ensembl.gene", fields="symbol",
                              species="human", as_dataframe=True, verbose=False)
        # mygene may return multi-column DF; extract symbol robustly
        if "symbol" in result.columns:
            mapping = result["symbol"].dropna().to_dict()
        else:
            # Fall back to notfound-aware dict
            mapping = {}
            for qid, row in result.iterrows():
                sym = row.get("symbol", None) if isinstance(row, dict) else getattr(row, "symbol", None)
                if sym and str(sym) not in ("nan", "None"):
                    mapping[qid] = sym
        print(f"  Mapped {len(mapping):,} genes to symbols")
        counts_df.index = [mapping.get(g, g) for g in counts_df.index]
        counts_df = counts_df.groupby(counts_df.index).mean()
    except Exception as e:
        print(f"  mygene mapping failed ({e}), keeping ENSEMBL IDs")
    return counts_df


def deconvolve(counts_df, ref_df):
    """NNLS deconvolution of bulk samples using scRNA-seq reference."""
    print("Running NNLS deconvolution ...")
    # Find common genes
    common = ref_df.index.intersection(counts_df.index)
    if len(common) < 100:
        print(f"  Only {len(common)} common genes — deconvolution may be unreliable")
    if len(common) < 20:
        return None

    ref = ref_df.loc[common].values.astype(float)
    bulk = counts_df.loc[common].values.astype(float)

    # Log-normalize bulk counts
    bulk = np.log1p(bulk / (bulk.sum(axis=0, keepdims=True) + 1) * 1e4)

    fracs = np.zeros((bulk.shape[1], ref.shape[1]))
    for i in range(bulk.shape[1]):
        coef, _ = nnls(ref, bulk[:, i])
        total = coef.sum()
        fracs[i] = coef / total if total > 0 else coef

    frac_df = pd.DataFrame(fracs, index=counts_df.columns,
                           columns=ref_df.columns)
    print(f"  Deconvolved {len(frac_df)} samples × {len(ref_df.columns)} cell types")
    return frac_df


def survival_in_tcga(frac_df, clin_df, hit_features):
    """Replicate top Korean cohort findings in TCGA."""
    # Match samples to clinical by submitter ID
    clin_df = clin_df.dropna(subset=["OS_days"])
    clin_df["sample_key"] = clin_df["submitter_id"].str[:12]

    frac_df.index = [i[:15] for i in frac_df.index]
    frac_df["sample_key"] = frac_df.index.str[:12]

    merged = frac_df.merge(clin_df[["sample_key", "OS_days", "os_event"]],
                           on="sample_key", how="inner")
    print(f"  Matched {len(merged)} TCGA samples with clinical data")
    if len(merged) < 20:
        print("  Too few matched samples — skipping TCGA survival")
        return

    ct_cols = [c for c in frac_df.columns if c != "sample_key"]
    fig, axes = plt.subplots(1, min(len(ct_cols), 4),
                             figsize=(4 * min(len(ct_cols), 4), 5))
    if len(ct_cols) == 1:
        axes = [axes]

    for ax, ct in zip(axes, ct_cols[:4]):
        med = merged[ct].median()
        low  = merged[merged[ct] <= med]
        high = merged[merged[ct] >  med]
        if len(low) < 5 or len(high) < 5:
            ax.set_visible(False)
            continue
        result = logrank_test(low["OS_days"], high["OS_days"],
                              event_observed_A=low["os_event"],
                              event_observed_B=high["os_event"])
        kmf = KaplanMeierFitter()
        kmf.fit(low["OS_days"],  event_observed=low["os_event"],  label=f"Low {ct}")
        kmf.plot_survival_function(ax=ax, ci_show=False, color="steelblue")
        kmf.fit(high["OS_days"], event_observed=high["os_event"], label=f"High {ct}")
        kmf.plot_survival_function(ax=ax, ci_show=False, color="firebrick")
        ax.set_title(f"{ct}\np={result.p_value:.3f}", fontsize=9)

    plt.suptitle("TCGA-STAD Validation: OS by deconvolved cell type fraction", y=1.02)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "tcga_km_os_validation.png"), dpi=150, bbox_inches="tight")
    plt.close()
    print("  saved: tcga_km_os_validation.png")

    merged.to_csv(os.path.join(OUT_DIR, "tcga_deconv_clinical.csv"), index=False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", default=os.path.join(
        PER_DIR, "korea_kim2022_annotated_scored.h5ad"))
    parser.add_argument("--max-samples", type=int, default=400)
    args = parser.parse_args()

    # Step 1: Build reference
    ref_cache = os.path.join(OUT_DIR, "scrnaseq_reference.csv")
    if os.path.exists(ref_cache):
        print("Loading cached scRNA-seq reference ...")
        ref_df = pd.read_csv(ref_cache, index_col=0)
    else:
        adata = sc.read_h5ad(args.ref)
        if "cell_type_coarse" not in adata.obs.columns:
            adata.obs["cell_type_coarse"] = "cluster_" + adata.obs["leiden"].astype(str)
        ref_df = build_reference(adata)

    # Step 2: Download TCGA
    hits = fetch_tcga_stad_files()
    counts_df = download_tcga_counts(hits, max_files=args.max_samples)
    if counts_df is None:
        print("Could not download TCGA counts. Exiting.")
        return

    counts_df = convert_ensembl_to_symbol(counts_df)

    # Step 3: Deconvolve
    frac_df = deconvolve(counts_df, ref_df)
    if frac_df is None:
        print("Deconvolution failed — insufficient gene overlap.")
        return
    frac_df.to_csv(os.path.join(OUT_DIR, "tcga_deconv_fractions.csv"))
    print("  saved: tcga_deconv_fractions.csv")

    # Step 4: Build UUID → TCGA barcode mapping from manifest
    uuid_to_barcode = {}
    for h in hits:
        fid = h.get("file_id", "")
        cases = h.get("cases", [{}])
        submitter = cases[0].get("submitter_id", "") if cases else ""
        if fid and submitter:
            uuid_to_barcode[fid] = submitter

    # Rename frac_df index: UUID sample keys → TCGA barcodes
    frac_df.index = [uuid_to_barcode.get(idx, idx) for idx in frac_df.index]

    # Step 5: Clinical survival
    clin_df = fetch_tcga_clinical()
    survival_in_tcga(frac_df, clin_df, hit_features=[])

    print(f"\n{'='*60}")
    print(f"TCGA validation complete. Outputs in: {OUT_DIR}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
