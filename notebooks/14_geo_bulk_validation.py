"""
14_geo_bulk_validation.py
Download two large gastric cancer bulk RNA-seq cohorts from GEO and validate
our TME gene set findings.

Cohorts:
  GSE84437  - 433 gastric cancer patients, OS available
  GSE62254  - 300 patients (ACRG), DFS + molecular subtype

Pipeline per cohort:
  1. Download series matrix / expression file
  2. Parse expression matrix + clinical metadata
  3. ssGSEA with same gene sets as script 13
  4. KM curves + Cox regression
  5. Forest plot

Outputs in outputs/geo_validation/
"""
import os, sys, gzip, io, warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import requests, gseapy as gp
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test

warnings.filterwarnings("ignore")

BASE    = os.path.expanduser("~/gastric_tme_project")
RAW_DIR = os.path.join(BASE, "data/external/geo_bulk")
OUT_DIR = os.path.join(BASE, "outputs/geo_validation")
os.makedirs(RAW_DIR, exist_ok=True)
os.makedirs(OUT_DIR, exist_ok=True)

GENE_SETS = {
    "SPP1_Macro":       ["SPP1", "APOC1", "APOE", "C1QB", "TREM2", "GPNMB", "FN1"],
    "FOLR2_Macro":      ["FOLR2", "LYVE1", "SELENOP", "MRC1", "CD163", "STAB1"],
    "Exhausted_CD8":    ["PDCD1", "LAG3", "HAVCR2", "TIGIT", "TOX", "ENTPD1", "CTLA4", "NR4A1"],
    "Effector_CD8":     ["GZMB", "PRF1", "IFNG", "CX3CR1", "KLRG1", "GNLY"],
    "Treg":             ["FOXP3", "IL2RA", "IKZF2", "TNFRSF18"],
    "Inflammatory_M1":  ["IL1B", "TNF", "CXCL9", "CXCL10", "IDO1"],
    "Immunosuppressive":["IL10", "TGM2", "CCL18", "ARG1"],
    "CAF":              ["ACTA2", "FAP", "PDGFRB", "POSTN", "MMP2"],
}

GEO_DATASETS = {
    "GSE84437": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE84nnn/GSE84437/matrix/GSE84437_series_matrix.txt.gz",
        "gpl": "GPL6947",  # Illumina HumanHT-12 V3.0
        # Series matrix has: "death: 0/1" and "duration overall survival: N" (months)
        "os_col":    "duration overall survival",   # exact column name after parsing
        "event_col": "death",
        "os_unit":   "months",
        "label":     "GSE84437 (n=433)",
    },
    "GSE62254": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE62nnn/GSE62254/matrix/GSE62254_series_matrix.txt.gz",
        "gpl": "GPL570",   # Affymetrix HG-U133 Plus 2.0 (or similar)
        # No survival in series matrix — ssGSEA only, no KM curves
        "os_col":    None,
        "event_col": None,
        "os_unit":   "months",
        "label":     "GSE62254 ACRG (n=300)",
    },
    "GSE26253": {
        "url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE26nnn/GSE26253/matrix/GSE26253_series_matrix.txt.gz",
        "gpl": "GPL8432",  # Illumina human-6 v2.0
        # Series matrix has: "recurrence free survival time (month): N" and
        # "status (0=non-recurrence, 1=recurrence): 0/1" -- this is RFS, not OS
        "os_col":    "recurrence free survival time (month)",
        "event_col": "status (0=non-recurrence, 1=recurrence)",
        "os_unit":   "months",
        "label":     "GSE26253 (n=432, RFS)",
    },
}


def download_geo_matrix(gse_id, url):
    """Download GEO series matrix file, cache locally."""
    cache_path = os.path.join(RAW_DIR, f"{gse_id}_series_matrix.txt.gz")
    if not os.path.exists(cache_path):
        print(f"  Downloading {gse_id} ...")
        r = requests.get(url, timeout=300, stream=True)
        if r.status_code != 200:
            print(f"  HTTP {r.status_code} for {url}")
            return None
        with open(cache_path, "wb") as f:
            for chunk in r.iter_content(1024 * 1024):
                f.write(chunk)
        print(f"  Downloaded to {cache_path}")
    else:
        print(f"  Using cached {gse_id}")
    return cache_path


def parse_geo_matrix(path):
    """Parse GEO series matrix into (expr_df, clin_df).

    Series matrix format:
      !Sample_* lines: clinical metadata
      !series_matrix_table_begin: start of expression matrix
      Rows = probes, Columns = samples
    """
    clin_rows = {}
    characteristics_lines = []  # !Sample_characteristics_ch1 repeats once per key -- collect all
    expr_lines = []
    in_matrix = False
    sample_ids = None

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("!Sample_geo_accession"):
                sample_ids = line.split("\t")[1:]
                sample_ids = [s.strip('"') for s in sample_ids]
            elif line.startswith("!Sample_characteristics_ch1"):
                vals = [v.strip('"') for v in line.split("\t")[1:]]
                characteristics_lines.append(vals)
            elif line.startswith("!Sample_"):
                key = line.split("\t")[0].replace("!Sample_", "").lower()
                vals = [v.strip('"') for v in line.split("\t")[1:]]
                clin_rows[key] = vals
            elif line == "!series_matrix_table_begin":
                in_matrix = True
            elif line == "!series_matrix_table_end":
                in_matrix = False
            elif in_matrix:
                expr_lines.append(line)

    if not expr_lines or sample_ids is None:
        return None, None

    # Expression matrix
    header = expr_lines[0].split("\t")
    header = [h.strip('"') for h in header]
    data = []
    index = []
    for line in expr_lines[1:]:
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        index.append(parts[0].strip('"'))
        data.append([float(v) if v not in ("", "null", "NA", "NaN") else np.nan
                     for v in parts[1:]])

    max_cols = len(sample_ids)
    data = [row[:max_cols] for row in data]
    expr_df = pd.DataFrame(data, index=index, columns=sample_ids[:len(data[0])])
    expr_df = expr_df.dropna(how="all")

    # Clinical
    clin_df = pd.DataFrame(clin_rows, index=sample_ids[:len(list(clin_rows.values())[0])])

    # Expand !Sample_characteristics_ch1 -- GEO repeats this tag once per
    # characteristic key (e.g. one line for "death: 0/1", another for
    # "duration overall survival: N"), each cell as "key: value". A naive
    # parse that stores them all under one "characteristics_ch1" dict key
    # silently overwrites every line but the last, losing most clinical
    # fields. Split each line's "key: value" cells into proper columns.
    for vals in characteristics_lines:
        vals = vals[:len(clin_df)]
        keys_seen = set()
        for v in vals:
            if ": " in v:
                keys_seen.add(v.split(": ", 1)[0].strip().lower())
        if len(keys_seen) != 1:
            continue  # mixed/malformed line, skip rather than guess
        col_name = next(iter(keys_seen))
        col_vals = []
        for v in vals:
            col_vals.append(v.split(": ", 1)[1].strip() if ": " in v else np.nan)
        if col_name not in clin_df.columns:
            clin_df[col_name] = col_vals

    return expr_df, clin_df


def load_gpl_mapping(gpl_id):
    """Download GPL platform soft file and return probe_id -> symbol dict."""
    cache_tsv = os.path.join(RAW_DIR, f"{gpl_id}_mapping.tsv")
    if os.path.exists(cache_tsv):
        df = pd.read_csv(cache_tsv, sep="\t")
        m = dict(zip(df["probe_id"].astype(str), df["symbol"].astype(str)))
        # Also add ILMN_ prefixed versions
        m.update({f"ILMN_{k}": v for k, v in m.items() if not k.startswith("ILMN_")})
        print(f"  {gpl_id}: loaded {len(df):,} probe mappings from cache")
        return m

    num = int(gpl_id[3:])
    prefix = f"GPL{(num // 1000)}nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{prefix}/{gpl_id}/soft/{gpl_id}_family.soft.gz"
    print(f"  Downloading GPL annotation: {url}")
    try:
        r = requests.get(url, timeout=180, stream=True)
        r.raise_for_status()
        content = b"".join(r.iter_content(1024 * 1024))
    except Exception as e:
        print(f"  GPL download failed: {e}")
        return {}

    mapping = {}
    in_table = False
    headers = None
    id_col = sym_col = None

    with gzip.open(io.BytesIO(content)) as gz:
        for raw in gz:
            line = raw.decode("utf-8", errors="replace").rstrip("\n")
            if "platform_table_begin" in line.lower():
                in_table = True; headers = None; continue
            if "platform_table_end" in line.lower():
                in_table = False; continue
            if not in_table:
                continue
            parts = line.split("\t")
            if headers is None:
                headers = [h.strip().upper() for h in parts]
                for i, h in enumerate(headers):
                    if h == "ID": id_col = i
                for i, h in enumerate(headers):
                    if h == "SYMBOL": sym_col = i; break
                if sym_col is None:
                    for i, h in enumerate(headers):
                        if "SYMBOL" in h: sym_col = i; break
                if sym_col is None:
                    for i, h in enumerate(headers):
                        if h == "ILMN_GENE": sym_col = i; break
                continue
            if id_col is None or sym_col is None or len(parts) <= max(id_col, sym_col):
                continue
            pid = parts[id_col].strip()
            sym = parts[sym_col].strip()
            if pid and sym and sym not in ("", "---", "NA", "N/A", "0"):
                mapping[pid] = sym
                mapping[f"ILMN_{pid}"] = sym

    n = len([k for k in mapping if not k.startswith("ILMN_")])
    print(f"  {gpl_id}: {n:,} probes mapped")
    rows = [(k, v) for k, v in mapping.items() if not k.startswith("ILMN_")]
    pd.DataFrame(rows, columns=["probe_id", "symbol"]).to_csv(cache_tsv, sep="\t", index=False)
    return mapping


def probe_to_symbol(expr_df, gse_id, gpl_id=None):
    """Map probe IDs to uppercase gene symbols. Uses GPL soft file for ILMN_ probes."""
    idx = expr_df.index

    # Detect Illumina probe IDs (ILMN_XXXXXXX)
    illumina_frac = sum(1 for i in idx[:100] if str(i).startswith("ILMN_")) / min(100, len(idx))
    # Detect numeric probe IDs
    numeric_frac  = sum(1 for i in idx[:100] if str(i) and str(i)[0].isdigit()) / min(100, len(idx))

    if illumina_frac > 0.3 and gpl_id:
        # ILMN_ probes: use GPL annotation file
        gpl_map = load_gpl_mapping(gpl_id)
        if gpl_map:
            expr_df.index = [str(i) for i in expr_df.index]
            expr_df.index = [gpl_map.get(str(i), gpl_map.get(f"ILMN_{str(i).lstrip('ILMN_')}", str(i)))
                             for i in expr_df.index]
            expr_df.index = [str(g).upper() for g in expr_df.index]
            expr_df = expr_df.groupby(expr_df.index).mean()
            expr_df = expr_df[~expr_df.index.str.match(r'^(ILMN_|[0-9])')]
            print(f"  GPL mapping: {len(expr_df):,} gene symbols")
            return expr_df

    if numeric_frac > 0.3 and not illumina_frac > 0.3:
        # Numeric/Affy probes: try mygene with entrezgene scope
        print(f"  Mapping {len(idx):,} numeric probe IDs via mygene ...")
        try:
            import mygene
            mg = mygene.MyGeneInfo()
            mapping = {}
            ids = idx.tolist()
            for i in range(0, len(ids), 10000):
                batch = ids[i:i+10000]
                res = mg.querymany(batch, scopes="reporter,entrezgene", fields="symbol",
                                   species="human", as_dataframe=True, verbose=False)
                if "symbol" in res.columns:
                    mapping.update(res["symbol"].dropna().to_dict())
            print(f"  Mapped {len(mapping):,} probes to gene symbols")
            expr_df.index = [mapping.get(str(g), str(g)) for g in expr_df.index]
        except Exception as e:
            print(f"  mygene mapping failed: {e}")
    else:
        print(f"  Index looks like gene symbols ({len(idx):,} rows) — uppercasing")

    # Uppercase, collapse, remove unmapped probes
    expr_df.index = [str(g).upper() for g in expr_df.index]
    expr_df = expr_df.groupby(expr_df.index).mean()
    expr_df = expr_df[~expr_df.index.str.match(r'^(ILMN_|[0-9])')]
    print(f"  Final: {len(expr_df):,} gene symbols")
    return expr_df


def extract_survival(clin_df, cfg):
    """Extract OS days and event from clinical dataframe."""
    # Print available columns for debugging
    print(f"  Clinical columns: {list(clin_df.columns[:20])}")

    # Try exact column names first (configured per dataset)
    os_col  = cfg.get("os_col")
    evt_col = cfg.get("event_col")

    # If None, search by keywords
    if os_col is None and cfg.get("os_keywords"):
        cols_l = [c.lower() for c in clin_df.columns]
        for kw in cfg["os_keywords"]:
            for i, col_l in enumerate(cols_l):
                if kw.lower() in col_l:
                    os_col = clin_df.columns[i]; break
            if os_col: break

    if evt_col is None and cfg.get("event_keywords"):
        cols_l = [c.lower() for c in clin_df.columns]
        for kw in cfg["event_keywords"]:
            for i, col_l in enumerate(cols_l):
                if kw.lower() in col_l:
                    evt_col = clin_df.columns[i]; break
            if evt_col: break

    print(f"  Survival columns: time={os_col}, event={evt_col}")

    if os_col is None or os_col not in clin_df.columns:
        return None

    surv = pd.DataFrame(index=clin_df.index)
    surv["OS_months"] = pd.to_numeric(clin_df[os_col], errors="coerce")
    surv["OS_days"]   = surv["OS_months"] * 30.44

    if evt_col is not None and evt_col in clin_df.columns:
        evt = pd.to_numeric(clin_df[evt_col], errors="coerce")
        if cfg.get("os_invert", False):
            evt = 1 - evt
        surv["os_event"] = evt
    else:
        surv["os_event"] = np.nan

    surv = surv.dropna(subset=["OS_days"])
    return surv


def km_plot(df, score_col, title, out_dir, min_events=5):
    valid = df[[score_col, "OS_days", "os_event"]].dropna()
    if valid["os_event"].sum() < min_events:
        print(f"  SKIP KM {score_col}: {valid['os_event'].sum()} events")
        return None

    med = valid[score_col].median()
    low  = valid[valid[score_col] <= med]
    high = valid[valid[score_col] >  med]

    res = logrank_test(low["OS_days"], high["OS_days"],
                       event_observed_A=low["os_event"],
                       event_observed_B=high["os_event"])
    p = res.p_value
    sig = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else "ns"))

    fig, ax = plt.subplots(figsize=(6, 5))
    kmf = KaplanMeierFitter()
    kmf.fit(low["OS_days"],  event_observed=low["os_event"],  label=f"Low (n={len(low)})")
    kmf.plot_survival_function(ax=ax, ci_show=True, color="steelblue")
    kmf.fit(high["OS_days"], event_observed=high["os_event"], label=f"High (n={len(high)})")
    kmf.plot_survival_function(ax=ax, ci_show=True, color="firebrick")
    ax.set_title(f"{title}\n{score_col}  p={p:.4f} {sig}", fontsize=9)
    ax.set_xlabel("Days"); ax.set_ylabel("Survival probability")
    plt.tight_layout()
    fname = f"km_{title.replace(' ', '_').replace('(', '').replace(')', '')}_{score_col}.png"
    plt.savefig(os.path.join(out_dir, fname), dpi=150)
    plt.close()
    print(f"    {score_col}: p={p:.4f} {sig}")
    return p


def run_geo_dataset(gse_id, cfg):
    print(f"\n{'='*60}")
    print(f"Processing {gse_id}: {cfg['label']}")
    print(f"{'='*60}")

    # Download
    path = download_geo_matrix(gse_id, cfg["url"])
    if path is None:
        return

    # Parse
    expr_df, clin_df = parse_geo_matrix(path)
    if expr_df is None:
        print(f"  Could not parse {gse_id}")
        return
    print(f"  Expression: {expr_df.shape[0]} probes x {expr_df.shape[1]} samples")

    # Map to gene symbols (pass GPL ID for proper Illumina probe annotation)
    expr_df = probe_to_symbol(expr_df, gse_id, gpl_id=cfg.get("gpl"))

    # Cache parsed expression
    expr_cache = os.path.join(RAW_DIR, f"{gse_id}_expr_symbols.parquet")
    expr_df.to_parquet(expr_cache)

    # ssGSEA — both expression index and gene sets must be uppercase
    expr_df.index = expr_df.index.str.upper()
    expr_df = expr_df.groupby(expr_df.index).mean()
    gene_sets_upper = {k: [g.upper() for g in v] for k, v in GENE_SETS.items()}
    print(f"  Running ssGSEA on {expr_df.shape[0]} genes x {expr_df.shape[1]} samples ...")
    try:
        result = gp.ssgsea(data=expr_df, gene_sets=gene_sets_upper, min_size=2,
                           outdir=None, sample_norm_method="rank",
                           no_plot=True, verbose=False)
        res = result.res2d
        # gseapy >= 1.0 returns long-format; older versions return wide (gene_sets x samples)
        if "Name" in res.columns and "Term" in res.columns:
            score_col = "NES" if "NES" in res.columns else "ES"
            scores = res.pivot(index="Name", columns="Term", values=score_col)
            scores.columns.name = None
        else:
            scores = res.T
        print(f"  ssGSEA: {scores.shape[0]} samples x {scores.shape[1]} gene sets")
    except Exception as e:
        print(f"  ssGSEA failed: {e}")
        return

    # Extract survival
    surv = extract_survival(clin_df, cfg)
    if surv is None:
        print(f"  No survival data found for {gse_id}")
        return
    print(f"  Survival: {len(surv)} patients, {surv['os_event'].sum():.0f} events")

    # Merge
    merged = scores.join(surv, how="inner")
    merged.to_csv(os.path.join(OUT_DIR, f"{gse_id}_ssgsea_clinical.csv"))

    # KM curves
    print(f"  KM curves:")
    results = {}
    for gs in GENE_SETS:
        if gs in merged.columns:
            p = km_plot(merged, gs, cfg["label"], OUT_DIR)
            if p is not None:
                results[gs] = p

    # Univariate Cox
    cox_rows = []
    for gs in GENE_SETS:
        df = merged[["OS_days", "os_event", gs]].dropna()
        if df["os_event"].sum() < 5:
            continue
        rng = df[gs].max() - df[gs].min()
        if rng > 0:
            df[gs] = (df[gs] - df[gs].min()) / rng
        try:
            cph = CoxPHFitter(penalizer=0.1)
            cph.fit(df, duration_col="OS_days", event_col="os_event")
            row = cph.summary.loc[gs].to_dict()
            row["feature"] = gs
            cox_rows.append(row)
        except Exception:
            pass

    if cox_rows:
        cox_df = pd.DataFrame(cox_rows).sort_values("exp(coef)")
        cox_df.to_csv(os.path.join(OUT_DIR, f"{gse_id}_cox.csv"), index=False)
        print(f"\n  Cox summary ({gse_id}):")
        print(cox_df[["feature", "exp(coef)", "p"]].to_string(index=False))

        # Forest plot
        fig, ax = plt.subplots(figsize=(8, max(4, len(cox_df) * 0.45)))
        y = range(len(cox_df))
        ax.errorbar(cox_df["exp(coef)"], y,
                    xerr=[cox_df["exp(coef)"] - cox_df["exp(coef) lower 95%"],
                          cox_df["exp(coef) upper 95%"] - cox_df["exp(coef)"]],
                    fmt="o", color="navy", ecolor="gray", capsize=4)
        ax.axvline(1, color="red", linestyle="--")
        ax.set_yticks(y); ax.set_yticklabels(cox_df["feature"], fontsize=9)
        ax.set_xlabel("Hazard Ratio (95% CI)")
        ax.set_title(f"{cfg['label']}: TME signatures -> OS")
        for i, row in cox_df.iterrows():
            p_str = f"p<0.001" if row['p'] < 0.001 else f"p={row['p']:.3f}"
            ax.text(cox_df["exp(coef) upper 95%"].max() * 1.05, list(y)[list(cox_df.index).index(i)],
                    p_str, va="center", fontsize=7)
        plt.tight_layout()
        plt.savefig(os.path.join(OUT_DIR, f"{gse_id}_cox_forest.png"), dpi=150, bbox_inches="tight")
        plt.close()
        print(f"  saved: {gse_id}_cox_forest.png")


def main():
    for gse_id, cfg in GEO_DATASETS.items():
        try:
            run_geo_dataset(gse_id, cfg)
        except Exception as e:
            print(f"  {gse_id} failed: {e}")

    print(f"\n{'='*60}")
    print(f"GEO validation complete. Outputs in: {OUT_DIR}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
