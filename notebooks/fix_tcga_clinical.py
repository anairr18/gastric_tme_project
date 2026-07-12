"""
fix_tcga_clinical.py
Re-fetch TCGA-STAD clinical data with correct GDC API field paths.
vital_status lives under demographic, not diagnoses.
"""
import os, json, requests
import numpy as np
import pandas as pd

BASE    = os.path.expanduser("~/gastric_tme_project")
RAW_DIR = os.path.join(BASE, "data/external/tcga_stad")
GDC_CASES = "https://api.gdc.cancer.gov/cases"

cache = os.path.join(RAW_DIR, "tcga_stad_clinical.csv")

print("Re-fetching TCGA-STAD clinical data with correct fields ...")
filters = {
    "op": "=",
    "content": {"field": "project.project_id", "value": "TCGA-STAD"}
}
params = {
    "filters": json.dumps(filters),
    "fields": (
        "submitter_id,case_id,"
        "demographic.vital_status,"
        "demographic.days_to_death,"       # days_to_death is in demographic, not diagnoses
        "diagnoses.days_to_last_follow_up"
    ),
    "format": "JSON",
    "size": "600"
}
r = requests.get(GDC_CASES, params=params, timeout=120)
r.raise_for_status()
hits = r.json()["data"]["hits"]
print(f"  Got {len(hits)} cases")

rows = []
for h in hits:
    demo  = h.get("demographic", {}) or {}
    diag  = (h.get("diagnoses") or [{}])[0]
    vital = demo.get("vital_status", "")
    dtd   = demo.get("days_to_death", None)          # demographic.days_to_death
    dtlf  = diag.get("days_to_last_follow_up", None)
    rows.append({
        "case_id":      h.get("case_id"),
        "submitter_id": h.get("submitter_id"),
        "vital_status": vital,
        "days_to_death":           dtd,
        "days_to_last_follow_up":  dtlf,
    })

df = pd.DataFrame(rows)
dtd_num  = pd.to_numeric(df["days_to_death"], errors="coerce")
dtlf_num = pd.to_numeric(df["days_to_last_follow_up"], errors="coerce")
# OS = days_to_death for dead; days_to_last_follow_up for alive; fallback to the other
df["OS_days"] = np.where(df["vital_status"] == "Dead", dtd_num, dtlf_num)
# If dead but days_to_death missing, use days_to_last_follow_up as minimum OS
mask_dead_no_dtd = (df["vital_status"] == "Dead") & df["OS_days"].isna()
df.loc[mask_dead_no_dtd, "OS_days"] = dtlf_num[mask_dead_no_dtd]
df["os_event"] = (df["vital_status"] == "Dead").astype(int)

print(f"  vital_status distribution:\n{df['vital_status'].value_counts()}")
print(f"  os_event: {df['os_event'].sum()} deaths / {len(df)} total")
print(f"  OS_days: {df['OS_days'].notna().sum()} non-null")

df.to_csv(cache, index=False)
print(f"  saved: {cache}")
