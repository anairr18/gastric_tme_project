"""
02_dataset_manifest.py
Canonical registry of all gastric cancer scRNA-seq datasets targeted for meta-analysis.
Run this to inspect the manifest, validate GEO accessions, and print a summary table.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", line_buffering=True)
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", line_buffering=True)

import os
import pandas as pd

BASE = os.path.expanduser("~/gastric_tme_project")

# ─── Dataset registry ───────────────────────────────────────────────────────
# Each entry describes one publicly available gastric cancer scRNA-seq dataset.
# Fields:
#   dataset_id   : short unique key used as file prefix throughout the pipeline
#   geo_accession: GEO/ArrayExpress accession (None = already local)
#   paper        : citation shorthand
#   journal      : publication venue
#   year         : publication year
#   n_cells_raw  : approximate raw cell count (pre-QC)
#   n_patients   : number of patients/donors
#   tissue       : tissue types profiled
#   has_treatment: whether the cohort has immunotherapy/treatment data
#   has_timepoints: longitudinal design
#   platform     : sequencing platform
#   species      : organism
#   notes        : any special considerations for download / preprocessing
DATASETS = [
    {
        "dataset_id":    "korea_kim2022",
        "geo_accession": None,
        "local_h5ad":    os.path.join(BASE, "data/raw/ge_korea_raw_data_count_matricies_raw_combined.h5ad"),
        "paper":         "Kim et al. 2022",
        "journal":       "unpublished / supplementary",
        "year":          2022,
        "n_cells_raw":   654770,
        "n_patients":    33,
        "tissue":        "tumor, adjacent normal, distal normal",
        "has_treatment": True,
        "has_timepoints":True,
        "platform":      "10x Chromium",
        "species":       "Human",
        "notes":         "Korean immunotherapy cohort; baseline + FU1 + FU2; clinical metadata in Table S1/S2. PRIMARY DATASET.",
    },
    {
        "dataset_id":    "kumar2022",
        "geo_accession": "GSE183904",
        "local_h5ad":    None,
        "paper":         "Kumar et al. 2022",
        "journal":       "Cancer Discovery",
        "year":          2022,
        "n_cells_raw":   200000,
        "n_patients":    31,
        "tissue":        "tumor, normal",
        "has_treatment": False,
        "has_timepoints":False,
        "platform":      "10x Chromium",
        "species":       "Human",
        "notes":         "Comprehensive atlas; 48 samples, 34 lineage states including rare populations. Paired with spatial transcriptomics.",
    },
    {
        "dataset_id":    "sathe2020",
        "geo_accession": "GSE150290",
        "local_h5ad":    None,
        "paper":         "Sathe et al. 2020",
        "journal":       "Clinical Cancer Research",
        "year":          2020,
        "n_cells_raw":   56167,
        "n_patients":    8,
        "tissue":        "tumor, normal, PBMC",
        "has_treatment": False,
        "has_timepoints":False,
        "platform":      "10x Chromium",
        "species":       "Human",
        "notes":         "7 GC patients + 1 IM patient; receptor-ligand TME network; includes paired PBMCs.",
    },
    {
        "dataset_id":    "zhang2021",
        "geo_accession": "GSE134520",
        "local_h5ad":    None,
        "paper":         "Zhang et al. 2019/2021",
        "journal":       "Gut",
        "year":          2021,
        "n_cells_raw":   30000,
        "n_patients":    12,
        "tissue":        "NAG, CAG, IM, early GC",
        "has_treatment": False,
        "has_timepoints":False,
        "platform":      "10x Chromium",
        "species":       "Human",
        "notes":         "Premalignant lesion series (NAG→CAG→IM→EGC). Key for understanding gastric carcinogenesis trajectory.",
    },
    {
        "dataset_id":    "diffuse_gc_2021",
        "geo_accession": "GSE167297",
        "local_h5ad":    None,
        "paper":         "Clin Cancer Res 2021",
        "journal":       "Clinical Cancer Research",
        "year":          2021,
        "n_cells_raw":   22464,
        "n_patients":    5,
        "tissue":        "tumor (superficial + deep layers), normal",
        "has_treatment": False,
        "has_timepoints":False,
        "platform":      "10x Chromium",
        "species":       "Human",
        "notes":         "Diffuse-type GC; spatial TME reprogramming by tumor invasion depth.",
    },
    {
        "dataset_id":    "tcell_exhaustion_2022",
        "geo_accession": "GSE206785",
        "local_h5ad":    None,
        "paper":         "Nat Commun 2022",
        "journal":       "Nature Communications",
        "year":          2022,
        "n_cells_raw":   166533,
        "n_patients":    10,
        "tissue":        "tumor, paratumor, blood",
        "has_treatment": False,
        "has_timepoints":False,
        "platform":      "10x Chromium",
        "species":       "Human",
        "notes":         "Alternative T cell exhaustion trajectory; ligand-receptor intercellular interaction network.",
    },
    {
        "dataset_id":    "helicobacter_2024",
        "geo_accession": "GSE212212",
        "local_h5ad":    None,
        "paper":         "eLife 2024",
        "journal":       "eLife",
        "year":          2024,
        "n_cells_raw":   83637,
        "n_patients":    None,
        "tissue":        "gastric tissue (healthy, H.pylori-, H.pylori+)",
        "has_treatment": False,
        "has_timepoints":False,
        "platform":      "10x Chromium",
        "species":       "Human",
        "notes":         "H. pylori-associated tumorigenesis atlas; healthy controls + GC with/without H.pylori.",
    },
]

# ─── Build manifest dataframe ────────────────────────────────────────────────
manifest = pd.DataFrame(DATASETS)
manifest_path = os.path.join(BASE, "data/external/dataset_manifest.csv")
os.makedirs(os.path.dirname(manifest_path), exist_ok=True)
manifest.to_csv(manifest_path, index=False)

print("=" * 80)
print("GASTRIC TME META-ANALYSIS — DATASET MANIFEST")
print("=" * 80)

display_cols = ["dataset_id", "geo_accession", "paper", "year", "n_cells_raw", "n_patients",
                "has_treatment", "has_timepoints"]
print(manifest[display_cols].to_string(index=False))

print(f"\nTotal datasets: {len(manifest)}")
print(f"Total cells (raw, approx): {manifest['n_cells_raw'].sum():,}")
print(f"Datasets with treatment data: {manifest['has_treatment'].sum()}")
print(f"Datasets with timepoints: {manifest['has_timepoints'].sum()}")
print(f"\nManifest saved to: {manifest_path}")

# ─── Validate local files ────────────────────────────────────────────────────
print("\n--- Local file status ---")
for _, row in manifest.iterrows():
    if row["local_h5ad"]:
        exists = os.path.exists(row["local_h5ad"])
        size_gb = os.path.getsize(row["local_h5ad"]) / 1e9 if exists else 0
        status = f"OK ({size_gb:.1f} GB)" if exists else "MISSING"
        print(f"  {row['dataset_id']:<25}  {status}")
    else:
        ext_dir = os.path.join(BASE, "data/external", row["dataset_id"])
        exists = os.path.isdir(ext_dir) and bool(os.listdir(ext_dir)) if os.path.isdir(ext_dir) else False
        status = "downloaded" if exists else "not yet downloaded"
        print(f"  {row['dataset_id']:<25}  [{row['geo_accession']}] {status}")

print("\nNext step: run 03_dataset_download.py to acquire GEO datasets.")
