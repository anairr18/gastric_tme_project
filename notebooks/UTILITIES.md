# Utility Scripts & One-Off Fixes

This document describes utility scripts that support the main analysis pipeline. These are **not** part of the core data analysis (01-18_*.py) but may be needed for data preparation or troubleshooting.

---

## Core Utilities

### 1. `fix_tcga_clinical.py`
**Purpose:** Correct GDC API field paths for TCGA-STAD clinical metadata retrieval.

**Context:** TCGA data access via GDC portal changes version-to-version. This script corrects field paths that may have shifted (e.g., `diagnoses.pathology_report_uuid` vs older formats).

**When to use:** If running `notebooks/13_tcga_ssgsea_cox.py` fails with "field not found" errors related to TCGA clinical data.

**What it does:**
- Fetches TCGA-STAD sample list from GDC
- Maps correct field paths for: diagnosis date, tumor grade, histology, stage, therapy
- Writes corrected `data/external/tcga_stad/tcga_stad_clinical.csv`

**Run:**
```bash
python fix_tcga_clinical.py
```

---

### 2. `fix_tcga_deconv_survival.py`
**Purpose:** Recalculate NNLS cell-type deconvolution with fixed clinical merging for TCGA-STAD survival analysis.

**Context:** NNLS deconvolution of TCGA bulk samples can fail if clinical metadata isn't correctly merged post-deconvolution.

**When to use:** If `notebooks/13_tcga_ssgsea_cox.py` produces inconsistent deconvolution results or missing survival data.

**What it does:**
- Runs NNLS deconvolution on TCGA-STAD gene counts (reference = Korean cohort signatures)
- Merges deconvolved fractions with clinical metadata
- Performs Kaplan-Meier + Cox regression on fibroblast fraction
- Outputs: `data/external/tcga_stad/tcga_nnls_deconv_with_survival.csv`

**Run:**
```bash
python fix_tcga_deconv_survival.py
```

---

### 3. `download_gpl_mappings.py`
**Purpose:** Fetch and cache GPL platform probe-to-gene mappings for GEO datasets.

**Context:** GEO datasets are often stored as probe IDs (e.g., GPL6947 for Illumina arrays). This script maps probes to genes for downstream analysis.

**When to use:** Automatically called by `notebooks/03_dataset_download.py` or `14_geo_bulk_validation.py` if mappings don't exist. Rarely needs manual execution.

**What it does:**
- For each GPL platform in dataset manifest (GPL6947, GPL8432, GPL20795, etc.):
  - Downloads `.annot.gz` from NCBI GEO (if available)
  - Parses probe ID → gene symbol mapping
  - Caches in `data/external/gpl_mappings.pkl` for reuse
- Falls back to direct family file parsing if `.annot.gz` unavailable (slower)

**Run:**
```bash
python download_gpl_mappings.py
```

---

## Supporting Notebooks (Non-Pipeline)

These are implemented but **not** part of the core 7-step Phase 1 or 6-step Phase 2 workflows. They generate additional analysis outputs but are optional for main conclusions.

### `08_survival_analysis.py`
Single-cohort Kaplan-Meier curves for Korean cohort, stratified by:
- TME composition (high vs low CAF, SPP1+ myeloid, exhausted T cells)
- Clinical features (MSI status, HER2, EBV, PDL1 CPS)
- Progression category

**Output:** `outputs/survival/km_*.png`, `cox_forest_*.csv`

---

### `09_tcga_validation.py`
TCGA-STAD Cox regression (standalone, predates integrated ssGSEA). Kept for legacy comparisons.

**Output:** `outputs/tcga_validation/cox_*.csv`, `forest_*.png`

---

### `10_subclustering.py`
High-resolution clustering of specific cell types:
- Macrophage subsets: SPP1+, APOE+, inflammatory, etc.
- T cell subsets: exhausted (PD1+, TIM3+, LAG3+), proliferating, naive

**Output:** `outputs/subclustering/macrophage_sumap.png`, `tcell_subtypes.csv`

---

### `11_cell_communication.py`
Cell-cell interaction analysis via LIANA (aggregates CellChat, Connectome, Tensor-cell2cell).
- Queries: SPP1+ macrophage → T cell, T cell → myeloid, etc.
- Outputs interaction networks (fast vs slow progressors)

**Output:** `outputs/cell_communication/liana_*.csv`, `*_dotplot.png`

---

### `12_trajectory.py`
T cell developmental trajectory (PAGA + pseudotime):
- Naive → effector → exhausted (PD1+, TIM3+ pathway)
- RNA velocity (optional, if splicing data available)

**Output:** `outputs/trajectory/paga_*.png`, `pseudotime_*.csv`

---

### `16_pathway_enrichment.py`
Gene ontology (GO) + KEGG enrichment on DEGs (differential expression by progression/cell type).
- Ranked by log2FC + adjusted p-value
- Outputs: heatmaps of top pathways per cluster

**Output:** `outputs/pathway_enrichment/go_heatmap.png`, `kegg_*.png`

---

## Troubleshooting Guide

### Issue: `notebooks/13_tcga_ssgsea_cox.py` fails with "field not found"
**Solution:** Run `fix_tcga_clinical.py`, then retry 13_*.py

### Issue: NNLS deconvolution has NaN survival data
**Solution:** Run `fix_tcga_deconv_survival.py` to re-merge clinical metadata correctly

### Issue: GEO bulk validation (14_*.py) times out downloading GPL mappings
**Solution:** Run `download_gpl_mappings.py` separately, then retry 14_*.py

### Issue: Insufficient memory during integration
**Solution:** `notebooks/05_integration.py --batch-size 32 --n-latent 15` (reduces memory footprint, slightly slower convergence)

---

## Pipeline Execution Order

### Recommended: Core pipeline (7 steps)
```
01_preprocessing.py
  → 02_dataset_manifest.py
  → 03_dataset_download.py (auto-calls download_gpl_mappings.py)
  → 04_standardized_qc.py
  → 05_integration.py
  → 06_cell_type_annotation.py
  → 07_meta_analysis.py
```

### Optional: Validation + figures
```
13_tcga_ssgsea_cox.py (may call fix_tcga_clinical.py if needed)
14_geo_bulk_validation.py
15_tme_predictor.py
17_figure_assembly.py
18_cross_dataset_replication.py
```

### Optional: Detailed analysis
```
08_survival_analysis.py (single-cohort KM)
09_tcga_validation.py (legacy TCGA)
10_subclustering.py
11_cell_communication.py
12_trajectory.py
16_pathway_enrichment.py
```

---

## Notes for Future Work

1. **GPT-4 mediated cell type refinement:** Consider fine-tuning CellTypist with external reference if custom cell type definitions needed.
2. **Splicing analysis:** 12_trajectory.py can incorporate RNA velocity if raw counts available.
3. **Spatial validation:** If spatial transcriptomics (Visium/MERFISH) data becomes available, create 19_spatial_validation.py.
4. **Patient-level immune profiling:** Aggregate cell-type fractions by patient_ID for prognostic biomarker development.

---

*Last updated: 2026-06-24*
