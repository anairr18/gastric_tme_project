# Gastric Cancer Tumor Microenvironment: scRNA-seq Analysis
### Advaith Nair | Korean Immunotherapy Cohort | Progress Report

---

## Slide 1 — Project Overview & Motivation

**Central question:** What features of the gastric cancer tumor microenvironment (TME) predict immunotherapy response, and how does TME composition evolve across treatment?

**Why this matters:**
- Gastric cancer is the 5th most common cancer worldwide; ~40% of patients receive immune checkpoint blockade (ICB)
- Response rates remain low (~15–20%), with no reliable TME biomarker for patient selection
- Longitudinal single-cell profiling offers a unique window into treatment-driven immune dynamics

**Dataset:**
- Source: Korean cohort (publicly available supplementary data)
- 33 patients receiving anti-PD-1 / combination immunotherapy
- 138 single-cell libraries (tumor, adjacent normal, distal normal) across 3 timepoints
- 654,770 raw cells profiled across 36,601 genes

---

## Slide 2 — Cohort Design

| Feature | Detail |
|---------|--------|
| Patients | 33 |
| Timepoints | Baseline → Follow-up 1 (FU1) → Follow-up 2 (FU2) |
| Sample types | Tumor (T), Adjacent Normal (AN), Distal Normal (DN) |
| Libraries | 138 total |
| Clinical outcomes | Best overall response (PD / SD / PR), PFS (days), OS (days) |
| Molecular subtypes | HER2, EBV, MSI, TCGA subtype (GS / CIN / MSI-High / EBV / Indeterminate) |
| PDL1 biomarker | PDL1 CPS score at baseline |
| Progression label | Slow Progressor (0) vs Fast Progressor (1) |

**Key clinical axis for downstream analysis:** Slow vs. Fast Progressor — binary label capturing differential treatment benefit, linked to PFS/OS outcomes.

---

## Slide 3 — Preprocessing Pipeline Overview

> Full pipeline implemented in `notebooks/01_preprocessing.py` (reproducible, checkpoint-aware)

```
Raw h5ad (654,770 cells × 36,601 genes)
        ↓
Step 1: Load + inspect
        ↓
Step 2: QC metrics (MT%, genes/cell, UMI/cell)
        ↓
Step 3: QC plots (violin, scatter)
        ↓
Step 4: Cell & gene filtering
        ↓
Step 5: Doublet detection (Scrublet, per-sample)
        ↓  ← checkpoint saved here
Step 6: Normalize (10k) → log1p → 3,000 HVGs
        ↓
Step 7: Scale → PCA (50 PCs) → kNN → UMAP → Leiden
        ↓
Step 8: Save processed AnnData
        ↓
Step 9: Merge clinical metadata (Tables S1 + S2)
        ↓
Final: gastric_processed.h5ad (429,867 cells × 3,000 HVGs)
```

---

## Slide 4 — QC Results: What We Kept and Why

### Pre-filtering snapshot
| Metric | Median | Mean |
|--------|--------|------|
| Genes per cell | 844 | 1,067 |
| UMI counts per cell | 2,217 | 3,998 |
| MT% | 7.88% | 20.53% |

The gap between mean and median MT% signals a population of low-quality / apoptotic cells inflating the average — exactly what the MT% filter is designed to remove.

### Cells removed
| Filter | Cells removed | Reason |
|--------|-------------|--------|
| < 200 genes | 92,664 | Empty droplets / ambient RNA |
| > 6,000 genes | 3,068 | Likely doublets (two cells captured as one) |
| MT% > 20% | 129,171 | Dying / apoptotic cells |
| **Total** | **224,903 (34.3%)** | |

### Genes removed
| Filter | Genes removed |
|--------|-------------|
| < 3 cells expressing | 3,910 (10.7%) |

### Final object
- **429,867 cells × 32,691 genes** (pre-HVG)
- **429,867 cells × 3,000 HVGs** (model-ready)

---

## Slide 5 — Doublet Detection

**Tool:** Scrublet, run independently per sample (138 samples)

**Why per-sample?** Scrublet simulates doublets by combining real cell profiles. Running globally would mix patient/sample identities, inflating doublet scores for biologically distinct populations — a known artifact.

| Metric | Value |
|--------|-------|
| Total predicted doublets | 4,161 / 429,867 (0.97%) |
| Per-sample range | 0.0% – 3.5% |
| Notable outliers (>3%) | E38_B (3.2%), E38_B_DN (3.5%), E38_F1 (3.3%) |

**Decision: doublets flagged but NOT removed.** Scores stored in `.obs['doublet_score']` and `.obs['predicted_doublet']`. This preserves flexibility — downstream analyses can choose to exclude them or use the score as a covariate.

---

## Slide 6 — Dimensionality Reduction & Clustering

| Step | Parameters | Rationale |
|------|-----------|-----------|
| HVG selection | 3,000 genes, seurat_v3, batch_key="sample" | Batch-aware selection deprioritizes patient/library-driven variance |
| PCA | 50 components, randomized SVD | Captures majority of transcriptional variance; standard for datasets of this scale |
| kNN graph | k=15, 50 PCs | Balances local resolution vs. noise for 430k cells |
| UMAP | Default | Visualization only — not used for clustering |
| Leiden | Resolution 0.5 | Conservative initial resolution; produces interpretable cluster count |

**Result: 45 Leiden clusters** at resolution 0.5

### Timepoint distribution of retained cells
| Timepoint | Cells | % |
|-----------|-------|---|
| Baseline | 180,678 | 42.0% |
| FU1 | 170,467 | 39.7% |
| FU2 | 78,722 | 18.3% |

FU2 underrepresentation is expected — some patients may have discontinued or progressed before the second follow-up timepoint.

---

## Slide 7 — Clinical Metadata Integration

All 10 clinical covariates from Table S1 (outcomes) and Table S2 (sample info) merged into `.obs` at the patient level and broadcast to all cells.

| Column | Source | Fill rate |
|--------|--------|-----------|
| `progression_category` | Table S1 | **100%** |
| `best_overall_response` | Table S1 | **100%** |
| `MSI_status` | Table S1 | **100%** |
| `EBV_status` | Table S1 | **100%** |
| `HER2_status` | Table S1 | **100%** |
| `TCGA_subtype` | Table S1 | **100%** |
| `PDL1_baseline_CPS` | Table S1 | **100%** |
| `PFS_days` | Table S1 | **100%** |
| `OS_days` | Table S1 | **100%** |
| `timepoint_label` | `.obs['timepoint']` | **100%** |

100% fill rate across all 429,867 cells — every cell can be stratified by any clinical covariate immediately.

---

## Slide 8 — Key Assumptions & Limitations

### QC thresholds
| Assumption | Threshold | Risk if wrong |
|------------|-----------|--------------|
| MT% < 20% = viable cell | 20% | Too low → loses metabolically active cells (e.g., cardiomyocytes, monocytes); 20% is appropriate for tumor tissue |
| > 6,000 genes = doublet | 6,000 | Could exclude highly active tumor cells; should validate against Scrublet scores |
| < 200 genes = empty | 200 | Standard; may retain a small fraction of low-complexity real cells |
| Gene in ≥ 3 cells = retained | 3 cells | Very permissive; keeps rare transcripts but adds sparsity |

### Batch correction
- **No explicit batch correction applied** (no Harmony, scVI, BBKNN)
- HVG selection used `batch_key="sample"` — this only adjusts which genes are selected, not the embedding
- With 33 patients and 138 libraries, **patient/batch effects are likely present in the PCA/UMAP**
- **This is the most significant outstanding assumption**: clusters at this stage may partially reflect technical rather than biological variation

### Doublets
- Doublets flagged but not removed — if they cluster together they could form spurious "doublet clusters" that contaminate cell type annotation

### Ambient RNA
- No ambient RNA decontamination (e.g., SoupX, DecontX) was applied
- In tumor samples this can introduce false positive expression of highly expressed genes (e.g., HBB, EPCAM) into non-expressing cell types

### Clinical metadata
- Progression category (Slow/Fast) is patient-level — all cells from a patient inherit the same label regardless of sample type or timepoint
- Assumes labels from Table S1 are correctly matched to sample identifiers

---

## Slide 9 — What Has Been Completed

| Component | Status |
|-----------|--------|
| Raw data loading and inspection | ✅ Complete |
| QC metric computation and visualization | ✅ Complete |
| Cell and gene filtering | ✅ Complete |
| Per-sample doublet detection (Scrublet) | ✅ Complete |
| Normalization, log1p transformation | ✅ Complete |
| Batch-aware HVG selection (3,000 genes) | ✅ Complete |
| Scaling, PCA (50 PCs) | ✅ Complete |
| kNN graph, UMAP, Leiden clustering (45 clusters) | ✅ Complete |
| Clinical metadata integration (100% fill) | ✅ Complete |
| Reproducible, checkpoint-aware pipeline | ✅ Complete |
| Cell type annotation | ❌ Not started |
| Batch correction | ❌ Not started |
| Differential expression analysis | ❌ Not started |
| TME composition analysis | ❌ Not started |
| Cell-cell communication | ❌ Not started |

---

## Slide 10 — Proposed Future Directions

### Priority 1: Batch Correction
Apply **Harmony** or **scVI** to remove patient/library-driven technical variance before cell type annotation. Without this, marker genes and cluster identities are unreliable. This is the single most important next step.

### Priority 2: Cell Type Annotation
Annotate the 45 Leiden clusters using:
- Known gastric TME markers (CD3/CD8 for T cells, CD68 for macrophages, EPCAM for epithelial, VIM/FAP for fibroblasts, CD19/MS4A1 for B cells)
- Automated tools: **CellTypist** (pretrained on immune atlases) or **scANVI** (semi-supervised)
- Manual review of top DE genes per cluster

### Priority 3: TME Composition Analysis
Once cell types are annotated:
- Quantify cell type fractions per patient per timepoint
- Compare compositions between **Slow vs. Fast Progressors** and **PR vs. PD**
- Test whether baseline TME composition predicts treatment outcome

### Priority 4: Differential Expression — Responders vs. Non-responders
- Pseudobulk DE (using **DESeq2 via pydeseq2**) per cell type
- Identify genes/pathways enriched in slow progressors at baseline
- Stratify by TCGA subtype (EBV+, MSI-High likely respond differently)

### Priority 5: Longitudinal Dynamics
- Track cell type abundance changes: Baseline → FU1 → FU2
- Focus on CD8+ T cell exhaustion trajectories and macrophage polarization (M1/M2)
- Use **scVelo** or **CellRank** for RNA velocity / fate inference in T cell and myeloid lineages

### Priority 6: Cell-Cell Communication
- Apply **CellChat** or **NicheNet** to infer ligand-receptor interactions between cell types
- Compare interaction networks in responders vs. non-responders
- Focus on PD-1/PD-L1, CXCL, TGF-β signaling axes

### Priority 7: Survival Correlation
- Correlate cell type abundance scores (from scRNA-seq) with PFS and OS using Cox regression
- Build a simple TME score (e.g., CD8:Treg ratio, M1:M2 ratio) and test as biomarker

---

## Slide 11 — Summary

> **429,867 high-quality cells across 33 patients and 3 treatment timepoints are ready for biological analysis.**

The preprocessing pipeline is complete, reproducible, and clinically annotated. The critical next step is batch correction followed by cell type annotation — without those two steps, all downstream biological conclusions are premature.

The dataset is uniquely powered for longitudinal TME analysis because it captures:
1. Pre-treatment immune landscape (what predicts response)
2. On-treatment dynamics (what changes in responders)
3. Multiple molecular subtypes (EBV, MSI-H, HER2) allowing stratified analysis

**The most scientifically impactful question this dataset can answer:** *Are there TME features at baseline that predict whether a gastric cancer patient will be a slow vs. fast progressor on immunotherapy, and which cell types drive that signal?*

---

*Report prepared by Advaith Nair | gastric_tme_project/notebooks/01_preprocessing.py*
