# scRNA-seq QC Summary 

## 1. Dataset Overview

| | |
|---|---|
| **Patients** | 33 |
| **Samples (libraries)** | 138 (tumor, adjacent normal, distal normal across timepoints) |
| **Timepoints** | Baseline, Follow-up 1 (FU1), Follow-up 2 (FU2) |
| **Raw cells** | 654,770 |
| **Raw genes** | 36,601 |

Sample types in `.obs['normal']`: tumor (NaN), adjacent normal (AN), distal normal (DN).

---

## 2. Pre-filtering QC Metrics

Mitochondrial genes detected: **13** (prefix `MT-`)

| Metric | Median | Mean |
|---|---|---|
| Genes per cell (`n_genes_by_counts`) | 844 | 1,067 |
| UMI counts per cell (`total_counts`) | 2,217 | 3,998 |
| % mitochondrial reads (`pct_counts_mt`) | 7.88% | 20.53% |

The elevated mean MT% relative to the median reflects a population of low-quality / dying cells with high MT content that are removed during filtering.

QC plots: `outputs/qc_plots/violin_qc_violin.png`, `scatter_mt_scatter.png`, `scatter_genes_scatter.png`

---

## 3. Cell and Gene Filtering

### Cell filters applied
| Filter | Cells removed |
|---|---|
| Fewer than 200 genes expressed | 92,664 |
| More than 6,000 genes expressed (multiplets / doublets) | 3,068 |
| MT% > 20% (low-quality / apoptotic) | 129,171 |
| **Total removed** | **224,903 (34.3%)** |

### Gene filter applied
| Filter | Genes removed |
|---|---|
| Expressed in fewer than 3 cells | 3,910 (10.7%) |

### Post-filtering
| | |
|---|---|
| **Cells retained** | **429,867** |
| **Genes retained** | **32,691** |

---

## 4. Doublet Detection (Scrublet)

Scrublet was run independently on each of the 138 samples to avoid cross-sample contamination of the doublet simulation. Per-sample doublet rates were generally low and biologically plausible.

| | |
|---|---|
| **Total predicted doublets** | 4,161 / 429,867 (0.97%) |
| Per-sample range | 0.0% – 3.5% |
| Notable outliers (>3%) | E38_B (3.2%), E38_B_DN (3.5%), E38_F1 (3.3%) |

Doublet scores and labels stored in `.obs['doublet_score']` and `.obs['predicted_doublet']`. Doublets were flagged but **not removed** at this stage to preserve flexibility for downstream analysis.

---

## 5. Normalization and Feature Selection

| Step | Parameters |
|---|---|
| Normalization | Total-count normalization to 10,000 counts per cell |
| Transformation | log1p |
| Highly variable genes (HVGs) | 3,000 genes, `seurat_v3` flavor, batch-corrected per sample |

HVG selection used per-sample batch correction (`batch_key="sample"`) to prioritize genes variable across biological conditions rather than technical batches.

---

## 6. Dimensionality Reduction and Clustering

| Step | Parameters |
|---|---|
| PCA | 50 components, randomized SVD |
| Nearest neighbours | k=15, using top 50 PCs |
| UMAP | Default parameters |
| Leiden clustering | Resolution 0.5 |

PCA variance ratio plot: `outputs/qc_plots/pca_variance_ratio_pca_variance.png`

**45 Leiden clusters** identified at resolution 0.5.

**Timepoint distribution of retained cells:**

| Timepoint | Cells |
|---|---|
| Baseline | 180,678 (42.0%) |
| FU1 | 170,467 (39.7%) |
| FU2 | 78,722 (18.3%) |

---

## 7. Clinical Metadata Integration

All clinical metadata from Table S1 (clinical outcomes) and Table S2 (sample contributions) were successfully appended to `.obs`. All 10 metadata columns achieved **100% fill rate** across all 429,867 cells.

| Column | Source | Fill rate |
|---|---|---|
| `progression_category` (Slow/Fast) | Table S1 | 100% |
| `best_overall_response` (PD/SD/PR) | Table S1 | 100% |
| `MSI_status` | Table S1 | 100% |
| `EBV_status` | Table S1 | 100% |
| `HER2_status` | Table S1 | 100% |
| `TCGA_subtype` | Table S1 | 100% |
| `PDL1_baseline_CPS` | Table S1 | 100% |
| `PFS_days` | Table S1 | 100% |
| `OS_days` | Table S1 | 100% |
| `timepoint_label` (Baseline/FU1/FU2) | `.obs['timepoint']` | 100% |

---

## 8. Output Files

| File | Description |
|---|---|
| `data/processed/gastric_processed.h5ad` | Final processed AnnData (429,867 cells × 3,000 HVGs) |
| `data/processed/checkpoint_post_scrublet.h5ad` | Intermediate checkpoint (post-QC, pre-normalization) |
| `outputs/qc_plots/violin_qc_violin.png` | Violin plots of QC metrics |
| `outputs/qc_plots/scatter_mt_scatter.png` | Total counts vs. MT% scatter |
| `outputs/qc_plots/scatter_genes_scatter.png` | Total counts vs. genes/cell scatter |
| `outputs/qc_plots/pca_variance_ratio_pca_variance.png` | PCA variance explained |
| `notebooks/01_preprocessing.py` | Reproducible preprocessing script |

---

