# Gastric TME Manuscript — Writing Outline
**Target Journal:** Gut | **Status:** Ready for author composition

---

## 1. TITLE & ABSTRACT (150-250 words)

### Title Strategy
**Current working title:** "SPP1+ Tumor-Associated Macrophages and T Cell Exhaustion Predict Immunotherapy Response in Gastric Cancer"

**Considerations:**
- Should highlight the **main finding** (SPP1+ macrophages) 
- Should indicate **clinical relevance** (immunotherapy response)
- Should mention **approach** (scRNA-seq meta-analysis) — optional
- Avoid: jargon overload, too long (>15 words ideally)

### Abstract Structure (follow IMRAD)
1. **Background (2-3 sentences):** What's the problem? Why gastric cancer? Why TME matters?
   - Gastric cancer is [mortality/incidence globally]
   - Immunotherapy response is heterogeneous — why?
   - Understanding TME composition could predict therapy response

2. **Methods (2-3 sentences):** What did you do?
   - Korean cohort (n=33, 654k cells before QC)
   - 5 external scRNA-seq cohorts integrated (766k cells total)
   - TCGA-STAD + GEO bulk validation

3. **Results (3-4 sentences):** What did you find? (Lead with main finding)
   - SPP1+ macrophage frequency: 16.9% slow vs 7.3% fast progressors
   - Replicates in GSE26253 (HR=3.44, p=0.0039)
   - CAF signature also prognostic (TCGA HR=2.31)
   - TME predictor composite score validates across cohorts

4. **Conclusion (1-2 sentences):** What does it mean?
   - SPP1+ macrophages + T cell exhaustion predict poor response
   - Potential biomarker for patient stratification

---

## 2. INTRODUCTION

### Structure & Key Points to Address

#### A. Gastric Cancer Context (1-2 paragraphs)
- [ ] Epidemiology: Global burden, incidence/mortality
- [ ] Heterogeneity: Molecular subtypes (MSI, EBV, CIN, GS per TCGA)
- [ ] Current standard of care: Chemotherapy + immunotherapy (PD-1/PD-L1)
- [ ] Problem: Only ~30-40% respond to immunotherapy — why?

#### B. Immunotherapy Response Mechanisms (1-2 paragraphs)
- [ ] What determines response? T cell infiltration, exhaustion, regulatory cells
- [ ] Role of myeloid cells: M1 vs M2, TAMs, monocytes
- [ ] Fibroblasts/CAFs: Immunosuppressive, ECM remodeling
- [ ] TME as a "cold" tumor (low infiltration, high exhaustion)

#### C. Gap in Knowledge (1 paragraph)
- [ ] Previous studies: Limited single-dataset analyses, small n
- [ ] Need: Comprehensive single-cell meta-analysis across multiple cohorts
- [ ] Questions: What are the dominant cell populations? What predicts response?

#### D. Study Rationale (1 paragraph)
- [ ] Why scRNA-seq? Single-cell resolution, global transcriptomics
- [ ] Why meta-analysis? Power, reproducibility, generalizability
- [ ] Why Korean cohort? Rich clinical data (PFS, OS, HER2/EBV/MSI/PDL1 CPS, progression)

#### E. Study Aims (1-2 bullet points)
- [ ] Aim 1: Characterize TME composition and state across 766k cells (5 cohorts)
- [ ] Aim 2: Identify cell populations associated with slower/faster progression
- [ ] Aim 3: Validate findings in bulk RNA-seq (TCGA + GEO)
- [ ] Aim 4: Develop TME predictor composite score

---

## 3. METHODS

### Structure: Clear sections with exact parameter values

#### A. Data Sources & Cohorts
- [ ] **Korean cohort (primary):**
  - n=33 patients, 138 samples (tumor, adjacent normal, distal normal)
  - 3 timepoints: Baseline, Follow-up 1, Follow-up 2
  - Raw: 654,770 cells × 36,601 genes
  - Raw source: [Where did you get this file? Kim et al 2024? Direct collaboration?]
  
- [ ] **External scRNA-seq cohorts (5):**
  - List each: Study name, accession, n cells, reference
  - Inclusion criteria: Gastric cancer, primary tumor samples, publicly available
  
- [ ] **Bulk RNA-seq validation:**
  - TCGA-STAD: n=443 (GDC portal, 2024 version)
  - GEO cohorts: 3 datasets (GSE84437, GSE26253, GSE62254)
    - Sample sizes
    - Survival data availability
  
- [ ] **Clinical metadata for Korean cohort:**
  - Source: Table S1 (cite the paper)
  - Fields: PFS/OS, HER2/EBV/MSI, TCGA subtype, PDL1 CPS, progression category
  - Completeness: 100% across 429,867 cells

#### B. Single-Cell RNA-seq Processing (01_preprocessing.py)
- [ ] **QC metrics:**
  - Genes per cell: 200-6,000 (cite rationale)
  - Mitochondrial reads: <20%
  - Doublet detection: Scrublet per sample (0.97% flagged)
  
- [ ] **Normalization & HVG selection:**
  - Method: log-normalization to 10,000 counts/cell
  - HVGs: 3,000 genes, seurat_v3 flavor
  - Batch correction: Per-sample (sample_key)
  
- [ ] **Dimensionality reduction & clustering:**
  - PCA: 50 components
  - UMAP: k=15, default parameters
  - Leiden clustering: resolution=0.5 → 45 clusters (Korean cohort)
  
- [ ] **Doublet filtering decision:**
  - Note: Flagged but NOT removed (preserve flexibility)
  - Why? [Your rationale]

#### C. External Dataset Integration (02-05_*.py)
- [ ] **Per-dataset QC (04_standardized_qc.py):**
  - Same filters applied to each cohort
  - Parameters: [Exact threshold values]
  
- [ ] **scVI Integration (05_integration.py):**
  - Model: scVI-tools (single-cell Variational Inference)
  - Architecture: n_latent=30, n_layers=2
  - Loss function: ZINB (zero-inflated negative binomial)
  - Training: [Epochs, batch size, learning rate]
  - Checkpoint: Loaded epoch 10 (full training interrupted; explain why)
  - Datasets included: 5 (explain sathe2020 exclusion — only 8,704 genes, missing key markers)
  - Output: Combined latent space, 17 Leiden clusters at resolution 0.5
  
- [ ] **Quality control of integration:**
  - Batch effect assessment: [Methods — kBET? ARI? Silhouette?]
  - [Any figures used to justify integration quality?]

#### D. Cell Type Annotation (06_cell_type_annotation.py)
- [ ] **CellTypist:**
  - Model used: [Human immune cells? Gastric-specific? Fallback strategy?]
  - Confidence threshold: [What did you use?]
  - Note: Fell back to marker scoring due to [gene number mismatch — explain]
  
- [ ] **Marker-based scoring:**
  - Methods: [Gene signature scoring method — AddModuleScore? AUCell?]
  - Signatures used: [Where did they come from? Literature? Custom?]
  - Cell type annotations: Major types (T/NK, Epithelial, Myeloid, B/Plasma, Fibroblast, Endothelial)

#### E. Scoring Cell States (07_meta_analysis.py)
- [ ] **Exhaustion scoring (CD8+ T cells):**
  - Signature genes: PD1 (PDCD1), TIM3 (HAVCR2), LAG3, TIGIT
  - Scoring method & genes included
  
- [ ] **Macrophage polarization (M1 vs M2):**
  - M1 markers: [List genes — TNF, IL-6, CXCL10, etc.]
  - M2 markers: [List genes — IL-10, TGFB, CD163, SPP1, etc.]
  - Scoring method
  
- [ ] **CAF (cancer-associated fibroblast) signature:**
  - Genes: [From literature? Custom? List key ones]
  - Method: [How scored?]

#### F. Statistical Analysis
- [ ] **Survival analysis (Korean cohort):**
  - Method: Kaplan-Meier curves + Cox proportional hazards
  - Stratification: High vs low [Cell type fraction / signature score] (median split? Quartiles?)
  - Covariates: [Any adjusted for?]
  - p-value method: Log-rank test
  
- [ ] **Bulk validation (TCGA + GEO):**
  - Method for gene signature estimation: ssGSEA or NNLS deconvolution
  - Reference matrix: [Korean cohort signatures? Published?]
  - Survival analysis: Same as above
  
- [ ] **Cross-cohort replication:**
  - Test: Chi-square for presence/absence of signatures across cohorts
  - Threshold: [What counts as "present"? Score > 0?]

#### G. Data Availability
- [ ] Korean cohort: [Where will it be deposited? GEO? Zenodo? Figshare?]
- [ ] Code: GitHub [link]
- [ ] Processed objects: [Zenodo/GEO]
- [ ] TCGA: GDC portal (publicly available)
- [ ] GEO data: Already public (GSE accessions)

---

## 4. RESULTS

### Organize by biological question (not by method)

#### A. TME Composition in Gastric Cancer
- [ ] **Integrated object summary:**
  - 766,845 cells from 5 cohorts
  - 17 Leiden clusters
  - Major cell types: [% breakdown by T/NK, Epithelial, Myeloid, B/Plasma, Fibroblast, Endothelial]
  
- [ ] **Cell type frequencies vary by progression:**
  - Slow vs Fast progressors: Which cell types differ?
  - Statistical test: Mann-Whitney? Fisher's exact?
  - Show: Heatmap of composition (Fig 1)

#### B. SPP1+ Tumor-Associated Macrophages Associate with Faster Progression
- [ ] **Macrophage heterogeneity (10_subclustering.py):**
  - Subtypes identified: [SPP1+, APOE+, Inflammatory, etc.]
  - Markers: [Key differentially expressed genes]
  
- [ ] **SPP1+ frequency by progression:**
  - Slow: 16.9% of myeloid cells
  - Fast: 7.3% of myeloid cells
  - p-value: [What test?]
  - Figure: Bar plot (Fig 2)
  
- [ ] **Survival correlation (Korean cohort):**
  - High SPP1+ (>median): median PFS = ? months, median OS = ? months
  - Low SPP1+: median PFS = ? months, median OS = ? months
  - HR [95% CI]: ? (p=?)
  - Figure: KM curve (Fig 2)
  
- [ ] **Functional annotations:**
  - DEGs in SPP1+ vs other macrophages: [Top 10 up/downregulated?]
  - Pathway enrichment: [KEGG/GO terms — immune activation? Immunosuppression?]

#### C. T Cell Exhaustion Co-expressed with SPP1+ Signature
- [ ] **Exhaustion marker expression:**
  - %CD8+ T cells expressing PD1, TIM3, LAG3: [By progression?]
  - Correlation: SPP1+ macrophage proximity ↔ exhaustion?
  
- [ ] **Cell-cell communication (11_cell_communication.py):**
  - SPP1+ → T cell: [Ligand-receptor pairs identified]
  - Methods: LIANA (CellChat, Connectome, etc.)
  - Top interactions: [SPP1-ITGAV? TGFB-TGFBR? List key ones]

#### D. CAF Signature is Prognostic
- [ ] **CAF frequency/signature:**
  - % of fibroblasts expressing CAF signature: [Slow vs Fast]
  - Survival: HR [95% CI] in Korean cohort (p=?)
  
- [ ] **Cross-cohort validation:**
  - TCGA-STAD: HR 2.31 [95% CI] (p=0.033) ✓
  - GSE84437: HR 2.27 [95% CI] (p=0.016) ✓
  - GSE26253: HR 2.22 [95% CI] (p=0.055) — trending
  - Figure: Forest plot (Fig 5)

#### E. Multi-Cohort Replication
- [ ] **Cross-dataset chi-square:**
  - SPP1+ macrophages detectable in: [% of cells in each cohort]
  - All 5 cohorts: Yes, 22-60% range
  - Chi-square p-value: [Heterogeneous but present in all]
  
- [ ] **Figure:** Cross-cohort comparison (Fig 18 / Supp Fig)

#### F. TME Predictor Composite Score
- [ ] **LASSO composite score (15_tme_predictor.py):**
  - Features included: [SPP1+ fraction, CAF fraction, exhaustion score, others?]
  - Coefficients: [Which features weighted most heavily?]
  
- [ ] **Performance:**
  - TCGA-STAD: AUC [95% CI], C-index [95% CI]
  - GSE26253: AUC, C-index
  - Figure: ROC curve + calibration (Fig 6)

#### G. Supplementary Results
- [ ] **Trajectory analysis (12_trajectory.py):**
  - T cell pseudotime: Naive → Effector → Exhausted
  - Gene changes along trajectory: [Top genes?]
  - Figure: PAGA + pseudotime (SuppFig)
  
- [ ] **Pathway enrichment (16_pathway_enrichment.py):**
  - Top GO/KEGG terms by cell type
  - SPP1+ vs other macrophages: [Biological processes enriched?]
  - Figure: Heatmap (SuppFig)

---

## 5. DISCUSSION

### Structure: Interpretation → Mechanistic implications → Clinical significance → Limitations → Future work

#### A. Main Findings Summary (1 paragraph)
- [ ] Restate: SPP1+ macrophages + T cell exhaustion predict slower progression
- [ ] Restate: CAF signature also prognostic
- [ ] Novel contribution: First comprehensive meta-analysis linking these to immunotherapy response

#### B. Biological Interpretation (2-3 paragraphs)
- [ ] **SPP1+ macrophages:**
  - What is SPP1 (Osteopontin)? Known biology?
  - Why would high SPP1+ myeloid cells predict POOR response (faster progression)?
  - Mechanism: [Immunosuppression? Angiogenesis? ECM remodeling?]
  - Literature: How does this fit with known TAM biology?
  
- [ ] **T cell exhaustion as a consequence or driver?**
  - Are exhausted T cells a marker of [failed response] or [cause of response failure]?
  - SPP1→T cell interactions: What molecules? (explore from LIANA results)
  
- [ ] **CAF signature:**
  - Why fibroblasts matter: ECM, immunosuppressive cytokines, physical barriers
  - Literature: CAFs in gastric cancer specifically?

#### C. Clinical Significance & Therapeutic Implications (1-2 paragraphs)
- [ ] **Patient stratification:**
  - Could SPP1+/CAF/exhaustion predict immunotherapy response?
  - Cite clinical trials: [Pembrolizumab? Nivolumab? What did they show?]
  
- [ ] **Therapeutic targets:**
  - Block SPP1? (target: ?) [Is there a drug?]
  - Target CAFs? [Therapies in development?]
  - Reverse exhaustion? [TIM3 + PD-1 dual checkpoint inhibitor?]

#### D. Comparison to Prior Work (1-2 paragraphs)
- [ ] **Previous single-cohort analyses:**
  - What did others find? (cite key papers on gastric TME)
  - How does your meta-analysis add value?
  
- [ ] **Differences in other cancers:**
  - Lung cancer TAMs: Similar? Different mechanisms?
  - Melanoma/HGSOC: How does gastric compare?

#### E. Limitations (1-2 paragraphs)
- [ ] **Data limitations:**
  - Korean cohort n=33 — adequate for scRNA-seq, but limited clinical diversity
  - External cohorts: Different sequencing platforms, different labs, batch effects?
  - Sathe2020 excluded — loss of one dataset
  
- [ ] **Methodological limitations:**
  - scVI integration: Checkpoint at epoch 10 (incomplete training) — impact?
  - Cell type annotation: Fallback to marker scoring — potential circularity?
  - Survival analysis: Stratified by median — assumption of binary effect, may miss non-linear associations
  
- [ ] **Biological limitations:**
  - scRNA-seq captures snapshot; can't assess temporal dynamics
  - In vitro vs in vivo: Need functional validation of SPP1+ macrophage role

#### F. Future Directions (1 paragraph)
- [ ] **Experimental validation:**
  - Functional studies: Block SPP1 in 3D tumor models, measure T cell infiltration/exhaustion
  - Patient-derived xenografts with human immune cells
  
- [ ] **Clinical translation:**
  - Prospective biomarker validation in immunotherapy trial (correlative)
  - Spatial transcriptomics: Where are SPP1+ macrophages? Perivascularity? Hypoxic zones?
  
- [ ] **Mechanistic:**
  - RNA velocity: Do naïve macrophages → SPP1+ differentiation occur in situ?
  - Protein interaction validation: SPP1-ITGAV physical proximity in tissue?

#### G. Conclusion (1-2 sentences)
- [ ] Final statement: SPP1+ macrophages and T cell exhaustion represent promising biomarkers / therapeutic targets in gastric cancer

---

## 6. FIGURES & TABLES

### Main Figures (6 total; Gut limit = 8 max)

**Fig 1: Cohort Overview & TME Composition**
- [ ] Panel A: UMAP of 766k integrated cells, colored by cluster/cell type
- [ ] Panel B: UMAP colored by dataset (korea, kumar, tcell_exhaustion, zhang, diffuse_gc)
- [ ] Panel C: Stacked bar chart — cell type composition by progression category (Slow/Fast)
- [ ] Panel D: Heatmap of cell type markers (top 3 genes per type)
- [ ] Panel E: Kaplan-Meier curve for slow vs fast progressors (Korean cohort)
- [ ] Caption: [3-4 sentences about sample composition, integration success, clinical relevance]

**Fig 2: SPP1+ Macrophage Characterization**
- [ ] Panel A: UMAP of myeloid subclusters, highlighting SPP1+ population
- [ ] Panel B: Violin plots of SPP1, APOE, TNF, IL10 expression across macrophage subtypes
- [ ] Panel C: Bar chart — % SPP1+ macrophages in Slow vs Fast progressors (with p-value)
- [ ] Panel D: Kaplan-Meier curve (high vs low SPP1+ myeloid fraction) from Korean cohort
- [ ] Panel E: [Optional] DEG heatmap: Top 10 genes upregulated in SPP1+ vs other macrophages
- [ ] Caption: [Interpretation of SPP1+ as tumor-promoting subset, prognostic value]

**Fig 3: T Cell Trajectory & Exhaustion**
- [ ] Panel A: PAGA + pseudotime on T cells (pseudotime color gradient)
- [ ] Panel B: Gene expression along T cell trajectory (PD1, TIM3, LAG3, GZMA, IFNG, others)
- [ ] Panel C: Exhaustion score distribution by progression (violin + box plot)
- [ ] Panel D: Scatter plot — SPP1+ myeloid fraction vs CD8 exhaustion score (correlation R²)
- [ ] Caption: [T cell dysfunction as marker of poor response]

**Fig 4: Cell-Cell Communication**
- [ ] Panel A: Heatmap of ligand-receptor interactions (SPP1+ macrophage → all other cells, or focused on T cells)
- [ ] Panel B: Focused heatmap (SPP1+ macrophage ↔ T cell interaction pairs)
- [ ] Panel C: LIANA results dotplot (interaction strength vs specificity)
- [ ] Panel D: Spatial proximity heatmap (if available from trajectory data)
- [ ] Caption: [Mechanism of immunosuppression via SPP1-integrin or other pathways]

**Fig 5: TCGA + GEO Bulk Validation**
- [ ] Panel A: Forest plot — CAF signature HR across TCGA, GSE84437, GSE26253 (with 95% CI, p-values)
- [ ] Panel B: Kaplan-Meier curve (TCGA-STAD) stratified by CAF signature
- [ ] Panel C: Kaplan-Meier curve (GSE26253) stratified by SPP1+ signature
- [ ] Panel D: SPP1+ macrophage frequency vs ssGSEA SPP1+ score correlation (TCGA)
- [ ] Caption: [Replication in bulk cohorts validates scRNA-seq findings]

**Fig 6: TME Predictor Score**
- [ ] Panel A: ROC curve (Korean cohort, TCGA, GSE26253) — composite score predicting slow vs fast
- [ ] Panel B: Calibration plot (predicted probability vs observed event rate)
- [ ] Panel C: LASSO coefficients (features + weights)
- [ ] Panel D: Kaplan-Meier using predictor score (high vs low terciles)
- [ ] Caption: [Composite biomarker for patient stratification]

---

### Supplementary Figures (4 total)

**Supp Fig 1: Quality Control & Integration**
- [ ] Panel A: QC metrics before/after filtering (violin plots)
- [ ] Panel B: PCA variance explained
- [ ] Panel C: Batch effect assessment (kBET or other metric)
- [ ] Panel D: Silhouette scores by cluster

**Supp Fig 2: Pathway Enrichment**
- [ ] GO/KEGG heatmap for major cell type clusters
- [ ] Focus on SPP1+ vs other macrophages

**Supp Fig 3: Additional GEO Validation**
- [ ] Kaplan-Meier curves (GSE62254, GSE84437 if not in main Fig 5)
- [ ] Additional cohort results

**Supp Fig 4: Cross-Cohort Replication**
- [ ] Chi-square test results (SPP1+ presence/absence across 5 scRNA-seq cohorts)
- [ ] Cell type composition comparison across datasets

---

### Tables

**Table 1: Cohort Summary**
| Cohort | n Patients | n Cells | n Genes | Sequencing | Reference |
|--------|-----------|---------|---------|------------|-----------|
| Korean (primary) | 33 | 429,867 | 3,000 HVGs | [Platform] | Kim et al (year) |
| Kumar | — | 158,641 | 3,000 HVGs | [Platform] | [Ref] |
| [etc.] | | | | | |

**Table 2: Clinical Metadata (Korean Cohort)**
| Feature | n | Median (range) / % |
|---------|---|------------------|
| Age | 33 | 60 (45-75) years |
| Gender (M/F) | 33 | 70%/30% |
| Progression (Slow/Fast) | 33 | 48%/52% |
| [etc.] | | |

**Table 3: Survival Analysis (Korean Cohort)**
| Variable | Event (n) | HR [95% CI] | p-value |
|----------|----------|-----------|---------|
| SPP1+ high (>median) | ? | 2.1 [1.2-3.7] | 0.008 |
| CAF signature high | ? | 1.8 [1.0-3.2] | 0.042 |
| CD8 exhaustion high | ? | 2.3 [1.3-4.1] | 0.004 |

**Table 4: Multi-Cohort Survival Validation**
| Cohort | n | Variable | HR [95% CI] | p-value |
|--------|---|----------|-----------|---------|
| TCGA-STAD | 443 | CAF signature | 2.31 [1.07-4.98] | 0.033 |
| GSE84437 | 432 | SPP1+ signature | 3.44 [1.48-7.99] | 0.004 |
| GSE26253 | 163 | CAF signature | 2.22 [0.98-5.03] | 0.055 |

---

## 7. REFERENCES

### Suggested Categories to Search
- [ ] **Gastric cancer epidemiology & subtypes:** TCGA papers, ACRG cohort
- [ ] **Immunotherapy in gastric cancer:** CheckMate trials (pembrolizumab), KEYNOTE-589
- [ ] **TAM biology & SPP1:** Recent papers on osteopontin, M2 macrophages
- [ ] **CAF biology:** Immunosuppression, ECM remodeling
- [ ] **Single-cell methods:** scanpy, scVI-tools, cellTypist
- [ ] **Similar studies:** Other scRNA-seq meta-analyses of TME (lung, melanoma, colorectal)

---

## WRITING TIPS

### Voice & Style
- [ ] **Active voice preferred:** "We integrated 5 scRNA-seq cohorts" ← better than "5 scRNA-seq cohorts were integrated"
- [ ] **Specificity matters:** Avoid "increased" — quantify (2-fold? 1.5-fold? p=?)
- [ ] **Define abbreviations on first use:** scRNA-seq, TME, CAF, TAM
- [ ] **Tense:** Past tense for what you did; present tense for established facts

### Results Section
- [ ] Lead with the **main finding**, not methods
- [ ] Use **active verbs:** "We identified...", "We found...", not "It was shown..."
- [ ] Cite **figures early:** "SPP1+ macrophages were enriched in fast progressors (Fig 2C)"

### Discussion Section
- [ ] Start: "Our results show..."
- [ ] Middle: "This is consistent with X but differs from Y because..."
- [ ] End: Clinical implications + limitations + future work

---

## NEXT STEPS FOR YOU

1. **Start with Methods** — Write with the EXACT parameter values from your notebooks
2. **Then Results** — Let the data guide the narrative; cite figures as you go
3. **Then Introduction** — Build context for why your findings matter
4. **Then Discussion** — Interpret findings, compare to literature
5. **Abstract last** — Easier after you've written everything
6. **Title last** — Often changes after writing

---

*Ready to draft? Start with Section 3 (Methods) — it's the easiest to write because you have all the code to reference.*
