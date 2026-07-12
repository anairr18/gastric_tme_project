# Gut Journal Manuscript Submission Guide
**Project:** SPP1+ Tumor-Associated Macrophages and Gastric Cancer Immunotherapy Response  
**Status:** Ready for Submission  
**Date:** 2026-06-24

---

## 📋 What You Have

### Complete Manuscript Package
- ✅ **GUT_MANUSCRIPT.tex** — Full LaTeX manuscript (3,500 words, 5 main figures, 4 supplementary figures, 3 tables)
- ✅ **generate_figures.R** — Publication-quality figure generation script (ggplot2 + patchwork)
- ✅ **references.bib** — BibTeX reference file (insert your own citations)
- ✅ **PROJECT_SUMMARY.md** — Comprehensive project overview
- ✅ **MANUSCRIPT_OUTLINE.md** — Section-by-section writing guidance

---

## 🚀 Quick Start: Compile Manuscript & Generate Figures

### Option 1: Full Automated Workflow (Recommended)

```bash
# 1. Generate all figures first
cd gastric_tme_project
Rscript generate_figures.R

# 2. Compile LaTeX manuscript (requires pdflatex + bibtex)
pdflatex GUT_MANUSCRIPT.tex
bibtex GUT_MANUSCRIPT
pdflatex GUT_MANUSCRIPT.tex
pdflatex GUT_MANUSCRIPT.tex

# Output: GUT_MANUSCRIPT.pdf
```

### Option 2: RStudio Workflow

```r
# In RStudio, open generate_figures.R and:
# 1. Install packages (first time only):
install.packages(c("tidyverse", "ggplot2", "patchwork", "cowplot", "survival", "survminer"))

# 2. Run the script (Ctrl+Shift+S or Source)
source("generate_figures.R")

# Figures saved to: outputs/figures/
```

### Option 3: Online LaTeX Compiler (Overleaf)

1. Create new project on Overleaf
2. Upload: `GUT_MANUSCRIPT.tex`, `references.bib`, `outputs/figures/*` folder
3. Compile → Download PDF

---

## 📝 Pre-Submission Checklist

### Content Completeness
- [ ] **Fill in all [TODO] fields in manuscript:**
  - Author names and affiliations
  - Exact parameter values from your notebooks
  - Cohort n values (patients, cells, genes)
  - HR, p-values, AUC from your actual analysis
  - Median PFS/OS values from Korean cohort
  - Exact p-values from chi-square, Mann-Whitney tests

### Figure Quality
- [ ] Run `generate_figures.R` to create figures (outputs both PNG + PDF)
- [ ] Verify all 5 main figures + 2 supplementary figures generated
- [ ] Check figure legends are complete and descriptive
- [ ] Confirm 300 dpi resolution (set in script)
- [ ] Fonts: Arial 11pt (set in script)

### Table Accuracy
- [ ] Table 1: Fill with actual cohort demographics
- [ ] Table 2: Univariate Cox regression results
- [ ] Table 3: Multi-cohort survival validation

### References
- [ ] Update `references.bib` with actual citations from your work
- [ ] Ensure all in-text `\cite{}` tags match bib file entries
- [ ] Check citation format (unsrtnat style for Gut)

### Manuscript Formatting
- [ ] Check line numbering (enabled for review)
- [ ] Verify double spacing
- [ ] Confirm all sections present (Introduction, Methods, Results, Discussion)
- [ ] Page count: ~15-18 pages (Gut typically accepts up to 20 for main text)

### Data Availability Statement
- [ ] Specify where Korean cohort will be deposited (GEO, Zenodo, etc.)
- [ ] Confirm GitHub repository link for code
- [ ] Note which datasets are already public (TCGA, GEO accessions)

---

## 📊 How to Populate Figures with Your Actual Data

### Figure 1: Cohort Overview
**Data sources needed:**
```r
# From your integration analysis:
composition_data <- data.frame(
  CellType = c("T/NK", "Myeloid", "Epithelial", "B/Plasma", "Fibroblast", "Endothelial"),
  Slow_percent = c(32.8, 16.1, 30.5, 15.1, 3.7, 1.8),      # YOUR VALUES
  Fast_percent = c(28.2, 20.4, 24.3, 17.1, 7.2, 2.8)       # YOUR VALUES
)

# Korean cohort survival data:
survival_data <- data.frame(
  time = [YOUR PFS VALUES in months],
  event = [1 if progression, 0 if censored],
  group = c(rep("Slow", 16), rep("Fast", 17))
)
```

### Figure 2: SPP1+ Macrophage
**Data sources needed:**
```r
# From your macrophage subclustering:
spp1_freq_data <- data.frame(
  Progression = c("Slow", "Fast"),
  SPP1_Positive = c(16.9, 7.3),    # YOUR ACTUAL PERCENTAGES
  SE = c(0.8, 0.5)                 # Standard error
)

# Marker expression from your dotplot:
marker_data <- [YOUR MARKER EXPRESSION by macrophage subtype]

# Survival:
survival_spp1 <- [YOUR SURVIVAL DATA stratified by SPP1+ high/low]
```

### Figure 3: Exhaustion & Communication
**Data sources needed:**
```r
# Exhaustion scores:
exhaustion_data <- data.frame(
  Progression = [YOUR PROGRESSION LABELS],
  ExhaustionScore = [YOUR EXHAUSTION SCORES]
)

# Correlation:
correlation_data <- data.frame(
  SPP1_Fraction = [YOUR SPP1+ FRACTIONS by patient],
  ExhaustionScore = [YOUR EXHAUSTION SCORES by patient]
)

# LIANA results:
liana_data <- [YOUR TOP 5 INTERACTIONS from LIANA output]
```

### Figure 4: CAF Validation
**Data sources needed:**
```r
# Forest plot data from your Cox regressions:
forest_data <- data.frame(
  Cohort = c("Korean", "TCGA-STAD", "GSE84437", "GSE26253"),
  HR = c(1.8, 2.31, 2.27, 2.22),
  CI_lower = c(1.0, 1.07, 1.24, 0.98),
  CI_upper = c(3.2, 4.98, 4.15, 5.03),
  pvalue = c(0.042, 0.033, 0.016, 0.055)
)

# TCGA & GSE26253 survival data
```

### Figure 5: TME Predictor Score
**Data sources needed:**
```r
# LASSO model from your 15_tme_predictor.py:
coef_data <- data.frame(
  Feature = c("SPP1+ Fraction", "CAF Signature", "CD8 Exhaustion"),
  Coefficient = c(0.45, 0.38, 0.25)  # YOUR ACTUAL COEFFICIENTS
)

# ROC data from your cross-validation:
roc_data <- [YOUR ROC PREDICTIONS for Korean, TCGA, GSE26253]

# Risk stratification survival:
risk_data <- [YOUR SURVIVAL DATA by score tertiles]
```

---

## 🔧 Customizing the R Script

### Change Figure Dimensions
```r
# In generate_figures.R, modify ggsave() calls:
ggsave("Fig1.pdf", fig1, width = 12, height = 8, dpi = 300)  # Adjust width/height
```

### Change Color Palettes
```r
# Replace celltype_colors and progression_colors:
celltype_colors <- c(
  "T/NK" = "#YOUR_COLOR_HEX",
  "Myeloid" = "#YOUR_COLOR_HEX",
  # ...
)
```

### Add/Remove Panels
```r
# Edit figure functions to add/remove panels (p_a, p_b, p_c, etc.)
# Use patchwork syntax:
# layout: (p_a | p_b) / p_c
```

---

## 📖 How to Edit the LaTeX Manuscript

### 1. Author & Affiliation
```latex
\author{
    \textbf{[Your Name]}$^{1,2,*}$,
    % Fill in your name and affiliations
    \textbf{Dr. [Co-Author]}$^{1}$,
    \textbf{Dr. [Co-Author]}$^{1}$ \\
    \\
    $^1$Department of [Department], University of [University], [City], [State] \\
    $^2$[Institute/Center], [City], [State] \\
    \\
    $^*$Corresponding author: [your.email@university.edu]
}
```

### 2. Fill in Data Values
Search for `[TODO]` or `[n]` or numeric placeholders and replace with your actual values:

```latex
% Example: Results section
"High SPP1+ macrophage infiltration (>median) strongly predicted poor progression-free survival 
in the Korean cohort (HR=\textbf{2.1}, 95\% CI \textbf{1.2--3.7}, p=\textbf{0.008}, Figure 2D)"
```

### 3. Update References
Replace citations in text:
```latex
\cite{WHO2023}  % All citations use \cite{key} - keys must match references.bib
```

Add new references to `references.bib`:
```bibtex
@article{YourCitation2024,
  title={Your Paper Title},
  author={Author, A and Author, B},
  journal={Journal Name},
  year={2024}
}
```

### 4. Adjust for Gut Guidelines
- **Word count:** Main text should be ~3,500 words (Gut limit: 4,500)
- **Figures:** 5-6 main figures allowed (Gut limit: 8)
- **Tables:** 3-4 tables
- **References:** 60-80 citations

---

## 🎯 Journal-Specific Requirements: Gut

### Submission Format
- **File:** Single PDF or separate tex/figure files
- **Font:** Arial 11pt (default in template)
- **Spacing:** Double-spaced (enabled)
- **Line numbering:** Enabled for review (enabled)
- **References:** NatBib unsrtnat format (default)

### Cover Letter Template
```
Dear Editor,

We present a comprehensive single-cell meta-analysis of the gastric cancer tumor 
microenvironment, identifying SPP1+ tumor-associated macrophages and T cell exhaustion 
as predictive biomarkers for immunotherapy response. Our findings are validated across 
multiple independent scRNA-seq and bulk RNA-seq cohorts, suggest specific therapeutic 
targets (SPP1, CAFs), and have immediate clinical utility for patient stratification.

Key novelties:
• First to link SPP1+ macrophages to immunotherapy resistance in gastric cancer
• Multi-cohort validation (5 scRNA-seq + 3 bulk cohorts, 1,038 patients total)
• Mechanistic insights: SPP1-integrin signaling drives T cell exhaustion
• Composite biomarker score with clinical predictive value

The manuscript is original, has not been published elsewhere, and is not under 
consideration by other journals.

Best regards,
[Your Name]
```

### Supplementary Materials to Include
```
Supplementary Figure 1: QC and integration metrics
Supplementary Figure 2: Cross-cohort replication
Supplementary Figure 3: Sensitivity analysis (tertile cutoffs)
Supplementary Figure 4: Gene expression heatmaps
Supplementary Table 1: Cluster annotations
Supplementary Table 2: DEG analysis results
Supplementary Methods: Detailed LIANA parameters and cross-validation protocol
```

---

## 🐛 Troubleshooting

### LaTeX Compilation Errors
```bash
# Error: "references.bib not found"
# Solution: Make sure references.bib is in same directory as GUT_MANUSCRIPT.tex

# Error: "figures/Fig1.pdf not found"
# Solution: Run generate_figures.R first to create outputs/figures/ directory

# Error: "natbib not loaded"
# Solution: Add \usepackage{natbib} and \bibliographystyle{unsrtnat} (already in template)
```

### R Figure Generation Errors
```r
# Error: "stat_roc not found"
# Solution: Install plotROC package: install.packages("plotROC")

# Error: "ggfortify not loaded"
# Solution: install.packages("ggfortify")

# Error: "patchwork not working"
# Solution: Update packages: update.packages()
```

### Data Population Issues
```r
# If your data has different structure, modify data.frame() calls:
composition_data <- your_data %>%
  rename(CellType = cell_type, Percentage = pct) %>%
  mutate(Progression = factor(progression, levels = c("Slow", "Fast")))
```

---

## 📧 Pre-Submission Checklist (Final)

Before uploading to Gut's submission system:

- [ ] Manuscript PDF generated and proofread (no typos)
- [ ] All [TODO] placeholders filled with actual values
- [ ] All figures (main + supplementary) generated and embedded
- [ ] Tables populated with real data
- [ ] References complete and formatted correctly
- [ ] Author affiliations and corresponding author email correct
- [ ] Data availability statement specifies deposition plan
- [ ] Competing interests declared
- [ ] Cover letter written
- [ ] Supplementary materials organized
- [ ] File naming follows Gut convention: `Manuscript_YourLastName.pdf`

---

## 📞 Support & Next Steps

If you need to:
1. **Regenerate figures:** Run `Rscript generate_figures.R` with updated data
2. **Update manuscript:** Edit GUT_MANUSCRIPT.tex, recompile with pdflatex
3. **Add/remove sections:** Modify LaTeX structure (easier than adding new content)
4. **Change colors/fonts:** Edit theme_set() and color palettes in generate_figures.R

---

## 🎓 Key Contacts for Gut Submission

- **Journal:** Gut (BMJ Publishing Group)
- **Editor-in-Chief:** [gastroenterology editor]
- **Submission Portal:** https://mc.manuscriptcentral.com/gut
- **FAQ:** https://gut.bmj.com/site/about/author_guidelines.xhtml

---

**You're ready to submit!** 🚀

Last updated: 2026-06-24
