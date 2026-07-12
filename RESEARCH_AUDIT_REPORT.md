# COMPREHENSIVE RESEARCH AUDIT REPORT
## Gastric TME Single-Cell Meta-Analysis Project

**Audit Date:** 2026-06-24  
**Auditor:** Senior Research Auditor  
**Project Status:** ⚠️ **CRITICAL ISSUES IDENTIFIED** — Ready for submission only after remediation

---

## EXECUTIVE SUMMARY

This project demonstrates **solid experimental design and multi-cohort validation**, but contains **one CRITICAL discrepancy** and several **MEDIUM-risk vulnerabilities** that must be addressed before publication or submission to Gut.

### Critical Issue (Must Fix)
- **Manuscript-Reality Mismatch on scVI Training:** Manuscript claims "epoch 10 checkpoint used" due to interrupted training, but actual evidence shows model trained to **epoch 120+** with final model generated 2026-06-16 16:17 (after all checkpoints complete).

### Risk Level Assessment
| Severity | Count | Impact |
|----------|-------|--------|
| **CRITICAL** | 1 | Methodological misrepresentation; affects Methods credibility |
| **HIGH** | 3 | Incomplete validation chains; missing robustness checks |
| **MEDIUM** | 4 | Statistical assumptions not fully justified; citation gaps |
| **LOW** | 5 | Documentation, clarity, and minor inconsistencies |

---

## SECTION 1: METHODOLOGICAL & DATA INTEGRITY

### 1.1 CRITICAL: scVI Training State Mismatch ⚠️

**Finding:**  
The manuscript states scVI was trained to epoch 10, but file system evidence contradicts this:

```
Actual Checkpoints Found:
- scvi_model_ckpt_ep010  (2026-06-15 23:42)
- scvi_model_ckpt_ep020  (2026-06-16 01:17)
- scvi_model_ckpt_ep030  (2026-06-16 02:26)
- ...
- scvi_model_ckpt_ep120  (2026-06-16 13:48)
- Final model.pt        (2026-06-16 16:17) ← AFTER all checkpoints

Integrated object output: 2026-06-17 00:46 (AFTER full training)
```

**What Manuscript Claims:**
> "Full training was interrupted at epoch 10 due to an unforeseen system shutdown. We assessed whether checkpoint loading at epoch 10 was adequate by examining the elbo curve, which had stabilized..."

**What Data Shows:**
Training ran to at least epoch 120, then final model was saved 16:17 on June 16, then integrated object was generated 2026-06-17 00:46. This is incompatible with "interrupted at epoch 10" narrative.

**Root Cause Analysis:**
- Either: Training actually completed, contradicting the interrupt claim
- Or: Model was trained to epoch 120, then REVERTED to epoch 10 for the published analysis (not documented)

**Impact:**
- **Credibility damage:** Reviewers will flag this as a "false limitation" designed to deflect criticism
- **Reproducibility question:** If epoch 10 was genuinely used, why are epoch 120 checkpoints present?
- **Trust issue:** If full training WAS done, claiming it wasn't undermines the Methods section

**How Hostile Reviewer Would Attack:**
> "Authors claim epoch 10 checkpoint was used due to interrupted training, yet file timestamps show the model trained to epoch 120 and the final integrated object was generated after this training completed. This is either a documentation error (actual model IS fully trained and claims are false) or a methodological error (training was performed but not disclosed). Either way, the Methods section cannot be trusted."

**Fix Required:**
Choose ONE and document:
- **Option A (Recommended):** "Model was trained for 120 epochs to convergence. Here are the ELBO curves/convergence metrics showing training quality."
- **Option B:** "Model was trained to epoch 10, then re-generated final outputs. Epoch 120 checkpoints are from an exploratory run and are not used in the paper."

---

### 1.2 HIGH: Incomplete Convergence Documentation

**Finding:**  
The manuscript claims "ELBO curve had stabilized" at epoch 10, but no convergence metrics are provided.

**Requirements NOT Met:**
- No ELBO loss trajectory plot in outputs/integration/
- No statement of "ELBO plateau threshold" (e.g., "decreased <0.1% per epoch")
- No comparison: epoch 10 ELBO vs. final ELBO vs. epoch 120 ELBO
- No discussion of divergence risk between epoch 10 and full training

**What's Missing:**
```python
# Should be documented in 05_integration.py output:
- ELBO at epoch 10: [VALUE]
- ELBO at epoch 120: [VALUE]
- Inference time (epoch 10 model): [TIME]
- Inference time (epoch 120 model): [TIME]
- Cluster stability (ARI): epoch 10 vs 120
```

**How Hostile Reviewer Would Attack:**
> "Authors claim ELBO 'stabilized' at epoch 10 but provide no convergence metrics, learning curves, or validation that epoch 10 yields equivalent results to full training. This is not a limitation; it's an incomplete analysis. We cannot assess whether the model quality is adequate."

**Fix Required:**
1. Generate ELBO curve plot (already computed during training)
2. Document: "ELBO decreased [X]% from epoch 10 to 120, [interpret stability]"
3. If epoch 10 ≠ full training: "Cluster stability metrics (ARI=[value]) confirm epoch 10 captures primary structure"

---

### 1.3 HIGH: SPP1+ Macrophage Frequency Reporting Inconsistency

**Finding:**  
The manuscript states two different numbers for SPP1+ macrophage frequency that are mathematically inconsistent:

```
Claim 1 (Abstract/Results):
"SPP1+ macrophages represented 7.3% of myeloid cells in fast progressors 
 versus 16.9% in slow progressors"

Claim 2 (Multi-cohort section):
"SPP1+ present in all 5 scRNA-seq cohorts: Korean 12.1%, Kumar 18.3%, 
 T cell exh 22.4%, Zhang 15.7%, diffuse GC 19.8%"
```

**The Discrepancy:**
- If Korean = 12.1%, how does it split into 7.3% (fast) vs. 16.9% (slow)?
- These numbers don't add up: Korean sample is n=33 (16 slow, 17 fast), so weighted average should be ~(16×16.9 + 17×7.3)/(33) = 11.8%, which is close to 12.1% ✓
- **BUT:** This means the manuscript is averaging across progression groups without clearly stating it in the multi-cohort section.

**How This Confuses Readers:**
The multi-cohort section presents raw percentages without noting these are mixed (slow + fast) populations, making cross-cohort comparisons appear uncontrolled.

**Fix Required:**
Clarify in multi-cohort section:
> "SPP1+ macrophage fraction (combined across progression groups) in each cohort: Korean 12.1% (composed of 16.9% in slow, 7.3% in fast groups), Kumar 18.3%, ..."

OR

> "SPP1+ stratified by progression in Korean cohort (primary): fast progressors 7.3%, slow progressors 16.9%. Frequencies in external cohorts (combined across cohorts, not stratified by progression): [values]"

---

### 1.4 MEDIUM: CAF Signature Definition Missing from Methods

**Finding:**  
Methods section states CAF signature is "45-gene signature from Öhlund et al. (Nature Reviews Cancer 2017)" but:
- **No gene list provided** in Methods or Supplementary
- **Not reproducible** without the exact 45 genes
- **Cannot verify** which genes drive the effect

**Location of Issue:**  
`MANUSCRIPT_FINAL.tex`, Methods, "Cell State Scoring":
> "CAF: 45-gene signature (Öhlund et al. 2017) including FAP, ACTA2, POSTN, COL1A1, etc."

The "etc." is a red flag — exact signature must be listed.

**How Hostile Reviewer Would Attack:**
> "Authors do not provide the CAF signature gene list. Reproducibility demands the exact 45 genes used. Did they use Öhlund's original list? Cherry-pick a subset? How can results be reproduced?"

**Fix Required:**
1. Add Supplementary Table listing all 45 CAF genes
2. State source explicitly: "CAF signature (Öhlund et al. 2017 Supplementary Table X)" or custom-derived
3. Include in Methods: "CAF genes: [list]"

---

### 1.5 MEDIUM: Doublet Flagging But Not Removal — Unresolved Choice

**Finding:**  
Scrublet detected 4,161 doublets (0.97% of cells), which were **flagged but retained** for analysis.

**Problem:**  
- No sensitivity analysis: What changes if doublets are removed?
- No justification: Why retain doublets if they can bias trajectory/communication analysis?
- Inconsistent with field practice: Most scRNA-seq papers remove flagged doublets

**Methods Statement:**  
> "Doublet detection via Scrublet identified 4,161 doublets (0.97%), which were flagged but retained for flexibility."

**Questions This Raises:**
- Flexibility for what? Not stated.
- Does retaining doublets affect LIANA cell-cell communication results (which can be sensitive to false interactions)?
- Does it affect trajectory/subclustering?

**Fix Required:**
Choose and justify:
- **Option A:** "Doublets were flagged and removed, retaining 425,706 cells (99.03%)"
- **Option B:** "Doublets were flagged but retained. Sensitivity analysis (Supplement Figure) shows [metric] is stable whether doublets included or excluded, confirming robustness."

---

## SECTION 2: STRUCTURAL & NARRATIVE COHESION

### 2.1 MEDIUM: Introduction Claims Not Addressed in Results

**Finding:**  
The Introduction sets up specific clinical questions that are not resolved in Results.

**Orphaned Argument:**
```
Introduction (Intro 3):
"Current approaches to predict immunotherapy response remain inadequate. 
PD-L1 CPS is used to guide treatment, yet its predictive value is modest 
(AUC 0.55-0.65). Many PD-L1 positive tumors fail to respond to checkpoint 
inhibition, while conversely, some PD-L1 negative tumors achieve durable 
responses..."

Results Section:
[Never compares composite score directly to PD-L1 CPS]
[Never shows: AUC(composite) vs AUC(PD-L1)]
[Never states whether composite score predicts response better than PD-L1]
```

**Impact:**  
Reader left hanging: "Okay, your score has AUC=0.71, but is that better than PD-L1's 0.55-0.65?"

**Fix Required:**
Add to Results:
> "Composite score achieved AUC=0.71 (95% CI 0.64-0.78) in external TCGA validation, compared to PD-L1 CPS's reported AUC=0.60-0.65, indicating improved predictive value."

OR explicitly address:
> "While PD-L1 CPS alone shows AUC~0.60, our composite score incorporating TME features achieves AUC=0.71, suggesting TME profiling adds independent prognostic value."

---

### 2.2 MEDIUM: "Multi-Cohort Validation" Language Overstatement

**Finding:**  
The paper uses "validation" and "replication" language inconsistently for bulk cohorts.

**Terminology Issue:**
- **scRNA-seq cohorts (5):** Truly "integration" — all profiled with same method
- **Bulk cohorts (TCGA, GEO):** NOT true replication — signature estimated via ssGSEA (not direct scRNA-seq)

**Misleading Statement (Results):**
> "This SPP1+ signature replicated in all five scRNA-seq cohorts... and validated in bulk RNA-seq (GSE26253: HR=3.44, p=0.004)"

**Problem:**  
The word "replicated" is correct for scRNA-seq (same platform, same metrics). But "validated in bulk" is weaker — bulk RNA-seq estimates the signature via ssGSEA, which is NOT the same as directly observing SPP1+ macrophages.

**How Reviewer Would Critique:**
> "Authors conflate true replication (5 scRNA-seq cohorts with direct cell-type observation) with signature-based projection onto bulk data (different modality, estimated signature). These are different types of evidence and should be labeled accordingly."

**Fix Required:**
Revise language:
> "SPP1+ signature was replicated across all five scRNA-seq cohorts (direct observation). To extend to bulk RNA-seq, we estimated the signature via ssGSEA in TCGA and GEO cohorts, which showed consistent prognostic value (GSE26253: HR=3.44, p=0.004)."

---

### 2.3 LOW: Discussion Mechanisms Not Tied to LIANA Results

**Finding:**  
Discussion proposes mechanisms (SPP1-integrin, TGF-β, IL-10) but doesn't clearly tie them to LIANA findings.

**Missing Link:**
```
Results show:
"LIANA identified SPP1-integrin (score 0.85) as top interaction"
"Other top interactions: TGFB1-TGFBR2 (0.78), IL10-IL10RA (0.72)"

Discussion claims:
"Mechanisms: (1) Direct integrin signaling via SPP1-ITGAV"
"(2) Paracrine immunosuppression via TGF-β and IL-10"

Missing: Explicit connection
"Our LIANA analysis identified [specific interactions], providing direct evidence for..."
```

**Fix Required:**
In Discussion Mechanisms, add:
> "Our cell-cell communication analysis (LIANA) identified SPP1-integrin as the top-scoring ligand-receptor pair (score 0.85, p<0.001) between SPP1+ macrophages and CD8+ T cells, providing direct evidence for mechanism (1). Additional enriched interactions included TGFB1-TGFBR2 (score 0.78) and IL10-IL10RA (0.72), supporting mechanisms (2) and (3)."

---

## SECTION 3: LITERATURE & CITATION GAPS

### 3.1 MEDIUM: Overstatement of "SPP1+ Discovery" Novelty

**Finding:**  
The manuscript claims SPP1+ macrophages as a novel finding, but does not cite prior literature on SPP1 in tumor-associated macrophages.

**Missing Context:**
- Is SPP1+ enrichment in TAMs known in other cancers? (Likely yes)
- Is SPP1's role in immunosuppression established? (Yes — cited: Rangaswami 2014)
- **What is NOVEL:** SPP1+ association with immunotherapy *resistance specifically in gastric cancer*

**Current Language:**
> "Our finding identifies SPP1+ macrophages as a key immunosuppressive population."

**How Reviewer Sees It:**
> "SPP1 in macrophages is well-established (Rangaswami et al., Cao et al.). The novelty here is specific to gastric cancer + immunotherapy response, but the manuscript over-claims novelty."

**Fix Required:**
Revise to be precise about novelty:
> "While SPP1 expression in tumor-associated macrophages has been documented across cancer types, its specific association with immunotherapy response in gastric cancer and co-expression with T cell exhaustion has not been previously characterized in a meta-analysis of this scale."

---

### 3.2 MEDIUM: Insufficient Citation of scRNA-seq Integration Methods

**Finding:**  
The Methods section cites scVI but does not cite key papers on scVI best practices and batch correction validation.

**Current Citation:**
> "scVI (Single-Cell Variational Inference, scvi-tools v0.20)"

**Missing:**
- Kang et al. (2017) or Tian et al. (2019) on batch correction validation
- No citation for why scVI was chosen over other methods (Seurat v3, Harmony, etc.)
- No justification for epoch 10 vs. full training (if epoch 10 truly used)

**Impact:**  
Readers cannot understand *why* this particular integration approach was selected or *how its quality was validated*.

**Fix Required:**
Add to Methods:
> "scVI (Single-Cell Variational Inference, scvi-tools v0.20; Lopez et al. 2018) was selected over alternative integration methods (Seurat v3, Harmony) due to its probabilistic framework, which is particularly suited to capturing gene-gene covariation across batches. Batch correction quality was assessed via [metric: LISI, Silhouette, or other]."

---

### 3.3 LOW: CAF/TAM Literature Slightly Dated

**Finding:**  
Some key references (Kalluri 2016, Ohlund 2017) are now >8 years old for a 2026 submission.

**Risk:**  
Competitors may cite more recent reviews (2022-2025) on CAF immunotherapy, making the background feel outdated.

**Recommendation:**  
Add 2-3 recent CAF papers (Google Scholar: "CAF immunotherapy 2024") to update context.

---

## SECTION 4: RISK & WEAKNESS ASSESSMENT

### 4.1 TOP 3 VULNERABILITIES (Hostile Reviewer Perspective)

#### **Vulnerability #1: scVI Training Misrepresentation (CRITICAL)**

**What Reviewer Will Say:**
> "The manuscript claims scVI training was interrupted at epoch 10 and uses a checkpoint from that point, but file evidence shows the model trained to epoch 120 and checkpoints were saved every 10 epochs. The final integrated object was generated AFTER this full training completed. This is a material misrepresentation of the methods. Either the authors' model IS fully trained (contradicting their limitation), or the authors performed training but chose to report using an incomplete checkpoint (unethical and undisclosed). Either way, the Methods section lacks credibility."

**Damage:** HIGH  
**Fix Priority:** IMMEDIATE

---

#### **Vulnerability #2: No Comparison to Standard Biomarkers (HIGH)**

**What Reviewer Will Say:**
> "The authors propose a composite score with AUC=0.71 as superior to current methods, but never compare it directly to PD-L1 CPS (the standard-of-care biomarker). The Introduction criticizes PD-L1's AUC=0.55-0.65, but Results never provide the head-to-head comparison: Is 0.71 materially better? In the same Korean cohort, what is PD-L1 CPS's AUC? Without this, the clinical value proposition is unsubstantiated."

**Damage:** MEDIUM  
**Fix Priority:** HIGH

---

#### **Vulnerability #3: Doublet Handling Unjustified (HIGH)**

**What Reviewer Will Say:**
> "The authors detected and flagged 4,161 doublets but retained them, citing 'flexibility' without definition. Doublets can introduce false cell-cell interactions (LIANA analysis), spurious subclusters, and trajectory artifacts. No sensitivity analysis demonstrates that results are robust to doublet inclusion. This is a methodological red flag."

**Damage:** MEDIUM  
**Fix Priority:** HIGH

---

### 4.2 ADDITIONAL VULNERABILITIES

| Vulnerability | Risk Level | Fix Effort | Impact |
|---|---|---|---|
| CAF signature gene list not provided | MEDIUM | Low (1 table) | Reproducibility blocking |
| LIANA mechanisms not explicitly tied to results | LOW | Low (1 paragraph) | Clarity only |
| SPP1+ novelty overstated relative to literature | MEDIUM | Low (rephrase) | Credibility |
| No PD-L1 CPS comparison in Korean cohort | HIGH | Medium (reanalysis) | Clinical relevance |
| Bulk cohort terminology ("validation" vs "projection") | MEDIUM | Low (rephrase) | Precision |
| Convergence metrics (ELBO) not documented | HIGH | Medium (generate plots) | Methods credibility |

---

## SECTION 5: ACTIONABLE REMEDIATION PLAN

### MUST-FIX BEFORE SUBMISSION (Blocking Issues)

**1. scVI Training Discrepancy (CRITICAL)**
- **Action:** Clarify whether full training (epoch 120) was performed
  - If YES: Update Methods: "Model was trained for 120 epochs to convergence"
  - If NO: Explain why epoch 120 checkpoints exist if training stopped at epoch 10
- **Timeline:** 1 hour (requires decision only)
- **Evidence to add:** ELBO plot, convergence metrics

**2. PD-L1 CPS Comparison (HIGH)**
- **Action:** Re-analyze Korean cohort to compute PD-L1 CPS AUC in same n=33 cohort
- **Steps:**
  1. Extract PD-L1 CPS from Korean cohort metadata
  2. Stratify patients by high/low CPS
  3. Run Cox regression: OS ~ PD-L1 CPS
  4. Calculate AUC using same Harrell's C-index approach
  5. Compare to composite score AUC
- **Timeline:** 2-3 hours
- **Add to Results:** "PD-L1 CPS alone achieved AUC=X (C-index=Y) vs. composite score AUC=0.71 (C-index=0.72)"

**3. Doublet Analysis (HIGH)**
- **Action:** Run sensitivity analysis
- **Steps:**
  1. Re-run key analyses excluding doublets
  2. Compare metrics: cluster stability (ARI), cell-type proportions, survival HRs
  3. Document: "SPP1+ HR with doublets=2.1 [1.2-3.7], without doublets=[value]"
- **Timeline:** 3-4 hours
- **Add to Supplement:** Sensitivity table showing results are robust

---

### SHOULD-FIX (High Value, Medium Effort)

**4. CAF Signature Documentation (MEDIUM)**
- **Action:** List all 45 CAF genes in Supplementary Table
- **Timeline:** 1 hour
- **Cite:** Öhlund et al. 2017 (or your source)

**5. Convergence Metrics (MEDIUM)**
- **Action:** Generate and include ELBO learning curve
- **Timeline:** 1 hour (re-run scVI or extract from logs)
- **Add to Supplement Figure:** ELBO vs. epoch plot

**6. LIANA-Mechanisms Link (MEDIUM)**
- **Action:** Revise Discussion to explicitly tie LIANA results to proposed mechanisms
- **Timeline:** 30 minutes
- **Add:** "Our LIANA analysis identified [interactions], providing direct evidence for..."

**7. PD-L1 Terminology Clarity (MEDIUM)**
- **Action:** Update Introduction to note that composite score includes PD-L1 CPS
- **Timeline:** 30 minutes
- **Add:** "Unlike PD-L1 CPS alone, our composite score incorporates TME features..."

---

### NICE-TO-HAVE (Polish, Low Effort)

**8. SPP1+ Novelty Precision (LOW)**
- **Action:** Rephrase novelty claims to acknowledge prior SPP1 literature
- **Timeline:** 15 minutes

**9. Bulk Cohort Terminology (LOW)**
- **Action:** Replace "validated in bulk" with "projected to bulk RNA-seq via ssGSEA"
- **Timeline:** 15 minutes

**10. Recent References (LOW)**
- **Action:** Add 2-3 recent CAF papers (2023-2025)
- **Timeline:** 30 minutes

---

## SECTION 6: PRIORITY ACTION ITEM LIST

### RANKED BY URGENCY & IMPACT

```
[BLOCKER] Fix scVI Training Documentation
├─ Decision Required: Full training (120 epochs) or interrupted (10 epochs)?
├─ Evidence: Checkpoints up to ep120, final model June 16 16:17
├─ Timeline: 1 hour (decision + documentation)
└─ Impact: Methods credibility — reviewers will flag this immediately

[BLOCKER] Add PD-L1 CPS Comparison to Korean Cohort
├─ Action: Calculate AUC(PD-L1) vs AUC(composite) in same n=33 cohort
├─ Timeline: 2-3 hours
└─ Impact: Justifies clinical value proposition

[BLOCKER] Doublet Sensitivity Analysis
├─ Action: Show SPP1+ HR, cluster stability, etc. stable if doublets removed
├─ Timeline: 3-4 hours
└─ Impact: Addresses major methodological concern

[HIGH] CAF Gene List in Supplement
├─ Action: List all 45 genes, cite Öhlund et al.
├─ Timeline: 1 hour
└─ Impact: Reproducibility

[HIGH] ELBO Convergence Plot
├─ Action: Generate learning curve, add to Supplement
├─ Timeline: 1 hour
└─ Impact: Supports convergence claims

[HIGH] Link LIANA to Mechanisms
├─ Action: Revise Discussion with "Our LIANA analysis identified..."
├─ Timeline: 30 minutes
└─ Impact: Logical flow

[MEDIUM] Terminology Precision (Validation vs. Projection)
├─ Action: Update bulk cohort language
├─ Timeline: 15 minutes
└─ Impact: Precision/clarity

[MEDIUM] Update References (Add 2023-2025 CAF papers)
├─ Action: Add 2-3 recent papers to Literature
├─ Timeline: 30 minutes
└─ Impact: Recency
```

---

## SECTION 7: SUMMARY & RECOMMENDATION

### Project Status Post-Audit

**Strengths:**
✅ Solid multi-cohort design (5 scRNA-seq + 3 bulk + TCGA)  
✅ Comprehensive analyses (integration, annotation, communication, trajectory)  
✅ Consistent findings across cohorts (SPP1+, CAF replication)  
✅ Strong external validation (HR=3.44 in GSE26253)  
✅ High-impact target (Gut journal submission)

**Critical Issues Requiring Immediate Fix:**
⚠️ scVI training state misrepresented (epoch 10 vs. actual epoch 120+ training)  
⚠️ No direct PD-L1 CPS comparison to justify biomarker superiority  
⚠️ Doublet methodology unjustified (flagged but not removed)

**Recommendation:**
**DO NOT SUBMIT** until:
1. scVI training discrepancy is resolved and documented
2. PD-L1 CPS comparison is added to Results
3. Doublet sensitivity analysis is performed and included in Supplement

**Estimated Time to Fix All Blockers:** 8-10 hours total work  
**Estimated Time to Polish (Nice-to-have):** 2-3 additional hours

**Post-Remediation Outlook:** Project will be **publication-ready** with strong methods transparency, direct clinical relevance comparison, and robust sensitivity analyses.

---

## APPENDIX: DETAILED TECHNICAL FINDINGS

### A1. scVI Checkpoint Timeline (Evidence)
```
Created:  2026-06-15 23:42 — scvi_model_ckpt_ep010
Created:  2026-06-16 01:17 — scvi_model_ckpt_ep020
Created:  2026-06-16 02:26 — scvi_model_ckpt_ep030
...
Created:  2026-06-16 13:48 — scvi_model_ckpt_ep120
Created:  2026-06-16 16:17 — scvi_model/model.pt (FINAL MODEL)
Generated: 2026-06-17 00:46 — gastric_meta_integrated.h5ad
```

### A2. Code Evidence (05_integration.py)
- Training configured for up to 400 epochs (line 351-352)
- Checkpoint callback saves every 10 epochs (line 300)
- Option to load from checkpoint exists (line use_checkpoint parameter)
- Main execution does NOT show --use-checkpoint flag in documented runs

### A3. Data Integrity Checks
```
✓ Korean cohort: 654,770 raw → 429,867 QC'd (34.3% filtered) ✓
✓ External cohorts: 5 datasets, 337,978 cells ✓
✓ Integrated: 766,845 cells, 17 Leiden clusters ✓
✓ Cell type fractions: T/NK 32.1%, Epithelial 28.7%, Myeloid 18.5% ✓
✓ SPP1+ frequencies across cohorts: 12.1%-22.4% (consistent) ✓
```

---

**END OF AUDIT REPORT**

---

## IMMEDIATE NEXT STEPS FOR PROJECT LEAD

1. **TODAY:** Schedule 30-min clarification meeting with data analysis team
   - Clarify: Was scVI trained to epoch 10 or epoch 120?
   - Reconcile: Why do epoch 120 checkpoints exist?
   - Decide: Update manuscript or provide full training metrics?

2. **THIS WEEK:** Execute blocker fixes (estimated 8-10 hours)
   - PD-L1 comparison analysis (coordinate with biostatistician)
   - Doublet sensitivity re-run (data analyst)
   - scVI documentation and ELBO plot (data analyst)

3. **NEXT WEEK:** Polish pass (estimated 2-3 hours)
   - CAF gene list curation
   - Terminology precision in Discussion/Methods
   - Recent reference additions

4. **FINAL REVIEW:** Before journal submission
   - Have uninvolved senior scientist read Methods section
   - Verify all claims are supported by data
   - Run hostile-reviewer exercise on Discussion

---

**Report prepared by:** Audit System  
**Verification status:** Grounded in notebook timestamps and data files  
**Confidence level:** HIGH for scVI training discrepancy, MEDIUM-HIGH for other findings
