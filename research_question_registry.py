#!/usr/bin/env python3
"""Create a data-linked, five-question registry for the manuscript analysis.

The registry is intentionally an analysis plan grounded in completed native
per-cohort discovery output. It does not convert provisional marker patterns
into biological claims. Every question names its eligible design, unit of
inference, required curation gate, and stopping rule.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


QUESTIONS = [
    {
        "question_id": "RQ1",
        "research_question": "Which reviewed cell states are reproducibly altered between tumour and non-tumour gastric tissue?",
        "rationale_from_current_data": "Native raw-count clustering produced candidate states in every major compartment across multiple cohorts; tumour-normal inference is possible only in cohorts with auditable tissue/sample metadata.",
        "eligible_design": "Paired tumour-normal where donor IDs are auditable; otherwise cohort-level descriptive effect only. Zhang is excluded.",
        "primary_compartments": "Myeloid;Fibroblast;Endothelial;T/NK;B/Plasma;Epithelial",
        "inferential_unit": "donor/sample",
        "primary_analysis": "Within-cohort state proportions, paired contrasts where available, random-effects pooling of matched contrast.",
        "required_gate": "Reviewed state name/rationale; >=3 samples per contributing condition; no excluded native cluster.",
        "advance_rule": "Report only states observed after harmonisation in >=3 cohorts with a defined matched contrast; otherwise descriptive/supplementary.",
        "claim_boundary": "State abundance differs by tissue context; no causal TME claim.",
    },
    {
        "question_id": "RQ2",
        "research_question": "Do recurrent myeloid macrophage/monocyte programs separate into portable versus context-specific tumour-associated states?",
        "rationale_from_current_data": "APOC1 is enriched in initially eligible myeloid clusters in Korea, Kumar, diffuse GC, and GSE206785; the prior broad SPP1 claim was heterogeneous.",
        "eligible_design": "Reviewed myeloid states in tumour-normal cohorts; donor-level pseudobulk programme scores.",
        "primary_compartments": "Myeloid",
        "inferential_unit": "donor/sample",
        "primary_analysis": "State abundance plus donor pseudobulk expression of C1QC/APOC1, SPP1/APOC1/TREM2/GPNMB, and IL1B/FCN1 programmes; random-effects and leave-one-cohort-out sensitivity.",
        "required_gate": "Manual approval of myeloid state labels and marker evidence in >=3 cohorts.",
        "advance_rule": "Call portable only with FDR <0.05, concordant direction in >=80% cohorts, leave-one-cohort-out stability, and prediction interval in one direction.",
        "claim_boundary": "Recurrent transcriptional state association, not universal macrophage programming or mechanism.",
    },
    {
        "question_id": "RQ3",
        "research_question": "Are matrix/inflammatory CAF and capillary/lymphatic endothelial states jointly altered in gastric tumour tissue?",
        "rationale_from_current_data": "Candidate fibroblast states occur in all six cohorts; candidate endothelial states occur in five. PDPN-positive fibroblast clusters remain manual-review candidates, while PLVAP-positive endothelial clusters are initially eligible in Kumar and Sathe.",
        "eligible_design": "Reviewed fibroblast/endothelial states in cohorts with auditable tumour-normal comparison.",
        "primary_compartments": "Fibroblast;Endothelial",
        "inferential_unit": "donor/sample",
        "primary_analysis": "Paired abundance contrasts and within-sample state correlation adjusted for cohort; spatial association only after single-cell state results are fixed.",
        "required_gate": "CAF and endothelial labels manually reviewed in >=3 cohorts; each state supported by >=3 samples per contributing condition.",
        "advance_rule": "Advance to spatial validation only for reproducible state changes, not for unreviewed ligand-receptor predictions.",
        "claim_boundary": "Co-occurrence/association, not physical interaction or causal signalling.",
    },
    {
        "question_id": "RQ4",
        "research_question": "How do reviewed T/NK states change across documented baseline and on-treatment Korean biopsies?",
        "rationale_from_current_data": "T/NK candidate clusters occur in all six cohorts. The Korean An/Mehta cohort is the only serial frontline chemoimmunotherapy cohort in the core dataset.",
        "eligible_design": "Korean within-patient baseline-to-FU1 comparison; FU2 secondary exploratory; response groups only after source RECIST audit.",
        "primary_compartments": "T/NK;Myeloid;Fibroblast;Endothelial",
        "inferential_unit": "patient",
        "primary_analysis": "Within-patient change in reviewed-state fraction and program score; paired permutation/Wilcoxon testing; FDR over prespecified states.",
        "required_gate": "Patient/timepoint mapping and response-label audit; no duplicated patient-timepoint-state rows.",
        "advance_rule": "Report longitudinal context regardless of response result; report response interaction only when patient-level labels and group sizes pass preflight.",
        "claim_boundary": "Treatment-associated change, not a treatment-response predictor or evidence of drug mechanism.",
    },
    {
        "question_id": "RQ5",
        "research_question": "Which epithelial states vary across the Zhang NAG-to-CAG-to-IM-to-early-GC continuum?",
        "rationale_from_current_data": "Zhang contains 18 initially eligible epithelial clusters and an audited disease continuum; its epithelial states should be analysed separately rather than pooled as tumour-normal.",
        "eligible_design": "Zhang GSE134520 only, after explicit stage/sample mapping; NAG, CAG, IM, EGC treated as ordered disease-context groups.",
        "primary_compartments": "Epithelial",
        "inferential_unit": "biopsy/sample",
        "primary_analysis": "Stage-stratified reviewed-state abundance and pre-specified trend tests; report individual biopsy values.",
        "required_gate": "Reviewed epithelial labels and verified stage metadata; no inference from sample names alone.",
        "advance_rule": "Exploratory panel only unless an independent premalignant cohort validates the same direction.",
        "claim_boundary": "Disease-continuum association within one cohort; no inferred lineage trajectory or causal progression claim.",
    },
]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--discovery-dir", type=Path, required=True)
    parser.add_argument("--audit-dir", type=Path, required=True)
    parser.add_argument("--context-registry", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    status = pd.read_csv(args.discovery_dir / "NATIVE_PER_COHORT_CLUSTER_STATUS.csv")
    audit = pd.read_csv(args.audit_dir / "NATIVE_CLUSTER_LINEAGE_AUDIT.csv")
    context = pd.read_csv(args.context_registry)
    questions = pd.DataFrame(QUESTIONS)
    questions.to_csv(args.output_dir / "MANUSCRIPT_RESEARCH_QUESTIONS.csv", index=False)

    eligibility = (
        audit.groupby(["broad_label", "curation_recommendation"], observed=True)
        .agg(n_clusters=("native_cluster", "size"), n_cells=("n_cells", "sum"), n_cohorts=("dataset_id", "nunique"))
        .reset_index()
    )
    eligibility.to_csv(args.output_dir / "CURRENT_NATIVE_STATE_EVIDENCE.csv", index=False)
    coverage = (
        status.groupby("dataset_id", observed=True)
        .agg(n_compartments=("broad_label", "nunique"), n_native_clusters=("n_native_clusters", "sum"), n_sampled_cells=("n_cells", "sum"))
        .reset_index()
        .merge(context, on="dataset_id", how="left", validate="one_to_one")
    )
    coverage.to_csv(args.output_dir / "COHORT_DESIGN_AND_DISCOVERY_COVERAGE.csv", index=False)

    gates = pd.DataFrame([
        {"gate": "Native state curation", "status": "pending_human_review", "requirement": "Reviewed state dictionary signed before harmonisation."},
        {"gate": "Cross-cohort harmonisation", "status": "blocked_by_curation", "requirement": "Same reviewed state in >=3 cohorts with documented marker rationale."},
        {"gate": "Donor/sample inference", "status": "blocked_by_cell_assignments_and_curation", "requirement": "Native cell assignments plus auditable donor/sample condition mapping."},
        {"gate": "Spatial validation", "status": "deferred", "requirement": "Only after a prespecified state association is reproducible at donor/sample level."},
        {"gate": "Manuscript figures", "status": "planned", "requirement": "Every panel generated from current tables and linked in claim ledger."},
    ])
    gates.to_csv(args.output_dir / "MANUSCRIPT_EXECUTION_GATES.csv", index=False)

    (args.output_dir / "README.md").write_text(
        "# Five-question manuscript registry\n\n"
        "These questions are derived from the completed six-cohort native state-discovery audit. They are not results. The analysis order is RQ1/RQ2/RQ3 after state curation, RQ4 after Korean clinical-label preflight, and RQ5 as a separate Zhang exploratory analysis.\n\n"
        "Do not open a spatial, interaction, survival, or manuscript-claim workstream until the relevant gate in `MANUSCRIPT_EXECUTION_GATES.csv` is satisfied.\n",
        encoding="utf-8",
    )
    (args.output_dir / "RUN_MANIFEST.json").write_text(json.dumps({
        "n_questions": len(QUESTIONS),
        "n_cohorts": int(status["dataset_id"].nunique()),
        "n_sampled_cells": int(status["n_cells"].sum()),
        "n_native_clusters": int(status["n_native_clusters"].sum()),
        "n_initially_eligible_clusters": int(audit["curation_recommendation"].eq("eligible_after_manual_state_annotation").sum()),
        "source_discovery_dir": str(args.discovery_dir),
        "source_audit_dir": str(args.audit_dir),
    }, indent=2) + "\n", encoding="utf-8")
    print(f"Research-question registry written to {args.output_dir}")


if __name__ == "__main__":
    main()
