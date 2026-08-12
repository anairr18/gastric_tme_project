#!/usr/bin/env python3
"""Run the mentor-requested systematic single-cell audit in Google Colab.

No integration is retrained here. The runner first writes durable Drive outputs
for cohort context, raw-input/method eligibility, and cohort-resolved labels.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def find_one(paths: list[Path], label: str) -> Path:
    found = next((path for path in paths if path.exists()), None)
    if found is None:
        raise FileNotFoundError(f"Could not find {label}. Checked:\n" + "\n".join(map(str, paths)))
    return found


def run(*command: str) -> None:
    print("+", " ".join(command), flush=True)
    subprocess.run(command, check=True)


def main() -> None:
    code = Path(__file__).resolve().parents[1]
    data_root = Path(os.environ.get("GASTRIC_TME_DATA_ROOT", "/content/drive/MyDrive/data"))
    project = Path(os.environ.get("GASTRIC_TME_PROJECT_ROOT", "/content/drive/MyDrive/gastric_tme_project"))
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    output = Path(os.environ.get("GASTRIC_TME_MENTOR_AUDIT_OUTPUT", str(project / "MENTOR_SYSTEMATIC_SINGLE_CELL_AUDIT" / timestamp)))
    output.mkdir(parents=True, exist_ok=True)
    manifest = find_one([data_root / "external" / "dataset_manifest.csv", code / "colab_transfer" / "dataset_manifest.csv"], "dataset manifest")
    atlas = find_one([
        Path(os.environ["GASTRIC_TME_ATLAS"]) if "GASTRIC_TME_ATLAS" in os.environ else project / "colab_rawcounts_sixcohort" / "gastric_meta_annotated_rawcounts_sixcohort.h5ad",
        project / "colab_rawcounts_sixcohort" / "gastric_meta_integrated_rawcounts_sixcohort.h5ad",
    ], "frozen annotated atlas")
    context_registry = code / "cohort_context_registry.csv"
    if not context_registry.exists():
        raise FileNotFoundError(f"Cohort context registry missing from code: {context_registry}")
    state_candidates = [
        Path(os.environ["GASTRIC_TME_STATE_LABELS"]) if "GASTRIC_TME_STATE_LABELS" in os.environ else project / "SUBMISSION_GRADE_STATE_META" / "CURATED_STATE_RESULTS" / "FROZEN_CURATED_STATE_CELL_LABELS.csv.gz",
        project / "SINGLE_CELL_ANALYSIS_AUTOMATED_COMPLETE" / "FROZEN_CURATED_STATE_CELL_LABELS.csv.gz",
    ]
    systematic_root = project / "SINGLE_CELL_ANALYSIS_AUTOMATED_COMPLETE"
    if systematic_root.exists():
        state_candidates.extend(sorted(systematic_root.rglob("FROZEN_CURATED_STATE_CELL_LABELS.csv.gz"), key=lambda path: path.stat().st_mtime, reverse=True))
    # The state audit is still useful at broad level when the frozen label table was not saved.
    state_labels = next((path for path in state_candidates if path.exists()), None)

    print(f"Code: {code}\nData: {data_root}\nAtlas: {atlas}\nOutput: {output}", flush=True)
    # Colab runtimes do not reliably retain AnnData after a restart. The audit
    # needs only this lightweight metadata/plotting stack; it does not retrain scVI.
    run(
        sys.executable,
        "-m",
        "pip",
        "install",
        "-q",
        "anndata>=0.10,<0.13",
        "h5py>=3.10",
        "scipy>=1.13",
        "matplotlib>=3.8",
        "seaborn>=0.13",
        "harmonypy>=0.0.10",
    )
    context_out = output / "01_COHORT_CONTEXT"
    raw_out = output / "02_RAW_INPUT_AND_METHOD_AUDIT"
    annotation_out = output / "03_PER_COHORT_ANNOTATION_AUDIT"
    run(sys.executable, str(code / "cohort_context_audit.py"), "--registry", str(context_registry), "--manifest", str(manifest), "--atlas", str(atlas), "--output-dir", str(context_out))
    run(sys.executable, str(code / "raw_integration_method_audit.py"), "--manifest", str(manifest), "--data-root", str(data_root), "--output-dir", str(raw_out))
    annotation_command = [sys.executable, str(code / "cohort_resolved_annotation_audit.py"), "--atlas", str(atlas), "--cohort-context", str(context_out / "COHORT_CONTEXT_FOR_SLIDES.csv"), "--output-dir", str(annotation_out)]
    if state_labels:
        annotation_command += ["--state-labels", str(state_labels)]
    run(*annotation_command)
    shutil.copy2(context_registry, output / "cohort_context_registry.csv")
    (output / "README.txt").write_text(
        "Mentor-requested systematic single-cell audit.\n\n"
        "01_COHORT_CONTEXT: exact cohort context and counts for slide 3.\n"
        "02_RAW_INPUT_AND_METHOD_AUDIT: raw-count, sample-ID, gene-overlap gates and scVI/Harmony decision.\n"
        "03_PER_COHORT_ANNOTATION_AUDIT: same frozen taxonomy in every cohort, marker support, coverage, and data-grounded questions.\n\n"
        "No integration was retrained and no new biological claims are made by this runner.\n",
        encoding="utf-8",
    )
    print(f"Completed mentor systematic audit: {output}", flush=True)


if __name__ == "__main__":
    main()
