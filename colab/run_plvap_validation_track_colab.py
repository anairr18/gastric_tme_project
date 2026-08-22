#!/usr/bin/env python3
"""Run the post-discovery validation track from a clean native-state run.

Stages: native-state reference profiles, marker-consistency review, optional
GSE251950 spatial mapping, and Korean patient-paired treatment context.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def run(*command: str) -> None:
    print("+", " ".join(command), flush=True)
    subprocess.run(command, check=True)


def latest_run(project: Path) -> Path:
    root = project / "REPRODUCIBLE_NATIVE_STATE_META"
    candidates = sorted((path for path in root.iterdir() if path.is_dir()), key=lambda path: path.stat().st_mtime, reverse=True)
    if not candidates:
        raise FileNotFoundError(f"No reproducible native-state run under {root}")
    return candidates[0]


def first_existing(candidates: list[Path], label: str) -> Path:
    found = next((path for path in candidates if path.exists()), None)
    if found is None:
        raise FileNotFoundError(f"Missing {label}. Checked:\n" + "\n".join(map(str, candidates)))
    return found


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, default=Path(os.environ.get("GASTRIC_TME_PROJECT", "/content/drive/MyDrive/gastric_tme_project")))
    parser.add_argument("--data-root", type=Path, default=Path(os.environ.get("GASTRIC_TME_DATA", "/content/drive/MyDrive/data")))
    parser.add_argument("--run", type=Path, help="Specific REPRODUCIBLE_NATIVE_STATE_META run; defaults to latest.")
    parser.add_argument("--download-spatial", action="store_true", help="Download and index GSE251950 before spatial mapping.")
    parser.add_argument("--spatial-permutations", type=int, default=499)
    parser.add_argument("--treatment-permutations", type=int, default=10000)
    args = parser.parse_args()
    base = Path(__file__).resolve().parents[1]
    run_root = args.run or latest_run(args.project)
    output = run_root / "07_POST_DISCOVERY_VALIDATION_TRACK"
    output.mkdir(parents=True, exist_ok=True)

    assignments = first_existing(list((run_root / "01_NATIVE_PER_COHORT_DISCOVERY").glob("NATIVE_PER_COHORT_CELL_ASSIGNMENTS.csv.gz")), "native cell assignments")
    ledger = next((path for path in (run_root / "04_CLUSTER_RESOLUTION").glob("ALL_*_NATIVE_CLUSTER_RESOLUTION_LEDGER.csv")), None)
    if ledger is None:
        raise FileNotFoundError("Missing native-cluster resolution ledger.")
    manifest = first_existing([args.data_root / "external" / "dataset_manifest.csv", run_root / "dataset_manifest_runtime.csv"], "dataset manifest")

    profiles = output / "01_STATE_REFERENCE_PROFILES"
    run(sys.executable, str(base / "export_native_state_reference_profiles.py"), "--assignments", str(assignments), "--resolution-ledger", str(ledger), "--manifest", str(manifest), "--data-root", str(args.data_root), "--output-dir", str(profiles))

    marker = output / "02_STATE_MARKER_CONSISTENCY"
    run(sys.executable, str(base / "state_marker_replication.py"), "--cohort-profiles", str(profiles / "COHORT_NATIVE_STATE_EXPRESSION_PROFILES.csv.gz"), "--output-dir", str(marker))

    spatial = output / "03_GSE251950_SPATIAL"
    spatial.mkdir(parents=True, exist_ok=True)
    extensions = args.data_root / "external" / "validation_extensions"
    spatial_index = extensions / "GSE251950" / "SPATIAL_INPUT_INDEX.csv"
    spatial_status: dict[str, str] = {"status": "not_run", "reason": "spatial index not available"}
    try:
        if args.download_spatial and not spatial_index.exists():
            run(sys.executable, str(base / "download_gse251950_visium.py"), "--output-root", str(extensions))
        if spatial_index.exists():
            run(sys.executable, str(base / "gse251950_spatial_validation.py"), "--spatial-root", str(extensions), "--reference-profiles", str(profiles / "FROZEN_NATIVE_STATE_REFERENCE_PROFILES.csv.gz"), "--output-dir", str(spatial), "--n-permutations", str(args.spatial_permutations))
            spatial_status = {"status": "completed", "reason": "spatial mapping completed"}
    except Exception as error:
        spatial.mkdir(parents=True, exist_ok=True)
        spatial_status = {"status": "excluded_or_failed", "reason": f"{type(error).__name__}: {error}"}
    (spatial / "SPATIAL_STAGE_STATUS.json").write_text(json.dumps(spatial_status, indent=2) + "\n", encoding="utf-8")

    treatment = output / "04_KOREAN_TREATMENT_CONTEXT"
    clinical = first_existing([
        args.project / "outputs" / "IF15_SUBMISSION_RESCUE" / "PATIENT_LEVEL_FEATURES_REBUILT.csv",
        args.project / "PATIENT_LEVEL_FEATURES_REBUILT.csv",
    ], "audited Korean patient-level clinical table")
    composition = run_root / "05_SAMPLE_LEVEL_COMPOSITION" / "PRIMARY_WORKING_STATE_SAMPLE_COMPOSITION.csv"
    run(sys.executable, str(base / "korean_native_state_treatment_context.py"), "--composition", str(composition), "--clinical", str(clinical), "--output-dir", str(treatment), "--permutations", str(args.treatment_permutations))

    (output / "VALIDATION_TRACK_MANIFEST.json").write_text(json.dumps({
        "generated_at_utc": datetime.now(timezone.utc).isoformat(), "native_run": str(run_root),
        "state_definition": str(marker / "STATE_MARKER_REPLICATION_SUMMARY.csv"),
        "spatial": spatial_status, "korean_treatment": str(treatment),
        "claim_boundary": "State-definition consistency is not independent validation; GSE251950 validates spatial localization in tumour slides, not tumour-normal depletion; Korean treatment analysis is patient-paired context, not prediction.",
    }, indent=2) + "\n", encoding="utf-8")
    print(f"Post-discovery validation track completed: {output}")


if __name__ == "__main__":
    main()
