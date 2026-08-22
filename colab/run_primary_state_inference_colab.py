#!/usr/bin/env python3
"""Resumable Colab driver for primary working-state inference.

The driver rebuilds native assignments one cohort at a time, exact-checks each
against the reviewed discovery partition, then calculates sample-level state
composition and cohort-wise/random-effects tumour-normal contrasts.  It never
uses excluded or secondary candidate clusters in primary testing.
"""

from __future__ import annotations

import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


CORE_IDS = [
    "korea_kim2022", "kumar2022", "sathe2020", "zhang2021",
    "diffuse_gc_2021", "tcell_exhaustion_2022",
]


def run(*command: str) -> None:
    print("+", " ".join(command), flush=True)
    subprocess.run(command, check=True)


def find_discovery(project_root: Path) -> Path:
    candidates = []
    for summary in project_root.rglob("NATIVE_PER_COHORT_CLUSTER_SUMMARY.csv"):
        candidate = summary.parent
        if (candidate / "NATIVE_PER_COHORT_CLUSTER_MARKERS.csv").exists():
            candidates.append(candidate)
    if not candidates:
        raise FileNotFoundError(
            "No completed native discovery directory found. Expected both "
            "NATIVE_PER_COHORT_CLUSTER_SUMMARY.csv and NATIVE_PER_COHORT_CLUSTER_MARKERS.csv under the project folder."
        )
    return max(candidates, key=lambda path: path.stat().st_mtime)


def main() -> None:
    code_root = Path(__file__).resolve().parents[1]
    project_root = Path(os.environ.get("GASTRIC_TME_PROJECT_ROOT", "/content/drive/MyDrive/gastric_tme_project"))
    data_root = Path(os.environ.get("GASTRIC_TME_DATA_ROOT", "/content/drive/MyDrive/data"))
    discovery = Path(os.environ["GASTRIC_TME_NATIVE_DISCOVERY"]) if os.environ.get("GASTRIC_TME_NATIVE_DISCOVERY") else find_discovery(project_root)
    manifest = data_root / "external" / "dataset_manifest.csv"
    reference_atlas = data_root / "processed" / "integrated" / "gastric_meta_annotated.h5ad"
    required = [manifest, reference_atlas, discovery / "NATIVE_PER_COHORT_CLUSTER_SUMMARY.csv", discovery / "NATIVE_PER_COHORT_CLUSTER_MARKERS.csv"]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required input(s):\n" + "\n".join(missing))

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_root = Path(os.environ.get("GASTRIC_TME_PRIMARY_STATE_RUN", str(project_root / "PRIMARY_WORKING_STATE_INFERENCE" / stamp)))
    run_root.mkdir(parents=True, exist_ok=True)
    audit_dir = run_root / "01_LINEAGE_AUDIT"
    provisional_dir = run_root / "02_PROVISIONAL_CURATION"
    resolution_dir = run_root / "03_ALL_CLUSTER_RESOLUTION"
    assignments_dir = run_root / "04_RECREATED_ASSIGNMENTS"
    combined_dir = run_root / "05_COMBINED_ASSIGNMENTS"
    composition_dir = run_root / "06_SAMPLE_LEVEL_COMPOSITION"
    meta_dir = run_root / "07_STATE_ABUNDANCE_META"

    print(f"Code: {code_root}")
    print(f"Project: {project_root}")
    print(f"Data: {data_root}")
    print(f"Discovery: {discovery}")
    print(f"Output: {run_root}")

    run(sys.executable, str(code_root / "native_cluster_lineage_audit.py"), "--discovery-dir", str(discovery), "--output-dir", str(audit_dir))
    run(sys.executable, str(code_root / "computational_provisional_native_curation.py"), "--audit", str(audit_dir / "NATIVE_CLUSTER_LINEAGE_AUDIT.csv"), "--markers", str(discovery / "NATIVE_PER_COHORT_CLUSTER_MARKERS.csv"), "--output-dir", str(provisional_dir))
    run(sys.executable, str(code_root / "resolve_all_native_cluster_dispositions.py"), "--dictionary", str(provisional_dir / "COMPUTATIONAL_PROVISIONAL_STATE_DICTIONARY.csv"), "--output-dir", str(resolution_dir))

    for cohort in CORE_IDS:
        checkpoint = assignments_dir / f"NATIVE_PER_COHORT_CELL_ASSIGNMENTS_{cohort}.csv.gz"
        if checkpoint.exists():
            print(f"Checkpoint exists; skipping recreation: {cohort}", flush=True)
            continue
        run(
            sys.executable, str(code_root / "recreate_native_cell_assignments.py"),
            "--manifest", str(manifest), "--data-root", str(data_root),
            "--reference-atlas", str(reference_atlas), "--discovery-dir", str(discovery),
            "--output-dir", str(assignments_dir), "--cohorts", cohort,
        )

    run(sys.executable, str(code_root / "combine_recreated_native_assignments.py"), "--assignment-dir", str(assignments_dir), "--discovery-dir", str(discovery), "--output-dir", str(combined_dir))
    run(sys.executable, str(code_root / "build_primary_working_state_composition.py"), "--assignments", str(combined_dir / "NATIVE_PER_COHORT_CELL_ASSIGNMENTS_ALL_COHORTS.csv.gz"), "--resolution-ledger", str(resolution_dir / "ALL_505_NATIVE_CLUSTER_RESOLUTION_LEDGER.csv"), "--output-dir", str(composition_dir))
    run(sys.executable, str(code_root / "state_abundance_meta_analysis.py"), "--curated-composition", str(composition_dir / "PRIMARY_WORKING_STATE_SAMPLE_COMPOSITION.csv"), "--output-dir", str(meta_dir))
    print(f"Completed primary working-state inference: {run_root}")


if __name__ == "__main__":
    main()
