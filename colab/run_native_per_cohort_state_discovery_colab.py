#!/usr/bin/env python3
"""Run native per-cohort state discovery and create its manual review package.

Designed for a fresh Colab runtime after ``git clone``/``git pull``. All
large inputs and outputs remain on Google Drive; only code lives in Git.
"""

from __future__ import annotations

import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


def run(*command: str) -> None:
    print("+", " ".join(command), flush=True)
    subprocess.run(command, check=True)


def latest(paths: list[Path], label: str) -> Path:
    existing = [path for path in paths if path.exists()]
    if not existing:
        raise FileNotFoundError(f"Could not find {label}. Checked:\n" + "\n".join(map(str, paths)))
    return max(existing, key=lambda path: path.stat().st_mtime)


def main() -> None:
    code = Path(__file__).resolve().parents[1]
    data = Path(os.environ.get("GASTRIC_TME_DATA_ROOT", "/content/drive/MyDrive/data"))
    project = Path(os.environ.get("GASTRIC_TME_PROJECT_ROOT", "/content/drive/MyDrive/gastric_tme_project"))
    recovery_root = project / "MENTOR_SYSTEMATIC_SINGLE_CELL_AUDIT"
    recovery_runs = list(recovery_root.glob("*KOREA_ID_AND_STATE_RECOVERY*")) if recovery_root.exists() else []
    recovery = latest(recovery_runs, "Korean identity-recovery run")
    manifest = recovery / "dataset_manifest_korea_id_repaired.csv"
    atlas = project / "colab_rawcounts_sixcohort" / "gastric_meta_annotated_rawcounts_sixcohort.h5ad"
    if not manifest.exists() or not atlas.exists():
        raise FileNotFoundError(
            "Native discovery requires the repaired manifest and frozen annotated atlas:\n"
            f"{manifest}\n{atlas}"
        )

    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    output = Path(os.environ.get(
        "GASTRIC_TME_NATIVE_DISCOVERY_OUTPUT",
        str(recovery / f"08_NATIVE_PER_COHORT_STATE_DISCOVERY_{timestamp}"),
    ))
    review_output = output / "MANUAL_CURATION_PACKAGE"
    output.mkdir(parents=True, exist_ok=True)
    print(f"Code: {code}\nData: {data}\nAtlas: {atlas}\nOutput: {output}", flush=True)

    run(
        sys.executable, "-m", "pip", "install", "-q",
        "scanpy>=1.11,<1.13", "anndata>=0.11,<0.13", "igraph", "leidenalg", "xlsxwriter",
    )
    run(
        sys.executable, str(code / "per_cohort_native_state_discovery.py"),
        "--manifest", str(manifest),
        "--data-root", str(data),
        "--reference-atlas", str(atlas),
        "--output-dir", str(output),
        "--max-cells-per-cohort", os.environ.get("GASTRIC_TME_MAX_CELLS_PER_COHORT", "30000"),
        "--min-cells-per-compartment", "100",
        "--resolution", "0.6",
        "--n-hvgs", "3000",
    )
    audit = output / "LINEAGE_AUDIT"
    run(
        sys.executable, str(code / "native_cluster_lineage_audit.py"),
        "--discovery-dir", str(output), "--output-dir", str(audit),
    )
    run(
        sys.executable, str(code / "native_cluster_review_workbook.py"),
        "--discovery-dir", str(output), "--audit-dir", str(audit), "--output-dir", str(review_output),
    )
    print("Completed native per-cohort discovery and curation package.", flush=True)
    print(f"Review workbook: {review_output / 'NATIVE_CLUSTER_REVIEW_WORKBOOK.xlsx'}", flush=True)


if __name__ == "__main__":
    main()
