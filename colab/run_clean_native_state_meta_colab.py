#!/usr/bin/env python3
"""Build a self-contained native-state discovery and sample-level meta-analysis.

This is the recovery path when an older native discovery lacks per-cell Leiden
assignments. It does not transfer labels from that older run. Instead it runs
one current, recorded environment from count-source audit through within-cohort
clustering, conservative computational curation, and sample-level inference.
All state labels remain provisional pending mentor review.
"""

from __future__ import annotations

import os
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def run(*command: str) -> None:
    print("+", " ".join(command), flush=True)
    subprocess.run(command, check=True)


def find_reference_atlas(project_root: Path, data_root: Path) -> Path:
    candidates = [
        project_root / "colab_rawcounts_sixcohort" / "gastric_meta_annotated_rawcounts_sixcohort.h5ad",
        project_root / "colab_rawcounts_sixcohort" / "gastric_meta_annotated.h5ad",
        data_root / "processed" / "integrated" / "gastric_meta_annotated.h5ad",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    raise FileNotFoundError("Annotated six-cohort atlas was not found in the supported Drive locations.")


def main() -> None:
    code_root = Path(__file__).resolve().parents[1]
    project_root = Path(os.environ.get("GASTRIC_TME_PROJECT_ROOT", "/content/drive/MyDrive/gastric_tme_project"))
    data_root = Path(os.environ.get("GASTRIC_TME_DATA_ROOT", "/content/drive/MyDrive/data"))
    manifest = data_root / "external" / "dataset_manifest.csv"
    reference_atlas = find_reference_atlas(project_root, data_root)
    repaired_korea = data_root / "processed" / "per_dataset" / "korea_kim2022_raw_counts_sample_qualified.h5ad"
    required = [manifest, reference_atlas, repaired_korea]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing required input(s):\n" + "\n".join(missing))

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    run_root = Path(os.environ.get(
        "GASTRIC_TME_CLEAN_NATIVE_RUN",
        str(project_root / "REPRODUCIBLE_NATIVE_STATE_META" / stamp),
    ))
    run_root.mkdir(parents=True, exist_ok=True)
    runtime_manifest = run_root / "dataset_manifest_runtime.csv"
    manifest_frame = pd.read_csv(manifest)
    manifest_frame["local_h5ad"] = manifest_frame.get(
        "local_h5ad", pd.Series(index=manifest_frame.index, dtype="object")
    ).astype("object")
    manifest_frame.loc[
        manifest_frame["dataset_id"].astype(str).eq("korea_kim2022"), "local_h5ad"
    ] = str(repaired_korea)
    manifest_frame.to_csv(runtime_manifest, index=False)

    raw_audit = run_root / "00_COUNT_SOURCE_AUDIT"
    native = run_root / "01_NATIVE_PER_COHORT_DISCOVERY"
    lineage = run_root / "02_LINEAGE_AUDIT"
    curation = run_root / "03_PROVISIONAL_CURATION"
    resolution = run_root / "04_CLUSTER_RESOLUTION"
    composition = run_root / "05_SAMPLE_LEVEL_COMPOSITION"
    meta = run_root / "06_STATE_ABUNDANCE_META"
    max_cells = os.environ.get("GASTRIC_TME_NATIVE_MAX_CELLS", "30000")

    print(f"Code: {code_root}")
    print(f"Project: {project_root}")
    print(f"Data: {data_root}")
    print(f"Korean count source: {repaired_korea}")
    print(f"Reference atlas: {reference_atlas}")
    print(f"Output: {run_root}")
    print(f"Native cell cap per cohort: {max_cells}")

    with (run_root / "environment_freeze.txt").open("w", encoding="utf-8") as handle:
        subprocess.run([sys.executable, "-m", "pip", "freeze"], stdout=handle, check=True)

    run(
        sys.executable, str(code_root / "raw_integration_method_audit.py"),
        "--manifest", str(runtime_manifest), "--data-root", str(data_root), "--output-dir", str(raw_audit),
    )
    run(
        sys.executable, str(code_root / "per_cohort_native_state_discovery.py"),
        "--manifest", str(runtime_manifest), "--data-root", str(data_root),
        "--reference-atlas", str(reference_atlas), "--output-dir", str(native),
        "--max-cells-per-cohort", str(max_cells), "--min-cells-per-compartment", "100",
        "--resolution", "0.6", "--n-hvgs", "3000",
    )
    run(
        sys.executable, str(code_root / "native_cluster_lineage_audit.py"),
        "--discovery-dir", str(native), "--output-dir", str(lineage),
    )
    run(
        sys.executable, str(code_root / "computational_provisional_native_curation.py"),
        "--audit", str(lineage / "NATIVE_CLUSTER_LINEAGE_AUDIT.csv"),
        "--markers", str(native / "NATIVE_PER_COHORT_CLUSTER_MARKERS.csv"),
        "--output-dir", str(curation),
    )
    run(
        sys.executable, str(code_root / "resolve_all_native_cluster_dispositions.py"),
        "--dictionary", str(curation / "COMPUTATIONAL_PROVISIONAL_STATE_DICTIONARY.csv"),
        "--output-dir", str(resolution),
    )
    run(
        sys.executable, str(code_root / "build_primary_working_state_composition.py"),
        "--assignments", str(native / "NATIVE_PER_COHORT_CELL_ASSIGNMENTS.csv.gz"),
        "--resolution-ledger", str(resolution / "ALL_505_NATIVE_CLUSTER_RESOLUTION_LEDGER.csv"),
        "--output-dir", str(composition),
    )
    run(
        sys.executable, str(code_root / "state_abundance_meta_analysis.py"),
        "--curated-composition", str(composition / "PRIMARY_WORKING_STATE_SAMPLE_COMPOSITION.csv"),
        "--output-dir", str(meta),
    )
    print(f"Completed clean native-state meta-analysis: {run_root}")


if __name__ == "__main__":
    main()
