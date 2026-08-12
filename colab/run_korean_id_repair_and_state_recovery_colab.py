#!/usr/bin/env python3
"""Repair Korean raw cell IDs, rerun the raw gate, and recover state coverage."""

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


def newest_label_table(project: Path) -> Path:
    roots = [project / "SINGLE_CELL_ANALYSIS_AUTOMATED_COMPLETE", project / "COMPARTMENT_STATE_DISCOVERY"]
    candidates: list[Path] = []
    for root in roots:
        if root.exists():
            candidates.extend(root.rglob("FROZEN_CURATED_STATE_CELL_LABELS.csv.gz"))
    if not candidates:
        raise FileNotFoundError("No FROZEN_CURATED_STATE_CELL_LABELS.csv.gz was found under the completed state-analysis outputs.")
    return max(candidates, key=lambda path: path.stat().st_mtime)


def main() -> None:
    code = Path(__file__).resolve().parents[1]
    data = Path(os.environ.get("GASTRIC_TME_DATA_ROOT", "/content/drive/MyDrive/data"))
    project = Path(os.environ.get("GASTRIC_TME_PROJECT_ROOT", "/content/drive/MyDrive/gastric_tme_project"))
    atlas = Path(os.environ.get("GASTRIC_TME_ATLAS", project / "colab_rawcounts_sixcohort" / "gastric_meta_annotated_rawcounts_sixcohort.h5ad"))
    manifest = Path(os.environ.get("GASTRIC_TME_MANIFEST", data / "external" / "dataset_manifest.csv"))
    raw = Path(os.environ.get("GASTRIC_TME_KOREA_RAW", data / "raw" / "ge_korea_raw_data_count_matricies_raw_combined.h5ad"))
    repaired = Path(os.environ.get("GASTRIC_TME_KOREA_REPAIRED_H5AD", data / "processed" / "per_dataset" / "korea_kim2022_raw_counts_sample_qualified.h5ad"))
    if not atlas.exists() or not manifest.exists() or not raw.exists():
        raise FileNotFoundError(f"Required input missing: atlas={atlas.exists()}, manifest={manifest.exists()}, korean_raw={raw.exists()}")
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    output = project / "MENTOR_SYSTEMATIC_SINGLE_CELL_AUDIT" / f"{timestamp}_KOREA_ID_AND_STATE_RECOVERY"
    output.mkdir(parents=True, exist_ok=True)
    print(f"Atlas: {atlas}\nKorean raw: {raw}\nRepaired Korean output: {repaired}\nResults: {output}")
    run(sys.executable, "-m", "pip", "install", "-q", "anndata>=0.10,<0.13", "h5py>=3.10", "scipy>=1.13", "matplotlib>=3.8", "seaborn>=0.13")
    repair = output / "01_KOREA_ID_REPAIR"
    run(sys.executable, str(code / "repair_korean_raw_cell_ids.py"), "--input-h5ad", str(raw), "--output-dir", str(repair), "--output-h5ad", str(repaired))

    repaired_manifest = pd.read_csv(manifest)
    if "dataset_id" not in repaired_manifest.columns or "korea_kim2022" not in set(repaired_manifest["dataset_id"].astype(str)):
        raise ValueError("Manifest must contain korea_kim2022.")
    repaired_manifest["local_h5ad"] = repaired_manifest.get("local_h5ad", pd.Series(index=repaired_manifest.index, dtype="object")).astype("object")
    repaired_manifest.loc[repaired_manifest["dataset_id"].astype(str).eq("korea_kim2022"), "local_h5ad"] = str(repaired)
    repaired_manifest_path = output / "dataset_manifest_korea_id_repaired.csv"
    repaired_manifest.to_csv(repaired_manifest_path, index=False)
    raw_gate = output / "02_RAW_GATE_AFTER_KOREA_ID_REPAIR"
    run(sys.executable, str(code / "raw_integration_method_audit.py"), "--manifest", str(repaired_manifest_path), "--data-root", str(data), "--output-dir", str(raw_gate))

    labels = Path(os.environ["GASTRIC_TME_STATE_LABELS"]) if "GASTRIC_TME_STATE_LABELS" in os.environ else newest_label_table(project)
    context = output / "03_COHORT_CONTEXT"
    run(sys.executable, str(code / "cohort_context_audit.py"), "--registry", str(code / "cohort_context_registry.csv"), "--manifest", str(repaired_manifest_path), "--atlas", str(atlas), "--output-dir", str(context))
    state_audit = output / "04_PER_COHORT_STATE_AUDIT"
    run(sys.executable, str(code / "cohort_resolved_annotation_audit.py"), "--atlas", str(atlas), "--cohort-context", str(context / "COHORT_CONTEXT_FOR_SLIDES.csv"), "--state-labels", str(labels), "--output-dir", str(state_audit))
    if os.environ.get("GASTRIC_TME_RUN_METHOD_BENCHMARK", "1") == "1":
        run(sys.executable, "-m", "pip", "install", "-q", "scanpy>=1.10,<1.13", "scvi-tools>=1.3,<1.5", "harmonypy>=0.0.10", "igraph>=0.11", "leidenalg>=0.10")
        benchmark = output / "05_SCVI_HARMONY_SENSITIVITY"
        run(
            sys.executable, str(code / "integration_method_sensitivity_benchmark.py"),
            "--manifest", str(repaired_manifest_path), "--data-root", str(data),
            "--preflight", str(raw_gate / "RAW_INPUT_PREFLIGHT.csv"), "--reference-atlas", str(atlas), "--output-dir", str(benchmark),
            "--max-cells-per-cohort", os.environ.get("GASTRIC_TME_BENCHMARK_CELLS_PER_COHORT", "1200"),
            "--max-epochs", os.environ.get("GASTRIC_TME_BENCHMARK_EPOCHS", "30"),
        )
    (output / "README.txt").write_text(
        "This run creates a new Korean raw-count H5AD with deterministic sample-qualified IDs, then reruns the raw integration gate, frozen-state coverage audit, and a controlled scVI/Harmony sensitivity benchmark. The original Korean H5AD is unchanged.\n",
        encoding="utf-8",
    )
    print(f"Completed Korean ID repair and state recovery: {output}")


if __name__ == "__main__":
    main()
