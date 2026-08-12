#!/usr/bin/env python3
"""Create the data-linked five-question registry without rerunning clustering."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


def main() -> None:
    code = Path(__file__).resolve().parents[1]
    project = Path(os.environ.get("GASTRIC_TME_PROJECT_ROOT", "/content/drive/MyDrive/gastric_tme_project"))
    runs = sorted(
        project.glob("MENTOR_SYSTEMATIC_SINGLE_CELL_AUDIT/*KOREA_ID_AND_STATE_RECOVERY*/08_NATIVE_PER_COHORT_STATE_DISCOVERY_*"),
        key=lambda path: path.stat().st_mtime,
        reverse=True,
    )
    runs = [path for path in runs if (path / "NATIVE_PER_COHORT_CLUSTER_MARKERS.csv").exists() and (path / "LINEAGE_AUDIT" / "NATIVE_CLUSTER_LINEAGE_AUDIT.csv").exists()]
    if not runs:
        raise FileNotFoundError("No completed native discovery and lineage audit found on Drive.")
    discovery = runs[0]
    output = discovery / "MANUSCRIPT_QUESTION_REGISTRY"
    command = [
        sys.executable, str(code / "research_question_registry.py"),
        "--discovery-dir", str(discovery),
        "--audit-dir", str(discovery / "LINEAGE_AUDIT"),
        "--context-registry", str(code / "cohort_context_registry.csv"),
        "--output-dir", str(output),
    ]
    print("+", " ".join(command), flush=True)
    subprocess.run(command, check=True)
    print(f"Question registry: {output}", flush=True)


if __name__ == "__main__":
    main()
