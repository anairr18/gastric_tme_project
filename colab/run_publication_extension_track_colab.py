#!/usr/bin/env python3
"""Run the evidence-building extension track without altering the frozen atlas.

Public tumour-only spatial cohorts are used solely for spatial ecology. GSE189926
is a public, CD45-selected pre/post chemoimmunotherapy cohort and is used only
for frozen immune-state context. Paired-spatial claims still require a truly
paired, auditable external spatial input.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import anndata as ad
import pandas as pd


GSE308624 = {
    "GSE308624_gastric_cancer_smi.h5ad": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE308624&format=file&file=GSE308624_gastric_cancer_smi.h5ad",
    "GSE308624_gastric_sample_clinical.xlsx": "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE308624&format=file&file=GSE308624_gastric_sample_clinical.xlsx",
}


def run(*command: str) -> None:
    print("+", " ".join(command), flush=True)
    subprocess.run(command, check=True)


def latest_native_run(project: Path) -> Path:
    root = project / "REPRODUCIBLE_NATIVE_STATE_META"
    candidates = sorted((item for item in root.iterdir() if item.is_dir()), key=lambda item: item.stat().st_mtime, reverse=True)
    if not candidates:
        raise FileNotFoundError(f"No reproducible native-state run under {root}")
    return candidates[0]


def first_existing(candidates: list[Path], description: str) -> Path:
    found = next((item for item in candidates if item.exists()), None)
    if found is None:
        raise FileNotFoundError(f"Missing {description}. Checked:\n" + "\n".join(str(item) for item in candidates))
    return found


def download(url: str, destination: Path, min_bytes: int = 1024) -> None:
    if destination.exists() and destination.stat().st_size >= min_bytes:
        print(f"[reuse] {destination}", flush=True)
        return
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".part")
    request = urllib.request.Request(url, headers={"User-Agent": "gastric-tme-publication-validation/1.0"})
    print(f"[download] {url}", flush=True)
    with urllib.request.urlopen(request, timeout=180) as response, temporary.open("wb") as handle:
        shutil.copyfileobj(response, handle)
    if temporary.stat().st_size < min_bytes:
        temporary.unlink(missing_ok=True)
        raise ValueError(f"Downloaded file was unexpectedly small: {destination.name}")
    temporary.replace(destination)


def inspect_optional_h5ad(path: Path, purpose: str) -> dict[str, object]:
    if not path.exists():
        return {"status": "missing", "purpose": purpose, "path": str(path)}
    atlas = ad.read_h5ad(path, backed="r")
    try:
        columns = {str(column).lower() for column in atlas.obs.columns}
        return {
            "status": "present", "purpose": purpose, "path": str(path),
            "n_cells": int(atlas.n_obs), "n_genes": int(atlas.n_vars),
            "has_counts": "counts" in atlas.layers,
            "has_patient": bool(columns & {"patient", "patient_id", "donor", "donor_id"}),
            "has_sample": bool(columns & {"sample", "sample_id", "orig.ident"}),
            "has_tissue_or_condition": bool(columns & {"tissue", "condition", "region", "pathology"}),
            "has_spatial_coordinates": "spatial" in atlas.obsm,
        }
    finally:
        atlas.file.close()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, default=Path("/content/drive/MyDrive/gastric_tme_project"))
    parser.add_argument("--data-root", type=Path, default=Path("/content/drive/MyDrive/data"))
    parser.add_argument("--run", type=Path)
    parser.add_argument("--download-gse308624", action="store_true")
    parser.add_argument("--download-gse189926", action="store_true")
    parser.add_argument("--download-tcga", action="store_true")
    parser.add_argument("--spatial-permutations", type=int, default=499)
    parser.add_argument("--treatment-permutations", type=int, default=10000)
    args = parser.parse_args()
    base = Path(__file__).resolve().parents[1]
    native_run = args.run or latest_native_run(args.project)
    output = native_run / "08_PUBLICATION_EXTENSION_TRACK"
    output.mkdir(parents=True, exist_ok=True)
    post = native_run / "07_POST_DISCOVERY_VALIDATION_TRACK"

    profiles = post / "01_STATE_REFERENCE_PROFILES" / "FROZEN_NATIVE_STATE_REFERENCE_PROFILES.csv.gz"
    composition = native_run / "05_SAMPLE_LEVEL_COMPOSITION" / "PRIMARY_WORKING_STATE_SAMPLE_COMPOSITION.csv"
    clinical = first_existing([
        args.project / "outputs" / "IF15_SUBMISSION_RESCUE" / "PATIENT_LEVEL_FEATURES_REBUILT.csv",
        args.project / "PATIENT_LEVEL_FEATURES_REBUILT.csv",
    ], "audited Korean clinical table")
    if not profiles.exists():
        raise FileNotFoundError("Missing frozen state reference profiles. Complete the post-discovery validation profile stage first.")
    if not composition.exists():
        raise FileNotFoundError(f"Missing zero-complete sample composition: {composition}")

    stages: dict[str, object] = {}
    korean = output / "01_KOREAN_LONGITUDINAL_AND_PATHOLOGY"
    run(sys.executable, str(base / "korean_native_state_treatment_context.py"), "--composition", str(composition), "--clinical", str(clinical), "--output-dir", str(korean / "treatment"), "--permutations", str(args.treatment_permutations))
    run(sys.executable, str(base / "korean_pathology_linkage.py"), "--composition", str(composition), "--clinical", str(clinical), "--output-dir", str(korean / "pathology"))
    stages["korean_longitudinal_and_pathology"] = "completed"

    extension_root = args.data_root / "external" / "validation_extensions"

    # This public single-cell cohort provides useful independent disease-site
    # context, but it has too few paired tumour-normal donors for confirmation.
    gse163558_root = args.data_root / "external" / "validation"
    gse163558_preflight = output / "00_GSE163558_PREFLIGHT"
    gse163558_programs = output / "00_GSE163558_FROZEN_PROGRAM_CONTEXT"
    try:
        run(sys.executable, str(base / "download_gse163558_validation.py"), "--data-root", str(args.data_root))
        run(
            sys.executable,
            str(base / "external_validation_preflight.py"),
            "--registry", str(base / "data" / "external" / "external_validation_registry.csv"),
            "--data-root", str(gse163558_root),
            "--output-dir", str(gse163558_preflight),
        )
        run(
            sys.executable,
            str(base / "frozen_external_program_validation.py"),
            "--preflight", str(gse163558_preflight / "EXTERNAL_VALIDATION_PREFLIGHT.csv"),
            "--output-dir", str(gse163558_programs),
        )
        run(
            sys.executable,
            str(base / "gse163558_context_summary.py"),
            "--scores", str(gse163558_programs / "EXTERNAL_FROZEN_PROGRAM_SAMPLE_SCORES.csv"),
            "--output-dir", str(gse163558_programs),
        )
        stages["gse163558_primary_adjacent_metastatic_context"] = (
            "completed; validation-only descriptive sample-pseudobulk context, not inferential tumour-normal replication"
        )
    except (subprocess.CalledProcessError, FileNotFoundError, ValueError) as error:
        gse163558_programs.mkdir(parents=True, exist_ok=True)
        (gse163558_programs / "GSE163558_CONTEXT_STATUS.txt").write_text(
            f"excluded_or_failed: {type(error).__name__}: {error}\n", encoding="utf-8"
        )
        stages["gse163558_primary_adjacent_metastatic_context"] = f"excluded_or_failed: {error}"

    gse308 = extension_root / "GSE308624"
    if args.download_gse308624:
        for name, url in GSE308624.items():
            download(url, gse308 / name, min_bytes=10000)
    spatial_h5ad = gse308 / "GSE308624_gastric_cancer_smi.h5ad"
    spatial_clinical = gse308 / "GSE308624_gastric_sample_clinical.xlsx"
    if spatial_h5ad.exists() and spatial_clinical.exists():
        spatial = output / "02_GSE308624_COSMX_SPATIAL_ECOLOGY"
        try:
            run(sys.executable, str(base / "spatial_gse308624_validation.py"), "--input-h5ad", str(spatial_h5ad), "--clinical", str(spatial_clinical), "--output-dir", str(spatial), "--n-permutations", str(args.spatial_permutations))
            stages["gse308624_tumour_spatial_ecology"] = "completed; tumour-only spatial association, not tumour-normal validation"
        except subprocess.CalledProcessError as error:
            spatial.mkdir(parents=True, exist_ok=True)
            (spatial / "SPATIAL_ECOLOGY_STATUS.txt").write_text(
                f"excluded_or_failed: {type(error).__name__}: {error}\n",
                encoding="utf-8",
            )
            stages["gse308624_tumour_spatial_ecology"] = f"excluded_or_failed: {error}"
    else:
        stages["gse308624_tumour_spatial_ecology"] = "not_run; download absent or incomplete"

    gse189 = extension_root / "GSE189926"
    gse189_expression = gse189 / "GSE189926_immune_pseudobulk_counts.csv"
    gse189_clinical = gse189 / "GSE189926_clinical_audited.csv"
    if args.download_gse189926:
        run(
            sys.executable,
            str(base / "download_gse189926_ici_pseudobulk.py"),
            "--output-dir", str(gse189), "--download-matrices",
        )
    if gse189_expression.exists() and gse189_clinical.exists():
        run(
            sys.executable,
            str(base / "gse189926_ici_immune_context.py"),
            "--expression", str(gse189_expression),
            "--clinical", str(gse189_clinical),
            "--reference-profiles", str(profiles),
            "--output-dir", str(output / "03_GSE189926_ICI_IMMUNE_CONTEXT"),
        )
        stages["gse189926_ici_immune_context"] = (
            "completed; public CD45-selected, frozen immune-state treatment context only"
        )
    else:
        stages["gse189926_ici_immune_context"] = (
            "not_run; pass --download-gse189926 to obtain public matrices and build auditable immune pseudobulk"
        )

    optional = {
        "paired_spatial_tumour_normal": inspect_optional_h5ad(extension_root / "PAIRED_GC_SPATIAL" / "paired_tumor_normal_spatial.h5ad", "requires >=3 paired tumour-normal patients and coordinates"),
        "external_ici_scrna": inspect_optional_h5ad(extension_root / "PRJEB60680_GSE315929" / "PRJEB60680_GSE315929.h5ad", "requires patient, sample, timepoint, response, and counts"),
    }
    ici_expression = extension_root / "PRJEB25780" / "PRJEB25780_expression.csv"
    ici_clinical = extension_root / "PRJEB25780" / "PRJEB25780_clinical.csv"
    optional["prjeb25780_bulk_ici"] = {"status": "present" if ici_expression.exists() and ici_clinical.exists() else "missing", "expression": str(ici_expression), "clinical": str(ici_clinical), "access": "https://www.ebi.ac.uk/ena/browser/view/PRJEB25780"}
    if ici_expression.exists() and ici_clinical.exists():
        run(sys.executable, str(base / "frozen_bulk_clinical_validation.py"), "--expression", str(ici_expression), "--clinical", str(ici_clinical), "--reference-profiles", str(profiles), "--output-dir", str(output / "03_PRJEB25780_ICI_RESPONSE"), "--mode", "ici_response")
        stages["prjeb25780_ici_response"] = "completed; frozen signature association only"
    else:
        stages["prjeb25780_ici_response"] = "not_run; processed expression and auditable clinical response table required"

    tcga_expression = extension_root / "TCGA_STAD_SUBTYPES" / "TCGA_STAD_expression.csv"
    tcga_clinical = extension_root / "TCGA_STAD_SUBTYPES" / "TCGA_STAD_clinical.csv"
    if args.download_tcga:
        run(
            sys.executable,
            str(base / "download_tcga_stad_2014.py"),
            "--output-dir", str(tcga_expression.parent),
        )
    optional["tcga_molecular_subtypes"] = {"status": "present" if tcga_expression.exists() and tcga_clinical.exists() else "missing", "expression": str(tcga_expression), "clinical": str(tcga_clinical)}
    if tcga_expression.exists() and tcga_clinical.exists():
        run(sys.executable, str(base / "frozen_bulk_clinical_validation.py"), "--expression", str(tcga_expression), "--clinical", str(tcga_clinical), "--reference-profiles", str(profiles), "--output-dir", str(output / "04_TCGA_MOLECULAR_SUBTYPES"), "--mode", "molecular_subtype")
        stages["tcga_molecular_subtypes"] = "completed; frozen signature association only"
    else:
        stages["tcga_molecular_subtypes"] = "not_run; matched expression and TCGA subtype table required"

    pd.DataFrame([{"extension": key, **value} for key, value in optional.items()]).to_csv(output / "OPTIONAL_EXTERNAL_INPUT_PREFLIGHT.csv", index=False)
    (output / "PUBLICATION_EXTENSION_MANIFEST.json").write_text(json.dumps({
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "native_run": str(native_run), "stages": stages, "optional_inputs": optional,
        "claim_boundary": "No external cohort enters scVI or state discovery. Spatial tumour-only cohorts evaluate local ecology only. Clinical associations use frozen signatures and are not predictive models.",
    }, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(stages, indent=2))
    print(f"Publication extension track completed: {output}")


if __name__ == "__main__":
    main()
