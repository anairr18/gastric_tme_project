#!/usr/bin/env python3
"""Audit raw inputs before deciding whether scVI or Harmony is appropriate.

The audit does not select an integration method from a UMAP. It first verifies
raw-count provenance, sample identifiers, gene overlap, and the biological
variables that must *not* be included in a technical batch key.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


CORE_IDS = ["korea_kim2022", "kumar2022", "sathe2020", "zhang2021", "diffuse_gc_2021", "tcell_exhaustion_2022"]
COUNT_LAYERS = ["counts", "raw_counts", "count", "counts_raw", "raw"]
SAMPLE_COLUMNS = ["sample_id", "sample", "orig.ident", "library_id", "analysis_unit"]
PATIENT_COLUMNS = ["patient_id", "patient", "donor_id", "donor"]
TISSUE_COLUMNS = ["tissue", "condition", "sample_type", "lesion_type", "stage"]
TREATMENT_COLUMNS = ["timepoint", "treatment", "best_overall_response", "response", "recist"]
BROAD_COLUMNS = ["cell_type_coarse", "cell_type", "broad_cell_type", "cell_type_coarse_fine"]


def col(frame: pd.DataFrame, candidates: list[str]) -> str | None:
    lookup = {str(name).lower(): str(name) for name in frame.columns}
    return next((lookup[name.lower()] for name in candidates if name.lower() in lookup), None)


def count_matrix(adata):
    for layer in COUNT_LAYERS:
        if layer in adata.layers:
            return adata.layers[layer], layer
    # Backed AnnData sparse datasets cannot be safely coerced by np.asarray.
    # Integer-count validation happens after chunked extraction in inspect_one.
    return adata.X, "X_unverified"


def matrix_sums(matrix, n_rows: int):
    if sparse.issparse(matrix):
        library = np.asarray(matrix.sum(axis=1)).ravel()
        genes = np.asarray((matrix > 0).sum(axis=1)).ravel()
        values = matrix.data
    elif isinstance(matrix, np.ndarray):
        array = np.asarray(matrix)
        library = array.sum(axis=1)
        genes = (array > 0).sum(axis=1)
        values = array.ravel()
    else:
        # AnnData backed sparse datasets do not expose .data or sum(axis=1) in
        # every supported version. Read small row blocks, never the full matrix.
        libraries, detected, value_parts = [], [], []
        for start in range(0, n_rows, 4096):
            block = matrix[start : min(n_rows, start + 4096)]
            if sparse.issparse(block):
                libraries.append(np.asarray(block.sum(axis=1)).ravel())
                detected.append(np.asarray((block > 0).sum(axis=1)).ravel())
                if len(value_parts) < 2:
                    value_parts.append(block.data[: max(0, 100000 - sum(len(part) for part in value_parts))])
            else:
                array = np.asarray(block)
                libraries.append(array.sum(axis=1))
                detected.append((array > 0).sum(axis=1))
                if len(value_parts) < 2:
                    value_parts.append(array.ravel()[: max(0, 100000 - sum(len(part) for part in value_parts))])
        library = np.concatenate(libraries) if libraries else np.array([])
        genes = np.concatenate(detected) if detected else np.array([])
        values = np.concatenate(value_parts) if value_parts else np.array([])
    return library, genes, values[: min(len(values), 100000)]


def candidate_paths(data_root: Path, dataset_id: str, manifest_row: pd.Series) -> list[Path]:
    paths: list[Path] = []
    supplied = manifest_row.get("local_h5ad")
    if isinstance(supplied, str) and supplied.strip():
        candidate = Path(supplied)
        paths += [candidate, data_root / candidate.name]
    paths += [
        data_root / "processed" / "per_dataset" / f"{dataset_id}_processed.h5ad",
        data_root / "processed" / "per_dataset" / f"{dataset_id}_metadata_only.h5ad",
    ]
    if dataset_id == "korea_kim2022":
        # A manifest-specified repaired source takes precedence. The original
        # raw file remains a fallback only when no repaired source is supplied.
        paths.append(data_root / "raw" / "ge_korea_raw_data_count_matricies_raw_combined.h5ad")
    unique: list[Path] = []
    for path in paths:
        if path not in unique:
            unique.append(path)
    return unique


def inspect_one(dataset_id: str, path: Path) -> tuple[dict[str, object], set[str], pd.DataFrame]:
    import anndata as ad

    data = ad.read_h5ad(path, backed="r")
    try:
        obs = data.obs.copy()
        var_names = pd.Index(data.var_names.astype(str).str.upper())
        matrix, layer = count_matrix(data)
        library, genes, values = matrix_sums(matrix, data.n_obs)
        sample_col, patient_col = col(obs, SAMPLE_COLUMNS), col(obs, PATIENT_COLUMNS)
        tissue_col, treatment_col = col(obs, TISSUE_COLUMNS), col(obs, TREATMENT_COLUMNS)
        raw_integer = bool(len(values) and np.allclose(values, np.rint(values), rtol=0, atol=1e-6) and np.nanmin(values) >= 0)
        if layer == "X_unverified" and raw_integer:
            layer = "X_integer_counts"
        record = {
            "dataset_id": dataset_id,
            "input_h5ad": str(path),
            "status": "eligible" if raw_integer and obs.index.is_unique else "blocked",
            "reason": "" if raw_integer and obs.index.is_unique else ("cell IDs are duplicated" if not obs.index.is_unique else "count layer is not non-negative integer counts"),
            "n_cells": int(data.n_obs), "n_genes": int(data.n_vars), "count_source": layer,
            "unique_cell_ids": bool(obs.index.is_unique), "integer_counts": raw_integer,
            "sample_column": sample_col or "", "patient_column": patient_col or "", "tissue_column": tissue_col or "", "treatment_column": treatment_col or "",
            "n_samples": int(obs[sample_col].dropna().astype(str).nunique()) if sample_col else 0,
            "n_patients": int(obs[patient_col].dropna().astype(str).nunique()) if patient_col else 0,
            "sample_metadata_complete": bool(sample_col and patient_col and tissue_col),
            "treatment_metadata_available": bool(treatment_col),
            "median_total_counts": float(np.nanmedian(library)), "median_detected_genes": float(np.nanmedian(genes)),
        }
        sample_frame = pd.DataFrame({"dataset_id": dataset_id, "library_size": library, "detected_genes": genes})
        sample_frame["sample_id"] = obs[sample_col].astype(str).to_numpy() if sample_col else "MISSING"
        sample_qc = sample_frame.groupby(["dataset_id", "sample_id"], observed=True).agg(n_cells=("sample_id", "size"), median_total_counts=("library_size", "median"), median_detected_genes=("detected_genes", "median")).reset_index()
        return record, set(var_names), sample_qc
    finally:
        data.file.close()


def write_decision(preflight: pd.DataFrame, overlaps: pd.DataFrame, destination: Path) -> dict[str, object]:
    eligible = preflight["status"].eq("eligible")
    all_raw = bool(eligible.all()) and len(preflight) == len(CORE_IDS)
    min_shared = int(overlaps.loc[overlaps["comparison"] == "all_core", "shared_genes"].iloc[0]) if (overlaps["comparison"] == "all_core").any() else 0
    if all_raw and min_shared >= 4000:
        status = "scvi_primary_harmony_sensitivity"
        rationale = (
            "All core datasets have auditable integer counts, unique cell IDs, and sufficient shared genes. "
            "Use scVI on raw counts with dataset_id as the technical batch key. Keep tissue, stage, timepoint, and response out of the batch key. "
            "Run Harmony only as a normalized-PC sensitivity analysis, then compare biological preservation and donor-level conclusions rather than choosing by visual mixing."
        )
    else:
        status = "blocked_pending_raw_input_repair"
        rationale = (
            "Do not retrain an integration model yet. At least one core input lacks auditable integer counts, unique IDs, or sufficient shared genes. "
            "Repair the indicated input and rerun this audit; neither scVI nor Harmony should be selected from incomplete provenance."
        )
    decision = {"status": status, "minimum_shared_genes": min_shared, "n_eligible_core_inputs": int(eligible.sum()), "rationale": rationale}
    destination.write_text("# Integration Method Decision\n\n" + f"**Status:** `{status}`\n\n{rationale}\n\n" + "## Required comparison\n\n- Unintegrated normalized PCA: diagnostic baseline.\n- scVI: primary only when raw-count gate passes.\n- Harmony: sensitivity analysis on the same normalized-PC input.\n- Judge methods with cohort mixing, label/marker preservation, donor-level state effects, and held-out mapping; never use UMAP appearance alone.\n", encoding="utf-8")
    return decision


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(args.manifest)
    if "dataset_id" not in manifest:
        raise ValueError("Dataset manifest needs dataset_id.")
    manifest = manifest.set_index("dataset_id", drop=False)

    rows, gene_sets, sample_qcs = [], {}, []
    for dataset_id in CORE_IDS:
        if dataset_id not in manifest.index:
            rows.append({"dataset_id": dataset_id, "status": "blocked", "reason": "missing from manifest"})
            continue
        candidates = candidate_paths(args.data_root, dataset_id, manifest.loc[dataset_id])
        path = next((candidate for candidate in candidates if candidate.exists()), None)
        if path is None:
            rows.append({"dataset_id": dataset_id, "status": "blocked", "reason": "no H5AD found", "searched": "; ".join(map(str, candidates))})
            continue
        try:
            row, genes, sample_qc = inspect_one(dataset_id, path)
            rows.append(row); gene_sets[dataset_id] = genes; sample_qcs.append(sample_qc)
        except Exception as error:
            rows.append({"dataset_id": dataset_id, "input_h5ad": str(path), "status": "blocked", "reason": f"{type(error).__name__}: {error}"})
    preflight = pd.DataFrame(rows)
    preflight.to_csv(args.output_dir / "RAW_INPUT_PREFLIGHT.csv", index=False)
    pd.concat(sample_qcs, ignore_index=True).to_csv(args.output_dir / "RAW_SAMPLE_QC_SUMMARY.csv", index=False) if sample_qcs else pd.DataFrame().to_csv(args.output_dir / "RAW_SAMPLE_QC_SUMMARY.csv", index=False)

    overlap_rows = []
    available = [genes for genes in gene_sets.values() if genes]
    if available:
        overlap_rows.append({"comparison": "all_core", "shared_genes": len(set.intersection(*available)), "n_datasets": len(available)})
    for dataset_id, genes in gene_sets.items():
        overlap_rows.append({"comparison": dataset_id, "shared_genes": len(genes), "n_datasets": 1})
    overlaps = pd.DataFrame(overlap_rows, columns=["comparison", "shared_genes", "n_datasets"])
    overlaps.to_csv(args.output_dir / "RAW_GENE_OVERLAP.csv", index=False)
    decision = write_decision(preflight, overlaps, args.output_dir / "INTEGRATION_METHOD_DECISION.md")
    (args.output_dir / "INTEGRATION_METHOD_DECISION.json").write_text(json.dumps(decision, indent=2) + "\n", encoding="utf-8")
    print(preflight[[column for column in ["dataset_id", "status", "n_cells", "n_genes", "count_source", "reason"] if column in preflight]].to_string(index=False))
    print(f"Integration-method audit written to {args.output_dir}")


if __name__ == "__main__":
    main()
