#!/usr/bin/env python3
"""Export count-derived cohort and frozen reference profiles for retained native states.

The exact source-row assignments from cohort-native clustering are reused. This
does not re-cluster cells or alter the working state dictionary.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

from raw_integration_method_audit import candidate_paths, count_matrix


KEY = ["dataset_id", "broad_label", "native_cluster"]


def log_cp10k(matrix):
    matrix = matrix.tocsr().astype(float) if sparse.issparse(matrix) else sparse.csr_matrix(matrix, dtype=float)
    totals = np.asarray(matrix.sum(axis=1)).ravel()
    if np.any(totals <= 0):
        raise ValueError(f"Selected native cells include {int(np.count_nonzero(totals <= 0))} zero-count cells.")
    matrix = sparse.diags(1e4 / totals) @ matrix
    matrix.data = np.log1p(matrix.data)
    return matrix


def collapse_duplicate_gene_symbols(matrix, genes: pd.Index):
    """Sum duplicate count columns before normalisation and marker profiling."""
    genes = pd.Index(genes.astype(str).str.upper())
    if not genes.duplicated().any():
        return matrix, genes, 0
    codes, unique = pd.factorize(genes, sort=True)
    source = matrix.tocsr() if sparse.issparse(matrix) else sparse.csr_matrix(matrix)
    projector = sparse.csr_matrix(
        (np.ones(len(codes), dtype=np.float32), (np.arange(len(codes)), codes)),
        shape=(len(codes), len(unique)),
    )
    return source @ projector, pd.Index(unique.astype(str)), int(genes.duplicated().sum())


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assignments", required=True, type=Path)
    parser.add_argument("--resolution-ledger", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--data-root", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    assignments = pd.read_csv(args.assignments, compression="infer")
    ledger = pd.read_csv(args.resolution_ledger)
    manifest = pd.read_csv(args.manifest).set_index("dataset_id", drop=False)
    for table in (assignments, ledger):
        table["native_cluster"] = table["native_cluster"].astype(str)
    retained = ledger.loc[ledger["analysis_tier"].eq("primary_working_state")].copy()
    if retained.empty:
        raise ValueError("No primary_working_state clusters in the resolution ledger.")
    states = retained[KEY + ["provisional_state_name"]].drop_duplicates()
    assigned = assignments.merge(states, on=KEY, how="inner", validate="many_to_one")
    if assigned.empty:
        raise ValueError("No assignment rows match the retained state dictionary.")
    if assigned.duplicated(["dataset_id", "source_row_index"]).any():
        raise ValueError("Duplicate source-row assignments prevent auditable profiling.")
    assigned["state_id"] = assigned["broad_label"].astype(str) + "::" + assigned["provisional_state_name"].astype(str)

    profile_parts: list[pd.DataFrame] = []
    audit_rows: list[dict[str, object]] = []
    import anndata as ad

    for dataset_id, cohort_rows in assigned.groupby("dataset_id", observed=True):
        if dataset_id not in manifest.index:
            raise ValueError(f"Manifest lacks {dataset_id}.")
        source_path = next(
            (path for path in candidate_paths(args.data_root, dataset_id, manifest.loc[dataset_id]) if path.exists()),
            None,
        )
        if source_path is None:
            raise FileNotFoundError(f"{dataset_id}: source H5AD was not found.")
        source = ad.read_h5ad(source_path, backed="r")
        try:
            rows = cohort_rows.sort_values("source_row_index").reset_index(drop=True)
            selected = rows["source_row_index"].to_numpy(dtype=int)
            if selected.min() < 0 or selected.max() >= source.n_obs:
                raise ValueError(f"{dataset_id}: source-row assignment lies outside the source H5AD.")
            subset = source[selected, :].to_memory()
            counts, layer = count_matrix(subset)
            counts, genes, duplicated = collapse_duplicate_gene_symbols(counts, subset.var_names)
            expression = log_cp10k(counts)
            for state_id, state_rows in rows.groupby("state_id", observed=True):
                indices = state_rows.index.to_numpy(dtype=int)
                mean_expression = np.asarray(expression[indices].mean(axis=0)).ravel()
                profile_parts.append(pd.DataFrame({
                    "dataset_id": dataset_id,
                    "state_id": state_id,
                    "gene": genes,
                    "mean_log_normalized_expression": mean_expression,
                    "n_state_cells": int(len(indices)),
                }))
            audit_rows.append({
                "dataset_id": dataset_id,
                "source_h5ad": str(source_path),
                "count_source": layer,
                "n_profiled_cells": int(len(rows)),
                "n_states": int(rows["state_id"].nunique()),
                "n_genes": int(len(genes)),
                "n_collapsed_duplicate_gene_symbols": duplicated,
            })
        finally:
            source.file.close()

    profiles = pd.concat(profile_parts, ignore_index=True)
    profiles.to_csv(args.output_dir / "COHORT_NATIVE_STATE_EXPRESSION_PROFILES.csv.gz", index=False, compression="gzip")
    weights = profiles["mean_log_normalized_expression"] * profiles["n_state_cells"]
    frozen = profiles.assign(weighted_expression=weights).groupby(["state_id", "gene"], observed=True).agg(
        weighted_expression=("weighted_expression", "sum"), n_reference_cells=("n_state_cells", "sum")
    ).reset_index()
    frozen["mean_log_normalized_expression"] = frozen["weighted_expression"] / frozen["n_reference_cells"]
    frozen = frozen.drop(columns="weighted_expression")
    frozen.to_csv(args.output_dir / "FROZEN_NATIVE_STATE_REFERENCE_PROFILES.csv.gz", index=False, compression="gzip")
    pd.DataFrame(audit_rows).to_csv(args.output_dir / "REFERENCE_PROFILE_INPUT_AUDIT.csv", index=False)
    (args.output_dir / "README_AND_CLAIM_BOUNDARIES.md").write_text(
        "# Native state reference profiles\n\n"
        "Profiles are CP10K-log1p means from the exact cells retained by native cohort clustering. "
        "They support state-definition consistency checks and external reference mapping; they do not independently "
        "validate labels that were initially curated with marker evidence.\n",
        encoding="utf-8",
    )
    (args.output_dir / "RUN_MANIFEST.json").write_text(json.dumps({
        "assignments": str(args.assignments), "resolution_ledger": str(args.resolution_ledger),
        "n_profile_rows": int(len(profiles)), "n_states": int(profiles["state_id"].nunique()),
        "status": "working_curation_pending_mentor_review",
    }, indent=2) + "\n", encoding="utf-8")
    print(f"Native state reference profiles written to {args.output_dir}")


if __name__ == "__main__":
    main()
