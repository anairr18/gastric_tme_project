#!/usr/bin/env python3
"""Discover native compartment clusters independently within each core cohort.

The frozen atlas supplies only an auditable broad-compartment correspondence.
Each cohort is normalized, reduced, neighboured, and clustered separately from
its own count matrix. Output clusters remain provisional until a reviewer
records a state name and marker rationale in the exported template.
"""

from __future__ import annotations

import argparse
import gc
import json
from pathlib import Path

import numpy as np
import pandas as pd

from integration_method_sensitivity_benchmark import (
    BROAD_COLUMNS,
    SAMPLE_COLUMNS,
    frozen_broad_label_lookup,
    label_keys,
    stratified_indices,
)
from raw_integration_method_audit import CORE_IDS, candidate_paths, col, count_matrix


def write_progress_snapshot(
    output_dir: Path,
    completed_dataset: str,
    input_rows: list[dict],
    marker_tables: list[pd.DataFrame],
    cluster_summaries: list[pd.DataFrame],
    templates: list[pd.DataFrame],
    status_rows: list[dict],
) -> None:
    """Persist a usable checkpoint after each cohort for interruption-safe review."""
    pd.DataFrame(input_rows).to_csv(output_dir / "NATIVE_PER_COHORT_INPUT_AUDIT_PARTIAL.csv", index=False)
    pd.DataFrame(status_rows).to_csv(output_dir / "NATIVE_PER_COHORT_CLUSTER_STATUS_PARTIAL.csv", index=False)
    if marker_tables:
        pd.concat(marker_tables, ignore_index=True).to_csv(
            output_dir / "NATIVE_PER_COHORT_CLUSTER_MARKERS_PARTIAL.csv", index=False
        )
    if cluster_summaries:
        pd.concat(cluster_summaries, ignore_index=True).to_csv(
            output_dir / "NATIVE_PER_COHORT_CLUSTER_SUMMARY_PARTIAL.csv", index=False
        )
    if templates:
        pd.concat(templates, ignore_index=True).to_csv(
            output_dir / "NATIVE_PER_COHORT_ANNOTATION_TEMPLATE_PARTIAL.csv", index=False
        )
    (output_dir / "STATUS.txt").write_text(
        f"Completed native clustering for {completed_dataset}. "
        f"Completed cohorts: {', '.join(sorted({row['dataset_id'] for row in input_rows}))}.\n",
        encoding="utf-8",
    )


def load_native_cohort(
    dataset_id: str,
    manifest_row: pd.Series,
    data_root: Path,
    reference_labels: dict[tuple[str, str], str],
    max_cells: int,
    seed: int,
):
    import anndata as ad

    source_path = next((path for path in candidate_paths(data_root, dataset_id, manifest_row) if path.exists()), None)
    if source_path is None:
        raise FileNotFoundError(f"{dataset_id}: no source H5AD found.")
    source = ad.read_h5ad(source_path, backed="r")
    try:
        _, layer = count_matrix(source)
        if layer == "X_unverified":
            layer = None
        source_labels = np.asarray([
            next((reference_labels[(dataset_id, key)] for key in candidates if (dataset_id, key) in reference_labels), "MISSING")
            for candidates in label_keys(dataset_id, source.obs)
        ])
        matched = np.flatnonzero(source_labels != "MISSING")
        if not len(matched):
            raise ValueError(f"{dataset_id}: no raw cells match the frozen broad-compartment reference.")
        selected = matched[stratified_indices(source.obs.iloc[matched], max_cells, seed)]
        data = source[selected, :].to_memory()
        if layer is not None:
            data.X = data.layers[layer].copy()
        data.obs["dataset_id"] = dataset_id
        data.obs["frozen_broad_label"] = source_labels[selected]
        # Raw IDs may be duplicated, notably in the Korean count source. The
        # source row index makes cluster assignments auditable and reusable.
        data.obs["source_row_index"] = selected.astype(np.int64)
        data.obs["source_cell_id"] = source.obs_names[selected].astype(str)
        sample_column = col(data.obs, SAMPLE_COLUMNS)
        data.obs["analysis_sample"] = data.obs[sample_column].astype(str) if sample_column else "MISSING_SAMPLE"
        return data, {
            "dataset_id": dataset_id,
            "input_h5ad": str(source_path),
            "source_cells": int(source.n_obs),
            "cells_matching_frozen_reference": int(len(matched)),
            "native_cells_sampled": int(data.n_obs),
            "sample_column": sample_column or "MISSING",
        }
    finally:
        source.file.close()


def native_cluster(
    data,
    resolution: float,
    n_hvgs: int,
    *,
    compute_umap: bool = True,
    compute_markers: bool = True,
):
    import scanpy as sc

    sc.pp.normalize_total(data, target_sum=1e4)
    sc.pp.log1p(data)
    sc.pp.highly_variable_genes(data, n_top_genes=min(n_hvgs, data.n_vars), flavor="seurat")
    data = data[:, data.var["highly_variable"].to_numpy()].copy()
    # Preserve log-normalized expression for interpretable marker effect sizes.
    data.raw = data.copy()
    sc.pp.scale(data, max_value=10)
    sc.tl.pca(data, n_comps=min(50, data.n_vars - 1), random_state=17)
    sc.pp.neighbors(data, n_neighbors=15, n_pcs=min(40, data.obsm["X_pca"].shape[1]), random_state=17)
    sc.tl.leiden(data, resolution=resolution, key_added="native_cluster", random_state=17)
    if compute_umap:
        sc.tl.umap(data, random_state=17)
    if compute_markers:
        sc.tl.rank_genes_groups(
            data, groupby="native_cluster", method="wilcoxon", n_genes=25, use_raw=True
        )
    return data


def export_compartment(data, dataset_id: str, broad_label: str, output_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    import matplotlib.pyplot as plt
    import scanpy as sc

    markers = []
    templates = []
    for cluster in sorted(data.obs["native_cluster"].astype(str).unique()):
        table = sc.get.rank_genes_groups_df(data, group=cluster)
        table["dataset_id"] = dataset_id
        table["broad_label"] = broad_label
        table["native_cluster"] = cluster
        markers.append(table)
        counts = data.obs.loc[data.obs["native_cluster"].astype(str).eq(cluster)]
        top_genes = ";".join(table["names"].astype(str).head(10))
        templates.append({
            "dataset_id": dataset_id,
            "broad_label": broad_label,
            "native_cluster": cluster,
            "n_cells": int(len(counts)),
            "n_samples": int(counts["analysis_sample"].nunique()),
            "top_marker_genes": top_genes,
            "provisional_state_name": "",
            "marker_rationale": "",
            "reviewer": "",
            "review_status": "pending_review",
        })
    marker_table = pd.concat(markers, ignore_index=True)
    summary = (
        data.obs.groupby("native_cluster", observed=True)
        .agg(n_cells=("native_cluster", "size"), n_samples=("analysis_sample", "nunique"))
        .reset_index()
        .assign(dataset_id=dataset_id, broad_label=broad_label)
    )
    figure = sc.pl.umap(
        data,
        color=["native_cluster", "analysis_sample"],
        show=False,
        return_fig=True,
        title=[f"{dataset_id}: {broad_label} native clusters", f"{dataset_id}: sample"],
    )
    stem = f"{dataset_id}__{broad_label.replace('/', '_').replace(' ', '_')}"
    figure.savefig(output_dir / f"{stem}_native_umap.png", dpi=300, bbox_inches="tight")
    figure.savefig(output_dir / f"{stem}_native_umap.pdf", bbox_inches="tight")
    plt.close(figure)
    assignments = data.obs.loc[:, [
        "dataset_id", "frozen_broad_label", "analysis_sample", "source_row_index", "source_cell_id", "native_cluster",
    ]].copy()
    assignments = assignments.rename(columns={"frozen_broad_label": "broad_label"})
    return marker_table, summary, pd.DataFrame(templates), assignments


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--reference-atlas", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--max-cells-per-cohort", type=int, default=30000)
    parser.add_argument("--min-cells-per-compartment", type=int, default=100)
    parser.add_argument("--resolution", type=float, default=0.6)
    parser.add_argument("--n-hvgs", type=int, default=3000)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    manifest = pd.read_csv(args.manifest).set_index("dataset_id", drop=False)
    missing = sorted(set(CORE_IDS) - set(manifest.index.astype(str)))
    if missing:
        raise ValueError(f"Manifest is missing core cohorts: {missing}")
    reference_labels = frozen_broad_label_lookup(args.reference_atlas)
    input_rows, marker_tables, cluster_summaries, templates, assignment_tables, status_rows = [], [], [], [], [], []

    for ordinal, dataset_id in enumerate(CORE_IDS):
        print(f"Loading native cohort: {dataset_id}", flush=True)
        data, input_row = load_native_cohort(
            dataset_id, manifest.loc[dataset_id], args.data_root, reference_labels,
            args.max_cells_per_cohort, 17 + ordinal,
        )
        input_rows.append(input_row)
        for broad_label in sorted(data.obs["frozen_broad_label"].astype(str).unique()):
            subset = data[data.obs["frozen_broad_label"].astype(str).eq(broad_label)].copy()
            if subset.n_obs < args.min_cells_per_compartment:
                status_rows.append({
                    "dataset_id": dataset_id, "broad_label": broad_label,
                    "n_cells": int(subset.n_obs), "status": "skipped_insufficient_cells",
                })
                continue
            print(f"  {broad_label}: {subset.n_obs} cells", flush=True)
            clustered = native_cluster(subset, args.resolution, args.n_hvgs)
            markers, summary, template, assignments = export_compartment(clustered, dataset_id, broad_label, args.output_dir)
            marker_tables.append(markers)
            cluster_summaries.append(summary)
            templates.append(template)
            assignment_tables.append(assignments)
            status_rows.append({
                "dataset_id": dataset_id,
                "broad_label": broad_label,
                "n_cells": int(subset.n_obs),
                "n_native_clusters": int(clustered.obs["native_cluster"].nunique()),
                "status": "clustered_native_raw_counts",
            })
            del clustered, subset
            gc.collect()
        del data
        gc.collect()
        write_progress_snapshot(
            args.output_dir, dataset_id, input_rows, marker_tables, cluster_summaries,
            templates, status_rows,
        )
        if assignment_tables:
            pd.concat(assignment_tables, ignore_index=True).to_csv(
                args.output_dir / "NATIVE_PER_COHORT_CELL_ASSIGNMENTS_PARTIAL.csv.gz", index=False, compression="gzip"
            )

    pd.DataFrame(input_rows).to_csv(args.output_dir / "NATIVE_PER_COHORT_INPUT_AUDIT.csv", index=False)
    pd.DataFrame(status_rows).to_csv(args.output_dir / "NATIVE_PER_COHORT_CLUSTER_STATUS.csv", index=False)
    pd.concat(marker_tables, ignore_index=True).to_csv(args.output_dir / "NATIVE_PER_COHORT_CLUSTER_MARKERS.csv", index=False)
    pd.concat(cluster_summaries, ignore_index=True).to_csv(args.output_dir / "NATIVE_PER_COHORT_CLUSTER_SUMMARY.csv", index=False)
    pd.concat(templates, ignore_index=True).to_csv(args.output_dir / "NATIVE_PER_COHORT_ANNOTATION_TEMPLATE.csv", index=False)
    pd.concat(assignment_tables, ignore_index=True).to_csv(
        args.output_dir / "NATIVE_PER_COHORT_CELL_ASSIGNMENTS.csv.gz", index=False, compression="gzip"
    )
    (args.output_dir / "README_AND_REVIEW_GATE.md").write_text(
        "# Native Per-Cohort State Discovery\n\n"
        "Each cohort and broad compartment was clustered independently from its own audited count matrix. The frozen atlas was used only to identify a broad-compartment correspondence and the retained Korean cell universe. Native clusters are not cross-cohort states. Review `NATIVE_PER_COHORT_ANNOTATION_TEMPLATE.csv` with marker evidence before harmonising labels or testing abundance.\n\n"
        "No tumour-normal, treatment, survival, or cell-cell interaction conclusion may be drawn from these discovery clusters alone. The unit of later inference must be the donor/sample, not the cell.\n",
        encoding="utf-8",
    )
    (args.output_dir / "RUN_MANIFEST.json").write_text(json.dumps({
        "scope": "six core cohorts; native within-cohort compartment clustering",
        "max_cells_per_cohort": args.max_cells_per_cohort,
        "min_cells_per_compartment": args.min_cells_per_compartment,
        "leiden_resolution": args.resolution,
        "n_hvgs": args.n_hvgs,
        "cell_assignment_table": "NATIVE_PER_COHORT_CELL_ASSIGNMENTS.csv.gz",
        "reference_atlas_used_only_for_broad_correspondence": str(args.reference_atlas),
    }, indent=2) + "\n", encoding="utf-8")
    print(f"Native per-cohort state discovery written to {args.output_dir}")


if __name__ == "__main__":
    main()
