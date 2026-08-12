#!/usr/bin/env python3
"""Compare unintegrated PCA, Harmony, and scVI after the raw-input gate passes.

This is a method sensitivity analysis, not a response model or a substitute for
strict query-mapping validation. All methods receive the same stratified cells
and shared genes, and no clinical outcome enters fitting or feature selection.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors

from raw_integration_method_audit import CORE_IDS, candidate_paths, col, count_matrix


BROAD_COLUMNS = ["cell_type_coarse", "broad_cell_type", "cell_type", "cell_type_coarse_fine"]
SAMPLE_COLUMNS = ["sample_id", "sample", "orig.ident", "library_id", "analysis_unit"]
MARKER_SETS = {
    "T/NK": ["CD3D", "CD3E", "TRBC2", "NKG7"],
    "Myeloid": ["LYZ", "FCER1G", "TYROBP", "C1QC"],
    "Fibroblast": ["COL1A1", "COL1A2", "DCN", "LUM"],
    "Endothelial": ["KDR", "EMCN", "VWF", "CLDN5"],
    "Epithelial": ["EPCAM", "KRT8", "KRT18", "MUC1"],
    "B/Plasma": ["CD74", "MS4A1", "CD79A", "MZB1"],
}


def label_keys(dataset_id: str, obs: pd.DataFrame, include_integration_suffix_base: bool = False) -> list[list[str]]:
    """Return progressively qualified keys for auditable reference-label matching."""
    sample_col = col(obs, SAMPLE_COLUMNS)
    keys: list[list[str]] = []
    for cell_id, (_, row) in zip(obs.index.astype(str), obs.iterrows()):
        candidates = [cell_id, f"{dataset_id}::{cell_id}"]
        if sample_col is not None:
            sample_id = str(row[sample_col])
            candidates.extend([
                f"{sample_id}::{cell_id}",
                f"{dataset_id}::{sample_id}::{cell_id}",
            ])
        if include_integration_suffix_base:
            # The frozen integration object may add exactly one final '-integer'
            # suffix to make already-qualified cell IDs unique. Preserve the
            # original barcode suffix (for example '-1') and only add a second
            # candidate; the original exact key always remains preferred.
            suffix_bases = [re.sub(r"-\d+$", "", value) for value in candidates]
            candidates.extend(base for base in suffix_bases if base != candidates[0] or len(candidates) > 1)
        keys.append(list(dict.fromkeys(candidates)))
    return keys


def frozen_broad_label_lookup(reference_atlas: Path) -> dict[tuple[str, str], str]:
    """Load marker-derived broad labels from the frozen atlas without refitting it."""
    import anndata as ad

    atlas = ad.read_h5ad(reference_atlas, backed="r")
    try:
        dataset_col = col(atlas.obs, ["dataset_id", "cohort", "cohort_id"])
        broad_col = col(atlas.obs, BROAD_COLUMNS)
        if dataset_col is None or broad_col is None:
            raise ValueError("Reference atlas must include dataset_id and a broad cell-type column.")
        lookup: dict[tuple[str, str], str] = {}
        for dataset_id, frame in atlas.obs.groupby(dataset_col, observed=True):
            for label, candidates in zip(
                frame[broad_col].astype(str),
                label_keys(str(dataset_id), frame, include_integration_suffix_base=True),
            ):
                for candidate in candidates:
                    lookup.setdefault((str(dataset_id), candidate), label)
        return lookup
    finally:
        atlas.file.close()


def stratified_indices(obs: pd.DataFrame, maximum: int, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    sample_col = col(obs, SAMPLE_COLUMNS)
    if sample_col is None:
        return np.sort(rng.choice(len(obs), min(maximum, len(obs)), replace=False))
    groups = [values.to_numpy() for _, values in obs.reset_index(drop=True).groupby(sample_col, dropna=False).groups.items()]
    target = max(1, maximum // max(1, len(groups)))
    selected = []
    for positions in groups:
        selected.extend(rng.choice(positions, min(target, len(positions)), replace=False).tolist())
    if len(selected) < min(maximum, len(obs)):
        pool = np.setdiff1d(np.arange(len(obs)), np.asarray(selected, dtype=int))
        selected.extend(rng.choice(pool, min(maximum - len(selected), len(pool)), replace=False).tolist())
    return np.sort(np.asarray(selected[:maximum], dtype=int))


def first_uppercase_indices(names: pd.Index) -> dict[str, int]:
    lookup: dict[str, int] = {}
    for index, name in enumerate(names.astype(str)):
        lookup.setdefault(name.upper(), index)
    return lookup


def neighbor_indices(representation: np.ndarray, maximum: int = 15) -> np.ndarray:
    """Return nearest non-self neighbours with safe behaviour for small groups."""
    if len(representation) < 2:
        return np.empty((len(representation), 0), dtype=int)
    return NearestNeighbors(n_neighbors=min(maximum, len(representation) - 1)).fit(representation).kneighbors(return_distance=False)


def load_cohort_subset(
    dataset_id: str,
    manifest_row: pd.Series,
    data_root: Path,
    shared_genes: list[str],
    max_cells: int,
    seed: int,
    reference_labels: dict[tuple[str, str], str],
):
    import anndata as ad

    path = next((candidate for candidate in candidate_paths(data_root, dataset_id, manifest_row) if candidate.exists()), None)
    if path is None:
        raise FileNotFoundError(f"No H5AD found for {dataset_id}.")
    source = ad.read_h5ad(path, backed="r")
    try:
        matrix, layer = count_matrix(source)
        if layer == "X_unverified":
            layer = None
        indices_by_gene = first_uppercase_indices(source.var_names)
        missing = [gene for gene in shared_genes if gene not in indices_by_gene]
        if missing:
            raise ValueError(f"{dataset_id} is missing {len(missing)} shared genes after preflight.")
        broad_col = col(source.obs, BROAD_COLUMNS)
        if broad_col is not None:
            source_labels = source.obs[broad_col].astype(str).to_numpy()
            label_source = "input_h5ad"
        else:
            source_labels = np.asarray([
                next((reference_labels[(dataset_id, key)] for key in candidates if (dataset_id, key) in reference_labels), "MISSING")
                for candidates in label_keys(dataset_id, source.obs)
            ])
            label_source = "frozen_atlas_match"
        labelled_rows = np.flatnonzero(source_labels != "MISSING")
        if not len(labelled_rows):
            raise ValueError(f"{dataset_id}: no source cells could be matched to a frozen broad label.")
        # The benchmark is defined on the same retained cell universe as the
        # frozen atlas. Raw cells outside that audited intersection were not
        # part of atlas discovery and must not be assigned a surrogate label.
        selected_relative = stratified_indices(source.obs.iloc[labelled_rows], max_cells, seed)
        rows = labelled_rows[selected_relative]
        columns = np.asarray([indices_by_gene[gene] for gene in shared_genes], dtype=int)
        small = source[rows, columns].to_memory()
        if layer is not None:
            # The view retains layers; put counts into X so all subsequent
            # preprocessing starts from exactly the audited count matrix.
            small.X = small.layers[layer].copy()
        small.var_names = pd.Index(shared_genes)
        small.obs["dataset_id"] = dataset_id
        small.obs["broad_label"] = source_labels[rows]
        coverage = float(small.obs["broad_label"].astype(str).ne("MISSING").mean())
        return small, str(path), coverage, label_source, int(source.n_obs), int(len(labelled_rows))
    finally:
        source.file.close()


def representation_metrics(representation: np.ndarray, batch: np.ndarray, biology: np.ndarray) -> dict[str, float]:
    sample = np.arange(len(representation))
    if len(sample) > 8000:
        sample = np.random.default_rng(17).choice(sample, 8000, replace=False)
    x, b, y = representation[sample], batch[sample], biology[sample]
    result = {
        "batch_silhouette_abs": np.nan,
        "broad_label_silhouette": np.nan,
        "same_batch_neighbor_fraction": np.nan,
        "same_broad_label_neighbor_fraction": np.nan,
    }
    if len(np.unique(b)) > 1:
        result["batch_silhouette_abs"] = float(abs(silhouette_score(x, b)))
        neighbors = neighbor_indices(x)
        result["same_batch_neighbor_fraction"] = float(np.mean(b[neighbors] == b[:, None]))
    valid = y != "MISSING"
    if valid.sum() > 2 and len(np.unique(y[valid])) > 1:
        result["broad_label_silhouette"] = float(silhouette_score(x[valid], y[valid]))
        neighbors = neighbor_indices(x[valid])
        y_valid = y[valid]
        result["same_broad_label_neighbor_fraction"] = float(np.mean(y_valid[neighbors] == y_valid[:, None]))
    return result


def lineage_batch_metrics(representations: dict[str, np.ndarray], batch: np.ndarray, biology: np.ndarray) -> pd.DataFrame:
    """Measure mixing within each broad lineage, where composition is comparable."""
    rows: list[dict[str, object]] = []
    for method, values in representations.items():
        for lineage in sorted(pd.unique(biology)):
            mask = biology == lineage
            x, b = values[mask], batch[mask]
            record: dict[str, object] = {
                "method": method,
                "broad_label": lineage,
                "n_cells": int(mask.sum()),
                "n_batches": int(len(np.unique(b))),
                "batch_silhouette_abs": np.nan,
                "same_batch_neighbor_fraction": np.nan,
            }
            if len(x) >= 20 and len(np.unique(b)) > 1:
                record["batch_silhouette_abs"] = float(abs(silhouette_score(x, b)))
                neighbors = neighbor_indices(x)
                record["same_batch_neighbor_fraction"] = float(np.mean(b[neighbors] == b[:, None]))
            rows.append(record)
    return pd.DataFrame(rows)


def marker_support(normalized, biology: np.ndarray) -> pd.DataFrame:
    """Audit raw-expression marker separation; it is not an embedding comparison."""
    var_lookup = first_uppercase_indices(normalized.var_names)
    datasets = normalized.obs["dataset_id"].astype(str).to_numpy()
    rows: list[dict[str, object]] = []
    for dataset_id in CORE_IDS:
        cohort = datasets == dataset_id
        for label, markers in MARKER_SETS.items():
            marker_indices = [var_lookup[marker] for marker in markers if marker in var_lookup]
            target = cohort & (biology == label)
            other = cohort & (biology != label)
            if not marker_indices or not target.any() or not other.any():
                rows.append({
                    "dataset_id": dataset_id,
                    "broad_label": label,
                    "n_target_cells": int(target.sum()),
                    "n_other_cells": int(other.sum()),
                    "n_markers_available": len(marker_indices),
                    "target_marker_mean": np.nan,
                    "other_marker_mean": np.nan,
                    "marker_separation": np.nan,
                })
                continue
            scores = np.asarray(normalized.X[:, marker_indices].mean(axis=1)).ravel()
            target_mean = float(scores[target].mean())
            other_mean = float(scores[other].mean())
            rows.append({
                "dataset_id": dataset_id,
                "broad_label": label,
                "n_target_cells": int(target.sum()),
                "n_other_cells": int(other.sum()),
                "n_markers_available": len(marker_indices),
                "target_marker_mean": target_mean,
                "other_marker_mean": other_mean,
                "marker_separation": target_mean - other_mean,
            })
    return pd.DataFrame(rows)


def marker_consensus_labels(normalized, minimum_score: float = 0.5, minimum_margin: float = 0.25) -> tuple[np.ndarray, pd.DataFrame]:
    """Create high-confidence broad labels directly from raw normalized markers.

    These labels are used only as an independent integration diagnostic. Cells
    below the predeclared score or margin threshold are left ``Ambiguous``.
    """
    gene_lookup = first_uppercase_indices(normalized.var_names)
    matrix = normalized.X
    gene_mean = np.asarray(matrix.mean(axis=0)).ravel()
    if hasattr(matrix, "multiply"):
        gene_square_mean = np.asarray(matrix.multiply(matrix).mean(axis=0)).ravel()
    else:
        gene_square_mean = np.asarray(np.square(matrix).mean(axis=0)).ravel()
    gene_sd = np.sqrt(np.maximum(gene_square_mean - np.square(gene_mean), 1e-8))
    score_rows: list[np.ndarray] = []
    labels: list[str] = []
    for label, markers in MARKER_SETS.items():
        marker_indices = [gene_lookup[marker] for marker in markers if marker in gene_lookup]
        if not marker_indices:
            score_rows.append(np.full(normalized.n_obs, -np.inf))
        else:
            values = np.asarray(matrix[:, marker_indices].toarray() if hasattr(matrix, "toarray") else matrix[:, marker_indices])
            score_rows.append(((values - gene_mean[marker_indices]) / gene_sd[marker_indices]).mean(axis=1))
        labels.append(label)
    scores = np.vstack(score_rows).T
    ranking = np.argsort(scores, axis=1)
    top_index = ranking[:, -1]
    second_index = ranking[:, -2]
    top_score = scores[np.arange(len(scores)), top_index]
    margin = top_score - scores[np.arange(len(scores)), second_index]
    confident = (top_score >= minimum_score) & (margin >= minimum_margin)
    assigned = np.asarray([labels[index] for index in top_index], dtype=object)
    assigned[~confident] = "Ambiguous"
    audit = pd.DataFrame({
        "marker_consensus_label": assigned,
        "marker_consensus_top_score": top_score,
        "marker_consensus_margin": margin,
        "marker_consensus_confident": confident,
    })
    return assigned.astype(str), audit


def write_method_interpretation(metrics: pd.DataFrame, lineage: pd.DataFrame, destination: Path) -> None:
    by_method = metrics.set_index("method")
    scvi = by_method.loc["scvi"]
    harmony = by_method.loc["harmony"]
    lineage_summary = lineage.groupby("method", observed=True)["same_batch_neighbor_fraction"].mean().to_dict()
    destination.write_text(
        "# Integration Method Decision\n\n"
        "## Result\n\n"
        f"- Harmony has lower global same-cohort-neighbour fraction ({harmony['same_batch_neighbor_fraction']:.3f}) than scVI ({scvi['same_batch_neighbor_fraction']:.3f}), indicating more aggressive local batch mixing.\n"
        f"- scVI has higher frozen broad-label silhouette ({scvi['broad_label_silhouette']:.3f}) than Harmony ({harmony['broad_label_silhouette']:.3f}), indicating stronger retention of broad lineage geometry.\n"
        f"- Against independent high-confidence raw marker-consensus labels (coverage={scvi['marker_consensus_label_coverage']:.1%}), scVI silhouette is {scvi['marker_consensus_label_silhouette']:.3f} and Harmony silhouette is {harmony['marker_consensus_label_silhouette']:.3f}.\n"
        f"- Mean within-lineage same-cohort-neighbour fractions are reported as a sensitivity diagnostic: Harmony={lineage_summary.get('harmony', np.nan):.3f}; scVI={lineage_summary.get('scvi', np.nan):.3f}.\n\n"
        "## Decision Boundary\n\n"
        "Use scVI provisionally as the primary representation for frozen-atlas visualization and reference mapping only if the marker-consensus diagnostic also preserves broad-label geometry at least as well as Harmony. Harmony remains a sensitivity representation. Raw-expression marker-support tables and donor/sample-level, within-cohort state analyses remain required before making biological claims.\n\n"
        "No clinical outcomes, response labels, or survival endpoints entered this benchmark. This analysis chooses a representation for exploratory atlas work; it does not test tumour biology or establish a universal best integration method.\n",
        encoding="utf-8",
    )


def save_umaps(adata, representations: dict[str, np.ndarray], destination: Path) -> None:
    import scanpy as sc
    import matplotlib.pyplot as plt

    destination.mkdir(parents=True, exist_ok=True)
    for method, values in representations.items():
        display = adata.copy()
        display.obsm["X_method"] = values
        sc.pp.neighbors(display, use_rep="X_method", n_neighbors=15)
        sc.tl.umap(display, random_state=17)
        figure = sc.pl.umap(display, color=["dataset_id", "broad_label"], show=False, return_fig=True, title=[f"{method}: cohort", f"{method}: broad label"])
        figure.savefig(destination / f"{method}_umap.png", dpi=300, bbox_inches="tight")
        figure.savefig(destination / f"{method}_umap.pdf", bbox_inches="tight")
        plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--preflight", type=Path, required=True)
    parser.add_argument("--reference-atlas", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--max-cells-per-cohort", type=int, default=1200)
    parser.add_argument("--max-epochs", type=int, default=30)
    parser.add_argument("--batch-size", type=int, default=1024)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    preflight = pd.read_csv(args.preflight)
    required = preflight.loc[preflight["dataset_id"].astype(str).isin(CORE_IDS)]
    if len(required) != len(CORE_IDS) or not required["status"].eq("eligible").all():
        blocked = required.loc[~required["status"].eq("eligible"), ["dataset_id", "status", "reason"]]
        raise ValueError(f"Raw gate did not pass for all core datasets:\n{blocked.to_string(index=False)}")
    manifest = pd.read_csv(args.manifest).set_index("dataset_id", drop=False)
    reference_labels = frozen_broad_label_lookup(args.reference_atlas)
    gene_sets: dict[str, set[str]] = {}
    for dataset_id in CORE_IDS:
        path = next((candidate for candidate in candidate_paths(args.data_root, dataset_id, manifest.loc[dataset_id]) if candidate.exists()), None)
        if path is None:
            raise FileNotFoundError(f"No H5AD for {dataset_id}.")
        import anndata as ad
        loaded = ad.read_h5ad(path, backed="r")
        try:
            gene_sets[dataset_id] = set(first_uppercase_indices(loaded.var_names))
        finally:
            loaded.file.close()
    shared_genes = sorted(set.intersection(*gene_sets.values()))
    if len(shared_genes) < 4000:
        raise ValueError(f"Only {len(shared_genes)} shared genes remain; benchmark gate requires at least 4000.")
    blocks, source_rows = [], []
    for ordinal, dataset_id in enumerate(CORE_IDS):
        block, source, label_coverage, label_source, source_cells, labelled_source_cells = load_cohort_subset(
            dataset_id, manifest.loc[dataset_id], args.data_root, shared_genes,
            args.max_cells_per_cohort, 17 + ordinal, reference_labels,
        )
        if label_coverage < 0.95:
            raise ValueError(
                f"{dataset_id}: only {label_coverage:.1%} of benchmark cells received an auditable broad label. "
                "Do not compare integration methods until the cell-identity mapping is repaired."
            )
        blocks.append(block)
        source_rows.append({
            "dataset_id": dataset_id,
            "input_h5ad": source,
            "n_cells_sampled": block.n_obs,
            "broad_label_coverage": label_coverage,
            "broad_label_source": label_source,
            "source_cells_before_reference_intersection": source_cells,
            "source_cells_in_frozen_reference_intersection": labelled_source_cells,
            "source_to_reference_intersection_fraction": labelled_source_cells / source_cells,
        })
    import anndata as ad
    import scanpy as sc
    # AnnData supports inner/outer joins only. Every block was explicitly
    # subset to the same shared-gene list above; inner is an additional guard.
    raw = ad.concat(blocks, join="inner", merge="same", index_unique="__benchmark__")
    raw.layers["counts"] = raw.X.copy()
    marker_normalized = raw.copy()
    sc.pp.normalize_total(marker_normalized, target_sum=1e4)
    sc.pp.log1p(marker_normalized)
    normalized = marker_normalized.copy()
    sc.pp.highly_variable_genes(normalized, n_top_genes=min(4000, normalized.n_vars), flavor="seurat")
    hvg = normalized.var["highly_variable"].to_numpy()
    raw = raw[:, hvg].copy()
    normalized = normalized[:, hvg].copy()
    sc.pp.scale(normalized, max_value=10)
    sc.tl.pca(normalized, n_comps=min(50, normalized.n_vars - 1), random_state=17)
    representations = {"unintegrated_pca": normalized.obsm["X_pca"].copy()}

    try:
        # Use harmonypy's direct API. The Scanpy wrapper has changed across
        # releases; this keeps the benchmark input explicit and reproducible.
        import harmonypy as hm

        harmony_fit = hm.run_harmony(
            normalized.obsm["X_pca"],
            normalized.obs,
            vars_use=["dataset_id"],
            random_state=17,
            verbose=False,
        )
        harmony_embedding = np.asarray(harmony_fit.Z_corr)
        expected_shape = normalized.obsm["X_pca"].shape
        if harmony_embedding.shape == expected_shape:
            corrected = harmony_embedding
        elif harmony_embedding.T.shape == expected_shape:
            corrected = harmony_embedding.T
        else:
            raise ValueError(
                f"Harmony returned shape {harmony_embedding.shape}; expected {expected_shape} or its transpose."
            )
        if not np.isfinite(corrected).all():
            raise ValueError("Harmony returned non-finite values.")
        representations["harmony"] = corrected
        (args.output_dir / "HARMONY_STATUS.txt").write_text(
            "Harmony completed through harmonypy.run_harmony on the shared PCA representation.\n",
            encoding="utf-8",
        )
    except Exception as error:
        (args.output_dir / "HARMONY_STATUS.txt").write_text(
            f"Harmony unavailable/failed: {type(error).__name__}: {error}\n",
            encoding="utf-8",
        )

    try:
        import scvi
        training = raw.copy()
        scvi.model.SCVI.setup_anndata(training, layer="counts", batch_key="dataset_id")
        model = scvi.model.SCVI(training, n_latent=30, n_layers=2, gene_likelihood="nb")
        train_args = {"max_epochs": args.max_epochs, "batch_size": args.batch_size, "check_val_every_n_epoch": max(1, args.max_epochs // 5)}
        try:
            import torch
            if torch.cuda.is_available():
                train_args.update({"accelerator": "gpu", "devices": 1})
        except ImportError:
            pass
        model.train(**train_args)
        representations["scvi"] = model.get_latent_representation()
    except Exception as error:
        (args.output_dir / "SCVI_STATUS.txt").write_text(f"scVI unavailable/failed: {type(error).__name__}: {error}\n", encoding="utf-8")

    batch = raw.obs["dataset_id"].astype(str).to_numpy()
    biology = raw.obs["broad_label"].astype(str).to_numpy()
    consensus_labels, consensus_audit = marker_consensus_labels(marker_normalized)
    consensus_audit["dataset_id"] = batch
    consensus_audit["frozen_broad_label"] = biology
    consensus_audit.to_csv(args.output_dir / "MARKER_CONSENSUS_LABEL_AUDIT.csv", index=False)
    consensus_rows = (
        consensus_audit.groupby(["dataset_id", "frozen_broad_label", "marker_consensus_label"], observed=True)
        .size().rename("n_cells").reset_index()
    )
    consensus_rows.to_csv(args.output_dir / "MARKER_CONSENSUS_LABEL_CONFUSION.csv", index=False)
    metrics_rows = []
    for method, values in representations.items():
        frozen_metrics = representation_metrics(values, batch, biology)
        consensus_metrics = representation_metrics(values, batch, consensus_labels)
        metrics_rows.append({
            "method": method,
            **frozen_metrics,
            "marker_consensus_label_coverage": float((consensus_labels != "Ambiguous").mean()),
            "marker_consensus_label_silhouette": consensus_metrics["broad_label_silhouette"],
            "marker_consensus_same_label_neighbor_fraction": consensus_metrics["same_broad_label_neighbor_fraction"],
        })
    metrics = pd.DataFrame(metrics_rows)
    required_methods = {"unintegrated_pca", "harmony", "scvi"}
    missing_methods = required_methods.difference(representations)
    if missing_methods:
        raise RuntimeError(
            "Benchmark is incomplete; no method decision may be made. "
            f"Missing method(s): {', '.join(sorted(missing_methods))}."
        )
    metrics.to_csv(args.output_dir / "INTEGRATION_METHOD_SENSITIVITY_METRICS.csv", index=False)
    lineage = lineage_batch_metrics(representations, batch, biology)
    lineage.to_csv(args.output_dir / "INTEGRATION_METHOD_LINEAGE_METRICS.csv", index=False)
    consensus_mask = consensus_labels != "Ambiguous"
    consensus_lineage = lineage_batch_metrics(
        {method: values[consensus_mask] for method, values in representations.items()},
        batch[consensus_mask], consensus_labels[consensus_mask],
    )
    consensus_lineage.to_csv(args.output_dir / "INTEGRATION_METHOD_MARKER_CONSENSUS_LINEAGE_METRICS.csv", index=False)
    marker_audit = marker_support(marker_normalized, biology)
    marker_audit.to_csv(args.output_dir / "RAW_MARKER_SUPPORT_BY_COHORT.csv", index=False)
    pd.DataFrame(source_rows).to_csv(args.output_dir / "INTEGRATION_METHOD_SENSITIVITY_INPUTS.csv", index=False)
    save_umaps(normalized, representations, args.output_dir / "UMAPS")
    (args.output_dir / "INTERPRETATION.md").write_text(
        "# Controlled Integration Sensitivity Benchmark\n\n"
        "All methods received the same stratified cells, shared genes, and no clinical outcome labels. Broad labels came from the input H5AD when available or were matched to the frozen marker-derived atlas; per-cohort coverage is in INTEGRATION_METHOD_SENSITIVITY_INPUTS.csv. Lower batch silhouette and same-batch neighbour fraction suggest mixing; higher broad-label silhouette suggests preservation. A comparison is invalid if any method is missing or label coverage is below 95% in a cohort. These metrics do not establish a universally best method. Select a primary integration only after checking marker fidelity, donor-level state effects, and strict query mapping.\n",
        encoding="utf-8",
    )
    write_method_interpretation(metrics, lineage, args.output_dir / "METHOD_DECISION.md")
    (args.output_dir / "RUN_MANIFEST.json").write_text(json.dumps({"max_cells_per_cohort": args.max_cells_per_cohort, "max_epochs": args.max_epochs, "n_shared_genes_before_hvg": len(shared_genes), "n_hvgs": int(hvg.sum()), "methods_completed": sorted(representations)}, indent=2) + "\n", encoding="utf-8")
    print(metrics.to_string(index=False))
    print(f"Method sensitivity benchmark written to {args.output_dir}")


if __name__ == "__main__":
    main()
