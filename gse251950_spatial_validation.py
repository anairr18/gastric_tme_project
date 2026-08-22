#!/usr/bin/env python3
"""Validation-only mapping of frozen curated states to GSE251950 Visium slides.

Spots are used for mapping, but every statistical test is performed at a
spatial-block level within a slide and then pooled at the slide level. The
external slides never contribute to the state reference or scVI atlas.
"""

from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, stats
from scipy.io import mmread
from scipy.spatial import cKDTree
from statsmodels.stats.multitest import multipletests


PRIORITY_STATES = {
    "PDPN_CTSK_CAF": "Fibroblast::PDPN_CTSK_inflammatory_CAF_candidate",
    "PLVAP_CA4_capillary_endothelial": "Endothelial::PLVAP_CA4_capillary_endothelial",
    "GZMK_effector_memory_T": "T/NK::GZMK_effector_memory_T",
}
MYELOID_HYPOXIA = ("HIF1A", "VEGFA", "CA9", "EGLN3", "ADM", "PDK1", "LDHA", "SLC2A1")
RELATIONSHIPS = {
    "CAF_to_capillary": ("PDPN_CTSK_CAF", "PLVAP_CA4_capillary_endothelial"),
    "CAF_to_GZMK_T": ("PDPN_CTSK_CAF", "GZMK_effector_memory_T"),
}


def bh(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.dropna()
    if not valid.empty:
        result.loc[valid.index] = multipletests(valid, method="fdr_bh")[1]
    return result


def read_visium_matrix(path: Path):
    if path.is_file() and path.suffix == ".h5":
        adata = sc.read_10x_h5(path)
    elif path.is_file() and "matrix.mtx" in path.name.lower():
        prefix = re.split(r"matrix\.mtx", path.name, maxsplit=1, flags=re.IGNORECASE)[0]
        def component(tokens: tuple[str, ...]) -> Path:
            matches = []
            for token in tokens:
                matches.extend(path.parent.glob(f"{prefix}*{token}*.tsv*"))
            if not matches:
                raise FileNotFoundError(f"No {'/'.join(tokens)} file found next to {path.name}")
            return sorted(set(matches))[0]
        features = pd.read_csv(component(("features", "genes")), sep="\t", header=None)
        barcodes = pd.read_csv(component(("barcodes",)), sep="\t", header=None).iloc[:, 0].astype(str)
        matrix = mmread(path).tocsr().T
        if matrix.shape != (len(barcodes), len(features)):
            raise ValueError(f"MTX dimension mismatch for {path.name}: {matrix.shape}, {len(barcodes)} barcodes, {len(features)} features")
        genes = features.iloc[:, 1] if features.shape[1] > 1 else features.iloc[:, 0]
        adata = sc.AnnData(X=matrix)
        adata.obs_names = pd.Index(barcodes)
        adata.var_names = pd.Index(genes.astype(str))
    elif path.is_dir():
        adata = sc.read_10x_mtx(path, var_names="gene_symbols", cache=False)
    else:
        raise FileNotFoundError(f"Unsupported Visium count input: {path}")
    adata.var_names = adata.var_names.astype(str).str.upper()
    adata.var_names_make_unique()
    adata.obs_names = adata.obs_names.astype(str)
    return adata


def read_coordinates(path: Path, barcodes: pd.Index) -> pd.DataFrame:
    raw = pd.read_csv(path, header=None)
    if raw.empty or raw.shape[1] < 3:
        raise ValueError(f"Spatial coordinate file is malformed: {path}")
    if str(raw.iloc[0, 0]).lower() in {"barcode", "barcodes"}:
        raw = raw.iloc[1:].copy()
    raw = raw.reset_index(drop=True)
    coordinates = pd.DataFrame(
        {
            "barcode": raw.iloc[:, 0].astype(str),
            "y": pd.to_numeric(raw.iloc[:, -2], errors="coerce"),
            "x": pd.to_numeric(raw.iloc[:, -1], errors="coerce"),
        }
    ).dropna()
    coordinates = coordinates.drop_duplicates("barcode").set_index("barcode")
    coordinates = coordinates.reindex(barcodes)
    if coordinates[["x", "y"]].isna().any().any():
        matched = int(coordinates.dropna().shape[0])
        raise ValueError(f"Only {matched}/{len(barcodes)} Visium barcodes have usable spatial coordinates.")
    return coordinates


def log_cp10k(matrix):
    matrix = matrix.copy().tocsr() if sparse.issparse(matrix) else np.asarray(matrix, dtype=float).copy()
    totals = np.asarray(matrix.sum(axis=1)).ravel() if sparse.issparse(matrix) else matrix.sum(axis=1)
    totals = np.maximum(totals, 1.0)
    if sparse.issparse(matrix):
        matrix = sparse.diags(1e4 / totals) @ matrix
        matrix.data = np.log1p(matrix.data)
    else:
        matrix = np.log1p(matrix * (1e4 / totals)[:, None])
    return matrix


def load_reference(path: Path) -> pd.DataFrame:
    reference = pd.read_csv(path, compression="infer")
    required = {"state_id", "gene", "mean_log_normalized_expression", "n_reference_cells"}
    missing = required - set(reference.columns)
    if missing:
        raise ValueError(f"Frozen reference profile missing columns: {sorted(missing)}")
    reference["gene"] = reference["gene"].astype(str).str.upper()
    pivot = reference.pivot_table(
        index="state_id", columns="gene", values="mean_log_normalized_expression", aggfunc="mean"
    ).fillna(0.0)
    missing_priority = set(PRIORITY_STATES.values()) - set(pivot.index)
    if missing_priority:
        raise ValueError(f"Frozen reference lacks prespecified states: {sorted(missing_priority)}")
    return pivot


def projected_nnls(spot_by_gene: np.ndarray, state_by_gene: np.ndarray, iterations: int = 250) -> np.ndarray:
    """Batched projected-gradient NNLS for many Visium spots."""
    profiles = state_by_gene.copy()
    gene_scale = np.maximum(profiles.std(axis=0), 0.05)
    profiles /= gene_scale
    observations = spot_by_gene / gene_scale
    gram = profiles @ profiles.T
    cross = observations @ profiles.T
    scale = float(np.linalg.eigvalsh(gram).max()) if gram.size else 1.0
    weights = np.maximum(cross / max(scale, 1e-6), 0.0)
    step = 1.0 / max(scale, 1e-6)
    for _ in range(iterations):
        weights = np.maximum(weights - step * (weights @ gram - cross), 0.0)
    totals = weights.sum(axis=1, keepdims=True)
    return np.divide(weights, totals, out=np.zeros_like(weights), where=totals > 0)


def map_slide(adata, coordinates: pd.DataFrame, reference: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, set[str]]:
    shared = sorted(set(adata.var_names) & set(reference.columns))
    variable = reference[shared].var(axis=0).sort_values(ascending=False)
    shared = variable.index[variable.gt(1e-5)].tolist()[:750]
    if len(shared) < 200:
        raise ValueError(f"Only {len(shared)} informative genes overlap the frozen state reference.")
    gene_positions = pd.Index(adata.var_names).get_indexer(shared)
    expression = log_cp10k(adata[:, gene_positions].X)
    expression = expression.toarray() if sparse.issparse(expression) else expression
    weights = projected_nnls(expression, reference.loc[:, shared].to_numpy(float))
    result = pd.DataFrame(weights, index=adata.obs_names, columns=reference.index)
    for label, state_id in PRIORITY_STATES.items():
        result[label] = result[state_id]
    available_hypoxia = [gene for gene in MYELOID_HYPOXIA if gene in adata.var_names]
    if len(available_hypoxia) < 4:
        raise ValueError(f"Hypoxia signature has only {len(available_hypoxia)}/8 genes on the slide.")
    hypoxia_positions = pd.Index(adata.var_names).get_indexer(available_hypoxia)
    hypoxia = log_cp10k(adata[:, hypoxia_positions].X)
    result["myeloid_hypoxia_score"] = np.asarray(hypoxia.mean(axis=1)).ravel()
    result = result.join(coordinates)
    coverage = pd.DataFrame(
        {
            "feature": [*PRIORITY_STATES, "myeloid_hypoxia_score"],
            "n_requested_genes": [len(shared)] * len(PRIORITY_STATES) + [len(MYELOID_HYPOXIA)],
            "n_available_genes": [len(shared)] * len(PRIORITY_STATES) + [len(available_hypoxia)],
        }
    )
    coverage["coverage_fraction"] = coverage.n_available_genes / coverage.n_requested_genes
    return result, coverage, set(adata.var_names.astype(str))


def spatial_blocks(frame: pd.DataFrame, bins: int = 6) -> np.ndarray:
    x = pd.qcut(frame["x"], q=min(bins, frame["x"].nunique()), labels=False, duplicates="drop")
    y = pd.qcut(frame["y"], q=min(bins, frame["y"].nunique()), labels=False, duplicates="drop")
    return (x.astype(str) + "_" + y.astype(str)).to_numpy()


def spatial_effect(frame: pd.DataFrame, source: str, target: str, seed: int, n_permutations: int) -> dict[str, float]:
    points = frame[["x", "y"]].to_numpy(float)
    tree = cKDTree(points)
    k = min(7, len(frame))
    _, neighbors = tree.query(points, k=k)
    neighbors = neighbors[:, 1:] if k > 1 else np.empty((len(frame), 0), dtype=int)
    if neighbors.shape[1] < 3:
        raise ValueError("Too few spatial neighbors for slide-level relationship test.")
    local_target = frame[target].to_numpy(float)[neighbors].mean(axis=1)
    blocks = spatial_blocks(frame)
    block_table = pd.DataFrame({"block": blocks, "source": frame[source].to_numpy(float), "local_target": local_target})
    block_table = block_table.groupby("block", observed=True).mean().dropna()
    if len(block_table) < 8:
        raise ValueError("Fewer than eight nonempty spatial blocks after coordinate binning.")
    if block_table.source.nunique() < 3 or block_table.local_target.nunique() < 3:
        raise ValueError("Spatial-block values are constant; the relationship is not estimable on this slide.")
    observed = float(stats.spearmanr(block_table.source, block_table.local_target).statistic)
    rng = np.random.default_rng(seed)
    permuted = np.asarray(
        [
            stats.spearmanr(block_table.source, rng.permutation(block_table.local_target)).statistic
            for _ in range(n_permutations)
        ]
    )
    p_value = float((1 + np.count_nonzero(np.abs(permuted) >= abs(observed))) / (n_permutations + 1))
    return {
        "spearman_rho": observed,
        "block_permutation_p_value": p_value,
        "n_spatial_blocks": int(len(block_table)),
        "n_spots": int(len(frame)),
    }


def random_effects(frame: pd.DataFrame) -> dict[str, float]:
    frame = frame.loc[np.isfinite(frame.spearman_rho)].copy()
    if len(frame) < 2:
        return {
            "n_slides": int(len(frame)), "pooled_rho": np.nan, "ci_low": np.nan, "ci_high": np.nan,
            "p_value": np.nan, "i2_percent": np.nan, "directional_concordance": np.nan,
        }
    values = frame.spearman_rho.to_numpy(float)
    variances = 1.0 / np.maximum(frame.n_spatial_blocks.to_numpy(float) - 3.0, 1.0)
    fixed_weights = 1 / variances
    fixed = float(np.sum(fixed_weights * values) / np.sum(fixed_weights))
    q = float(np.sum(fixed_weights * (values - fixed) ** 2))
    df = len(values) - 1
    c = float(np.sum(fixed_weights) - np.sum(fixed_weights**2) / np.sum(fixed_weights))
    tau2 = max(0.0, (q - df) / c) if c > 0 else 0.0
    weights = 1 / (variances + tau2)
    pooled = float(np.sum(weights * values) / np.sum(weights))
    se = float(math.sqrt(1 / np.sum(weights)))
    return {
        "n_slides": int(len(frame)),
        "pooled_rho": pooled,
        "ci_low": pooled - 1.96 * se,
        "ci_high": pooled + 1.96 * se,
        "p_value": float(2 * stats.norm.sf(abs(pooled / se))),
        "i2_percent": max(0.0, (q - df) / q) * 100 if q > 0 else 0.0,
        "directional_concordance": float(max((values > 0).mean(), (values < 0).mean())),
    }


def plot_maps(slide_id: str, frame: pd.DataFrame, output: Path) -> None:
    features = [*PRIORITY_STATES, "myeloid_hypoxia_score"]
    fig, axes = plt.subplots(1, len(features), figsize=(4.2 * len(features), 3.8))
    for axis, feature in zip(axes, features):
        image = axis.scatter(frame.x, frame.y, c=frame[feature], s=1.2, cmap="magma", linewidths=0)
        axis.set_title(feature.replace("_", " "), fontsize=9)
        axis.set_axis_off()
        fig.colorbar(image, ax=axis, fraction=0.045, pad=0.02)
    fig.suptitle(f"GSE251950 {slide_id}: frozen-state spatial mapping", y=1.02)
    fig.tight_layout()
    fig.savefig(output / f"{slide_id}_FROZEN_STATE_MAPS.png", dpi=300, bbox_inches="tight")
    fig.savefig(output / f"{slide_id}_FROZEN_STATE_MAPS.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_effects(meta: pd.DataFrame, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, max(3.2, 0.9 * len(meta) + 1.2)))
    values = meta.pooled_rho.to_numpy(float)
    y = np.arange(len(meta))
    ax.errorbar(values, y, xerr=[values - meta.ci_low, meta.ci_high - values], fmt="o", color="#1f4e79", capsize=3)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_yticks(y, meta.relationship.str.replace("_", " "))
    ax.set_xlabel("Pooled slide-level spatial-block Spearman correlation")
    ax.set_title("GSE251950 frozen-state spatial associations")
    fig.tight_layout()
    fig.savefig(output / "GSE251950_SPATIAL_ASSOCIATION_FOREST.png", dpi=300, bbox_inches="tight")
    fig.savefig(output / "GSE251950_SPATIAL_ASSOCIATION_FOREST.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--spatial-root", required=True, type=Path)
    parser.add_argument("--reference-profiles", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--n-permutations", type=int, default=499)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    index_path = args.spatial_root / "GSE251950" / "SPATIAL_INPUT_INDEX.csv"
    index = pd.read_csv(index_path)
    reference = load_reference(args.reference_profiles)
    qc_rows, coverage_rows, effects, exclusions, all_weights, available_gene_sets = [], [], [], [], [], []
    for ordinal, row in index.iterrows():
        if row.get("status") != "ready":
            qc_rows.append({"slide_id": row.slide_id, "status": row.status, "reason": "input index blocked"})
            continue
        try:
            adata = read_visium_matrix(Path(row.matrix_path))
            coordinates = read_coordinates(Path(row.coordinates_path), adata.obs_names)
            mapped, coverage, available_genes = map_slide(adata, coordinates, reference)
            mapped.insert(0, "slide_id", row.slide_id)
            mapped.to_csv(args.output_dir / f"{row.slide_id}_FROZEN_STATE_WEIGHTS.csv.gz", compression="gzip")
            plot_maps(str(row.slide_id), mapped, args.output_dir)
            qc_rows.append({"slide_id": row.slide_id, "status": "ready", "reason": "", "n_spots": adata.n_obs, "n_genes": adata.n_vars, "n_reference_genes": int(coverage.n_available_genes.max())})
            coverage.insert(0, "slide_id", row.slide_id)
            coverage_rows.append(coverage)
            all_weights.append(mapped)
            available_gene_sets.append(available_genes)
            for relationship, (source, target) in RELATIONSHIPS.items():
                try:
                    effect = spatial_effect(mapped, source, target, seed=481 + ordinal, n_permutations=args.n_permutations)
                    effects.append({"slide_id": row.slide_id, "relationship": relationship, **effect})
                except ValueError as exc:
                    exclusions.append({"slide_id": row.slide_id, "relationship": relationship, "reason": str(exc)})
        except Exception as exc:
            qc_rows.append({"slide_id": row.slide_id, "status": "blocked", "reason": f"{type(exc).__name__}: {exc}"})
    qc = pd.DataFrame(qc_rows)
    qc.to_csv(args.output_dir / "GSE251950_SPATIAL_SLIDE_QC.csv", index=False)
    index = index.drop(columns=[column for column in ["n_spots", "n_genes", "status"] if column in index]).merge(
        qc[[column for column in ["slide_id", "status", "n_spots", "n_genes"] if column in qc]], on="slide_id", how="left"
    )
    index.to_csv(index_path, index=False)
    if not all_weights:
        raise RuntimeError("No GSE251950 slides passed frozen-state mapping and spatial-coordinate validation.")
    pd.concat(coverage_rows, ignore_index=True).to_csv(args.output_dir / "GSE251950_SPATIAL_GENE_COVERAGE.csv", index=False)
    common_genes = sorted(set.intersection(*available_gene_sets))
    (args.output_dir / "GSE251950_COMMON_SPATIAL_GENES.txt").write_text("\n".join(common_genes) + "\n", encoding="utf-8")
    pd.concat(all_weights, ignore_index=True).to_csv(args.output_dir / "GSE251950_FROZEN_STATE_SPOT_WEIGHTS.csv.gz", index=False, compression="gzip")
    slide_effects = pd.DataFrame(effects)
    if not slide_effects.empty:
        slide_effects["block_permutation_fdr_bh"] = bh(slide_effects.block_permutation_p_value)
    slide_effects.to_csv(args.output_dir / "GSE251950_SLIDE_SPATIAL_EFFECTS.csv", index=False)
    pd.DataFrame(exclusions).to_csv(args.output_dir / "GSE251950_SPATIAL_RELATIONSHIP_EXCLUSIONS.csv", index=False)
    meta_rows = []
    for name, group in slide_effects.groupby("relationship", observed=True):
        meta_rows.append({"relationship": name, **random_effects(group)})
    for name in RELATIONSHIPS:
        if name not in {row["relationship"] for row in meta_rows}:
            meta_rows.append({"relationship": name, "n_slides": 0, "pooled_rho": np.nan, "ci_low": np.nan, "ci_high": np.nan, "p_value": np.nan, "i2_percent": np.nan, "directional_concordance": np.nan})
    meta = pd.DataFrame(meta_rows)
    meta["fdr_bh"] = bh(meta.p_value)
    meta["claim_boundary"] = "Slide-level spatial association after block permutation; not physical interaction or causal signalling."
    meta.to_csv(args.output_dir / "GSE251950_SPATIAL_RANDOM_EFFECTS.csv", index=False)
    plot_effects(meta, args.output_dir)
    print(meta.to_string(index=False))


if __name__ == "__main__":
    main()
