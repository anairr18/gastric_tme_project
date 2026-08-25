#!/usr/bin/env python3
"""Frozen, section-level spatial ecology analysis in GSE308624 CosMx data."""

from __future__ import annotations

import math
import argparse
import os
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.spatial import cKDTree
from statsmodels.stats.multitest import multipletests


BASE = Path(__file__).resolve().parent
DATA_ROOT = Path(os.environ.get("GASTRIC_TME_DATA_ROOT", BASE / "data"))
INPUT = DATA_ROOT / "external" / "spatial" / "GSE308624" / "GSE308624_gastric_cancer_smi.h5ad"
CLINICAL = DATA_ROOT / "external" / "spatial" / "GSE308624" / "GSE308624_gastric_sample_clinical.xlsx"
OUT_DIR = BASE / "outputs" / "SPATIAL_VALIDATION"

SIGNATURES = {
    "SPP1_myeloid": ["SPP1", "APOC1", "APOE", "C1QA", "C1QB", "TREM2", "GPNMB"],
    "C1Q_myeloid": ["C1QA", "C1QB", "C1QC", "APOE"],
    "CAF": ["COL1A1", "COL1A2", "COL3A1", "FAP", "PDGFRA", "LUM", "DCN", "THBS2", "POSTN"],
    "Exhausted_T": ["PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "TOX", "CXCL13"],
    "Cytotoxic_T": ["NKG7", "GNLY", "PRF1", "GZMB", "IFNG"],
    "PDPN_CTSK_CAF": ["PDPN", "CTSK", "COL6A3", "CXCL14"],
    "PLVAP_CA4_capillary_endothelial": ["PLVAP", "CA4", "FLT1", "RBP7"],
    "GZMK_effector_memory_T": ["GZMK", "GZMA", "NKG7", "CCL5"],
}

# Optional, prespecified extension. It is kept outside the default result set
# so existing SPP1 analyses remain exactly reproducible until this panel-coverage
# check and the supporting spatial run are explicitly requested.
HYPOXIA_SIGNATURE = {
    "Hypoxia_myeloid": ["HIF1A", "VEGFA", "CA9", "LDHA", "PDK1", "SLC2A1"],
}

RELATIONSHIPS = {
    "SPP1-myeloid to CAF field": ("SPP1_myeloid", "Myeloid cell", "CAF", "Fibroblast"),
    "C1Q-myeloid to CAF field": ("C1Q_myeloid", "Myeloid cell", "CAF", "Fibroblast"),
    "SPP1-myeloid to exhausted-T field": (
        "SPP1_myeloid",
        "Myeloid cell",
        "Exhausted_T",
        "T_cell",
    ),
    "C1Q-myeloid to exhausted-T field": (
        "C1Q_myeloid",
        "Myeloid cell",
        "Exhausted_T",
        "T_cell",
    ),
    "SPP1-myeloid to cytotoxic-T field": (
        "SPP1_myeloid",
        "Myeloid cell",
        "Cytotoxic_T",
        "T_cell",
    ),
    "CAF to capillary-endothelial field": (
        "PDPN_CTSK_CAF",
        "Fibroblast",
        "PLVAP_CA4_capillary_endothelial",
        "Endothelial",
    ),
    "CAF to GZMK-T field": (
        "PDPN_CTSK_CAF",
        "Fibroblast",
        "GZMK_effector_memory_T",
        "T_cell",
    ),
}

HYPOXIA_RELATIONSHIPS = {
    "Hypoxia-myeloid to CAF field": ("Hypoxia_myeloid", "Myeloid cell", "CAF", "Fibroblast"),
}

MIN_SOURCE_CELLS = 20
MIN_TARGET_CELLS = 20
N_NEIGHBORS = 5
N_PERMUTATIONS = 250
PERMUTATION_BATCH_SIZE = 512


def load_scores(include_hypoxia: bool = False) -> tuple[pd.DataFrame, pd.DataFrame]:
    atlas = ad.read_h5ad(INPUT, backed="r")
    signatures = SIGNATURES | (HYPOXIA_SIGNATURE if include_hypoxia else {})
    genes = list(atlas.var_names.astype(str))
    available = set(genes)
    coverage_rows = []
    used_genes = []
    for signature, markers in signatures.items():
        present = [gene for gene in markers if gene in available]
        coverage_rows.append(
            {
                "signature": signature,
                "n_requested": len(markers),
                "n_available": len(present),
                "coverage_fraction": len(present) / len(markers),
                "available_genes": ";".join(present),
                "missing_genes": ";".join(gene for gene in markers if gene not in available),
            }
        )
        used_genes.extend(present)
    used_genes = sorted(set(used_genes))
    if any(row["coverage_fraction"] < 0.6 for row in coverage_rows):
        atlas.file.close()
        raise ValueError("At least one frozen signature has less than 60% panel coverage.")

    obs = atlas.obs[["sample", "cell_type"]].copy().reset_index(drop=True)
    coords = np.asarray(atlas.obsm["spatial"], dtype=float)
    matrix = atlas[:, used_genes].X
    matrix = matrix.to_memory() if hasattr(matrix, "to_memory") else matrix
    matrix = matrix.copy() if sparse.issparse(matrix) else np.asarray(matrix).copy()
    atlas.file.close()

    gene_index = {gene: index for index, gene in enumerate(used_genes)}
    transformed = matrix.copy()
    if sparse.issparse(transformed):
        transformed.data = np.log1p(transformed.data)
    else:
        transformed = np.log1p(transformed)
    for signature, markers in signatures.items():
        indices = [gene_index[gene] for gene in markers if gene in gene_index]
        values = transformed[:, indices].mean(axis=1)
        obs[signature] = np.asarray(values).ravel()
    obs["x"] = coords[:, 0]
    obs["y"] = coords[:, 1]
    obs["sample"] = pd.to_numeric(obs["sample"], errors="raise").astype(int)
    clinical = pd.read_excel(CLINICAL)
    clinical["sample"] = pd.to_numeric(clinical["sample"], errors="raise").astype(int)
    obs = obs.merge(clinical, on="sample", how="left", validate="many_to_one")
    if obs["Stage"].isna().any():
        raise ValueError("Clinical metadata did not map to every spatial cell.")
    return obs, pd.DataFrame(coverage_rows)


def batched_spearman_permutation_pvalue(
    source_values: np.ndarray,
    target_values: np.ndarray,
    rng: np.random.Generator,
    n_permutations: int,
    batch_size: int = PERMUTATION_BATCH_SIZE,
) -> tuple[float, float]:
    """Compute an exact Monte Carlo Spearman permutation p-value in batches.

    Permuting source values and then ranking them is equivalent to permuting
    their precomputed ranks. This avoids repeated ``scipy.stats.spearmanr``
    calls while preserving the original two-sided permutation null.
    """
    source_ranks = stats.rankdata(source_values).astype(np.float64, copy=False)
    target_ranks = stats.rankdata(target_values).astype(np.float64, copy=False)
    source_ranks -= source_ranks.mean()
    target_ranks -= target_ranks.mean()
    source_norm = np.linalg.norm(source_ranks)
    target_norm = np.linalg.norm(target_ranks)
    if source_norm == 0 or target_norm == 0:
        return float("nan"), float("nan")
    source_ranks /= source_norm
    target_ranks /= target_norm
    rho = float(source_ranks @ target_ranks)
    extreme = 0
    for start in range(0, n_permutations, batch_size):
        current = min(batch_size, n_permutations - start)
        # Generator.permuted shuffles every row independently along axis 1.
        permuted = rng.permuted(np.broadcast_to(source_ranks, (current, len(source_ranks))), axis=1)
        extreme += int(np.count_nonzero(np.abs(permuted @ target_ranks) >= abs(rho)))
    return rho, float((1 + extreme) / (n_permutations + 1))


def section_spatial_effects(
    obs: pd.DataFrame,
    seed: int = 481,
    relationships: dict | None = None,
    n_permutations: int = N_PERMUTATIONS,
    permutation_batch_size: int = PERMUTATION_BATCH_SIZE,
) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    relationships = relationships or RELATIONSHIPS
    for section, section_data in obs.groupby("sample", sort=True):
        coordinates = section_data[["x", "y"]].to_numpy(float)
        labels = section_data["cell_type"].astype(str).to_numpy()
        for relationship, (source_score, source_type, target_score, target_type) in relationships.items():
            source_positions = np.flatnonzero(cell_type_mask(labels, source_type))
            target_positions = np.flatnonzero(cell_type_mask(labels, target_type))
            if len(source_positions) < MIN_SOURCE_CELLS or len(target_positions) < MIN_TARGET_CELLS:
                continue
            source_values = section_data.iloc[source_positions][source_score].to_numpy(float)
            target_values = section_data.iloc[target_positions][target_score].to_numpy(float)
            if np.unique(source_values).size < 3 or np.unique(target_values).size < 3:
                continue
            tree = cKDTree(coordinates[target_positions])
            k = min(N_NEIGHBORS, len(target_positions))
            distances, neighbors = tree.query(coordinates[source_positions], k=k)
            if k == 1:
                neighbors = neighbors[:, None]
                distances = distances[:, None]
            local_target = target_values[neighbors].mean(axis=1)
            rho, permutation_p = batched_spearman_permutation_pvalue(
                source_values,
                local_target,
                rng,
                n_permutations,
                batch_size=permutation_batch_size,
            )
            _, p_value = stats.spearmanr(source_values, local_target)
            rows.append(
                {
                    "dataset_id": "GSE308624",
                    "section": int(section),
                    "relationship": relationship,
                    "n_source_cells": len(source_positions),
                    "n_target_cells": len(target_positions),
                    "spearman_rho": float(rho),
                    "analytic_p_value": float(p_value),
                    "within_section_permutation_p_value": permutation_p,
                    "mean_neighbor_distance_mm": float(np.mean(distances)),
                    "stage": section_data["Stage"].iloc[0],
                    "grade": section_data["Grade"].iloc[0],
                    "age": section_data["Age"].iloc[0],
                    "sex": section_data["Sex"].iloc[0],
                }
            )
    effects = pd.DataFrame(rows)
    effects["permutation_fdr_bh"] = multipletests(
        effects["within_section_permutation_p_value"], method="fdr_bh"
    )[1]
    return effects


def cell_type_mask(labels: np.ndarray, requested: str) -> np.ndarray:
    """Match documented coarse-label spelling variants without crossing lineages."""
    normalized = np.char.lower(np.char.replace(np.char.replace(labels.astype(str), "_", ""), " ", ""))
    token = requested.lower().replace("_", "").replace(" ", "")
    aliases = {
        "fibroblast": ("fibroblast",),
        "endothelial": ("endothelial",),
        "tcell": ("tcell", "tlymphocyte"),
        "myeloidcell": ("myeloid", "macrophage", "monocyte"),
    }
    accepted = aliases.get(token, (token,))
    return np.asarray([any(value == alias or value.startswith(alias) for alias in accepted) for value in normalized])


def random_effects_fisher(effects: pd.DataFrame) -> dict:
    valid = effects[(effects["n_source_cells"] >= 4) & effects["spearman_rho"].notna()].copy()
    z = np.arctanh(np.clip(valid["spearman_rho"].to_numpy(float), -0.999999, 0.999999))
    variance = 1.0 / (valid["n_source_cells"].to_numpy(float) - 3.0)
    fixed_weights = 1.0 / variance
    fixed_mean = float(np.sum(fixed_weights * z) / np.sum(fixed_weights))
    q = float(np.sum(fixed_weights * (z - fixed_mean) ** 2))
    df = len(z) - 1
    c_term = float(np.sum(fixed_weights) - np.sum(fixed_weights**2) / np.sum(fixed_weights))
    tau2 = max(0.0, (q - df) / c_term) if df > 0 and c_term > 0 else 0.0
    weights = 1.0 / (variance + tau2)
    pooled = float(np.sum(weights * z) / np.sum(weights))
    standard_error = math.sqrt(1.0 / np.sum(weights))
    positive_count = int((valid["spearman_rho"] > 0).sum())
    return {
        "n_sections": len(valid),
        "n_source_cells": int(valid["n_source_cells"].sum()),
        "pooled_rho": float(np.tanh(pooled)),
        "ci_low": float(np.tanh(pooled - 1.96 * standard_error)),
        "ci_high": float(np.tanh(pooled + 1.96 * standard_error)),
        "p_value": float(2 * stats.norm.sf(abs(pooled / standard_error))),
        "tau2_fisher_z": tau2,
        "i2_percent": max(0.0, (q - df) / q) * 100 if q > 0 else 0.0,
        "q_p_value": float(stats.chi2.sf(q, df)) if df > 0 else np.nan,
        "median_section_rho": float(valid["spearman_rho"].median()),
        "positive_section_fraction": positive_count / len(valid),
        "section_wilcoxon_p_value": float(stats.wilcoxon(valid["spearman_rho"]).pvalue),
        "section_sign_test_p_value": float(
            stats.binomtest(positive_count, len(valid), 0.5).pvalue
        ),
        "sections_permutation_fdr_lt_0_05": int((valid["permutation_fdr_bh"] < 0.05).sum()),
    }


def pool_effects(effects: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for relationship, group in effects.groupby("relationship", sort=False):
        rows.append({"relationship": relationship, **random_effects_fisher(group)})
    pooled = pd.DataFrame(rows)
    pooled["meta_fdr_bh"] = multipletests(pooled["p_value"], method="fdr_bh")[1]
    pooled["section_wilcoxon_fdr_bh"] = multipletests(
        pooled["section_wilcoxon_p_value"], method="fdr_bh"
    )[1]
    pooled["section_sign_test_fdr_bh"] = multipletests(
        pooled["section_sign_test_p_value"], method="fdr_bh"
    )[1]
    return pooled


def plot_results(effects: pd.DataFrame, pooled: pd.DataFrame, path: Path) -> None:
    relationships = list(pooled["relationship"])
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
    for position, relationship in enumerate(relationships):
        values = effects.loc[effects["relationship"].eq(relationship), "spearman_rho"]
        jitter = np.linspace(-0.15, 0.15, len(values)) if len(values) else []
        axes[0].scatter(np.repeat(position, len(values)) + jitter, values, s=16, alpha=0.55)
    axes[0].axhline(0, color="#6B7280", linestyle="--", linewidth=1)
    axes[0].set_xticks(range(len(relationships)), relationships, rotation=35, ha="right")
    axes[0].set_ylabel("Within-section spatial Spearman rho")
    axes[0].set_title("Section-level spatial effects")

    y = np.arange(len(pooled))
    effect = pooled["pooled_rho"].to_numpy()
    axes[1].errorbar(
        effect,
        y,
        xerr=[effect - pooled["ci_low"].to_numpy(), pooled["ci_high"].to_numpy() - effect],
        fmt="o",
        color="#176B87",
        capsize=3,
    )
    axes[1].axvline(0, color="#6B7280", linestyle="--", linewidth=1)
    axes[1].set_yticks(y, relationships)
    axes[1].set_xlabel("Random-effects pooled rho")
    axes[1].set_title("Across independent tumour sections")
    axes[1].grid(axis="x", alpha=0.2)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", type=Path, help="CosMx H5AD; defaults to the repository data path")
    parser.add_argument("--clinical", type=Path, help="Clinical workbook; defaults to the repository data path")
    parser.add_argument("--output-dir", type=Path, help="Output directory; defaults to the repository output path")
    parser.add_argument("--include-hypoxia", action="store_true", help="run the preregistered hypoxia-to-CAF supporting association")
    parser.add_argument("--n-permutations", type=int, default=N_PERMUTATIONS, help="within-section block-permutation count; use 10000 for final figures")
    parser.add_argument("--permutation-batch-size", type=int, default=PERMUTATION_BATCH_SIZE, help="permutation rows evaluated at once")
    args = parser.parse_args()
    global INPUT, CLINICAL, OUT_DIR
    if args.input_h5ad:
        INPUT = args.input_h5ad
    if args.clinical:
        CLINICAL = args.clinical
    if args.output_dir:
        OUT_DIR = args.output_dir
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    obs, coverage = load_scores(include_hypoxia=args.include_hypoxia)
    relationships = RELATIONSHIPS | (HYPOXIA_RELATIONSHIPS if args.include_hypoxia else {})
    if args.n_permutations < 250:
        raise ValueError("At least 250 permutations are required for a spatial result.")
    effects = section_spatial_effects(
        obs,
        relationships=relationships,
        n_permutations=args.n_permutations,
        permutation_batch_size=args.permutation_batch_size,
    )
    pooled = pool_effects(effects)
    coverage.to_csv(OUT_DIR / "GSE308624_SIGNATURE_COVERAGE.csv", index=False)
    effects.to_csv(OUT_DIR / "GSE308624_SECTION_SPATIAL_EFFECTS.csv", index=False)
    pooled.to_csv(OUT_DIR / "SPATIAL_VALIDATION_RESULTS.csv", index=False)
    plot_results(effects, pooled, OUT_DIR / "Figure_GSE308624_Spatial_Validation.png")
    supported = pooled[pooled["section_wilcoxon_fdr_bh"] < 0.05]
    text = f"""# GSE308624 Frozen Spatial Validation

Frozen project signatures were scored in CosMx cells from independent gastric tumour sections.
For each tumour section, source-program scores were related to mean target-program expression among
the five nearest lineage-matched cells. Labels were permuted within each section, and section-specific
correlations were combined with random effects.

- Prespecified spatial relationships tested: {len(pooled)}.
- Relationships with section-level Wilcoxon FDR <0.05: {len(supported)}.
- Signature coverage range: {coverage.coverage_fraction.min():.2f}-{coverage.coverage_fraction.max():.2f}.

The assay has a targeted panel. Spatial association supports local ecology, not direct ligand-receptor
signalling or causality. Each tumour section contributes one effect; this is not a tumour-normal comparison.
Published fine labels were not used; only the supplied coarse lineage labels and frozen gene programs
entered the primary analysis.
"""
    (OUT_DIR / "SPATIAL_VALIDATION_SUMMARY.md").write_text(text, encoding="utf-8")
    print(coverage.to_string(index=False))
    print(pooled.to_string(index=False))


if __name__ == "__main__":
    main()
