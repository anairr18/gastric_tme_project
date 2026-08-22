"""Condition-aware, donor/sample-level state-abundance meta-analysis.

Input is the curated compartment-discovery output.  Sample conventions are
explicitly encoded from the audited source cohorts; unresolved samples are not
silently admitted to tumour-normal meta-analysis.
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
from scipy import stats


KUMAR_NORMAL = {f"sample{i}" for i in [1, 4, 6, 9, 11, 21, 23, 25, 31, 35, 37]}
PAIRED_COHORTS = {"korea_kim2022", "sathe2020", "diffuse_gc_2021", "tcell_exhaustion_2022"}


def bh_fdr(values: pd.Series) -> pd.Series:
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.dropna()
    if valid.empty:
        return result
    order = valid.sort_values().index
    raw = valid.loc[order].to_numpy() * len(order) / np.arange(1, len(order) + 1)
    result.loc[order] = np.minimum(np.minimum.accumulate(raw[::-1])[::-1], 1.0)
    return result


def sample_metadata(dataset_id: str, sample: str) -> dict[str, str | None]:
    """Metadata rules documented in the project audit and source manifests."""
    if dataset_id == "korea_kim2022":
        patient = re.match(r"^(E\d+)", sample)
        if "_AN" in sample:
            return {"patient": patient.group(1) if patient else None, "tissue": "normal", "condition": "adjacent_normal", "mapping_rule": "korea_sample_suffix_AN"}
        if "_DN" in sample:
            return {"patient": patient.group(1) if patient else None, "tissue": "normal", "condition": "distal_normal", "mapping_rule": "korea_sample_suffix_DN"}
        if sample.endswith("_B"):
            return {"patient": patient.group(1) if patient else None, "tissue": "tumor", "condition": "tumor_BL", "mapping_rule": "korea_sample_suffix_B"}
        if "_F1" in sample:
            return {"patient": patient.group(1) if patient else None, "tissue": "tumor", "condition": "tumor_FU1", "mapping_rule": "korea_sample_suffix_F1"}
        if "_F2" in sample:
            return {"patient": patient.group(1) if patient else None, "tissue": "tumor", "condition": "tumor_FU2", "mapping_rule": "korea_sample_suffix_F2"}
    if dataset_id == "kumar2022":
        return {"patient": None, "tissue": "normal" if sample in KUMAR_NORMAL else "tumor", "condition": "tumor_normal", "mapping_rule": "kumar_audited_sample_set"}
    if dataset_id == "sathe2020":
        match = re.match(r"^Pat(\d+)-([AB])$", sample)
        # The audited project mapping documents the first 24 patients as paired.
        if match and int(match.group(1)) <= 24:
            return {"patient": f"Pat{int(match.group(1)):02d}", "tissue": "normal" if match.group(2) == "A" else "tumor", "condition": "tumor_normal", "mapping_rule": "cho_sample_pattern_PatAB"}
    if dataset_id == "diffuse_gc_2021":
        match = re.match(r"^(Pt\d+)_(Normal|Superficial|Deep)_CountMatrix$", sample)
        if match:
            return {"patient": match.group(1), "tissue": "normal" if match.group(2) == "Normal" else "tumor", "condition": match.group(2).lower(), "mapping_rule": "jeong_sample_pattern_layer"}
    if dataset_id == "tcell_exhaustion_2022" and re.match(r"^\d+[NT]$", sample):
        return {"patient": sample[:-1], "tissue": "normal" if sample.endswith("N") else "tumor", "condition": "tumor_normal", "mapping_rule": "kang_sample_suffix_NT"}
    if dataset_id == "zhang2021":
        for stage in ["NAG", "CAG", "IMS", "IMW", "EGC"]:
            if stage in sample:
                canonical = "IM" if stage in {"IMS", "IMW"} else stage
                tissue = "normal" if canonical == "NAG" else ("tumor" if canonical == "EGC" else "premalignant")
                return {"patient": None, "tissue": tissue, "condition": canonical, "mapping_rule": "zhang_stage_in_sample_id"}
    return {"patient": None, "tissue": None, "condition": None, "mapping_rule": "unresolved"}


def attach_metadata(table: pd.DataFrame) -> pd.DataFrame:
    metadata = [sample_metadata(dataset, sample) for dataset, sample in zip(table["dataset_id"], table["sample"])]
    result = table.copy()
    result[["source_patient", "source_tissue", "source_condition", "mapping_rule"]] = pd.DataFrame(metadata, index=result.index)
    return result


def cohort_effect(frame: pd.DataFrame, cohort: str, state: str) -> dict | None:
    subset = frame.loc[(frame["dataset_id"] == cohort) & (frame["state_id"] == state)]
    if cohort == "korea_kim2022":
        subset = subset.loc[subset["source_condition"].isin(["tumor_BL", "adjacent_normal"])]
    else:
        subset = subset.loc[subset["source_tissue"].isin(["tumor", "normal"])]
    if subset.empty:
        return None
    paired = cohort in PAIRED_COHORTS and subset["source_patient"].notna().any()
    if paired:
        pivot = subset.pivot_table(index="source_patient", columns="source_tissue", values="state_fraction", aggfunc="mean", observed=True)
        if "tumor" not in pivot.columns or "normal" not in pivot.columns:
            return None
        pivot = pivot.dropna(subset=["tumor", "normal"])
        if len(pivot) < 3:
            return None
        delta = pivot["tumor"] - pivot["normal"]
        effect = float(delta.mean())
        se = float(delta.std(ddof=1) / math.sqrt(len(delta)))
        return {"cohort": cohort, "state_id": state, "n_tumor": len(delta), "n_normal": len(delta), "n_pairs": len(delta),
                "design": "paired", "effect": effect, "standard_error": se}
    tumor = subset.loc[subset["source_tissue"] == "tumor", "state_fraction"]
    normal = subset.loc[subset["source_tissue"] == "normal", "state_fraction"]
    if len(tumor) < 3 or len(normal) < 3:
        return None
    effect = float(tumor.mean() - normal.mean())
    se = float(math.sqrt(tumor.var(ddof=1) / len(tumor) + normal.var(ddof=1) / len(normal)))
    return {"cohort": cohort, "state_id": state, "n_tumor": len(tumor), "n_normal": len(normal), "n_pairs": 0,
            "design": "unpaired", "effect": effect, "standard_error": se}


def random_effects(effects: pd.DataFrame) -> dict | None:
    values = effects.loc[(effects["standard_error"] > 0) & np.isfinite(effects["standard_error"])].copy()
    if len(values) < 3:
        return None
    variances = values["standard_error"].to_numpy() ** 2
    observed = values["effect"].to_numpy()
    fixed_weights = 1 / variances
    fixed_effect = np.sum(fixed_weights * observed) / np.sum(fixed_weights)
    q = float(np.sum(fixed_weights * (observed - fixed_effect) ** 2))
    df = len(values) - 1
    c = float(np.sum(fixed_weights) - np.sum(fixed_weights ** 2) / np.sum(fixed_weights))
    tau2 = max(0.0, (q - df) / c) if c > 0 else 0.0
    weights = 1 / (variances + tau2)
    pooled = float(np.sum(weights * observed) / np.sum(weights))
    se = float(math.sqrt(1 / np.sum(weights)))
    z = pooled / se if se else np.nan
    p_value = float(2 * stats.norm.sf(abs(z))) if np.isfinite(z) else np.nan
    ci_low, ci_high = pooled - 1.96 * se, pooled + 1.96 * se
    i2 = max(0.0, (q - df) / q) * 100 if q > 0 else 0.0
    pi_sd = math.sqrt(tau2 + se ** 2)
    return {"n_cohorts": len(values), "pooled_effect": pooled, "standard_error": se, "ci_low": ci_low, "ci_high": ci_high,
            "p_value": p_value, "tau2": tau2, "i2_percent": i2, "prediction_low": pooled - 1.96 * pi_sd,
            "prediction_high": pooled + 1.96 * pi_sd, "directional_concordance": float(max((observed > 0).mean(), (observed < 0).mean()))}


def loco_direction(effects: pd.DataFrame) -> bool:
    full = random_effects(effects)
    if full is None or full["pooled_effect"] == 0:
        return False
    directions = []
    for cohort in effects["cohort"].unique():
        result = random_effects(effects.loc[effects["cohort"] != cohort])
        if result is not None:
            directions.append(np.sign(result["pooled_effect"]) == np.sign(full["pooled_effect"]))
    return bool(directions and all(directions))


def zhang_trends(frame: pd.DataFrame) -> pd.DataFrame:
    stage_order = {"NAG": 0, "CAG": 1, "IM": 2, "EGC": 3}
    data = frame.loc[(frame["dataset_id"] == "zhang2021") & frame["source_condition"].isin(stage_order)].copy()
    rows = []
    for state, group in data.groupby("state_id", observed=True):
        sample_table = group.groupby(["sample", "source_condition"], observed=True)["state_fraction"].mean().reset_index()
        if sample_table["source_condition"].nunique() < 3:
            continue
        x = sample_table["source_condition"].map(stage_order).to_numpy()
        rho, p_value = stats.spearmanr(x, sample_table["state_fraction"].to_numpy())
        rows.append({"state_id": state, "n_samples": len(sample_table), "rho_stage": rho, "p_value": p_value,
                     "mean_NAG": sample_table.loc[sample_table["source_condition"] == "NAG", "state_fraction"].mean(),
                     "mean_CAG": sample_table.loc[sample_table["source_condition"] == "CAG", "state_fraction"].mean(),
                     "mean_IM": sample_table.loc[sample_table["source_condition"] == "IM", "state_fraction"].mean(),
                     "mean_EGC": sample_table.loc[sample_table["source_condition"] == "EGC", "state_fraction"].mean()})
    result = pd.DataFrame(rows)
    if not result.empty:
        result["fdr"] = bh_fdr(result["p_value"])
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--curated-composition", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--marker-summary", type=Path, help="Optional leave-one-cohort-out marker replication summary.")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    composition = pd.read_csv(args.curated_composition)
    if "zero_filled_absent_state" not in composition.columns:
        raise ValueError(
            "Curated composition is not zero-complete. Re-run "
            "build_primary_working_state_composition.py before meta-analysis; "
            "otherwise samples with an absent state are dropped from the contrast."
        )
    if composition["zero_filled_absent_state"].isna().any():
        raise ValueError("zero_filled_absent_state contains missing values; stop before inference.")
    # The stage-test workflow writes the reviewed label as ``reviewed_state``;
    # the standalone finalizer uses ``reviewed_state_label``. Accept either
    # explicit, semantically identical column without guessing a state label.
    if "reviewed_state_label" not in composition.columns:
        if "reviewed_state" not in composition.columns:
            raise ValueError(
                "Curated composition needs reviewed_state_label or reviewed_state. "
                f"Available columns: {sorted(composition.columns)}"
            )
        composition = composition.rename(columns={"reviewed_state": "reviewed_state_label"})
    composition = attach_metadata(composition)
    if "state_test_tier" in composition.columns:
        composition = composition.loc[composition["state_test_tier"].eq("candidate_confirmatory")].copy()
    composition["state_id"] = composition["compartment"].astype(str) + "::" + composition["reviewed_state_label"].astype(str)
    composition.to_csv(args.output_dir / "STATE_ABUNDANCE_METADATA_HARMONIZED.csv", index=False)
    mapping = composition[["dataset_id", "sample", "source_patient", "source_tissue", "source_condition", "mapping_rule"]].drop_duplicates()
    mapping["eligible_for_condition_analysis"] = mapping["source_condition"].notna() & mapping["source_tissue"].notna()
    mapping.to_csv(args.output_dir / "SAMPLE_CONDITION_MAPPING_AUDIT.csv", index=False)
    mapping.loc[~mapping["eligible_for_condition_analysis"]].to_csv(
        args.output_dir / "SAMPLE_CONDITION_MAPPING_EXCLUSIONS.csv", index=False
    )

    effect_rows = []
    for state in sorted(composition["state_id"].unique()):
        for cohort in sorted(composition["dataset_id"].unique()):
            result = cohort_effect(composition, cohort, state)
            if result is not None:
                effect_rows.append(result)
    effects = pd.DataFrame(effect_rows, columns=[
        "cohort", "state_id", "n_tumor", "n_normal", "n_pairs", "design", "effect", "standard_error",
    ])
    effects.to_csv(args.output_dir / "TUMOR_NORMAL_STATE_COHORT_EFFECTS.csv", index=False)
    meta_rows = []
    for state, group in effects.groupby("state_id", observed=True):
        result = random_effects(group)
        if result is not None:
            meta_rows.append({"state_id": state, **result, "loco_stable_direction": loco_direction(group),
                              "cohorts": ";".join(sorted(group["cohort"].unique()))})
    meta = pd.DataFrame(meta_rows, columns=[
        "state_id", "n_cohorts", "pooled_effect", "standard_error", "ci_low", "ci_high", "p_value",
        "tau2", "i2_percent", "prediction_low", "prediction_high", "directional_concordance",
        "loco_stable_direction", "cohorts",
    ])
    if not meta.empty:
        if args.marker_summary:
            marker_summary = pd.read_csv(args.marker_summary)
            required_marker_columns = {"state_id", "state_marker_gate"}
            missing_marker_columns = required_marker_columns - set(marker_summary.columns)
            if missing_marker_columns:
                raise ValueError(f"Marker summary is missing: {sorted(missing_marker_columns)}")
            meta = meta.merge(
                marker_summary[["state_id", "state_marker_gate"]], on="state_id", how="left", validate="one_to_one"
            )
        else:
            meta["state_marker_gate"] = "not_run"
        meta["fdr"] = bh_fdr(meta["p_value"])
        meta["publication_class"] = np.where(
            (meta["fdr"] < 0.05) & (meta["directional_concordance"] >= 0.8) & meta["loco_stable_direction"]
            & (meta["i2_percent"] < 50) & ((meta["prediction_low"] > 0) | (meta["prediction_high"] < 0))
            & meta["state_marker_gate"].eq("pass"),
            "portable_strict", "exploratory_or_heterogeneous"
        )
    meta = meta.reindex(columns=[
        "state_id", "n_cohorts", "pooled_effect", "standard_error", "ci_low", "ci_high", "p_value",
        "tau2", "i2_percent", "prediction_low", "prediction_high", "directional_concordance",
        "loco_stable_direction", "cohorts", "state_marker_gate", "fdr", "publication_class",
    ])
    meta.to_csv(args.output_dir / "TUMOR_NORMAL_STATE_RANDOM_EFFECTS_META.csv", index=False)
    trends = zhang_trends(composition)
    trends.to_csv(args.output_dir / "ZHANG_PREMALIGNANT_STATE_TRENDS.csv", index=False)

    if not meta.empty:
        plot = meta.sort_values("fdr").head(15).iloc[::-1]
        fig, ax = plt.subplots(figsize=(10, max(5, 0.38 * len(plot))))
        ax.errorbar(plot["pooled_effect"], range(len(plot)), xerr=[plot["pooled_effect"] - plot["ci_low"], plot["ci_high"] - plot["pooled_effect"]], fmt="o", color="#1f4e79")
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_yticks(range(len(plot)), plot["state_id"])
        ax.set_xlabel("Tumour - normal state fraction (random-effects pooled difference)")
        ax.set_title("Curated state-abundance meta-analysis")
        fig.tight_layout()
        fig.savefig(args.output_dir / "STATE_ABUNDANCE_META_FOREST.png", dpi=300, bbox_inches="tight")
        fig.savefig(args.output_dir / "STATE_ABUNDANCE_META_FOREST.pdf", bbox_inches="tight")
        plt.close(fig)
    if not trends.empty:
        top = trends.sort_values("fdr").head(2)
        stages = ["NAG", "CAG", "IM", "EGC"]
        fig, ax = plt.subplots(figsize=(7.5, 4.5))
        for _, row in top.iterrows():
            values = [row.get(f"mean_{stage}", np.nan) for stage in stages]
            ax.plot(stages, values, marker="o", linewidth=2, label=row["state_id"])
        ax.set_ylabel("Mean curated state fraction per biopsy")
        ax.set_title("Zhang gastric disease continuum: top exploratory state trends")
        ax.legend(frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(args.output_dir / "ZHANG_STAGE_STATE_TRAJECTORY.png", dpi=300, bbox_inches="tight")
        fig.savefig(args.output_dir / "ZHANG_STAGE_STATE_TRAJECTORY.pdf", bbox_inches="tight")
        plt.close(fig)
    print(f"State-abundance meta-analysis written to {args.output_dir}")


if __name__ == "__main__":
    main()
