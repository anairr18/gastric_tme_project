from pathlib import Path

import pandas as pd

from spatial_gse308624_validation import RELATIONSHIPS, load_scores, section_spatial_effects, pool_effects


BASE = Path(__file__).resolve().parents[1]


def test_spatial_signature_coverage_and_section_mapping():
    obs, coverage = load_scores()
    assert obs["sample"].nunique() == 71
    assert len(obs) == 200182
    assert (coverage["coverage_fraction"] >= 0.6).all()
    assert obs["Stage"].notna().all()


def test_spatial_outputs_are_section_level_and_prespecified():
    output = BASE / "outputs" / "SPATIAL_VALIDATION"
    effects = pd.read_csv(output / "GSE308624_SECTION_SPATIAL_EFFECTS.csv")
    pooled = pd.read_csv(output / "SPATIAL_VALIDATION_RESULTS.csv")
    assert not effects.duplicated(["section", "relationship"]).any()
    assert set(pooled["relationship"]) == set(RELATIONSHIPS)
    assert pooled["n_sections"].between(20, 71).all()
    assert pooled["section_wilcoxon_fdr_bh"].between(0, 1).all()
    assert effects["within_section_permutation_p_value"].between(0, 1).all()
