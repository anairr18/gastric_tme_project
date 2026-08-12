#!/usr/bin/env python3
"""Create an auditable sample-qualified cell-ID repair for Korean raw counts.

The original H5AD is never modified. The repaired copy preserves its original
barcode in ``obs['original_obs_name']`` and writes a row-level mapping table.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from datetime import datetime, timezone

import pandas as pd


SAMPLE_COLUMNS = ["sample_id", "sample", "orig.ident", "library_id", "analysis_unit"]


def find_column(frame: pd.DataFrame, candidates: list[str]) -> str | None:
    lookup = {str(name).lower(): str(name) for name in frame.columns}
    return next((lookup[name.lower()] for name in candidates if name.lower() in lookup), None)


def repair_mapping(obs: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, object]]:
    sample_col = find_column(obs, SAMPLE_COLUMNS)
    if sample_col is None:
        raise ValueError("Korean raw metadata lacks a sample identifier; cannot create sample-qualified cell IDs.")
    samples = obs[sample_col].astype("string")
    if samples.isna().any() or samples.str.strip().eq("").any():
        raise ValueError("Korean raw metadata has missing/blank sample identifiers; repair stopped.")
    original = pd.Index(obs.index.astype(str), name="original_obs_name")
    base = samples.astype(str) + "::" + original.astype(str)
    rank = base.groupby(base, sort=False).cumcount()
    repaired = base.where(rank.eq(0), base + "::duplicate_" + rank.astype(str))
    mapping = pd.DataFrame(
        {
            "source_row": range(len(obs)),
            "sample": samples.astype(str).to_numpy(),
            "original_obs_name": original.to_numpy(),
            "sample_qualified_base": base.to_numpy(),
            "within_base_duplicate_rank": rank.to_numpy(),
            "repaired_obs_name": repaired.to_numpy(),
        }
    )
    if mapping["repaired_obs_name"].duplicated().any():
        raise RuntimeError("Internal repair failure: repaired IDs are still duplicated.")
    report = {
        "n_cells": int(len(mapping)),
        "sample_column": sample_col,
        "original_obs_names_unique": bool(original.is_unique),
        "sample_qualified_before_rank_unique": bool(pd.Index(base).is_unique),
        "n_duplicate_original_obs_name_rows": int(mapping["original_obs_name"].duplicated(keep=False).sum()),
        "n_duplicate_sample_qualified_rows": int(mapping["sample_qualified_base"].duplicated(keep=False).sum()),
        "n_rows_with_deterministic_rank_suffix": int(rank.gt(0).sum()),
        "repaired_obs_names_unique": bool(pd.Index(repaired).is_unique),
    }
    return mapping, report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--output-h5ad", type=Path, help="Write a new H5AD with repaired obs names; source file is retained unchanged.")
    parser.add_argument("--skip-write", action="store_true", help="Only create the mapping/audit, even if --output-h5ad was supplied.")
    args = parser.parse_args()
    if not args.input_h5ad.exists():
        raise FileNotFoundError(f"Korean raw H5AD not found: {args.input_h5ad}")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    import anndata as ad

    # First pass reads only obs metadata, so the audit is saved before a large rewrite.
    backed = ad.read_h5ad(args.input_h5ad, backed="r")
    try:
        mapping, report = repair_mapping(backed.obs.copy())
    finally:
        backed.file.close()
    mapping.to_csv(args.output_dir / "KOREA_RAW_CELL_ID_REPAIR_MAPPING.csv.gz", index=False)
    (args.output_dir / "KOREA_RAW_CELL_ID_REPAIR_AUDIT.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    pd.DataFrame([report]).to_csv(args.output_dir / "KOREA_RAW_CELL_ID_REPAIR_AUDIT.csv", index=False)
    if args.skip_write or args.output_h5ad is None:
        print(json.dumps(report, indent=2))
        print("Audit complete; no H5AD was written.")
        return

    if args.output_h5ad.exists():
        existing = ad.read_h5ad(args.output_h5ad, backed="r")
        try:
            valid_existing = existing.obs_names.is_unique and len(existing.obs) == len(mapping)
        finally:
            existing.file.close()
        if valid_existing:
            print(f"Reusing validated repaired H5AD: {args.output_h5ad}")
            return
        archived = args.output_h5ad.with_name(
            args.output_h5ad.stem + "_invalid_" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + args.output_h5ad.suffix
        )
        args.output_h5ad.rename(archived)
        print(f"Existing repaired H5AD failed uniqueness validation; preserved as {archived}")
    args.output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    print(f"Loading source counts to create repaired H5AD: {args.input_h5ad}", flush=True)
    data = ad.read_h5ad(args.input_h5ad)
    if len(data.obs) != len(mapping):
        raise RuntimeError("Source changed between audit and rewrite; stop rather than assigning IDs to different rows.")
    data.obs["original_obs_name"] = mapping["original_obs_name"].to_numpy()
    data.obs["sample_qualified_cell_id"] = mapping["repaired_obs_name"].to_numpy()
    data.obs_names = pd.Index(mapping["repaired_obs_name"].to_numpy())
    if not data.obs_names.is_unique:
        raise RuntimeError("Repaired AnnData obs names are not unique.")
    print(f"Writing repaired raw-count H5AD: {args.output_h5ad}", flush=True)
    data.write_h5ad(args.output_h5ad, compression="lzf")
    (args.output_dir / "KOREA_RAW_CELL_ID_REPAIR_MANIFEST.json").write_text(
        json.dumps({**report, "input_h5ad": str(args.input_h5ad), "output_h5ad": str(args.output_h5ad), "source_modified": False}, indent=2) + "\n",
        encoding="utf-8",
    )
    print("Korean raw-cell-ID repair complete.")


if __name__ == "__main__":
    main()
