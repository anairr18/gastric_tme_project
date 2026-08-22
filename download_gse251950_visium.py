#!/usr/bin/env python3
"""Download, unpack, and index GSE251950 Visium inputs without using them for atlas fitting."""

from __future__ import annotations

import argparse
import re
import shutil
import tarfile
import urllib.request
from pathlib import Path

import pandas as pd


GEO_TAR = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE251950&format=file&file=GSE251950_RAW.tar"
GSM_PATTERN = re.compile(r"GSM7990\d+")


def download(url: str, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() and target.stat().st_size > 100_000_000:
        print(f"Using existing archive: {target}")
        return
    print(f"Downloading {url}")
    urllib.request.urlretrieve(url, target)


def unpack(archive: Path, destination: Path) -> None:
    marker = destination / ".extracted"
    if not marker.exists():
        destination.mkdir(parents=True, exist_ok=True)
        with tarfile.open(archive, "r") as tar:
            for member in tar.getmembers():
                resolved = (destination / member.name).resolve()
                if not str(resolved).startswith(str(destination.resolve())):
                    raise ValueError(f"Unsafe archive member: {member.name}")
            tar.extractall(destination)
        marker.write_text("GSE251950 GEO supplementary archive extracted\n", encoding="utf-8")
    # GEO RAW archives often contain one compressed archive per Visium slide.
    # Expand those once, while leaving ordinary matrix .gz files untouched.
    for nested in sorted(destination.rglob("*")):
        if not nested.is_file() or nested == archive or not tarfile.is_tarfile(nested):
            continue
        nested_dir = nested.parent / f"{nested.stem}_unpacked"
        if nested_dir.exists():
            continue
        nested_dir.mkdir(parents=True, exist_ok=True)
        with tarfile.open(nested, "r:*") as tar:
            for member in tar.getmembers():
                resolved = (nested_dir / member.name).resolve()
                if not str(resolved).startswith(str(nested_dir.resolve())):
                    raise ValueError(f"Unsafe nested archive member: {member.name}")
            tar.extractall(nested_dir)


def slide_id(path: Path) -> str:
    match = GSM_PATTERN.search(str(path))
    if match:
        return match.group(0)
    for part in reversed(path.parts):
        if part.lower().startswith("gc"):
            return part
    return path.parent.name


def is_matrix_file(path: Path) -> bool:
    lower = path.name.lower()
    return lower.endswith(".h5") or "matrix.mtx" in lower


def is_coordinate_file(path: Path) -> bool:
    lower = path.name.lower()
    return (lower.endswith(".csv") or lower.endswith(".csv.gz")) and ("position" in lower or "spatial" in lower)


def index_inputs(root: Path) -> pd.DataFrame:
    raw = root / "raw"
    records: dict[str, dict[str, object]] = {}
    for file in raw.rglob("*"):
        if not file.is_file():
            continue
        identifier = slide_id(file)
        record = records.setdefault(identifier, {"slide_id": identifier, "source_dir": str(file.parent)})
        if is_matrix_file(file):
            current = Path(str(record["matrix_path"])) if record.get("matrix_path") else None
            # Prefer a standard filtered 10x H5 when both H5 and MTX are present.
            if current is None or (file.suffix.lower() == ".h5" and current.suffix.lower() != ".h5"):
                record["matrix_path"] = str(file)
        if is_coordinate_file(file):
            record["coordinates_path"] = str(file)
    rows = []
    for identifier, record in sorted(records.items()):
        matrix = Path(str(record.get("matrix_path", ""))) if record.get("matrix_path") else None
        coordinates = Path(str(record.get("coordinates_path", ""))) if record.get("coordinates_path") else None
        if matrix is None and coordinates is None:
            continue
        # Dimensions are populated by the validation runner after Scanpy reads the matrix.
        rows.append(
            {
                **record,
                "matrix_path": str(matrix) if matrix else "",
                "coordinates_path": str(coordinates) if coordinates else "",
                "status": "ready" if matrix and coordinates else "blocked_missing_matrix_or_coordinates",
                "n_genes": 0,
                "n_spots": 0,
            }
        )
    return pd.DataFrame(rows).drop_duplicates(["slide_id", "matrix_path", "coordinates_path"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--url", default=GEO_TAR)
    parser.add_argument("--skip-download", action="store_true")
    args = parser.parse_args()
    root = args.output_root / "GSE251950"
    archive = root / "GSE251950_RAW.tar"
    if not args.skip_download:
        download(args.url, archive)
    if not archive.exists():
        raise FileNotFoundError(f"GSE251950 archive is missing: {archive}")
    unpack(archive, root / "raw")
    index = index_inputs(root)
    if index.empty:
        raise ValueError("No GSE251950 spatial inputs were discovered after extraction.")
    index.to_csv(root / "SPATIAL_INPUT_INDEX.csv", index=False)
    print(index[["slide_id", "status", "matrix_path", "coordinates_path"]].to_string(index=False))


if __name__ == "__main__":
    main()
