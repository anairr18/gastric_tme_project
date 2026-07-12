"""
03_dataset_download.py
Download all GEO datasets listed in the manifest.

Usage:
    python 03_dataset_download.py                   # download all pending
    python 03_dataset_download.py GSE183904         # download one accession
    python 03_dataset_download.py --list            # show status only

Strategy per dataset:
  1. Use GEOparse to fetch the GEO SOFT record and locate supplementary files.
  2. Download all supplementary files to data/external/<dataset_id>/raw/.
  3. If a single .h5ad exists, use it directly; otherwise try to convert from
     10x MEX (barcodes.tsv.gz / features.tsv.gz / matrix.mtx.gz) or
     from a dense CSV/TSV count matrix.
  4. Write a standardized raw h5ad to data/external/<dataset_id>/<dataset_id>_raw.h5ad
     so that 04_standardized_qc.py has a uniform input regardless of source format.

Note: GEO download sizes are large (1–20 GB per dataset). Run on a fast connection.
"""
import sys, io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", line_buffering=True)
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8", line_buffering=True)

import os
import gzip
import shutil
import argparse
import warnings
import urllib.request
import pandas as pd
import numpy as np
import scipy.sparse as sp
import scanpy as sc
import GEOparse

warnings.filterwarnings("ignore")

BASE     = os.path.expanduser("~/gastric_tme_project")
EXT_DIR  = os.path.join(BASE, "data/external")
MANIFEST = os.path.join(EXT_DIR, "dataset_manifest.csv")


# ─── helpers ─────────────────────────────────────────────────────────────────

def dataset_dir(dataset_id):
    return os.path.join(EXT_DIR, dataset_id)


def raw_dir(dataset_id):
    d = os.path.join(dataset_dir(dataset_id), "raw")
    os.makedirs(d, exist_ok=True)
    return d


def out_h5ad(dataset_id):
    return os.path.join(dataset_dir(dataset_id), f"{dataset_id}_raw.h5ad")


def download_file(url, dest_path, label=""):
    if os.path.exists(dest_path):
        print(f"    already exists: {os.path.basename(dest_path)}")
        return
    print(f"    downloading {label or os.path.basename(dest_path)} ...")
    try:
        urllib.request.urlretrieve(url, dest_path)
        size_mb = os.path.getsize(dest_path) / 1e6
        print(f"    -> {size_mb:.1f} MB")
    except Exception as e:
        print(f"    FAILED: {e}")


def gunzip_file(gz_path, out_path):
    if os.path.exists(out_path):
        return
    with gzip.open(gz_path, "rb") as f_in, open(out_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)


# ─── format converters ────────────────────────────────────────────────────────

def load_10x_mex(mex_dir):
    """Load 10x MEX directory into AnnData."""
    return sc.read_10x_mtx(mex_dir, var_names="gene_symbols", cache=False)


def load_h5(h5_path):
    """Load 10x HDF5 (.h5) file into AnnData."""
    return sc.read_10x_h5(h5_path)


def load_h5ad(h5ad_path):
    return sc.read_h5ad(h5ad_path)


def detect_and_load(raw_dir_path, dataset_id):
    """
    Auto-detect the format of files in raw_dir_path and return an AnnData.
    Priority: .h5ad > .h5 (10x) > 10x MEX > dense matrix.
    """
    files = os.listdir(raw_dir_path)

    # .h5ad (already anndata)
    h5ads = [f for f in files if f.endswith(".h5ad")]
    if h5ads:
        print(f"  detected format: .h5ad")
        adatas = [load_h5ad(os.path.join(raw_dir_path, f)) for f in sorted(h5ads)]
        if len(adatas) == 1:
            return adatas[0]
        # multiple samples → concatenate
        for i, (f, a) in enumerate(zip(sorted(h5ads), adatas)):
            a.obs["sample_file"] = os.path.splitext(f)[0]
        combined = sc.concat(adatas, label="sample_file",
                             keys=[os.path.splitext(f)[0] for f in sorted(h5ads)],
                             merge="same")
        return combined

    # 10x HDF5
    h5s = [f for f in files if f.endswith(".h5") and not f.endswith(".h5ad")]
    if h5s:
        print(f"  detected format: 10x .h5")
        adatas = []
        for f in sorted(h5s):
            a = load_h5(os.path.join(raw_dir_path, f))
            a.obs["sample_file"] = os.path.splitext(f)[0]
            adatas.append(a)
        if len(adatas) == 1:
            return adatas[0]
        return sc.concat(adatas, label="sample_file",
                         keys=[os.path.splitext(f)[0] for f in sorted(h5s)],
                         merge="same")

    # 10x MEX (matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz)
    mtx_files = [f for f in files if f.endswith("matrix.mtx.gz") or f == "matrix.mtx"]
    if mtx_files:
        print(f"  detected format: 10x MEX")
        return load_10x_mex(raw_dir_path)

    # Subdirectories — each is a sample in 10x MEX format
    subdirs = [f for f in files if os.path.isdir(os.path.join(raw_dir_path, f))]
    mex_subdirs = [s for s in subdirs
                   if any(x in os.listdir(os.path.join(raw_dir_path, s))
                          for x in ["matrix.mtx.gz", "matrix.mtx", "barcodes.tsv.gz"])]
    if mex_subdirs:
        print(f"  detected format: 10x MEX (multi-sample subdirs, {len(mex_subdirs)} samples)")
        adatas = []
        for s in sorted(mex_subdirs):
            a = load_10x_mex(os.path.join(raw_dir_path, s))
            a.obs["sample"] = s
            a.obs_names = [f"{s}_{bc}" for bc in a.obs_names]
            adatas.append(a)
        return sc.concat(adatas, merge="same")

    raise ValueError(
        f"Could not detect count matrix format in {raw_dir_path}. "
        f"Files present: {files}"
    )


# ─── GEO download ─────────────────────────────────────────────────────────────

def download_geo(dataset_id, geo_accession):
    print(f"\n{'='*60}")
    print(f"Downloading {dataset_id} ({geo_accession})")
    print(f"{'='*60}")

    target_h5ad = out_h5ad(dataset_id)
    if os.path.exists(target_h5ad):
        print(f"  Raw h5ad already exists: {target_h5ad}")
        return True

    rdir = raw_dir(dataset_id)

    # Fetch GEO metadata
    print(f"  fetching GEO metadata for {geo_accession} ...")
    try:
        gse = GEOparse.get_GEO(geo=geo_accession, destdir=rdir, silent=True)
    except Exception as e:
        print(f"  GEOparse fetch failed: {e}")
        return False

    # Collect supplementary file URLs from GSE and all GSMs
    supp_urls = []

    # Series-level supplementary files
    if hasattr(gse, "metadata") and "supplementary_file" in gse.metadata:
        for url in gse.metadata["supplementary_file"]:
            if url and url.lower() != "none":
                supp_urls.append(url.strip())

    # Sample-level supplementary files
    for gsm_name, gsm in gse.gsms.items():
        if "supplementary_file_1" in gsm.metadata:
            for url in gsm.metadata["supplementary_file_1"]:
                if url and url.lower() != "none":
                    supp_urls.append(url.strip())

    supp_urls = list(dict.fromkeys(supp_urls))  # deduplicate, preserve order
    print(f"  found {len(supp_urls)} supplementary file(s)")

    if not supp_urls:
        print("  WARNING: no supplementary files found via GEOparse.")
        print("  Attempting direct FTP listing ...")
        # Fallback: try the standard GEO FTP path
        accession_prefix = geo_accession[:6] + "nnn"
        ftp_url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{accession_prefix}/{geo_accession}/suppl/"
        print(f"  FTP URL (manual check needed): {ftp_url}")
        print(f"  Please manually download to: {rdir}")
        return False

    # Download each supplementary file
    for url in supp_urls:
        fname = os.path.basename(url.split("?")[0])
        dest  = os.path.join(rdir, fname)
        download_file(url, dest, fname)

    # Decompress tar archives if present
    tar_files = [f for f in os.listdir(rdir) if f.endswith(".tar") or f.endswith(".tar.gz")]
    for tf in tar_files:
        print(f"  extracting {tf} ...")
        import tarfile
        with tarfile.open(os.path.join(rdir, tf)) as tar:
            tar.extractall(rdir)

    # Decompress .gz files that aren't part of 10x MEX
    for gz_file in os.listdir(rdir):
        if gz_file.endswith(".gz") and not any(
            gz_file.endswith(x) for x in ["matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz"]
        ):
            out_name = gz_file[:-3]
            out_path = os.path.join(rdir, out_name)
            if not os.path.exists(out_path):
                print(f"  decompressing {gz_file} ...")
                gunzip_file(os.path.join(rdir, gz_file), out_path)

    # Convert to h5ad
    print(f"  converting to AnnData ...")
    try:
        adata = detect_and_load(rdir, dataset_id)
    except Exception as e:
        print(f"  conversion failed: {e}")
        return False

    # Tag with source metadata
    adata.obs["geo_accession"] = geo_accession
    adata.obs["dataset_id"]    = dataset_id

    print(f"  shape: {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print(f"  saving raw h5ad to {target_h5ad} ...")
    adata.write_h5ad(target_h5ad)
    print(f"  done.")
    return True


# ─── main ─────────────────────────────────────────────────────────────────────

def print_status(manifest):
    print("\n--- Download status ---")
    for _, row in manifest.iterrows():
        did  = row["dataset_id"]
        acc  = row["geo_accession"]
        if pd.isna(acc) or acc is None:
            print(f"  {did:<25}  (local — no download needed)")
            continue
        h5ad = out_h5ad(did)
        rdir = raw_dir(did)
        if os.path.exists(h5ad):
            size_gb = os.path.getsize(h5ad) / 1e9
            print(f"  {did:<25}  [{acc}]  h5ad ready ({size_gb:.1f} GB)")
        elif os.path.isdir(rdir) and os.listdir(rdir):
            n = len(os.listdir(rdir))
            print(f"  {did:<25}  [{acc}]  raw files present ({n} files), h5ad not yet built")
        else:
            print(f"  {did:<25}  [{acc}]  NOT downloaded")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("accession", nargs="?", help="specific GEO accession to download")
    parser.add_argument("--list", action="store_true", help="show status and exit")
    args = parser.parse_args()

    if not os.path.exists(MANIFEST):
        print("Manifest not found. Run 02_dataset_manifest.py first.")
        sys.exit(1)

    manifest = pd.read_csv(MANIFEST)
    geo_rows = manifest[manifest["geo_accession"].notna()]

    if args.list:
        print_status(manifest)
        return

    if args.accession:
        rows = geo_rows[geo_rows["geo_accession"] == args.accession]
        if rows.empty:
            print(f"Accession {args.accession} not in manifest.")
            sys.exit(1)
    else:
        rows = geo_rows

    results = {}
    for _, row in rows.iterrows():
        did = row["dataset_id"]
        acc = row["geo_accession"]
        ok  = download_geo(did, acc)
        results[did] = "OK" if ok else "FAILED"

    print("\n--- Summary ---")
    for did, status in results.items():
        print(f"  {did:<25}  {status}")

    print_status(manifest)


if __name__ == "__main__":
    main()
