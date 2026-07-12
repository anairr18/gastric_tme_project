"""Download GPL probe->gene mappings for GEO datasets.

IMPORTANT: GEO's "{GPL}_family.soft.gz" bundles every GSM/GSE record ever
linked to that platform, not just the platform's own annotation table.
For a heavily-used platform like GPL6947 that file is 11.5 GB compressed
(confirmed via HEAD request) -- downloading and parsing it in memory ran
a python process up to 12+ GB RAM with no end in sight. Two fixes:
  1. GPL6947 has a dedicated, curated "annot.gz" file (~7 MB) containing
     just the ID -> Gene symbol table. Use that instead of the family file.
  2. GPL8432 has no annot.gz (404), so we fall back to streaming the
     family.soft.gz and STOP as soon as we hit the first
     "platform_table_end" marker -- the platform's own data table is
     the very first record in the file (confirmed: fully contained
     within the first 20 MB of an 896 MB compressed file), so we never
     download or decompress the rest.
"""
import requests, gzip, os
import pandas as pd

RAW_DIR = os.path.expanduser("~/gastric_tme_project/data/external/geo_bulk")


def _parse_table_stream(line_iter, gpl_id):
    mapping = {}
    in_table = False
    headers = None
    id_col = sym_col = None

    for line in line_iter:
        lo = line.lower()
        if "platform_table_begin" in lo:
            in_table = True
            headers = None
            continue
        if "platform_table_end" in lo:
            break  # done -- this is the platform's own table, stop reading
        if not in_table:
            continue
        parts = line.split("\t")
        if headers is None:
            headers = [h.strip().upper() for h in parts]
            for i, h in enumerate(headers):
                if h == "ID":
                    id_col = i
            for i, h in enumerate(headers):
                if h == "GENE SYMBOL" or h == "SYMBOL":
                    sym_col = i
                    break
            if sym_col is None:
                for i, h in enumerate(headers):
                    if "SYMBOL" in h:
                        sym_col = i
                        break
            if sym_col is None:
                for i, h in enumerate(headers):
                    if h == "ILMN_GENE":
                        sym_col = i
                        break
            print(f"  Columns: {headers[:6]} | id={id_col}, sym={sym_col}")
            continue
        if id_col is None or sym_col is None or len(parts) <= max(id_col, sym_col):
            continue
        pid = parts[id_col].strip()
        sym = parts[sym_col].strip()
        if pid and sym and sym not in ("", "---", "NA", "N/A", "0"):
            mapping[pid] = sym

    return mapping


def download_via_annot(gpl_id, prefix):
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{prefix}/{gpl_id}/annot/{gpl_id}.annot.gz"
    print(f"  Trying curated annot file: {url}")
    r = requests.get(url, timeout=120, stream=True)
    if r.status_code == 404:
        print("  No annot.gz available for this platform")
        return None
    r.raise_for_status()
    content = b"".join(r.iter_content(1024 * 1024))
    print(f"  Got {len(content)/1e6:.1f} MB (annot)")

    def line_iter():
        with gzip.open(__import__("io").BytesIO(content)) as gz:
            for raw in gz:
                yield raw.decode("utf-8", errors="replace").rstrip("\n")

    return _parse_table_stream(line_iter(), gpl_id)


def download_via_family_stream(gpl_id, prefix):
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/platforms/{prefix}/{gpl_id}/soft/{gpl_id}_family.soft.gz"
    print(f"  Streaming family file (stop after first platform table): {url}")
    r = requests.get(url, timeout=300, stream=True)
    r.raise_for_status()

    def line_iter():
        with gzip.open(r.raw) as gz:
            for raw in gz:
                yield raw.decode("utf-8", errors="replace").rstrip("\n")

    mapping = _parse_table_stream(line_iter(), gpl_id)
    r.close()
    return mapping


def download_gpl_mapping(gpl_id):
    cache_tsv = os.path.join(RAW_DIR, f"{gpl_id}_mapping.tsv")
    if os.path.exists(cache_tsv):
        df = pd.read_csv(cache_tsv, sep="\t")
        m = dict(zip(df["probe_id"].astype(str), df["symbol"].astype(str)))
        m.update({f"ILMN_{k}": v for k, v in m.items() if not k.startswith("ILMN_")})
        print(f"  {gpl_id}: loaded {len(df):,} probe mappings from cache")
        return m

    num = int(gpl_id[3:])
    prefix = f"GPL{num // 1000}nnn"

    mapping = download_via_annot(gpl_id, prefix)
    if not mapping:
        mapping = download_via_family_stream(gpl_id, prefix)

    n = len(mapping)
    print(f"  {gpl_id}: {n:,} probes mapped")
    rows = list(mapping.items())
    pd.DataFrame(rows, columns=["probe_id", "symbol"]).to_csv(cache_tsv, sep="\t", index=False)

    mapping.update({f"ILMN_{k}": v for k, v in mapping.items() if not k.startswith("ILMN_")})
    return mapping


if __name__ == "__main__":
    for gpl in ["GPL6947", "GPL8432"]:
        m = download_gpl_mapping(gpl)
        test = ["ILMN_1343048", "ILMN_1651209", "ILMN_2344887"]
        for t in test:
            print(f"  {t} -> {m.get(t, 'NOT FOUND')}")
        print()

    print("Done.")
