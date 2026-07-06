#!/usr/bin/env python3
"""Reconnaissance: pull ONLY the metadata header lines of candidate GEO series
matrices to confirm platform + sample composition (no expression download)."""
import gzip, io, sys, urllib.request

CANDIDATES = [
    "GSE30528",  # Woroniecka glomeruli
    "GSE30529",  # Woroniecka tubulointerstitium
    "GSE1009",   # glomeruli DN (old)
    "GSE47183",  # ERCB glomeruli (multi-disease)
    "GSE47184",  # ERCB tubulointerstitium
    "GSE104948", # ERCB glomerular (multi CKD)
    "GSE104954", # ERCB tubulointerstitial (multi CKD)
    "GSE96804",  # glomeruli RNA-seq DN
    "GSE142025", # kidney RNA-seq DN
]

def matrix_url(gse):
    stub = gse[:-3] + "nnn"
    return f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/matrix/{gse}_series_matrix.txt.gz"

def probe(gse):
    url = matrix_url(gse)
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "recon/1.0"})
        raw = urllib.request.urlopen(req, timeout=60).read()
    except Exception as e:
        print(f"\n### {gse}: ERROR {e}")
        return
    txt = gzip.decompress(raw).decode("utf-8", "replace")
    keep = {}
    n_samples = None
    for line in txt.splitlines():
        if line.startswith("!series_matrix_table_begin"):
            break
        if line.startswith("!Series_title") or line.startswith("!Series_platform_id") \
           or line.startswith("!Sample_title") or line.startswith("!Sample_source_name") \
           or line.startswith("!Sample_characteristics") or line.startswith("!Platform"):
            key = line.split("\t", 1)[0].strip("!")
            vals = [v.strip('"') for v in line.split("\t")[1:]]
            keep.setdefault(key, []).append(vals)
            if key == "Sample_title":
                n_samples = len(vals)
    print(f"\n### {gse}  (n_samples={n_samples})  {matrix_url(gse)}")
    if "Series_title" in keep:
        print("  TITLE:", keep["Series_title"][0][0])
    if "Series_platform_id" in keep:
        print("  PLATFORM:", ";".join(v for row in keep["Series_platform_id"] for v in row))
    # summarize source names & characteristics (unique value counts)
    for key in ("Sample_source_name_ch1", "Sample_source_name"):
        if key in keep:
            from collections import Counter
            c = Counter(keep[key][0])
            print("  SOURCE:", dict(c))
    for i, row in enumerate(keep.get("Sample_characteristics_ch1", [])):
        from collections import Counter
        c = Counter(row)
        # only print if it looks discriminative (<= 12 unique values)
        if 1 < len(c) <= 12:
            print(f"  CHAR[{i}]:", dict(c))

if __name__ == "__main__":
    for gse in CANDIDATES:
        probe(gse)
