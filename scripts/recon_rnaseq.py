#!/usr/bin/env python3
"""Recon candidate human DN kidney RNA-seq studies: groups + supplementary counts."""
import gzip, urllib.request, re
from collections import Counter

CANDIDATES = ["GSE142025", "GSE175759", "GSE162830", "GSE185011", "GSE199437", "GSE204880"]

def stub(g): return g[:-3] + "nnn"

def get(url, timeout=60):
    return urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "r/1"}), timeout=timeout).read()

def header(gse):
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub(gse)}/{gse}/matrix/{gse}_series_matrix.txt.gz"
    try:
        txt = gzip.decompress(get(url)).decode("utf-8", "replace")
    except Exception as e:
        return f"  (matrix header error: {e})"
    keep = {}
    for line in txt.splitlines():
        if line.startswith("!series_matrix_table_begin"): break
        if line.startswith(("!Series_title", "!Sample_organism_ch1", "!Sample_source_name_ch1",
                            "!Sample_characteristics_ch1", "!Sample_title", "!Series_platform_id")):
            k = line.split("\t", 1)[0].strip("!"); v = [x.strip('"') for x in line.split("\t")[1:]]
            keep.setdefault(k, []).append(v)
    out = []
    if "Series_title" in keep: out.append("  TITLE: " + keep["Series_title"][0][0])
    if "Series_platform_id" in keep: out.append("  PLAT: " + ";".join(keep["Series_platform_id"][0]))
    if "Sample_organism_ch1" in keep: out.append("  ORG: " + str(dict(Counter(keep["Sample_organism_ch1"][0]))))
    n = len(keep.get("Sample_title", [[]])[0]); out.append(f"  N samples: {n}")
    if "Sample_source_name_ch1" in keep:
        out.append("  SRC: " + str(dict(Counter(keep["Sample_source_name_ch1"][0]))))
    for i, row in enumerate(keep.get("Sample_characteristics_ch1", [])):
        c = Counter(row)
        if 1 < len(c) <= 12: out.append(f"  CHAR[{i}]: {dict(c)}")
    return "\n".join(out)

def suppl(gse):
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub(gse)}/{gse}/suppl/"
    try:
        html = get(url, 40).decode("utf-8", "replace")
    except Exception as e:
        return [f"(suppl error: {e})"]
    files = sorted(set(re.findall(r'href="([^"]+)"', html)))
    return [f for f in files if not f.startswith(("/", "?")) and f not in ("", "..")]

for gse in CANDIDATES:
    print(f"\n===== {gse} =====")
    print(header(gse))
    print("  SUPPL:", suppl(gse))
