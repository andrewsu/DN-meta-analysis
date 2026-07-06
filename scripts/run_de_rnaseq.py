#!/usr/bin/env python3
"""RNA-seq Stage 1: uniform processing + per-study DE for human whole-kidney
DN cohorts (GSE142025, GSE162830, GSE166239). Study-matched design; each study
is normalized independently and only standardized effects are pooled later."""
import os, io, gzip, tarfile, json, urllib.request, warnings, re
import numpy as np, pandas as pd
import lib
warnings.filterwarnings("ignore")
os.makedirs("data/geo", exist_ok=True); os.makedirs("results", exist_ok=True)
UA = {"User-Agent": "dn-meta/1.0"}


def _get(url, timeout=180):
    return urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout).read()


def gsm_group_map(gse, char_index):
    """GSM -> raw group string, from the series-matrix metadata header."""
    stub = gse[:-3] + "nnn"
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/matrix/{gse}_series_matrix.txt.gz"
    txt = gzip.decompress(_get(url)).decode("utf-8", "replace")
    acc, chars = None, []
    for line in txt.splitlines():
        if line.startswith("!series_matrix_table_begin"): break
        if line.startswith("!Sample_geo_accession"):
            acc = [v.strip('"') for v in line.split("\t")[1:]]
        if line.startswith("!Sample_characteristics_ch1"):
            chars.append([v.strip('"') for v in line.split("\t")[1:]])
    return dict(zip(acc, chars[char_index]))


def ensembl_to_symbol(ids):
    """Map Ensembl gene IDs -> symbols via mygene.info (batched POST)."""
    clean = [i.split(".")[0] for i in ids]
    out = {}
    for k in range(0, len(clean), 900):
        chunk = clean[k:k + 900]
        body = ("q=" + ",".join(chunk) + "&scopes=ensembl.gene&fields=symbol&species=human").encode()
        req = urllib.request.Request("https://mygene.info/v3/query", data=body,
                                     headers={**UA, "Content-Type": "application/x-www-form-urlencoded"})
        for d in json.loads(urllib.request.urlopen(req, timeout=90).read()):
            if "symbol" in d and d.get("query"):
                out[d["query"]] = d["symbol"]
    return {orig: out[c] for orig, c in zip(ids, clean) if c in out}


# ---------------------------------------------------------------- loaders ----
def load_GSE142025():
    """Per-sample symbol-keyed normalized-log files inside RAW.tar."""
    tf = tarfile.open("data/geo/GSE142025_RAW.tar")
    series = {}
    for m in tf.getnames():
        if not m.endswith(".gz"): continue
        gsm = m.split("_")[0]
        s = pd.read_csv(io.BytesIO(gzip.decompress(tf.extractfile(m).read())),
                        sep="\t", index_col=0)
        series[gsm] = s.iloc[:, 0]
    expr = pd.DataFrame(series)
    grp = gsm_group_map("GSE142025", 0)   # 'group: Advanced_DN|Early_DN|Control'
    lab = {g: ("DN" if "DN" in grp[g] else ("CTRL" if "Control" in grp[g] else None)) for g in expr.columns}
    return expr, lab


def load_GSE162830():
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE162nnn/GSE162830/suppl/GSE162830_ING_quantile_normalized_final.csv.gz"
    lib.download(url, "data/geo/GSE162830.csv.gz")
    expr = pd.read_csv("data/geo/GSE162830.csv.gz", index_col=0)
    lab = {c: ("DN" if c.startswith("DIA") else ("CTRL" if c.startswith("REF") else None)) for c in expr.columns}
    return expr, lab


def load_GSE166239():
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE166nnn/GSE166239/suppl/GSE166239_Nordbo_et_al_counts.txt.gz"
    lib.download(url, "data/geo/GSE166239_counts.txt.gz")
    counts = pd.read_csv("data/geo/GSE166239_counts.txt.gz", sep="\t", index_col=0)
    logcpm = lib.counts_to_log2cpm(counts)
    m = ensembl_to_symbol(list(logcpm.index))
    logcpm = logcpm.loc[[i for i in logcpm.index if i in m]]
    logcpm.index = [m[i] for i in logcpm.index]
    lab = {c: ("DN" if c.startswith("T2DN") else ("CTRL" if c.startswith("Ctrl") else None)) for c in logcpm.columns}
    return logcpm, lab


LOADERS = {
    "GSE142025": (load_GSE142025, "Whole kidney biopsy (Fan 2019)"),
    "GSE162830": (load_GSE162830, "Whole kidney, DN vs reference (2021)"),
    "GSE166239": (load_GSE166239, "Whole FFPE biopsy (Nordbo 2021)"),
}

summary = []
for gse, (loader, label) in LOADERS.items():
    print(f"\n===== {gse} : {label} =====")
    expr, lab = loader()
    keep = [c for c in expr.columns if lab.get(c) in ("DN", "CTRL")]
    expr = expr[keep]
    is_dn = np.array([lab[c] == "DN" for c in keep])
    print(f"  samples: DN={int(is_dn.sum())} CTRL={int((~is_dn).sum())} (dropped {len(lab)-len(keep)})")
    genes = lib.collapse_by_symbol(expr).dropna(axis=0)
    genes, logged = lib.maybe_log2(genes)
    print(f"  genes={genes.shape[0]}  log2-applied-now={logged}")
    de = lib.moderated_ttest(genes, is_dn)
    eff = lib.hedges_g(genes, is_dn)
    de.to_csv(f"results/rnaseq_de_{gse}.csv")
    eff.to_csv(f"results/rnaseq_effect_{gse}.csv")
    nsig = int((de["adj.P.Val"] < 0.05).sum())
    print(f"  moderated-t FDR<0.05: {nsig}")
    summary.append(dict(id=gse, label=label, n_dn=int(is_dn.sum()),
                        n_ctrl=int((~is_dn).sum()), n_genes=int(genes.shape[0]), n_sig=nsig))

pd.DataFrame(summary).to_csv("results/rnaseq_dataset_summary.csv", index=False)
print("\n=== RNA-seq dataset summary ===")
print(pd.DataFrame(summary).to_string(index=False))
