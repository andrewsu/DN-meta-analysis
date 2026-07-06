#!/usr/bin/env python3
"""Stage 1: per-study differential expression + effect sizes (glomerular DN)."""
import os, warnings, json
import numpy as np, pandas as pd
import GEOparse
import lib
warnings.filterwarnings("ignore")
os.makedirs("data/geo", exist_ok=True); os.makedirs("results", exist_ok=True)

FTP = "https://ftp.ncbi.nlm.nih.gov/geo/series"

def murl(gse, sub=None):
    stub = gse[:-3] + "nnn"
    fn = f"{gse}-{sub}" if sub else gse
    return f"{FTP}/{stub}/{gse}/matrix/{fn}_series_matrix.txt.gz"

# label classifiers: source_name(lower) -> 'DN' | 'CTRL' | None
def cls_30528(s):
    return "DN" if "dkd" in s else ("CTRL" if "control" in s else None)
def cls_96804(s):
    return "DN" if "diabetic nephropathy" in s else ("CTRL" if "nephrectom" in s else None)
def cls_104948(s):
    return "DN" if "glom-dn" in s else ("CTRL" if "glom-ld" in s else None)

DATASETS = [
    dict(id="GSE30528", gse="GSE30528", sub=None, gpl="GPL571",
         kind="gene_symbol_col", classify=cls_30528,
         label="Glomeruli (Woroniecka 2011)"),
    dict(id="GSE96804", gse="GSE96804", sub=None, gpl="GPL17586",
         kind="gene_assignment", classify=cls_96804,
         label="Glomeruli (Pan 2018)"),
    dict(id="GSE104948", gse="GSE104948", sub="GPL22945", gpl="GPL22945",
         kind="symbol_col", classify=cls_104948,
         label="Glomeruli, ERCB (Grayson/Ju 2018)"),
]

summary = []
for d in DATASETS:
    print(f"\n===== {d['id']} : {d['label']} =====")
    mpath = f"data/geo/{d['id']}_matrix.txt.gz"
    lib.download(murl(d["gse"], d["sub"]), mpath)
    expr, meta = lib.parse_series_matrix(mpath)
    gsms = list(expr.columns)
    src = meta.get("Sample_source_name_ch1") or meta.get("Sample_source_name")
    src_map = dict(zip(meta["Sample_geo_accession"], src))
    labels = np.array([d["classify"](str(src_map[g]).lower()) for g in gsms], dtype=object)
    keep = np.array([l in ("DN", "CTRL") for l in labels])
    expr = expr.loc[:, np.array(gsms)[keep]]
    labels = labels[keep]
    is_dn = labels == "DN"
    n_dn, n_ctrl = int(is_dn.sum()), int((~is_dn).sum())
    print(f"  samples kept: DN={n_dn}  CTRL={n_ctrl}  (dropped {int((~keep).sum())})")

    expr, logged = lib.maybe_log2(expr)
    print(f"  log2-transformed on load: {logged}; probes={expr.shape[0]}")

    gpl = GEOparse.get_GEO(geo=d["gpl"], destdir="data/geo", silent=True)
    p2s = lib.build_probe2symbol(gpl.table, d["kind"])
    genes = lib.collapse_to_genes(expr.dropna(how="all"), p2s)
    # drop genes with any NaN across kept samples
    genes = genes.dropna(axis=0)
    print(f"  genes after collapse: {genes.shape[0]}")

    de = lib.moderated_ttest(genes, is_dn)
    eff = lib.hedges_g(genes, is_dn)
    de.to_csv(f"results/de_{d['id']}.csv")
    eff.to_csv(f"results/effect_{d['id']}.csv")
    n_sig = int((de["adj.P.Val"] < 0.05).sum())
    print(f"  moderated-t DE genes FDR<0.05: {n_sig}")
    summary.append(dict(id=d["id"], label=d["label"], gpl=d["gpl"],
                        n_dn=n_dn, n_ctrl=n_ctrl, n_genes=int(genes.shape[0]),
                        n_sig_fdr05=n_sig, logged=bool(logged)))

pd.DataFrame(summary).to_csv("results/dataset_summary.csv", index=False)
print("\n=== dataset summary ===")
print(pd.DataFrame(summary).to_string(index=False))
