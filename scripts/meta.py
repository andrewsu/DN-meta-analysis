#!/usr/bin/env python3
"""Stage 2: random-effects meta-analysis across studies (+ Stouffer cross-check)."""
import numpy as np, pandas as pd
from scipy import stats
import lib

IDS = ["GSE30528", "GSE96804", "GSE104948"]

# ---- primary: random-effects meta of Hedges' g -----------------------------
effects = {i: pd.read_csv(f"results/effect_{i}.csv", index_col=0) for i in IDS}
meta = lib.meta_analyze(effects, min_studies=2)

# per-study g / logFC columns for the top-gene forest plots + transparency
for i in IDS:
    meta[f"g_{i}"] = effects[i]["g"].reindex(meta.index)
    meta[f"logFC_{i}"] = effects[i]["logFC"].reindex(meta.index)

# ---- secondary: weighted Stouffer of per-study moderated-t p-values --------
de = {i: pd.read_csv(f"results/de_{i}.csv", index_col=0) for i in IDS}
all_genes = sorted(set().union(*[set(d.index) for d in de.values()]))
signed_z = pd.DataFrame(index=all_genes)
weights = []
for i in IDS:
    d = de[i].reindex(all_genes)
    p = d["P.Value"].clip(lower=1e-300)
    z = np.sign(d["logFC"]) * stats.norm.isf(p / 2.0)
    signed_z[i] = z
    weights.append(np.sqrt(d["n_dn"].dropna().iloc[0] + d["n_ctrl"].dropna().iloc[0]))
stouf = lib.stouffer(signed_z, weights)
meta = meta.join(stouf)

meta.to_csv("results/meta_results.csv")

# ---- summaries --------------------------------------------------------------
sig = meta[(meta["adj.P.Val"] < 0.05) & (meta["dir_consistent"])]
sig_all3 = sig[sig["k"] == 3]
print(f"genes meta-analyzed (k>=2): {len(meta)}")
print(f"  present in all 3 studies (k=3): {(meta['k']==3).sum()}")
print(f"significant FDR<0.05 & direction-consistent: {len(sig)}")
print(f"  of those, present in all 3 studies: {len(sig_all3)}")
print(f"  up in DN: {(sig_all3['pooled_g']>0).sum()}   down in DN: {(sig_all3['pooled_g']<0).sum()}")
conc = (np.sign(meta['pooled_g']) == np.sign(meta['stouffer_z']))
print(f"sign concordance RE vs Stouffer: {conc.mean():.3f}")

cols = ["k", "pooled_g", "se", "P.Value", "adj.P.Val", "I2", "mean_logFC",
        "dir_consistent", "stouffer_p", "g_GSE30528", "g_GSE96804", "g_GSE104948"]
up = sig_all3.sort_values("pooled_g", ascending=False).head(25)
dn = sig_all3.sort_values("pooled_g").head(25)
print("\n=== TOP 25 UP IN DN (k=3, FDR<0.05, consistent) ===")
print(up[cols].round(3).to_string())
print("\n=== TOP 25 DOWN IN DN (k=3, FDR<0.05, consistent) ===")
print(dn[cols].round(3).to_string())

up.to_csv("results/top_up.csv"); dn.to_csv("results/top_down.csv")
