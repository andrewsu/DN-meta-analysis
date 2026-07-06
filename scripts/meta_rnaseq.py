#!/usr/bin/env python3
"""RNA-seq Stage 2: random-effects meta-analysis + immediate-early-gene (IEG)
artifact handling. Mirrors scripts/meta.py but on the whole-kidney RNA-seq set."""
import numpy as np, pandas as pd
from scipy import stats
import lib

IDS = ["GSE142025", "GSE162830", "GSE166239"]

effects = {i: pd.read_csv(f"results/rnaseq_effect_{i}.csv", index_col=0) for i in IDS}
meta = lib.meta_analyze(effects, min_studies=2)
for i in IDS:
    meta[f"g_{i}"] = effects[i]["g"].reindex(meta.index)

# Stouffer cross-check
de = {i: pd.read_csv(f"results/rnaseq_de_{i}.csv", index_col=0) for i in IDS}
allg = sorted(set().union(*[set(d.index) for d in de.values()]))
sz = pd.DataFrame(index=allg); w = []
for i in IDS:
    d = de[i].reindex(allg)
    sz[i] = np.sign(d["logFC"]) * stats.norm.isf(d["P.Value"].clip(lower=1e-300) / 2)
    w.append(np.sqrt(d["n_dn"].dropna().iloc[0] + d["n_ctrl"].dropna().iloc[0]))
meta = meta.join(lib.stouffer(sz, w))

meta["is_IEG"] = meta.index.isin(lib.IEG_GENES)
meta.to_csv("results/rnaseq_meta_results.csv")

# primary high-confidence set: all 3 studies, consistent, FDR<0.05, NOT an IEG artifact
robust = meta[(meta.k == 3) & meta.dir_consistent & (meta["adj.P.Val"] < 0.05) & (~meta.is_IEG)]
robust.sort_values("P.Value").to_csv("results/rnaseq_robust_DEGs.csv")

print(f"genes meta-analyzed (k>=2): {len(meta)};  in all 3: {(meta.k==3).sum()}")
print(f"high-confidence DEGs (k=3, FDR<0.05, consistent, non-IEG): {len(robust)}  "
      f"[{(robust.pooled_g>0).sum()} up / {(robust.pooled_g<0).sum()} down]")
conc = (np.sign(meta.pooled_g) == np.sign(meta.stouffer_z))
print(f"sign concordance RE vs Stouffer: {conc.mean():.3f}")

cols = ["k", "pooled_g", "adj.P.Val", "I2", "mean_logFC", "g_GSE142025", "g_GSE162830", "g_GSE166239"]
print("\n=== TOP 20 UP in DN (non-IEG, by meta p) ===")
print(robust[robust.pooled_g > 0].sort_values("P.Value").head(20)[cols].round(3).to_string())
print("\n=== TOP 20 DOWN in DN (non-IEG, by meta p) ===")
print(robust[robust.pooled_g < 0].sort_values("P.Value").head(20)[cols].round(3).to_string())

# ---- IEG artifact quantification --------------------------------------------
ieg = meta[meta.is_IEG].sort_values("pooled_g")
print("\n=== IMMEDIATE-EARLY / STRESS GENES (flagged & excluded from primary) ===")
print(ieg[["k", "pooled_g", "adj.P.Val", "I2"]].round(3).to_string())
n_ieg_sig = int((ieg["adj.P.Val"] < 0.05).sum())
print(f"IEGs significant (FDR<0.05): {n_ieg_sig}/{len(ieg)}; "
      f"down: {int((ieg.pooled_g<0).sum())} up: {int((ieg.pooled_g>0).sum())}")

# ---- where do the SuLab/ARCHS4 top genes land here? -------------------------
sulab = ["FOS", "FOSB", "EGR1", "NR4A1", "DUSP1", "OAS2", "RSAD2"]
have = [g for g in sulab if g in meta.index]
print("\n=== SuLab/ARCHS4 reported genes in THIS RNA-seq meta ===")
print(meta.loc[have, ["k", "pooled_g", "adj.P.Val", "I2", "is_IEG"]].round(3).to_string())
