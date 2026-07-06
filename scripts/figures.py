#!/usr/bin/env python3
"""Stage 3: figures + significance-ranked tables."""
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

plt.rcParams.update({"figure.dpi": 130, "font.size": 10, "axes.grid": True,
                     "grid.alpha": .25})
UP, DOWN, NS = "#c0392b", "#2471a3", "#b0b0b0"
IDS = ["GSE30528", "GSE96804", "GSE104948"]

meta = pd.read_csv("results/meta_results.csv", index_col=0)
eff = {i: pd.read_csv(f"results/effect_{i}.csv", index_col=0) for i in IDS}

# significance-ranked robust set (present in all 3, consistent direction)
robust = meta[(meta.k == 3) & (meta.dir_consistent) & (meta["adj.P.Val"] < 0.05)].copy()
robust["neglog10fdr"] = -np.log10(robust["adj.P.Val"].clip(lower=1e-300))
robust = robust.sort_values(["adj.P.Val", "pooled_g"], key=lambda s: s if s.name!="pooled_g" else -s.abs())
robust.sort_values("P.Value").to_csv("results/robust_meta_DEGs.csv")

# ----- 1. Volcano ------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 6.5))
m = meta[meta.k == 3].copy()
y = -np.log10(m["adj.P.Val"].clip(lower=1e-300))
sig = (m["adj.P.Val"] < 0.05) & m["dir_consistent"]
ax.scatter(m.loc[~sig, "pooled_g"], y[~sig], s=6, c=NS, alpha=.4, linewidths=0, label="ns")
ax.scatter(m.loc[sig & (m.pooled_g > 0), "pooled_g"], y[sig & (m.pooled_g > 0)], s=8, c=UP, alpha=.6, linewidths=0, label="up in DN")
ax.scatter(m.loc[sig & (m.pooled_g < 0), "pooled_g"], y[sig & (m.pooled_g < 0)], s=8, c=DOWN, alpha=.6, linewidths=0, label="down in DN")
ax.axhline(-np.log10(0.05), ls="--", c="k", lw=.8)
lab = pd.concat([robust.sort_values("pooled_g").head(9),
                 robust.sort_values("pooled_g").tail(9),
                 robust.loc[robust.index.intersection(["NPHS1", "TGFBI", "C1QA", "COL1A1", "MMP2"])]])
for g, r in lab.iterrows():
    ax.annotate(g, (r["pooled_g"], -np.log10(max(r["adj.P.Val"], 1e-300))),
                fontsize=7, ha="center", va="bottom")
ax.set_xlabel("Pooled effect size (random-effects Hedges' g)  —  DN vs control")
ax.set_ylabel(r"$-\log_{10}$ FDR")
ax.set_title("Diabetic-nephropathy glomerular meta-analysis (3 studies)")
ax.legend(loc="upper center", ncol=3, frameon=False)
fig.tight_layout(); fig.savefig("results/fig1_volcano.png"); plt.close(fig)

# ----- 2. Forest plots for hallmark / robust genes --------------------------
genes = ["NPHS1", "TGFBI", "C1QA", "COL1A1", "MMP2", "EZH2", "TJP1", "TSPYL5"]
genes = [g for g in genes if g in robust.index or g in meta.index]
fig, axes = plt.subplots(2, 4, figsize=(15, 7), sharex=False)
for ax, g in zip(axes.ravel(), genes):
    gs, los, his = [], [], []
    for i in IDS:
        if g in eff[i].index:
            gi, vi = eff[i].at[g, "g"], eff[i].at[g, "var_g"]
            se = np.sqrt(vi); gs.append(gi); los.append(gi - 1.96 * se); his.append(gi + 1.96 * se)
        else:
            gs.append(np.nan); los.append(np.nan); his.append(np.nan)
    ypos = np.arange(len(IDS), 0, -1)
    col = [UP if v > 0 else DOWN for v in gs]
    for yv, gv, lo, hi, c in zip(ypos, gs, los, his, col):
        ax.plot([lo, hi], [yv, yv], c=c, lw=1.5)
        ax.plot(gv, yv, "s", c=c, ms=6)
    # pooled diamond
    pg, se = meta.at[g, "pooled_g"], meta.at[g, "se"]
    lo, hi = pg - 1.96 * se, pg + 1.96 * se
    ax.fill([lo, pg, hi, pg], [0, 0.25, 0, -0.25], c="k")
    ax.axvline(0, c="k", lw=.6, ls="-")
    ax.set_yticks(list(ypos) + [0]); ax.set_yticklabels(IDS + ["POOLED"], fontsize=8)
    ax.set_ylim(-0.8, len(IDS) + 0.8)
    ax.set_title(f"{g}  (FDR={meta.at[g,'adj.P.Val']:.1e}, I²={meta.at[g,'I2']:.0f}%)", fontsize=9)
    ax.set_xlabel("Hedges' g")
fig.suptitle("Per-study effect sizes and random-effects pooled estimate (DN vs control)", y=1.02)
fig.tight_layout(); fig.savefig("results/fig2_forest.png", bbox_inches="tight"); plt.close(fig)

# ----- 3. Top genes bar chart (by significance) -----------------------------
topN = 20
up = robust[robust.pooled_g > 0].sort_values("P.Value").head(topN)
dn = robust[robust.pooled_g < 0].sort_values("P.Value").head(topN)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 7))
a1.barh(range(len(up))[::-1], up["pooled_g"], color=UP); a1.set_yticks(range(len(up))[::-1]); a1.set_yticklabels(up.index, fontsize=8)
a1.set_title(f"Top {topN} UP in DN (by meta p-value)"); a1.set_xlabel("pooled Hedges' g")
a2.barh(range(len(dn))[::-1], dn["pooled_g"], color=DOWN); a2.set_yticks(range(len(dn))[::-1]); a2.set_yticklabels(dn.index, fontsize=8)
a2.set_title(f"Top {topN} DOWN in DN (by meta p-value)"); a2.set_xlabel("pooled Hedges' g")
fig.tight_layout(); fig.savefig("results/fig3_topgenes.png"); plt.close(fig)

# ----- 4. Cross-study concordance -------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(14, 4.6))
pairs = [("GSE30528", "GSE96804"), ("GSE30528", "GSE104948"), ("GSE96804", "GSE104948")]
for ax, (a, b) in zip(axes, pairs):
    common = eff[a].index.intersection(eff[b].index)
    x, yv = eff[a].loc[common, "logFC"], eff[b].loc[common, "logFC"]
    ok = np.isfinite(x) & np.isfinite(yv)
    r = np.corrcoef(x[ok], yv[ok])[0, 1]
    ax.scatter(x, yv, s=4, alpha=.2, c="#444")
    ax.axhline(0, c="k", lw=.5); ax.axvline(0, c="k", lw=.5)
    ax.set_xlabel(f"log2FC {a}"); ax.set_ylabel(f"log2FC {b}")
    ax.set_title(f"r = {r:.2f}")
fig.suptitle("Cross-study concordance of per-gene log2 fold-changes (DN vs control)")
fig.tight_layout(); fig.savefig("results/fig4_concordance.png"); plt.close(fig)

print("figures written to results/  (fig1..fig4)")
print(f"robust DEGs (k=3, consistent, FDR<0.05): {len(robust)}")
print("\nMost significant UP:", ", ".join(up.index[:15]))
print("Most significant DOWN:", ", ".join(dn.index[:15]))
# also emit compact tables for the report
for name, tab in [("up", up), ("dn", dn)]:
    tab[["k","pooled_g","se","adj.P.Val","I2","mean_logFC","g_GSE30528","g_GSE96804","g_GSE104948"]]\
        .round(3).to_csv(f"results/report_top_{name}.csv")
