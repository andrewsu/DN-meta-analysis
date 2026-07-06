#!/usr/bin/env python3
"""RNA-seq Stage 3: compare RNA-seq meta vs glomerular microarray meta vs the
SuLab/ARCHS4 result; make figures; emit comparison tables."""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
import lib

plt.rcParams.update({"figure.dpi": 130, "font.size": 10, "axes.grid": True, "grid.alpha": .25})
UP, DOWN, NS, IEGC = "#c0392b", "#2471a3", "#c9c9c9", "#e67e22"

marr = pd.read_csv("results/meta_results.csv", index_col=0)          # glomerular microarray
rseq = pd.read_csv("results/rnaseq_meta_results.csv", index_col=0)   # whole-kidney RNA-seq

def robust(df):
    return df[(df.k == 3) & df.dir_consistent & (df["adj.P.Val"] < 0.05) &
              (~df.index.isin(lib.IEG_GENES))]

rob_m, rob_r = robust(marr), robust(rseq)

# ---------- cross-platform concordance ----------
shared = marr.index.intersection(rseq.index)
x = marr.loc[shared, "pooled_g"]; y = rseq.loc[shared, "pooled_g"]
ok = np.isfinite(x) & np.isfinite(y)
r_all = np.corrcoef(x[ok], y[ok])[0, 1]
both_sig = rob_m.index.intersection(rob_r.index)
agree = (np.sign(marr.loc[both_sig, "pooled_g"]) == np.sign(rseq.loc[both_sig, "pooled_g"]))
print(f"shared genes: {ok.sum()}  | Pearson r(pooled_g) = {r_all:.3f}")
print(f"robust DEGs microarray={len(rob_m)}  RNA-seq={len(rob_r)}  overlap={len(both_sig)} "
      f"({agree.mean()*100:.0f}% same direction)")
# hypergeometric-ish enrichment of overlap
N = len(shared); a, b = len(rob_m.index.intersection(shared)), len(rob_r.index.intersection(shared))
exp = a * b / N
print(f"overlap expected by chance ~{exp:.0f}; observed {len(both_sig)} (={len(both_sig)/exp:.1f}x)")

pd.DataFrame({"microarray_g": marr.loc[both_sig, "pooled_g"],
              "rnaseq_g": rseq.loc[both_sig, "pooled_g"]}).to_csv("results/cross_platform_overlap.csv")

# ---------- hallmark genes in both ----------
hall = ["NPHS1", "NPHS2", "PTPRO", "MAGI2", "PODXL", "TGFBI", "COL1A1", "COL1A2",
        "C1QA", "MMP2", "MMP11", "LUM", "VCAM1", "CCL2"]
def look(df, g):
    if g in df.index:
        return f"{df.at[g,'pooled_g']:+.2f} (FDR {df.at[g,'adj.P.Val']:.1e})"
    return "n/a"
print("\n=== hallmark genes: microarray-glom | RNA-seq-whole ===")
for g in hall:
    print(f"  {g:8s}  micro {look(marr,g):26s}  rnaseq {look(rseq,g)}")

# ============================ FIGURES ============================
# Fig A: the IEG-artifact story (RNA-seq) --------------------------
m = rseq[rseq.k == 3].copy()
isieg = m.index.isin(lib.IEG_GENES)
fig, (a1, a2) = plt.subplots(1, 2, figsize=(14, 6))
# volcano
y1 = -np.log10(m["adj.P.Val"].clip(lower=1e-300))
a1.scatter(m.loc[~isieg, "pooled_g"], y1[~isieg], s=6, c=NS, alpha=.5, linewidths=0)
a1.scatter(m.loc[isieg, "pooled_g"], y1[isieg], s=40, c=IEGC, edgecolor="k", linewidths=.4, label="immediate-early / stress gene", zorder=5)
a1.axhline(-np.log10(0.05), ls="--", c="k", lw=.8)
for g in ["FOS", "FOSB", "EGR1", "NR4A1", "DUSP1"]:
    if g in m.index:
        a1.annotate(g, (m.at[g, "pooled_g"], -np.log10(max(m.at[g, "adj.P.Val"], 1e-300))), fontsize=8, color=IEGC, weight="bold")
a1.set_xlabel("pooled Hedges' g (DN vs control)"); a1.set_ylabel(r"$-\log_{10}$ FDR")
a1.set_title("RNA-seq meta: ARCHS4 top genes have large g but fail FDR"); a1.legend(frameon=False, loc="upper left")
# |g| vs I2
a2.scatter(np.abs(m.loc[~isieg, "pooled_g"]), m.loc[~isieg, "I2"], s=6, c=NS, alpha=.5, linewidths=0)
a2.scatter(np.abs(m.loc[isieg, "pooled_g"]), m.loc[isieg, "I2"], s=40, c=IEGC, edgecolor="k", linewidths=.4, zorder=5, label="IEG / stress")
a2.axhline(50, ls="--", c="k", lw=.6)
for g in ["FOS", "NR4A1", "DUSP1", "FOSB", "EGR1", "ATF3"]:
    if g in m.index: a2.annotate(g, (abs(m.at[g, "pooled_g"]), m.at[g, "I2"]), fontsize=8, color=IEGC, weight="bold")
a2.set_xlabel("|pooled Hedges' g|"); a2.set_ylabel("I²  between-study heterogeneity (%)")
a2.set_title("IEGs: big effect, extreme heterogeneity (I²≈95%) → not reproducible"); a2.legend(frameon=False, loc="lower right")
fig.tight_layout(); fig.savefig("results/fig5_ieg_artifact.png"); plt.close(fig)

# Fig B: cross-platform concordance -------------------------------
fig, ax = plt.subplots(figsize=(6.8, 6.4))
ax.scatter(x[ok], y[ok], s=5, c=NS, alpha=.25, linewidths=0)
sig_both = shared.intersection(both_sig)
ax.scatter(marr.loc[sig_both, "pooled_g"], rseq.loc[sig_both, "pooled_g"], s=14,
           c=["#c0392b" if v > 0 else "#2471a3" for v in marr.loc[sig_both, "pooled_g"]], alpha=.8, linewidths=0)
ax.axhline(0, c="k", lw=.5); ax.axvline(0, c="k", lw=.5)
ax.plot([-3, 3], [-3, 3], ls=":", c="k", lw=.6)
ax.set_xlabel("pooled g — glomerular microarray meta"); ax.set_ylabel("pooled g — whole-kidney RNA-seq meta")
ax.set_title(f"Cross-modality concordance (r={r_all:.2f}; {len(both_sig)} shared DEGs)")
for g in ["NPHS1", "PTPRO", "MAGI2", "TGFBI", "MMP11", "COL1A1", "C1QA", "FBXO21"]:
    if g in sig_both or (g in shared and g in rob_r.index):
        if g in shared:
            ax.annotate(g, (marr.at[g, "pooled_g"], rseq.at[g, "pooled_g"]), fontsize=7)
fig.tight_layout(); fig.savefig("results/fig6_cross_platform.png"); plt.close(fig)

# Fig C: forest — artifact genes vs reproduced genes --------------
eff = {i: pd.read_csv(f"results/rnaseq_effect_{i}.csv", index_col=0) for i in ["GSE142025", "GSE162830", "GSE166239"]}
genes = ["NR4A1", "DUSP1", "FOS", "PTPRO", "MMP11", "TACC2"]
fig, axes = plt.subplots(2, 3, figsize=(13, 6.5))
for ax, g in zip(axes.ravel(), genes):
    ys = np.arange(3, 0, -1)
    for yv, i in zip(ys, eff):
        if g in eff[i].index:
            gi = eff[i].at[g, "g"]; se = np.sqrt(eff[i].at[g, "var_g"])
            c = UP if gi > 0 else DOWN
            ax.plot([gi - 1.96 * se, gi + 1.96 * se], [yv, yv], c=c, lw=1.4); ax.plot(gi, yv, "s", c=c, ms=6)
    pg, se = rseq.at[g, "pooled_g"], rseq.at[g, "se"]
    ax.fill([pg - 1.96 * se, pg, pg + 1.96 * se, pg], [0, .22, 0, -.22], c="k")
    ax.axvline(0, c="k", lw=.6)
    ax.set_yticks(list(ys) + [0]); ax.set_yticklabels(list(eff) + ["POOLED"], fontsize=7)
    tag = "IEG artifact" if g in lib.IEG_GENES else "reproducible"
    ax.set_title(f"{g} — {tag}\nFDR={rseq.at[g,'adj.P.Val']:.2f}, I²={rseq.at[g,'I2']:.0f}%", fontsize=9)
    ax.set_xlabel("Hedges' g")
fig.suptitle("Why IEGs fail: consistent-looking but heterogeneous across studies (I²≈95%)", y=1.02)
fig.tight_layout(); fig.savefig("results/fig7_rnaseq_forest.png", bbox_inches="tight"); plt.close(fig)

print("\nfigures: fig5_ieg_artifact, fig6_cross_platform, fig7_rnaseq_forest")
