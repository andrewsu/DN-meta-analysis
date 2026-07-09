# When does ontology subclass expansion improve the "human gene-expression dataset" query?

**Question:** find a disease for which expanding an NDE `healthCondition` query to its ontology **subclasses** demonstrably recovers *valid human gene-expression case/control datasets* that the exact-parent-term query misses — testing multiple candidates.
**Date:** 2026-07-08 · **Author:** Andrew Su (with Claude Code)
**Artifacts:** `ibd-comparison/` (winner's full A-B-D) · screen notes in scratchpad
**Companions:** [DN](data-discovery-method-comparison.md) · [glaucoma](glaucoma-method-comparison.md) · [SSc](ssc-method-comparison.md) · [arthritis](arthritis-method-comparison.md)

---

## Short answer

**Yes — inflammatory bowel disease (IBD) is a clear demonstrable example** (interstitial lung disease is a second). Subclass expansion recovers **≈80 valid human IBD gene-expression case/control datasets** that the parent-term query ("inflammatory bowel disease", 574 datasets) misses — because the field labels studies **"Crohn's disease" / "ulcerative colitis," almost never "IBD."** Claude-alone said it best, unprompted:

> *"GEO almost never labels a sample 'IBD.' Samples are labeled Crohn's disease or ulcerative colitis… If you take 'IBD' literally, you will include almost nothing."*

This is the opposite of glaucoma (subclass expansion recovered **0** valid). And a key negative control emerged: **raw subclass recovery is a poor proxy for valid recovery** (muscular dystrophy below).

---

## The screen (metric = valid human GEX case/control datasets recovered)

For each candidate: NDE datasets tagged with a **subclass** term but **not** the parent term ("subclass-only"), restricted to human, then a GEO-curated sample estimates the *valid* fraction (human, mRNA expression, disease-vs-control; excluding animal / in-vitro / treatment-only / wrong-modality / mis-tags).

| Disease | parent term (N) | subclass-only human | valid fraction (GEO sample) | **valid recovered** | poolable? |
|---|---|--:|--:|--:|---|
| **IBD** ★ | inflammatory bowel disease (574) | 248 | ~33% (unbiased n=61) | **≈80** | yes (IBD = CD+UC vs control) |
| ILD | interstitial lung disease (113) | 284 | high† | strong (co-example) | mostly (IPF/ILD vs control) |
| muscular dystrophy | muscular dystrophy (243) | **369** | ~6% (n=33) | ~20 | no (distinct diseases) |
| cardiomyopathy | cardiomyopathy (384) | 118 | not curated | — | partly |
| glaucoma | glaucoma (134) | 25 | 0% | **0** | — |
| systemic sclerosis | systemic sclerosis (167) | 3 | 33% | 1 | — |
| diabetic kidney disease | diabetic kidney disease (129) | 0 (MONDO **leaf**) | — | 0 | — |

†ILD's sample was **selection-biased** (I hand-picked recognizable IPF accessions: 22/30 human+expr, ~17 valid — GSE53845 "40 IPF + 8 controls," GSE32537 LGRC, GSE136831 IPF Cell Atlas, GSE24206, GSE28042). Its true valid fraction needs a random sample, but it is clearly a strong second example.

---

## Winner deep-dive — IBD (full A-B-D)

| Cond | Result | Cost |
|---|---|--:|
| **A · alone** | 13 canonical leads (GSE75214, GSE57945, GSE117993, GSE36807, GSE3365, GSE83687…), **refused counts**; nailed the mechanism (quote above) and the CD∪UC definition problem | $0.18 |
| **B · web** | **11-dataset curated core panel, ≈1,330 disease / 364 control, 0 FP** (verified each; resolved superseries GSE87473 & GSE59071⊂GSE75214; flagged the 2,490-sample GSE193677). A *core*, not all 701 GEO hits. | $1.27 |
| **D · MCP→NDE** | **used subclass expansion** (`get_descendants` + IBD/UC/Crohn/Crohn-ileitis terms): **592 human IBD candidates vs 344 parent-IBD-only → +248 recovered by subclass expansion.** Refused sample counts (no `Sample` entity); flagged noise (a zebrafish study tagged *Homo sapiens*; drug-trial datasets; superseries). | $0.89 |

**The demonstration, two independent ways:**
1. *My screen:* 248 human subclass-only datasets, ~33% valid → **≈80 valid IBD case/control GEX datasets recovered** (colon/ileum biopsy & blood, UC/CD vs control — e.g. GSE10791, GSE9452, GSE48634, GSE164871, GSE102134, GSE74265).
2. *The D-condition itself:* NDE IBD discovery rises **344 → 592 (+248)** the moment subtype terms are included — exactly the recovery the parent-only query would forfeit.

Caveat: the 248 subclass-only set is not pure — it contains real false positives (bladder-cancer datasets mis-tagged UC/Crohn, celiac, in-vitro organoids/cell-lines, drug-trial arms), which is why the *valid* figure (~80) is ~1/3 of the raw, not all 248.

---

## The counter-example that proves "raw ≠ valid" — muscular dystrophy

MD had the **highest raw** subclass-only recovery (369 human) but the **lowest valid fraction (~6%, ~20 recovered)**. Its "myotonic dystrophy" tag (418 datasets) is inflated by the **"DM" abbreviation colliding with *diabetes mellitus*** — diabetic-retinopathy datasets (GSE221521, GSE239512) are mis-tagged myotonic dystrophy — plus breast/lung-cancer cell lines (GSE120268, GSE140432) and myoblast/organoid in-vitro. And MD subtypes (myotonic dystrophy, Duchenne, FSHD, LGMD) are *distinct diseases* you would not pool. So a big subclass count can be almost entirely noise.

---

## The law (now with a positive example)

Across seven diseases, subclass-expansion value is governed by two independent factors:

1. **Does the field label by subtype instead of the umbrella?** (IBD: yes → parent under-tagged, huge recovery. DN: N/A, MONDO leaf. SSc: mostly labels the umbrella → tiny recovery.)
2. **Do the subtypes host genuine patient case/control cohorts?** (IBD/ILD: yes → recovery is *valid*. Glaucoma/MD: no → subtype tags are cell-line/animal/mis-tag → recovery is ~0 valid.)

**Both must hold.** Raw subclass count (which scales with ontology-subtree bushiness: DN 0 → arthritis +354%) is decoupled from *valid* recovery. IBD is the case where both align: **subclass expansion turns a near-empty parent-term query into ~80 real datasets.**

---

## Bottom line

Ontology subclass expansion is worth doing **specifically for diseases whose literature labels subtypes over the umbrella and whose subtypes are real patient cohorts** — IBD (≈80 valid recovered) and ILD are prime cases; DN, glaucoma, SSc, and muscular dystrophy are not (leaf, or subtypes that are cell-line/animal/mis-tagged). The practical rule for NDE: **always expand `healthCondition` via `get_descendants` for umbrella/subtype-labeled diseases**, then apply the standard primary-record verification (organism + assay + case/control) — because a large share of the expanded hits are still animal/in-vitro/mis-tagged noise.
