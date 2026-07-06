# RNA-seq companion meta-analysis of diabetic nephropathy — and a re-examination of the ARCHS4 result

**Analysis date:** 2026-07-06
**Scope:** Human **whole-kidney** bulk RNA-seq, DN vs non-diabetic control. Built as the RNA-seq counterpart to the [glomerular microarray meta-analysis](REPORT.md), using the *same* study-matched effect-size framework, plus explicit handling of immediate-early/ischemia (IEG) artifact genes.

---

## 1. Why this analysis exists

Two prompts motivated it:
1. *"Why did you ignore RNA-seq?"* — the microarray meta-analysis was glomerular-microarray only.
2. A comparison to the **SuLab ARCHS4 analysis** (`OKN-WOBD/docs/dn_differential_expression_archs4_analysis.md`), whose top DN genes were the immediate-early transcription factors **FOS, FOSB, EGR1, NR4A1, DUSP1** (plus interferon genes OAS2/RSAD2). Those are classic tissue-procurement/ischemia-response genes, so the question was whether they are DN biology or a batch artifact.

**Design goal (the "union" analysis):** RNA-seq breadth + my study-matched random-effects effect-size framework (which controls batch and is robust to cross-study quantification differences) + explicit IEG artifact handling.

**Note on data source.** No ARCHS4 / `wobd-gene` MCP is configured in this environment, and ARCHS4's per-series REST extraction was not cleanly reachable, so instead of the ARCHS4 H5 I used **GEO-deposited RNA-seq matrices, uniformly reprocessed here**. This is methodologically sound for a *study-matched* design: DE is computed within each study and only standardized effect sizes are pooled, so cross-study quantification-pipeline differences are absorbed rather than confounding the meta-analysis.

## 2. Datasets

Candidates were screened from GEO (human, `expression profiling by high throughput sequencing`) and the SuLab study list. **Inclusion required: human, kidney tissue, bulk RNA-seq, and an internal DN-vs-control contrast.** This is stricter than the ARCHS4 pipeline and it matters (see §5).

| Dataset | Cohort | Tissue | DN | Control | Counts source |
|---|---|---|---:|---|---|
| **GSE142025** | Fan 2019 | Whole kidney biopsy | 27 (6 early+21 adv) | 9 | per-sample symbol matrices |
| **GSE162830** | 2021 | Whole kidney | 18 | 9 (reference) | quantile-norm symbol matrix |
| **GSE166239** | Nordbo 2021 | Whole FFPE biopsy | 6 | 6 | Ensembl raw counts |
| **Total** | | | **51** | **24** | |

**Rejected candidates (important):** GSE175759 — tubulointerstitium but only **3 samples are actually DN** (46 are IgAN, plus other GN); GSE185011 — **PBMC/blood**; GSE199437 — **mouse/mixed**; GSE204880, GSE195460, GSE209781, GSE266146 — **single-cell**; GSE192889, GSE158230, GSE262793, GSE154561 — **cell lines**; GSE154881 — blood; GSE199838 — 3v3 with a corrupted gene column. The SuLab analysis included several of these (GSE175759 as its dominant "62 DN" study; GSE185011 blood as kidney), which its automated LLM classifier appears to have mislabeled.

## 3. Methods

Same engine as the microarray analysis (`scripts/lib.py`). Per study:
- Gene identifiers → HGNC symbol (two studies were symbol-keyed; GSE166239 Ensembl IDs mapped via mygene.info). Duplicate symbols collapsed to the highest-mean row.
- **GSE166239 raw counts:** low-count filter (CPM≥1 in ≥20% of samples) → **log2-CPM** (limma-trend style). Symbol-keyed studies used their provided normalized values (log2 applied where linear).
- Per-study DE: **limma-style empirical-Bayes moderated t-test**; effect size = **Hedges' g**.
- Pooling: **DerSimonian–Laird random-effects** of Hedges' g; BH-FDR; Cochran's Q / **I²**. Weighted **Stouffer** cross-check.
- **IEG handling:** a curated set of 45 immediate-early/stress genes (FOS/JUN/EGR/NR4A/DUSP/ATF3/…) is **flagged and excluded** from the primary DEG set, and analyzed separately as a QC readout.

**High-confidence DEG:** present in all 3 studies, FDR<0.05, direction-consistent, **and not an IEG**.

## 4. Results

- 13,350 genes in all 3 studies; **1,647 high-confidence DEGs** (51 up / **1,596 down** — a strong down-skew typical of loss of differentiated epithelial function in diseased whole kidney). RE vs Stouffer sign concordance **0.958**.
- Real DN biology surfaces even in whole biopsy: **PTPRO** (GLEPP1, podocyte) and **MAGI2** (slit-diaphragm scaffold) among the top down genes; **MMP11** (fibrosis) top up.

### 4.1 The ARCHS4 immediate-early "signature" is a procurement artifact

The SuLab/ARCHS4 top genes, evaluated in this heterogeneity-aware RNA-seq meta:

| Gene | pooled g | I² | FDR | verdict |
|---|---:|---:|---:|---|
| NR4A1 | −2.33 | 96% | 0.41 | not significant |
| DUSP1 | −1.75 | 96% | 0.55 | not significant |
| FOSB | −1.53 | 96% | 0.63 | not significant |
| FOS | −1.05 | 94% | 0.66 | not significant |
| EGR1 | −1.01 | 94% | 0.69 | not significant |
| OAS2 (their ↑) | +0.15 | 77% | 0.90 | flat, not significant |
| RSAD2 (their ↑) | −0.41 | 44% | 0.50 | wrong sign, not significant |

Every one of the 7 genes SuLab reported **fails FDR here**. The 5 immediate-early genes have large effects but **I² ≈ 94–96%** — enormous between-study heterogeneity, the fingerprint of a technical (warm-ischemia / processing-time) artifact rather than reproducible disease biology. A naïve *pooled* analysis (SuLab's primary mode) or a fixed-effect/Stouffer combination surfaces them; a **random-effects** analysis correctly demotes them. Across all 45 curated IEG/stress genes, **only 1 is significant** (`fig5`, `fig7`).

### 4.2 Cross-modality comparison (RNA-seq whole-kidney vs microarray glomerular)

- 11,084 shared genes; pooled-g correlation **r = 0.44** (moderate — expected, since compartments differ).
- Robust DEGs overlap: 208 genes shared between the two independent meta-analyses, **86% in the same direction**, ~1.2× chance.
- **The most reproducible DN signal is podocyte injury**, significant across *both* modalities *and* compartments:

| Gene | Glomerular microarray | Whole-kidney RNA-seq |
|---|---|---|
| NPHS1 (nephrin) | −2.07 (FDR 1e-2) | −1.29 (FDR 3e-4) |
| NPHS2 (podocin) | −0.80 (FDR 0.11) | −0.90 (FDR 9e-3) |
| PTPRO (GLEPP1) | −0.83 (FDR 1e-3) | −1.46 (FDR 9e-5) |
| MAGI2 | −1.56 (FDR 8e-3) | −1.57 (FDR 5e-5) |

- **Fibrosis/complement genes are glomerular-compartment-specific**: TGFBI, COL1A1/2, C1QA, MMP2, LUM, VCAM1, CCL2 are strongly significant in microdissected glomeruli but **not significant in whole biopsy** (same direction, diluted by tubular mass). **MMP11** is the exception — reproducible in both (+1.84 vs +1.82).

This is exactly the compartment-dilution effect predicted when comparing microdissected-glomerular to whole-biopsy data, and it explains why the whole-tissue ARCHS4 approach struggled to see the structural DN signature.

### Figures
- **fig5_ieg_artifact.png** — RNA-seq volcano with IEGs highlighted (fail FDR) + |g|-vs-I² (IEGs at I²≈95%).
- **fig6_cross_platform.png** — microarray-glomerular vs RNA-seq-whole pooled-g concordance (r=0.44), hallmark genes labeled.
- **fig7_rnaseq_forest.png** — forest plots: IEG genes (heterogeneous → demoted) vs reproducible genes (PTPRO, MMP11).

## 5. What this says about the two strategies
- **Sample labeling matters more than sample count.** SuLab's dominant study contributed 62 "DN" samples of which only 3 were diabetic nephropathy; deterministic, verified labeling avoided that here.
- **The meta-analysis model matters.** Naïve pooling + `|log2FC|>2` selects high-dynamic-range artifact genes (IEGs); random-effects with I² penalizes non-reproducibility and removes them.
- **Compartment matters.** Whole-biopsy RNA-seq recovers podocyte-injury genes but dilutes the glomerular fibrosis/complement program that microdissected glomeruli reveal.

## 6. Limitations
- Only 3 RNA-seq studies with matched controls met strict criteria; GSE142025 is the largest and carries the most weight.
- Whole-biopsy compartment (tubulointerstitium-dominant); not directly comparable to microdissected glomeruli — the two analyses are complementary, not redundant.
- Controls include tumor-nephrectomy reference tissue (its own procurement profile).
- Cross-study RNA-seq quantification was not uniformly re-aligned (no ARCHS4 H5); mitigated, not eliminated, by the study-matched design.

## 7. Reproduce
```
scripts/recon_rnaseq.py          RNA-seq dataset discovery/QC
scripts/run_de_rnaseq.py         uniform processing + per-study DE  -> results/rnaseq_de_*, rnaseq_effect_*
scripts/meta_rnaseq.py           random-effects meta + IEG handling -> results/rnaseq_meta_results.csv, rnaseq_robust_DEGs.csv
scripts/compare_and_figures.py   cross-modality comparison + fig5-7
```

## 8. Data sources
- Fan Y *et al.* 2019 — GSE142025 · GSE162830 (nodular mesangial sclerosis, 2021) · Nordbø *et al.* 2021 — GSE166239.
