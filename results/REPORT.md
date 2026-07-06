# Pooled (meta-analysis) identification of differentially expressed genes in diabetic nephropathy

**Analysis date:** 2026-07-06
**Scope:** Human **glomerular** transcriptome, diabetic nephropathy / diabetic kidney disease (DN/DKD) vs. non-diabetic control, cross-platform microarray meta-analysis.

---

## 1. Objective

Identify genes robustly differentially expressed in the glomerular compartment of the kidney in diabetic nephropathy by pooling multiple independent public expression datasets, so that findings are not driven by any single cohort or platform.

## 2. Datasets

I searched GEO (NCBI eutils, `db=gds`) for human DN kidney expression studies and screened candidates by tissue compartment, platform, and group composition. To keep the pooled comparison biologically coherent I restricted to a **single compartment (glomeruli)** and to **microarray** platforms with normalized series-matrix data and a clean disease-vs-control contrast. Three independent studies on three different Affymetrix platforms met these criteria:

| Dataset | Cohort (ref.) | Compartment | DN | Control (type) | Platform | Genes tested |
|---|---|---|---:|---|---|---:|
| **GSE30528** | Woroniecka *et al.* 2011 | Glomeruli | 9 | 13 (nephrectomy/donor) | GPL571 · HG-U133A_2 | 13,237 |
| **GSE96804** | Pan *et al.* 2018 | Glomeruli | 41 | 20 (unaffected tumor-nephrectomy) | GPL17586 · HTA-2.0 | 30,905 |
| **GSE104948** (GPL22945 subseries) | ERCB, Ju/Grayson *et al.* | Glomeruli | 7 | 18 (living donor) | GPL22945 · HG-U133 Plus 2 ENTREZG CDF | 12,074 |
| **Total** | | | **57** | **51** | 3 platforms | |

Rejected/observed candidates: GSE30529 & GSE47184 (tubulointerstitium — different compartment); GSE47183 glomerular (only tumor-nephrectomy "controls", no healthy donors); GSE142025 (bulk-biopsy RNA-seq — different modality/compartment mix); GSE1009 (3 vs 3, obsolete HuGeneFL platform — underpowered). The chosen three give three biological cohorts, three platforms, and both nephrectomy and living-donor controls.

## 3. Methods

All analysis is in Python (numpy/scipy/pandas); no R/limma dependency. Statistics were re-implemented from the primary literature so every step is transparent (`scripts/lib.py`).

**Data processing (per study)**
1. Downloaded normalized `series_matrix` from GEO; parsed expression + sample metadata directly.
2. Assigned each sample DN vs control from `Sample_source_name`; samples of any other diagnosis (e.g. FSGS, MCD, RPGN in the ERCB series) were **excluded**.
3. Log2 scale verified/enforced (all three were already log2).
4. Probe → **HGNC gene symbol** mapping per platform (GPL571 `Gene Symbol`; HTA-2.0 `gene_assignment` parsed; ENTREZG CDF `Symbol`). Multiple probes per gene collapsed to the **highest-mean probe** (max-mean collapse). Genes with missing values dropped.

**Per-study differential expression** — limma-style **empirical-Bayes moderated t-test** (Smyth 2004): per-gene variances shrunk toward a fitted scaled inverse-χ² prior (`fitFDist` via trigamma inversion), moderated *t*, and Benjamini–Hochberg FDR. Positive log2FC = up in DN. (`results/de_<GSE>.csv`)

**Effect size (per study)** — **Hedges' g** standardized mean difference with small-sample correction *J* and its variance. (`results/effect_<GSE>.csv`)

**Pooling (primary)** — per gene, **DerSimonian–Laird random-effects** meta-analysis of Hedges' g across studies: pooled g, standard error, z, p, between-study variance τ², Cochran's Q and **I² heterogeneity**. Random effects (not fixed) was chosen because studies differ in platform and control type. FDR by BH over all genes present in ≥2 studies. (`results/meta_results.csv`)

**Pooling (secondary / cross-check)** — sample-size-weighted **Stouffer** combination of the per-study moderated-t p-values (signed by fold-change direction). Used only to confirm the random-effects result.

**High-confidence DEG definition:** present in **all 3 studies**, meta **FDR < 0.05**, and **same direction in all three** studies.

## 4. Results

- **11,472** genes were measured in all three studies; **12,792** in ≥2.
- **1,714 high-confidence DEGs** (all-3, FDR<0.05, direction-consistent): **1,022 up** and **692 down** in DN. (`results/robust_meta_DEGs.csv`)
- The two independent pooling methods agree: **91.5%** sign concordance between random-effects pooled g and Stouffer z. Pairwise per-gene log2FC correlations between studies are all positive (see `fig4`).

**Biological validation (direction is correct):** the single most down-regulated hallmark, **NPHS1 (nephrin)** — the podocyte slit-diaphragm gene lost in DN — is down in all three cohorts, alongside podocyte/junction genes **TJP1 (ZO-1), MPP5, GJA1**. Up-regulated genes recapitulate the known DN glomerulosclerosis + inflammation program: ECM/fibrosis (**TGFBI, COL1A1, COL1A2, COL16A1, LUM, LAMC3, MFAP2, LTBP1, MMP2, PCOLCE, FMOD**), complement/macrophage/immune (**C1QA, VSIG4, MS4A4A, MS4A6A, CLEC10A, STAB1, ZAP70**), and growth-factor/Wnt signaling (**PDGFB, WNT2B, DACT1**). Down-regulated metabolic genes (**LPL, ALB, CYP3A7, G6PC**) fit metabolic dysregulation of the diabetic glomerulus.

### Top 20 UP in DN (ranked by meta p-value; all FDR≈0, I²≈0 → most reproducible)

| Gene | pooled g | FDR | I²% | g GSE30528 | g GSE96804 | g GSE104948 |
|---|---:|---:|---:|---:|---:|---:|
| COL16A1 | 1.88 | 0.0 | 0 | 1.73 | 1.74 | 2.52 |
| ISLR | 1.73 | 0.0 | 0 | 1.86 | 1.61 | 1.94 |
| CLEC2D | 1.69 | 0.0 | 0 | 2.01 | 1.60 | 1.65 |
| OSBPL3 | 1.60 | 0.0 | 0 | 1.46 | 1.49 | 2.08 |
| TMEM243 | 1.59 | 0.0 | 0 | 1.43 | 1.69 | 1.52 |
| LTBP1 | 1.59 | 0.0 | 0 | 1.19 | 1.66 | 1.89 |
| FADS2 | 1.58 | 0.0 | 0 | 1.08 | 1.69 | 1.96 |
| DACT1 | 1.56 | 0.0 | 0 | 1.46 | 1.39 | 2.24 |
| LUM | 1.56 | 0.0 | 0 | 1.97 | 1.53 | 1.29 |
| PDGFB | 1.52 | 0.0 | 0 | 1.21 | 1.53 | 1.85 |
| PTEN | 1.47 | 0.0 | 0 | 1.40 | 1.32 | 2.04 |
| WNT2B | 1.46 | 0.0 | 0 | 1.14 | 1.45 | 1.93 |
| LAMC3 | 1.47 | 0.0 | 0 | 1.63 | 1.24 | 2.01 |
| MFAP2 | 1.43 | 0.0 | 0 | 1.54 | 1.33 | 1.61 |

### Top 20 DOWN in DN (ranked by meta p-value)

| Gene | pooled g | FDR | I²% | g GSE30528 | g GSE96804 | g GSE104948 |
|---|---:|---:|---:|---:|---:|---:|
| HSPA1L | -2.01 | 0.0 | 0 | -1.74 | -2.11 | -2.06 |
| ERVMER34-1 | -1.95 | 0.0 | 0 | -1.85 | -2.18 | -1.56 |
| LPL | -1.91 | 0.0 | 0 | -1.95 | -1.83 | -2.10 |
| ZBTB10 | -1.83 | 0.0 | 0 | -2.18 | -1.69 | -1.89 |
| CCDC6 | -1.74 | 0.0 | 0 | -1.49 | -1.73 | -2.06 |
| MMP28 | -1.73 | 0.0 | 0 | -2.20 | -1.60 | -1.68 |
| SUCO | -1.70 | 0.0 | 0 | -1.57 | -1.79 | -1.62 |
| CYP3A7 | -1.60 | 0.0 | 0 | -1.34 | -1.60 | -1.94 |
| SPTB | -1.60 | 0.0 | 0 | -1.68 | -1.38 | -2.19 |
| ALB | -1.58 | 0.0 | 0 | -1.12 | -1.70 | -1.86 |
| NLK | -1.53 | 0.0 | 0 | -2.05 | -1.43 | -1.34 |
| FZD2 | -1.53 | 0.0 | 0 | -1.87 | -1.63 | -1.05 |

*Hallmark podocyte gene **NPHS1** (pooled g = -2.07, FDR = 0.012) and fibrosis driver **TGFBI** (pooled g = +2.62, FDR = 0.022) are direction-consistent across all three studies; their FDR is higher only because their effect sizes are more heterogeneous between cohorts (see forest plot).*

### Figures (`results/`)
- **fig1_volcano.png** — pooled effect size vs −log10 FDR; up (red) / down (blue).
- **fig2_forest.png** — per-study g ±95% CI + random-effects diamond for hallmark genes (NPHS1, TGFBI, C1QA, COL1A1, MMP2, EZH2, TJP1, TSPYL5).
- **fig3_topgenes.png** — top 20 up / down genes by meta p-value.
- **fig4_concordance.png** — cross-study log2FC concordance (all pairwise r > 0).

## 5. Limitations
- **Compartment-specific:** glomerular only; tubulointerstitial DN (e.g. GSE30529/GSE47184) would give a different, complementary signature.
- **Controls are heterogeneous** (living donors vs unaffected tumor-nephrectomy tissue), a recognized issue in kidney transcriptomics.
- **Between-study heterogeneity** is high for the largest-effect genes (I² up to ~90%), driven mainly by the small custom-CDF ERCB set showing larger standardized effects; random-effects pooling handles this conservatively (wider CIs, larger p-values), so the largest-|g| genes are *not* necessarily the most significant. The most reproducible genes are the low-I² ones in the tables above.
- Microarray probe→gene collapse (max-mean) and single-symbol assignment for multi-mapping probes are standard but lossy choices.
- No independent RNA-seq validation cohort included; effect sizes are relative, not absolute.

## 6. Reproducibility
```
scripts/recon.py     GEO reconnaissance (dataset discovery/QC)
scripts/lib.py       parsing, probe→gene, moderated t-test, Hedges' g, DerSimonian-Laird, BH-FDR
scripts/run_de.py    Stage 1: download + per-study DE + effect sizes  -> results/de_*, effect_*
scripts/meta.py      Stage 2: random-effects meta + Stouffer cross-check -> results/meta_results.csv
scripts/figures.py   Stage 3: figures + significance-ranked tables
```
Key outputs: `results/meta_results.csv` (all genes), `results/robust_meta_DEGs.csv` (1,714 high-confidence DEGs), `results/top_up.csv` / `top_down.csv`, `results/dataset_summary.csv`.

## 7. References (data)
- Woroniecka KI *et al.* Transcriptome analysis of human diabetic kidney disease. *Diabetes* 2011 — GSE30528.
- Pan Y *et al.* Exon-level expression profiling of the diabetic nephropathy glomerulus, 2018 — GSE96804.
- European Renal cDNA Bank (ERCB); Ju W / Grayson PC *et al.* — GSE104948 (GPL22945 glomerular subseries).
