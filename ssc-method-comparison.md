# Data-discovery methods compared — Systemic sclerosis (generalization test #2, on-mission)

**Task:** *"Find datasets to run a pooled analysis to find differentially expressed genes in systemic sclerosis (scleroderma)."*
**Date:** 2026-07-07 · **Author:** Andrew Su (with Claude Code)
**Artifacts:** `ssc-comparison/` · **Harness:** `run_dn_condition_comparison.sh` (`DISEASE="systemic sclerosis (scleroderma)"`)
**Companions:** [data-discovery-method-comparison.md](data-discovery-method-comparison.md) (DN) · [glaucoma-method-comparison.md](glaucoma-method-comparison.md)

> **Why this disease:** systemic sclerosis is squarely **on NDE's NIAID mission** (immune-mediated / autoimmune) — the test is whether NDE performs *better* in its focus area than off-mission (DN, glaucoma). It also has a larger footprint: **167** NDE `healthCondition` datasets (146 human), ~30% above the DN/glaucoma benchmark.

---

## TL;DR

- **On-mission did NOT make NDE cleaner — it made the problem bigger.** NDE returned a **174-record raw pool it could not filter**; a primary-record scan shows **~27% are false positives** (non-human or non-gene-expression) before any sample-level work.
- **Web (B) handled the large disease well** — a rigorous, verified 16-dataset curation with **zero false positives** — but only on a re-run: its **first attempt timed out** (900s), so web's runtime is **high-variance at large N**, not a categorical wall.
- **Web has a freshness edge:** it found **four 2025 datasets that are absent from NDE's mirror** (confirmed) — NDE lags the primary archive.
- Both prior false-positive mechanisms **recur and compound** in NDE: the DN species-mislabel (mouse studies tagged "Homo sapiens") *and* the glaucoma wrong-modality (non-coding RNA as gene expression).
- **A (alone) recall recovered to 8** and is a strict **subset of B** — SSc is well-catalogued, unlike glaucoma.

---

## Conditions & token cost

| Condition | Result | Tokens | Turns | Cost |
|---|---|--:|--:|--:|
| **A · alone** | ✅ 8 datasets (all ⊂ B) | 27,512 | 1 | $0.33 |
| **B · web** | ✅ 16 datasets, 0 FP — **on re-run**; 1st attempt timed out at 900s | ~(re-run) | 29 | $1.87 |
| **D · MCP→NDE** | ✅ but a **174-record raw pool, not a curated set** | 760,230 | 21 | $2.25 |

B's failed first attempt took a 63-WebFetch path and blew the 900s cap; the successful re-run took a leaner 26-fetch path and finished in 29 turns. Same task, very different runtime — **high variance**, not an inability to scale.

---

## Primary result — finding the relevant datasets

- **A (alone): 8 datasets**, all of which B also found (**A ⊂ B**). Recall recovered vs glaucoma's 1 — SSc has a well-known DE literature.
- **B (web): 16 datasets** (14 bulk case-control + 2 single-cell), verified, **0 false positives**. It excluded morphea, IPF, IPAH, drug-treatment timepoints, cultured-fibroblast arms, and a superseries — and **self-corrected a mis-recalled accession** (GSE181818, "not SSc"). Totals: **≈366 disease / ≈175 control** (verified core) to **≈441 / ≈225** (full curated).
- **D (MCP→NDE): a 174-record candidate pool, unfiltered.** Steer held (12 `nde` queries), but D explicitly declined to curate: *"nde cannot answer per-sample N, tissue, assay, or the sample-level exclusions… no sample-level records."*

### B's included bulk case-control set (verified)

| Accession | Tissue | Assay | N dis / ctrl |
|---|---|---|---|
| [GSE9285](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE9285) | Skin | 2-color array | 24 / 6 |
| [GSE32413](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32413) | Skin | array | 22 / 9 |
| [GSE58095](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58095) | Skin | array | 66 / 36 |
| [GSE76807](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE76807) | Skin | array | 10 / 5 |
| [GSE95065](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95065) | Skin | array | ~17 / ~15 |
| [GSE130955](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130955) | Skin | RNA-seq | 48 / 33 |
| [GSE181549](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181549) | Skin | array | 113 / 44 |
| [GSE45485](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45485) | Skin | array | MMF trial (subseries of GSE59787) |
| [GSE59785](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE59785) | Skin | array | MMF trial (subseries of GSE59787) |
| [GSE308096](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE308096) 🆕 | (blood/skin) | RNA-seq + miRNA | 5 / 5 |
| [GSE19617](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE19617) | Blood (PBMC) | array | 36 / 10 |
| [GSE33463](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33463) | Blood (PBMC) | array | 69 / 41 |
| [GSE48149](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE48149) | Lung | array | ~18 / ~9 |
| [GSE317056](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE317056) 🆕 | Lung | RNA-seq | ~12 / ~12 |

Plus single-cell (separate modality): [GSE320020](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE320020) 🆕, [GSE334710](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE334710) 🆕. 🆕 = 2025 dataset **absent from NDE's mirror**.

---

## Method overlap

No Venn is drawn: D returned an unfiltered 174-record pool rather than a curated set, so a three-way overlap of *candidate sets* would be misleading. The relationships:

- **A ⊂ B** — all 8 of alone's datasets are in web's curated set.
- **B ∩ D-pool = 12** of B's 16 — the shared, older datasets.
- **4 of B's datasets are absent from NDE** (GSE308096, GSE317056, GSE320020, GSE334710 — all 2025): NDE's GEO mirror lags. Confirmed by direct NDE query.
- **D's pool holds ~160 datasets beyond B's picks** — mostly unfiltered material (see below), not additional verified candidates.

---

## GEO false-positive verification — D's raw pool (192 records)

Every accession D surfaced was downloaded from GEO and classified by organism + assay modality:

| Category | Count | Verdict |
|---|--:|---|
| Human + expression-profiling | **140 (73%)** | eligible — still need sample-level triage D couldn't do |
| Non-human-only (mouse/rat models) | **31** | ❌ FP — NDE mislabels ~half as "Homo sapiens" |
| Human, non-coding-RNA only (miRNA/lncRNA) | **4** | ❌ FP for a gene/mRNA DE analysis |
| Human, other-assay only (methylation/ATAC/genotyping) | **17** | ❌ FP |

**≈52 of 192 (27%) are false positives at the pool level.**

**Species-mislabel confirmed at scale (the DN/GSE33744 mechanism):** of 7 sampled non-human-only series, **4 carry a spurious "Homo sapiens" tag in NDE** (GSE100212, GSE226375, GSE61728, GSE71991 — all `Mus musculus` in GEO). This explains why D reported excluding only "17 mouse/rat" while GEO shows **31** — NDE's text-mined species tags let ~half the mouse studies through a human filter.

**Web (B) introduced zero false positives** even at this scale — it verified each candidate against the primary GEO record itself.

---

## Stretch result — per-sample counts

- **B (web):** ≈366 disease / ≈175 control (verified core), with explicit caveats about **patient overlap** across the Whitfield/Dartmouth-Northwestern skin cohorts and the Hopkins PAH blood cohorts (unique-patient count materially lower than the arithmetic sum; not resolvable from metadata).
- **D (NDE):** cannot — no sample records.
- **A (alone):** memory only.

---

## Generalization — three diseases

| | DN (off-mission) | Glaucoma (off-mission) | **SSc (on-mission)** |
|---|--|--|--|
| NDE `healthCondition` count | 129 | 134 | **167** |
| A (alone) recall | 11 | 1 | 8 (⊂ B) |
| B (web) included / FPs | 9 / 0 | 6 / 0 | **16 / 0** (but 1st run timed out) |
| D (NDE) included / FPs | 21 / 2 | 9 / 3 | **174 raw pool; ~52 (27%) pool-level FP; couldn't filter** |
| FP mechanism(s) | species-mislabel, cell-line | wrong-modality (miRNA) | **both, compounded** |

**Held across all three:** every false positive traces to NDE; web introduces none; NDE's FP mechanism is metadata/description-driven; NDE has no sample-level data; primary-record verification is mandatory.

**On-mission result (vs the "NDE should be better here" prior):** being in NDE's focus area bought **breadth, not precision** — more coverage, a *larger* absolute FP burden, and (at ~170 datasets) no filtered list at all. NDE also **lagged the archive** (missed 4 recent datasets web found).

**New regime finding:** the best method is **disease-size-dependent**. Web wins at small N (DN, glaucoma) *and* stays clean at large N (SSc) — but its runtime is high-variance and it can time out. NDE gives breadth but can't filter and lags on freshness. At large N the effective add-on is a **programmatic GEO-metadata filter** (organism + assay type) — the ~2-minute step that cut D's 192-pool to 140 human-expression candidates.

---

## Bottom line

Across DN, glaucoma, and SSc the *method-level* findings are robust: **web produces clean, verified, current candidate sets (zero false positives in all three)**; **NDE offers broad discovery but description-driven false positives, no sample counts, and a lagging mirror**; **primary-record verification is non-negotiable**. On-mission focus (SSc) gave NDE more coverage but not more precision, and at ~170 datasets the practical pipeline is **web or a GEO query for candidates → programmatic organism+assay filter → sample-level curation from GEO SOFT files**, with patient-overlap de-duplication that no conversational method resolves from metadata alone.

---

## Addendum — subclass-aware recall (ontology expansion)

Ontology subclass expansion is an NDE/MCP-only capability (`get_descendants`), and the SSc **D run actually used it** — the stream shows it looked up SSc's URI, enumerated subtypes, and excluded morphea. Expanding systemic sclerosis (MONDO:0005100) to its **7 MONDO subclasses** (diffuse/limited cutaneous SSc, diffuse scleroderma, CREST, …) changes the NDE count only **167 → 170 (+1.8%)** — SSc studies almost always carry the umbrella "systemic sclerosis" tag, so there is little subclass-only recall to recover. This is the low-gain end of a strongly disease-dependent effect (glaucoma +19%, DN 0%).
