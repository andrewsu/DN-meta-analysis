# Data-discovery methods compared — Arthritis (generalization test #3, broad umbrella)

**Task:** *"Find datasets to run a pooled analysis to find differentially expressed genes in arthritis."*
**Date:** 2026-07-08 · **Author:** Andrew Su (with Claude Code)
**Artifacts:** `arthritis-comparison/` · **Harness:** `run_dn_condition_comparison.sh` (`DISEASE="arthritis"`)
**Companions:** [DN](data-discovery-method-comparison.md) · [glaucoma](glaucoma-method-comparison.md) · [SSc](ssc-method-comparison.md)

> **Why this disease:** "arthritic joint disease" (MONDO:0005578) is the **broad-umbrella extreme** — the opposite of DN (a MONDO leaf). It has **53 MONDO descendants** and its subclasses are individually huge (RA 2,006, OA 424, PsA 110, JIA 82 datasets in NDE). It stress-tests what the methods and subclass expansion do when the query term is *too general*.

---

## TL;DR

- **All three methods independently recognized the task is ill-posed** at this scope — "arthritis" is not one disease (pooling RA + OA + JIA + septic arthritis is biologically incoherent) — and each responded characteristically rather than producing a bogus pooled set.
- **Subclass expansion explodes recall** here: exact "arthritic joint disease" **603 → 2,737 subclass-aware** (all species; **603 → 1,559** human GEO). But **valid *poolable* recall is undefined** for the umbrella — you must first scope to one subtype.
- **A (alone)** gave the sharpest conceptual pushback + 16 leads but no verified data; **B (web)** scoped to a **verified 9-dataset RA/OA/JIA core (423 disease / 141 control, 0 FP)**; **D (NDE)** delivered a **discovery inventory (1,559 human GEO; 451 arthritis-only; 60 superseries)** and refused sample-level claims.

---

## Token cost & effort

| Condition | Total tokens | Turns | Cost |
|---|--:|--:|--:|
| **A · alone** | 28,388 | 1 | $0.35 |
| **B · web** | 286,546 | 17 | $1.14 |
| **D · MCP→NDE** | 288,808 | 17 | $1.16 |

(B did *not* time out this time — it scoped rather than attempting to screen all 868 candidates.)

---

## Primary result — how each method handled the umbrella

- **A · alone — conceptual pushback + leads.** Opened with *"'Arthritis' is not one disease group… biologically incoherent unless that's an explicit design choice,"* refused to invent counts, listed **16 candidate GSEs as "leads only, counts unverified,"** and specified the GEO+GEOmetadb pipeline to do it properly. Best methodological framing; zero verified data (no tools).
- **B · web — verified curated core.** A GEO query returned **868 human arthritis expression series**; B did not screen all 868 but built a **verified 9-dataset core** from the RA/OA/JIA literature, each checked against its GEO record:
  - Synovium: GSE1919, GSE55235, GSE55457, GSE55584, GSE77298, GSE82107, GSE89408
  - Blood/PBMC: GSE15573, GSE13501
  - **Totals: 423 disease / 141 control** (RA 224, OA 63, JIA 136); adult RA/OA-only = 287/82.
  - Excluded a drug-response longitudinal set (GSE93272) and a dual-platform set (GSE12021); correctly ruled GSE55235/57/84 a 3-center study (not a duplicating superseries). **0 false positives.** Flagged JIA as a scope judgment call.
- **D · MCP→NDE — discovery inventory.** Used subclass expansion (label contains "arthritis") to enumerate **1,559 human GEO arthritis datasets** (RA 1,221, OA 272, PsA 76, JIA 70). It **refused sample-level claims** ("no sample records… control is a sample attribute, invisible here"), and surfaced two scoping facts web didn't quantify: **only 451 of 1,559 are arthritis-*only*** (1,108 are multi-disease), and **60 contain "SuperSeries"** (double-count hazard). Offered the accession list at any chosen scope.

No Venn / union table: the task is ill-posed at the umbrella and D returned an inventory rather than a curated set.

---

## Subclass-aware recall (the whole story here)

| | exact "arthritic joint disease" | subclass-aware |
|---|--:|--:|
| All species | 603 | **2,737** (+354%) |
| Human GEO (D's count) | — | **1,559** |

Subclass expansion is **essential** for discovery at an umbrella term (without it you'd retrieve 603 and miss RA's 2,006, etc.) — but it inflates a pool you then can't pool. The "+354%" recall buys candidates across incoherent diseases; **valid poolable recall requires first scoping to a subtype**, which turns arthritis back into a tractable ~130–170-scale problem (like the other three diseases).

---

## Generalization across all four diseases — the subclass-expansion spectrum

| Disease | MONDO structure | exact | subclass-aware (raw) | **valid gain from subclasses** |
|---|---|--:|--:|--:|
| **DN** | leaf (0 subclasses) | 129 | 129 (0%) | 0 |
| **SSc** | 7 subclasses | 167 | 170 (+1.8%) | **1** (GSE9285, a real cohort) |
| **Glaucoma** | 50 subclasses | 134 | 159 (+19%) | 0 (all in-vitro/treatment/mis-tag) |
| **Arthritis** | 53 subclasses, umbrella | 603 | 2,737 (+354%) | undefined — task ill-posed |

**The law that emerges:** *raw* recall gain from subclass expansion scales with how bushy the ontology subtree is (leaf → umbrella spans 0% → +354%), but **valid gain is decoupled from it** — it depends on whether the subtypes host genuine patient case-control cohorts (SSc: one does; glaucoma: none; arthritis: yes, but across incoherent diseases you shouldn't pool). Subclass expansion maximizes recall, not usefulness.

---

## Bottom line

At a broad umbrella term, the "find datasets for a pooled DE analysis" task is **ill-posed**, and the methods' value is in *flagging* that rather than answering it: **A** gives the crispest conceptual warning, **B** the most useful pragmatic move (scope to a verified subtype core), **D** the best landscape inventory (with the multi-disease/superseries hazards quantified). Subclass expansion is required for discovery here but produces a pool spanning distinct diseases; the correct next step for any of the methods is to **scope to a specific arthritis subtype** and then run the DN/glaucoma/SSc-style pipeline (discover → programmatic organism+assay filter → sample-level curation from GEO).
