# Data-discovery methods compared: Claude alone vs. web vs. OKN knowledge graphs

**Task:** *"Find datasets to run a pooled analysis to find differentially expressed genes in diabetic nephropathy."*
**Date:** 2026-07-07 · **Author:** Andrew Su (with Claude Code)
**Artifacts:** `dn-comparison-run1/` · **Harness:** `run_dn_condition_comparison.sh` · **Registry fix:** [frink-okn/okn-registry#441](https://github.com/frink-okn/okn-registry/pull/441)

---

## TL;DR

- We ran the identical dataset-finding task under four isolated conditions and measured token cost, dataset recall, and whether each could deliver **verified per-sample disease/control counts** under strict inclusion criteria.
- **Claude + web search (D)** produced the best answer for this task — a verified, de-duplicated dataset table with **138 DN / 114 control** samples segregated by compartment — at **$1.19**, the second-cheapest run.
- **Claude + OKN MCP** was the most expensive and, *left to route itself, answered the least*: it selected the wrong graph (Gene Expression Atlas, a 243-study curated subset) and returned **1 dataset** for **$1.97**.
- **Explicitly steering the MCP to the NDE graph** fixed discovery dramatically (~20 human datasets, correct exclusions, double-count catches) for **$2.27** — but NDE **structurally cannot** return per-sample counts (no `Sample`/`CellType` entity), so it stopped at *47/53 verifiable + 12 datasets with counts absent* and told the user to finish at GEO.
- Root cause of the mis-routing is **registry-description scope**, not model reasoning — addressed in PR #441.

---

## 1. Conditions

| ID | Condition | Tools available | Isolation |
|----|-----------|-----------------|-----------|
| **A** | Claude **alone** | none (no web, no MCP, no file/exec, no ToolSearch) | fresh headless session, clean CWD |
| **D** | Claude **+ web** | `WebSearch`, `WebFetch` | fresh headless session, clean CWD |
| **C (unsteered)** | Claude **+ OKN MCP** | `mcp__okn-mcp__*` (all 35 graphs), self-routed | fresh headless session, clean CWD |
| **C (steered)** | Claude **+ OKN MCP → NDE** | same, but system prompt names the `nde` graph | fresh headless session, clean CWD |

**Why separate sessions:** to guarantee **zero cross-condition leakage**, each condition is a distinct `claude -p` invocation (no `--continue`/`--resume`), launched from a clean `mktemp -d` working directory so this repo's `CLAUDE.md`/auto-memory (which already contains DN findings) is not loaded as a hint. Tool boundaries are enforced with `--allowedTools`/`--disallowedTools` + `--strict-mcp-config` — never `--dangerously-skip-permissions`. Model pinned to Opus 4.8 for all runs.

---

## 2. Token cost & effort

Summed across **all** models per run (top-level `.usage` under-counts and hides a secondary Haiku model — see §7).

| Condition | Input | Output | Cache-create | Cache-read | **Total tokens** | Turns | **Cost (USD)** |
|---|--:|--:|--:|--:|--:|--:|--:|
| **A · alone** | 3,429 | 7,445 | 19,369 | 0 | **30,243** | 1 | **$0.39** |
| **D · web** | 127,931 | 31,485 | 32,950 | 125,302 | **317,668** | 21 | **$1.19** |
| **C · unsteered → GxA** | 5,014 | 23,297 | 92,813 | 871,538 | **992,662** | 35 | **$1.97** |
| **C · steered → NDE** | 4,756 | 42,698 | 89,601 | 566,634 | **703,689** | 16 | **$2.27** |

Note the **cost/quality inversion**: the two MCP runs are the most expensive yet (unsteered) the least useful. MCP cost is dominated by **cache-read of large tool payloads** (`list_graphs`/`route_query` ≈90 KB each, NDE query results 10–17 KB each) re-read across many turns.

---

## 3. What each condition delivered

| Condition | Datasets found | Sample counts (strict) | Verifiable? |
|---|--:|---|---|
| **A · alone** | ~11 (from memory) | refused firm total; soft ~70–90 DN / ~55–60 ctrl | ✗ self-flagged unverifiable |
| **D · web** | 9 (curated) | **138 DN / 114 control** (verified floor 116/70) | ✅ from GEO records |
| **C · unsteered → GxA** | **1** (GSE1009) | none (GxA stores no sample counts) | n/a |
| **C · steered → NDE** | ~20 human | **47 DN / 53 control verifiable** + 12 datasets "N not reported" | ⚠️ partial (structural limit) |

**A (alone)** — Opened honestly ("no tools, these are training-data recollections, verify before use"), listed the canonical series (GSE30528/30529/30122/96804/1009/131882/142153 + ERCB), and correctly flagged the GSE30122 superseries and ERCB patient-reuse traps. Accurate recall, but cutoff-bound and unverifiable by itself.

**D (web)** — Best answer. Verified, de-duplicated table with compartment segregation: Glomerular 65/57, Tubulointerstitial 20/35, Whole-kidney bulk (GSE142025) 27/9, snRNA-seq (GSE131882) 3/3 subjects, PBMC (GSE142153) 23/10. Correctly excluded animal/in-vitro/treatment arms and resolved ERCB reprocessing overlaps. Uniquely surfaced **prior published pooled analyses** (panels + batch methods).

**C (unsteered → GxA)** — Self-routed to `gene-expression-atlas-okn` because the query says "differentially expressed genes." That graph is a **243-study curated subset**, so it found only **GSE1009** and (correctly) concluded "not enough for a meaningful pooled analysis," pointing the user elsewhere. Never queried NDE.

**C (steered → NDE)** — With the graph named, ran 12 SPARQL queries all on `nde`, resolved DKD to `MONDO_0005016`, screened 95 human-tagged datasets, and produced an excellent triage: ~20 human DN-vs-control expression series, correct exclusion of mouse/in-vitro/treatment/other-disease records, and real curation smarts (caught **GSE195799** = OVE26 mice mislabeled human; **GSE195460** reusing **GSE131882**'s snRNA samples; the GSE30122 superseries). But it hit NDE's ceiling on counts (§5).

---

## 4. Finding 1 — unsteered MCP routes to the wrong graph

Given the whole OKN server, Claude picked **Gene Expression Atlas over NDE** for a *dataset-discovery* query. Two causes, both in the registry metadata, not the model:

1. **Lexical/semantic pull** — the query's "differentially expressed genes / pooled analysis" matches GxA's name and its description ("differential gene expression… log2 fold changes… **cross-study meta-analyses**") far better than NDE's "NIAID Data Ecosystem."
2. **NDE self-excludes** — its description scopes it to "infectious and immune-mediated disease," so a *diabetic nephropathy* (non-infectious) query reads as out of scope — even though NDE actually holds **128+ DKD datasets** (verified: 129 under `healthCondition = MONDO_0005016`).

The built-in `route_query` keyword matcher was **no help** — for this query it ranked non-biomedical graphs (soil, aging, hydrology) on top and surfaced neither biomedical graph. The model's own reading of `list_graphs`/`get_description` summaries is what drove selection, which is why **description text is the fix**.

---

## 5. Finding 2 — steering fixes discovery, but NDE has a structural ceiling

Naming the graph flipped C from 1 dataset to a full discovery+triage run. But NDE is a **dataset-metadata index** (Schema.org: `Dataset`, `DefinedTerm`, `species`, `healthCondition`, `includedInDataCatalog`…). It has:

- **no `Sample` entity** → cannot count DN vs control samples,
- **no `CellType` entity** → cannot segregate by cell type,
- **empty `abstract`** and **no `measurementTechnique`/sample-count predicate** → per-group N is only knowable when a free-text `description` happens to state it.

So even steered, NDE could verify only **47 DN / 53 control** (8 datasets that state N in text) and had to mark 12 qualifying datasets "NR" and defer to GEO for the rest. **This is the discovery-vs-granularity boundary:** NDE is excellent for *which datasets exist*; the *sample-level curation* still requires GEO — exactly the step condition D performed.

---

## 6. Cross-validation across methods

The methods **corroborate each other**, which raises confidence:

- **GSE1009** — the only dataset the GxA graph contained — is confirmed real and appears in A, D, and NDE.
- Every flagship accession A produced from memory (GSE30528, GSE30529, GSE30122, GSE96804, GSE131882, GSE99339, GSE142153) is independently confirmed by web (D) and/or NDE.
- D's ~138/114 (headless) is consistent with the earlier interactive curation (~142/109) to within a few ERCB samples.

Consensus **human DN-vs-control backbone** for the pool: **glomerular** GSE30528, GSE96804, GSE1009, GSE104948; **tubulointerstitial** GSE30529, GSE104954; **whole-kidney bulk** GSE142025; **snRNA-seq** GSE131882; **PBMC** GSE142153 — with GSE30122/GSE47183/GSE99339/GSE99325 excluded as superseries/ERCB duplicates.

---

## 7. What each resource is actually for

| Resource | Best question | Returns | Not for |
|---|---|---|---|
| **Claude alone** | rapid scoping / "what's the landscape?" | recalled candidate list + design traps | anything requiring verification or post-cutoff data |
| **Claude + web** | "find datasets **and** verify sample counts" | verified table, counts, precedent, methods | structured/queryable reuse |
| **NDE graph** | "**discover** which datasets exist for disease X" | queryable dataset metadata + GEO links | per-sample counts, cell-type splits, expression values |
| **Gene Expression Atlas graph** | "give me the **computed DE values** for disease/tissue X" | genes + log2FC + adjusted p for curated studies | comprehensive dataset discovery (it's a 243-study subset) |

The last two are **complementary, not competing** — NDE finds the datasets; GxA gives DE values for its curated subset. The mis-routing happened because both descriptions blurred this line.

---

## 8. The registry fix (PR #441)

[frink-okn/okn-registry#441](https://github.com/frink-okn/okn-registry/pull/441) — *"Better define the scope of the nde and Gene Expression Atlas graphs"* (draft) — sharpens both descriptions so routing is correct from both sides:

- **`nde`**: lead with dataset discovery + gene-expression/transcriptomics/RNA-seq + pooled/meta-analysis; state coverage is broad (deepest for IID but mirrors GEO/70+ repos); add a guardrail that it holds metadata/links, not per-sample or DE values; add `keywords`/`example_queries`.
- **`gene-expression-atlas-okn`**: lead with its real purpose (retrieve computed DE values for a curated subset); add a "subset, not comprehensive" caveat; reframe "cross-study meta-analyses"; point dataset-discovery queries to NDE.

**Validation after deploy:** re-run `route_query` on this query (NDE should rank above GxA) and confirm an **unsteered** agent now selects NDE for discovery without a hint.

---

## 9. Operational lessons (reproducibility)

Non-obvious issues surfaced while building the harness — worth recording:

1. **Token accounting** — sum `.modelUsage` across models, not top-level `.usage` (which zeroed `cache_creation` in one test while `modelUsage` showed 5,975). Runs also silently invoke a secondary **Haiku** model; only the `modelUsage` sum captures it.
2. **Streaming** — use `--output-format stream-json --verbose`. Non-streaming `--output-format json` buffers the whole response and, on long generations, exposed a hang under WSL2; streaming keeps bytes flowing and lets you inspect partial progress.
3. **Timeouts** — NDE queries are heavy (256K datasets, 10–17 KB results); steered C needed **~900s**, whereas the tiny GxA graph finished well under 420s. Always run each condition under a `timeout` with retries so a stall fails fast instead of blocking.
4. **"Claude alone" must also disallow `ToolSearch`/`Monitor`** — otherwise the model burns the whole budget *hunting* for database tools and thinking, never answering. Pair the deny-list with an `--append-system-prompt` telling it it has no tools.
5. **`jq` `// false` pitfall** — `jq -r '.is_error // true'` misreads a successful run (`is_error:false`) as failure, because jq's `//` treats `false` as empty; with retries enabled this silently **doubles spend**. Use bare `.is_error`.
6. **`gh pr edit` vs projects-classic** — editing title/body via `gh pr edit` can be aborted by a projects-classic GraphQL deprecation error; use the REST API (`gh api -X PATCH repos/{o}/{r}/pulls/{n}`) instead.

---

## 10. Bottom line

For *"find datasets **and** sample counts for a pooled DN DE analysis,"* **Claude + web (D)** wins on quality-per-dollar. **NDE (steered)** is the strongest *discovery/triage* front-end but must hand off to GEO for sample-level counts. **Claude alone** is a cheap, honest scoping sketch. **Unsteered MCP** spent the most to answer the least — a routing-metadata problem, now addressed in PR #441.

**Recommended workflow:** NDE (discover candidate datasets) → GEO/web (verify per-sample DN/control counts, de-duplicate ERCB/superseries) → GxA (pull computed DE values where a curated study exists) → pool within compartment.

---

### Appendix — reproduction

```bash
# all three conditions, isolated sessions, token table + CSV
./run_dn_condition_comparison.sh
# single condition, larger budget (NDE needs it):
RUN_ONLY="C_claude_mcp_nde" TIMEOUT=900 ./run_dn_condition_comparison.sh
```

Outputs in `dn-comparison-run1/`: `*.answer.md` (each answer), `*.json` (result events), `*.stream.jsonl` (full event stream), `token_usage.csv`. Unsteered C preserved as `C_UNSTEERED.*`.
