# Backend Conformance Audit

Audit date: 2026-08-20

## Scope and verdict

This is a source-level audit of Shennong's public analysis wrappers, optional R
backends, managed Python runners, parameter forwarding, and stored-result
boundaries. It asks whether the repository currently proves that a Shennong
call agrees with its declared upstream call or executable reference workflow.

The package already has strong runtime, result-schema, and AutoZyme
acceleration checks. At the start of this audit it did not have a package-wide
backend-conformance gate. This change set adds the initial static admission gate
and six C1 pilots, but most backend paths remain a legacy-unverified backlog.
Most older tests establish that a workflow runs and returns the expected shape;
they do not compare the wrapper against an independently executed upstream
reference recipe. The findings below are audit findings only. None of the
listed implementation, registry, provenance, or result-validation defects is
fixed in this change set.

## Follow-up: clustering and enrichment, 2026-08-21

A focused direct-oracle pass added three further subjects: the complete
no-batch Seurat clustering policy pipeline, MSigDB-style ORA, and MSigDB-style
GSEA. The ORA and GSEA pilots use a fixed local term table so the direct
clusterProfiler comparison has no network dependency. Their version envelope
binds clusterProfiler 4.20.0, enrichit 0.2.1, and msigdbr 26.1.1; enrichit is
material because its overlap-size and significance filtering changed between
releases.

Confirmed and repaired in that focused change:

- `sn_run_cluster()` formerly reused normalization/HVG/PCA/graph state after
  same-shaped count values changed. Stage and checkpoint signatures now bind a
  blockwise selected-layer digest and relevant metadata values.
- A pre-existing lexically earlier `_snn` graph could be selected after
  `FindNeighbors()` overwrote `RNA_snn` in place. The current call's expected
  graph name now wins, and ambiguous fallback fails closed.
- Seurat CCA/RPCA accepted `integration_control$new.reduction` but returned and
  recorded a hard-coded reduction name; it now propagates the actual reduction
  and restores the original split count layers.
- Ordinary clustering no longer requires HGNChelper; only custom blocked gene
  symbols cross that dependency boundary. Log workflows now retain the removed
  blocked-HVG set, and resolution is excluded from the HVG cache key when rare
  feature selection is disabled.
- `sn_enrich()` formerly forced permissive upstream cutoffs and kept rows by raw
  p-value alone. It now delegates the effective p/adjusted-p/q rules to the
  versioned upstream, supports an explicit ORA universe and common size/
  adjustment controls, rejects silent duplicate-rank collapse, and requires an
  explicit mode for a numeric formula RHS.
- clusterProfiler 4.20 GSEA executes enrichit, not the fgsea namespace. The
  workflow therefore no longer activates or advertises the unrelated fgsea
  patch, while retaining the clusterProfiler cache scope through result
  storage provenance.

These remain pilot, not release-conformant, claims. Grouped GSEA
(`gene | score ~ group`), GO/KEGG resource and ID-mapping parity, actual msigdbr
resource fixtures, batch/Harmony/CCA/RPCA stage comparators, fresh-process
evidence, and real/OOD inputs remain open C2--C4 work. The earlier Milo,
spatial-registry, Python-runner, importer, and generic provenance findings below
also remain open.

## Follow-up: real workflow acceleration and observability, 2026-08-21

A fresh-process three-arm runner now compares direct upstream, unaccelerated
Shennong, and accelerated Shennong on the same call shape. The governed local
fixture samples 10,000 cells from each of three independent public PBMC
captures, retains 20,453 common features, and has 30,000 cells total. Its
derived SHA-256 is
`20b884fb6c8f20c932ae0eb3e56484e40089386208653156b447439032945013`.

With the exact pinned AutoZyme revision and three rotating repetitions, all
declared comparators passed. Shennong-off/Shennong-on median ratios were 9.11x
for LISI, 33.26x for UCell, 4.49x for Seurat NormalizeData, 3.04x for Assay5
merge, 1.69x for JoinLayers, and 1.75x for scDblFinder. UCell also reduced
whole-worker peak RSS by about 1.46 GiB. The scDblFinder path is a different
tradeoff: peak worker RSS increased by about 290 MiB. The GO annotation-cache
scope was only 1.04x in its focused
one-repetition test and is not evidence for a faster ORA/GSEA statistical
kernel.

The audit initially ran against an installed official AutoZyme 0.3.1 build with
no `RemoteSha`. That exposed a gate defect: `strict = FALSE` admitted both
upstream-version and AutoZyme-source drift. The gate now always requires the
pinned revision or an exact Shennong-vendored patch; relaxed mode applies only
to upstream version labels. The final main report uses an isolated library with
RemoteSha `8fc2e9c3a7f70302f97589aaa9b0395dcf86f9bc`. Seurat 5.5.1 and UCell
2.17.0 remain explicit upstream-version-drift evidence and are not promoted to
strict admission by this host audit.

The first 30,000-cell normalization report also exposed avoidable wrapper work:
the standard single-count-layer path still captured full temporary layer
snapshots. A no-temporary-layer branch reduced Shennong-off from 3.45 to 3.01
seconds (direct 2.69 seconds), reduced worker RSS from about 3.25 GiB to 2.27
GiB, and moved the pinned accelerated arm to 0.671 seconds. Focused layer,
preprocessing, and direct-normalization contracts pass after the change.

The benchmark records eligibility, static target intersection, scoped
activation, parity, elapsed time, RSS, source identity, and rollback. It leaves
`fast_path_hit` unknown because AutoZyme currently provides no per-call hit
counter. This is intentionally more conservative than equating an active patch
with executed accelerated code.

The resulting automatic policy is the guarded subset `cellchat`,
`clusterprofiler`, `lisi`, `nichenetr`, `scdblfinder`, `seurat`,
`seurat_merge`, `soupx`, and `ucell`. Coralysis, standalone decontX, broad
Seurat targets, JoinLayers, tradeSeq, and WGCNA remain explicit-only even when
they appear in the wider patch catalog or benchmark inventory.

The Coralysis patch has a target intersection but no current Shennong default
fast-path intersection. Its released guard requires the literal
`batch.label = "benchmark_batch"`, `build.train.set = FALSE`, an empty
`build.train.params`, and a fixed `k = 4`, `L = 30`, `C = 1`, single-thread
signature. `sn_run_cluster(integration_method = "coralysis")` passes the real
metadata column name and generated `build.train.params` (`nhvg` and `p`), in
addition to its workflow defaults. The patch can therefore be active while
always delegating this workflow call to upstream. The existing 8,274-cell PBMC
grid identifies Coralysis integration (roughly 38--47 seconds) as the slowest
clustering stage, but the current patch does not optimize it. Generalizing and
revalidating that guard is a higher priority than another activation benchmark.

An opt-in SQLite layer now complements scheduled benchmarks with real usage
frequency and runtime evidence. After explicit enablement it instruments 257 of
267 function exports: analysis, plot/get/list/store, IO, validation,
administration, backend and rio adapters, excluding ten usage-control APIs to
avoid recursion. It records explicit development/production/test/benchmark
mode and sanitized parameter hashes, and summarizes root calls without
double-counting nested time
(`R/usage_tracking.R:405-471,931-970,1918-1990`). Managed DBI delivery remains
local-outbox-first and requires explicit versioned `remote_research` consent
before a connection or flush. The current enable-time namespace replacement
does not intercept a function reference saved or imported earlier; the status
API exposes that limitation (`R/usage_tracking.R:1642-1658`). Usage rows are
not conformance evidence by themselves; they identify which admitted contracts
deserve optimization effort and which parameter envelopes occur in practice.

The generated public-API inventory plus clustering and enrichment method
matrices add static completeness gates for formals, selectors, pairwise axes,
dispatch cells, and high-risk cases.
They do not execute those cases: `runtime_classified` is a classification bit,
not a passed runtime cell. A passed claim still requires the direct upstream
oracle and comparator defined by an executable backend contract
(`inst/conformance/public-api-parameters.json`,
`inst/conformance/matrices/clustering-methods.json`, and
`inst/conformance/matrices/enrichment-methods.json`).

## Current surface and evidence

- `NAMESPACE` exports 267 function symbols: 253 `sn_*` functions and 14 `rio`
  adapter functions. The repository runtime matrix classifies all `sn_*`
  exports across `scripts/real-data/coverage.csv` and
  `scripts/real-data/coverage-exclusions.csv`; the completeness assertion is in
  `tests/testthat/test-real-data-runtime-coverage.R:103-129`.
- The shipped method registry contains 64 method entries: 40 R, 13 pixi, 10
  external-adapter, and one CLI entry under `inst/methods/`. Its required schema
  is descriptive and does not require an implementation symbol, upstream
  version, reference call, equivalence class, or conformance evidence
  (`R/analysis_registry.R:14-18`).
- The common analysis-result envelope requires input, parameters, tables,
  models, diagnostics, warnings, and provenance
  (`R/analysis_result.R:1-7,116-147`). Registered legacy analytical collections
  and runtime artifacts are separated in `R/interpretation.R:1-97`; generic
  analytical results are stored under `object@misc$analysis_results`
  (`R/analysis_result.R:463-506`).
- AutoZyme already supplies a useful conformance vocabulary: exact scoped,
  numeric-tolerance scoped, and approximate patches, with exact tested upstream
  versions (`R/acceleration.R:24-123`). The vendored merge contract also records
  its upstream commit, supported call shape, fallback surface, and exact
  comparators (`inst/autozyme/patches/seurat_merge/manifest.yml:13-75`).
- Existing exact tests cover SeuratObject merge and JoinLayers snapshots,
  input immutability, fallback conditions, and lazy arguments
  (`tests/testthat/test-seurat-merge-patch.R:157-334` and
  `tests/testthat/test-seurat-joinlayers-patch.R:78-155`). SoupX stochastic
  integer rounding is compared directly with upstream
  (`tests/testthat/test_clustering.R:3029-3061`).

## Why existing coverage is not conformance evidence

The real-data coverage runner traces whether a mapped Shennong export was
called while an article rendered (`scripts/real-data/run-runtime-coverage.R:188-219,433-500`).
It does not execute an independent upstream arm or compare numerical results.
External backends are disabled during that run
(`scripts/real-data/run-runtime-coverage.R:166-178`).

The general AutoZyme benchmarks compare the same Shennong wrapper with an
acceleration patch disabled and enabled. CellChat and WGCNA are examples
(`scripts/real-data/benchmark-autozyme.R:97-164,167-231`); the single-cell
benchmark does the same for `sn_find_doublets()`
(`scripts/real-data/benchmark-single-cell-autozyme.R:194-231,354-357`). This
proves patch equivalence, not wrapper-to-upstream conformance. Workflow-hook
tests primarily inspect that acceleration scopes occur before target calls
(`tests/testthat/test-autozyme-workflow-hooks.R:66-238`).

Ordinary CI also cannot serve as the complete backend gate. `R CMD check`
normally runs with `--no-tests` (`.github/workflows/R-CMD-check.yaml:82-84`),
and the coverage job installs only a selected optional dependency set
(`.github/workflows/test-coverage.yaml:45-76`). Backend tests guarded by
`skip_if_not_installed()` can therefore remain unexecuted in the baseline job.

## Confirmed audit findings

### Milo annotation retrieval uses a different field from storage

`sn_store_milo()` writes the annotation column name as `annotation_by`
(`R/analysis_metrics.R:2870-2882`), while `sn_get_milo_result()` reads
`stored$annotation_col` (`R/analysis_metrics.R:2931-2934`). Annotation-only
retrieval therefore cannot apply the requested filter. The existing test also
applies a SpatialFDR cutoff, which happens to select the expected row and masks
the field mismatch (`tests/testthat/test_milo.R:121-140`).

### Spatial registry claims have drifted from executable paths

- The registry describes Moran's I as a Squidpy/pixi method
  (`inst/methods/spatial.yml:3-10`), but `sn_find_spatial_features()` executes a
  Shennong R k-nearest-neighbor and permutation implementation
  (`R/analysis_spatial.R:121-143`). This should be treated as a native Shennong
  method with an independent oracle, not as Squidpy parity.
- The BANKSY entry names package `BANKSY`
  (`inst/methods/spatial.yml:23-30`), while DESCRIPTION and executable code use
  `Banksy` (`DESCRIPTION:45-47` and `R/analysis_spatial.R:189-205`). Availability
  checks use the registry package string with `requireNamespace()`
  (`R/analysis_registry.R:72-92`), so the case difference is material.
- `stlearn` and `squidpy` are marked implemented pixi methods
  (`inst/methods/spatial.yml:33-50`), but the canonical domain and neighborhood
  functions require an explicit `backend_control$runner` or precomputed result
  for those methods (`R/analysis_spatial.R:226-245,317-342`). They currently
  have adapter-schema boundaries rather than direct upstream fitting contracts.
- The separate runtime coverage inventory has already drifted: it labels
  `sn_deconvolve_bulk()` as BisqueRNA even though the public implementation
  supports BayesPrism and CIBERSORTx
  (`scripts/real-data/coverage.csv:84` and
  `R/analysis_deconvolution.R:1-15,90-174`). This demonstrates that duplicated
  method inventories are not a reliable truth source.

### scArches and stLearn names exceed the algorithms run by their default scripts

The default `scarches_run.py` path imports scArches but performs Scanpy total
normalization, log transformation, variable-feature selection, scaling, and
PCA. Its manifest explicitly calls the output a PCA latent rather than trained
reference mapping (`inst/pixi/scarches/scripts/scarches_run.py:59-100`). The R
wrapper also selects an existing `data` layer before `counts` when no layer is
specified (`R/package_tools.R:1205-1208,1272-1277`), so an already normalized
layer can enter the Python normalize/log path again.

Likewise, the stLearn runner imports stlearn but its executable work is Scanpy
normalization, log transformation, and PCA
(`inst/pixi/stlearn/scripts/stlearn_run.py:37-64`). These public names require
either a narrower declared adapter contract or a real method-specific upstream
reference workflow before they can claim algorithm-level conformance.

### The generic Python importer can omit method-defining outputs

The object-level Python boundary exports the selected sparse matrix and
metadata (`R/package_tools.R:1292-1332`). On import, it adds `obs.csv` and only
CSV files whose names match latent, PCA, UMAP, or embedding
(`R/package_tools.R:1394-1435`), then stores a top-level runtime manifest
(`R/package_tools.R:1437-1467`).

Consequently, outputs such as Tangram `mapping.csv`
(`inst/pixi/tangram/scripts/tangram_run.py:54-66`) and Squidpy graph/enrichment
state retained only in H5AD (`inst/pixi/squidpy/scripts/squidpy_run.py:46-54`)
are not imported as canonical Shennong tables or graphs by this generic path.
A successful process and manifest are transport evidence, not proof that the
method-defining result survived the Python-to-R boundary.

## Provenance, parameter, and result-validation gaps

- Default provenance records Shennong and R versions, a seed, and a timestamp,
  but not the executed backend package version or source revision
  (`R/analysis_result.R:47-73`). Most result inputs are counts and labels rather
  than immutable input digests.
- Stored single-cell DE omits `features`, `only_pos`, `subset_levels`,
  `min_cells_per_sample`, forwarded dots, and backend controls from its replay
  metadata (`R/analysis_de.R:602-624`).
- Bulk DE records the design, contrast, and high-level shrink flags but omits
  controls such as edgeR robustness, DESeq2 minimum-count filtering, and the
  actual shrink implementation/fallback from `parameters`
  (`R/analysis_bulk.R:211-248,337-345`).
- Program scoring records only `min_genes` and level; backend controls and a
  signature digest are absent (`R/program_scoring.R:265-290`). Communication
  stores only sender and receiver in `parameters`, despite a much wider backend
  surface (`R/analysis_communication.R:553-580,737-765`).
- Regulatory storage does not retain `minsize` or forwarded arguments
  (`R/analysis_regulatory.R:73-127,153-163`). Milo storage does not retain the
  graph, neighborhood, refinement, normalization, and FDR controls used by
  `sn_run_milo()` (`R/analysis_metrics.R:2645-2781,2870-2882`).
- Type-specific primary-column validation currently exists only for annotation,
  program scoring, trajectory, state priority, Scissor, and bulk survival
  (`R/analysis_result.R:30-40`). DE, enrichment, communication, CNV, spatial,
  Milo, and several other families can pass the generic validator with any data
  frame in `tables$primary`.
- Method availability checks package/executable/pixi presence but do not report
  conformance status, tested versions, or an evidence artifact
  (`R/analysis_registry.R:72-113`). `implemented = true` is therefore not an
  admission gate for backend equivalence.

## Recommended first migration sequence

1. **Establish one contract registry and fix static truth drift.** Add a
   `contract_id`, implementation kind (`transparent_wrapper`,
   `reference_workflow`, `adapter_schema`, `native_shennong`, or
   `transport_only`), upstream version/revision, equivalence class, fixture, and
   evidence path. Generate coverage metadata from that source instead of
   maintaining backend names independently. Correctness fixes identified in
   this audit, including the Milo field mismatch, should be separate reviewed
   changes.
2. **Migrate transparent deterministic wrappers first.** Start with
   `sn_normalize_data(method = "seurat")` (`R/preprocessing.R:406-431`),
   `sn_find_de()` marker/contrast paths (`R/analysis_de.R:32-139,515-639`),
   UCell scoring (`R/program_scoring.R:72-118`), and LISI
   (`R/analysis_metrics.R:43-78`). Compare every forwarded argument, input
   immutability, names/order, values, conditions, and stored-result round trips.
3. **Add composed reference recipes for the core workflows.** Cover the
   no-batch Seurat clustering path first, then Harmony and CCA/RPCA. Compare
   normalization, HVGs, scaling, PCA, graph, clusters, and projection stage by
   stage. Shennong's default blocked-feature workflow is intentional and must be
   reproduced by the reference arm rather than compared with raw Seurat
   defaults (`R/analysis_clustering.R:2882-2918,4206-4567`).
4. **Migrate statistical pipelines.** Add independent edgeR, DESeq2, and limma
   recipes for pseudobulk and standalone bulk DE, followed by Propeller and
   Milo. Compare the aggregated inputs, design and contrast direction, filter
   mask, raw backend result, standardized result, and retrieval behavior
   (`R/analysis_de.R:141-275`, `R/analysis_bulk.R:211-345`, and
   `R/analysis_abundance.R:124-411`).
5. **Extend already strong preprocessing evidence.** Convert the existing
   scDblFinder and SoupX acceleration comparisons into three-arm tests: direct
   upstream reference, Shennong with AutoZyme disabled, and Shennong with it
   enabled. Add non-default layers, grouping, filtered cells, zero-count policy,
   and decontX/decontPro paths (`R/preprocessing.R:1119-1448,1739-1927`).
6. **Gate heavy R and Python backends separately.** SingleR, Slingshot,
   tradeSeq, CellChat, WGCNA, nnSVG, BANKSY, CopyKAT, and similar direct
   backends should run in versioned optional jobs. Adapter-only and transport
   Python methods should first prove configuration parity, cell/feature order,
   artifact checksums, and lossless result import. Neural/GPU methods require
   invariant or distributional comparators rather than byte identity.

Each executable backend contract should use three arms where applicable:
independent upstream reference, unaccelerated Shennong, and accelerated
Shennong. Exact, numeric-tolerance, invariant, distributional, adapter-schema,
and native-oracle equivalence must remain distinct. A method should not become
`implemented = true` as a direct backend until its declared contract passes for
the admitted upstream version.
