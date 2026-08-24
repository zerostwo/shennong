# Shennong Maintainer Status

Last updated: 2026-08-22

CellTypist Seurat export now preserves sparsity end to end. The selected layer
is written as MatrixMarket plus exact gene/cell sidecars; the default retains
gene-by-cell orientation with `--transpose-input`, while the explicit false
branch writes a sparse cell-by-gene transpose. A 20,000-by-10,000 matrix with
three nonzero values round-trips as sparse with an MTX file below 1 KB, and the
helper contains no dense matrix/data-frame conversion. The focused clustering
suite passes 531 assertions with the pre-existing BBKNN warning and 20 optional
dependency skips in the current dependency profile. BPCells still crosses into an in-memory sparse Matrix object
at O(nnz); a zero-materialization BPCells writer remains separate work.

Pkgdown now starts from a research-workflow map rather than a flat article
list. Four real-data narratives connect immune-cohort, tumour/clinical,
cell-dynamics and spatial questions to module articles and explicit visual
checkpoints. Five article groups, a local Bootstrap 5 theme, accessible
light/dark colors, responsive cards, print/reduced-motion styles and a local
sidebar replace the external tidytemplate/font/analytics runtime. The final
core real-data build at `site/dev` contains 24 article HTML pages and 65
audited figure assets across the 15 mapped core articles. Runtime tracing observed all 94 declared core
analysis/visualization functions across 15 mapped executable articles, with
zero missing functions, zero failed articles and zero attempted downloads; 58
optional-backend cells remain explicitly extended. The same-build manifest
records 678 static-site files. Local-link validation passes and the site has no
external runtime scripts or stylesheets.

The repository redundancy audit found no exact duplicate tracked source or
duplicate R function definitions. A guarded cleanup script removed 20.71 MiB
of reproducible check/build products, logs, test plots and Python bytecode plus
one stale Git worktree record. The stale 9.01 MiB pkgdown tree was then removed
and replaced by the current real-data site. Cleanup deliberately retained 2.50 GiB of research
outputs, 1.30 GiB of real fixtures, and 1.86 GiB of Coralysis benchmark inputs
whose exact base object is not yet governed. The classification and future
canonical-cache work are recorded in `RedundancyAudit-2026-08-21.md`.

Usage observability now instruments 257 of 267 exports after explicit
`sn_enable_usage_tracking()`: every analysis, plot/get/list/store, IO,
validation, administration and rio adapter except ten tracking-control APIs.
The SQLite v2 contract records session mode and
material package versions plus per-run invocation/configuration ordinals,
nested parentage, sanitized parameter JSON/hash, elapsed and CPU time,
warnings/errors, entrypoint/component scope, consent/outbox state and
activation-only AutoZyme evidence. Focused usage tests pass 127
assertions. A 100-call microbenchmark measured about 6.14 ms of incremental
wall overhead per tiny call; after excluding the initial database insert from
the stored analytical timer, its persisted median matched the approximately
2 ms function body. Scientific return values and RNG state were unchanged.
Managed remote DBI delivery now requires a tamper-checked, versioned
`remote_research` consent, uses the SQLite file as an outbox, and is explicit or
on-disable rather than a network dependency of scientific calls. Delivery
selects only sessions with the exact consent receipt, applies field-level data
category gates, omits PID/error text, and keeps differently consented rows
pending. A temporary second SQLite database exercises the generic DBI contract,
unique remote IDs, idempotent retry and remote query path. The current
namespace-binding model explicitly reports that function references saved
before enablement are outside interception.

The generated public API matrix inventories all 253 `sn_*` functions and every
formal/selector digest; its gate passes 682 assertions. The clustering matrix
accounts for all three normalization and eleven integration methods with
selector-exhaustive, pairwise and high-risk case requirements. A separate
24-cell enrichment matrix inventories ORA/GSEA, GO/KEGG/MSigDB and
grouped/ungrouped dispatch, marking grouped GSEA unsupported and only the two
existing MSigDB pilots as pilot evidence. These artifacts
are static coverage/admission requirements: `runtime_classified` means the
parameter was classified, not that its case executed, and neither matrix is a
claim that every pending backend has passed a real execution contract.

A governed 30,000-cell PBMC benchmark fixture now combines three independent
public captures, 10,000 cells each, over 20,453 common features; its SHA-256 is
`20b884fb6c8f20c932ae0eb3e56484e40089386208653156b447439032945013`.
The fresh-process, three-arm benchmark uses the exact pinned AutoZyme revision
and rotating order. Three-repetition medians with passing output comparators
are LISI 9.11x, UCell 33.26x, Seurat NormalizeData 4.49x, Assay5 merge 3.04x,
JoinLayers 1.69x, and scDblFinder 1.75x. The repeated clusterProfiler GO
annotation scope was only
1.04x in its focused one-repetition test; fgsea has no current `sn_enrich()`
call intersection. scDblFinder is retained as a separate time/memory tradeoff
result because its speedup raised peak worker RSS by about 290 MiB.
The standard-counts preparation fast branch reduced the unaccelerated
NormalizeData wrapper from 3.45 to 3.01 seconds against a 2.69-second direct
call and reduced its worker RSS by nearly 1 GiB relative to the pre-fix report.
The automatic patch set is now the explicit nine-patch subset `cellchat`,
`clusterprofiler`, `lisi`, `nichenetr`, `scdblfinder`, `seurat`,
`seurat_merge`, `soupx`, and `ucell`, each behind its operation/input guard; it
is narrower than the manual catalog. Coralysis and
WGCNA guards do not match Shennong's scientific defaults; the pinned build has
no standalone decontX provider; broad Seurat includes unvalidated/unsafe
targets (especially RPCA anchors); default JoinLayers does not match its
counts-only guard; and tradeSeq lacks an mgcv dependency guard. Those paths are
explicit-only. NormalizeData is the sole automatic broad-Seurat operation,
RPCA/CCA and default JoinLayers fail closed, all-default scDblFinder is strict,
and owned label-transfer/simulation merges now scope the exact Assay5 patch.
Profile-guided backend cleanup moves only composition/ROE grouping to a shared
base-R contingency kernel. Randomized 10k/20k differential tests plus the
existing composition suite pass 74 assertions, preserving factors, NA rules,
multi-column grouping, `min_cells`, proportions and RO/E. Repo-wide tidyverse
replacement is rejected: PBMC30k wrapper overhead is already negligible in
most paths, grouped QC was faster in dplyr than data.table, and tibble result
classes are public contracts. A zero-materialization BPCells-to-MatrixMarket
writer remains separate work; the current sparse handoff is O(nnz).
The complete local testthat suite finishes with
`FAIL 0 | WARN 2 | SKIP 66 | PASS 4277`. The warnings are the pre-existing
BBKNN command-log warning and sandbox `timedatectl` discovery; skips are
optional backend packages or the unavailable local figure fixture. Focused
usage/matrix/composition tests pass 883 assertions; the AutoZyme hooks plus
clustering/enrichment conformance group passes 410 with one BPCells skip.
Focused clustering passes 531 with the known warning and 20 optional skips,
and program scoring passes 49 with one optional GSVA skip.
`R CMD build .` completes. After the full local suite, the rebuilt tarball
passes `_R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual --no-tests` with
`Status: OK`, including examples and vignette reconstruction; all three
installed conformance matrix paths were then opened and parsed explicitly.
The real-data pkgdown build under `site/dev` contains the new usage references,
`runtime-observability.html`, the revised `autozyme-benchmarks.html`, and the
research-workflow map; its same-build manifest and publisher dry audit are
regenerated after final source changes. Runtime coverage observes 94/94 core
functions across 15/15 successful articles with no download attempts. The
publisher dry audit accepts 678 static files and 65 mapped core figure assets,
and a separate scan resolves 292 HTML pages with zero broken relative links.
The strict backend-conformance profile passes 381 assertions without warnings
or skips.

Backend conformance v1 now defines the repository's AutoZyme-inspired
scientific admission contract: direct upstream oracle versus Shennong with
acceleration disabled, followed by an optional accelerated arm. Six
machine-readable pilot contracts now cover Seurat log-normalization, silhouette
widths, edgeR quasi-likelihood bulk DE, the full no-batch Seurat clustering
policy pipeline, and MSigDB-style ORA/GSEA. A frozen inventory
accounts for all 64 methods that already had `implemented: true`; a new method
cannot enter that historical backlog and must instead supply an admitted
contract. The dedicated strict CI profile installs all pilot dependencies,
forces single-thread execution, and permits no dependency skip. The focused
local run passes 363 assertions with no warnings or skips. These pilots are C0
and C1 evidence only; none is yet represented as fresh-process C2 or release
admission evidence.

BPCells-backed scran normalization now crosses the unsupported
BPCells-to-scran boundary through direct sparse `dgCMatrix` materialization.
The original BPCells counts layer remains attached to the returned Seurat
object. Focused preprocessing tests pass 119 assertions, and an end-to-end
120-cell `sn_run_cluster()` scran plus Harmony regression reproducer completes
with a Harmony reduction and cluster metadata. The complete local suite passes
with `FAIL 0 | WARN 1 | SKIP 3 | PASS 3269`; the warning is the existing BBKNN
command-log warning. Source build, structural `R CMD check`, and the complete
pkgdown site rebuild also pass.

Parameter-grid clustering now writes optional atomic RDS checkpoints after each
completed combination and resumes matching calls without re-entering completed
runs. The comparison manifest schema is 2.1.0 and includes a run-level
performance table plus integration-specific elapsed time, R heap peak, and
Linux pixi backend process-tree RSS. `sn_integration_control_template()` exposes
the complete supported control surface for all eleven method names. Focused
clustering validation passes 554 assertions with the existing BBKNN command-log
warning, including a two-resolution checkpoint/resume test that mocks the
compute entry point on recovery to prove no recomputation occurs.
The complete local suite passes with `FAIL 0 | WARN 1 | SKIP 3 | PASS 3265`.
The real PBMC1K/PBMC3K/PBMC4K validation also passes on 8,274 cells and
19,721 features across 12 grid runs and six unique embeddings: raw counts are
unchanged, no default t-SNE is created, every run has performance columns, and
a second identical invocation resumes all 12 runs from the persisted
checkpoint before rerunning the scIB evaluation. The source package builds,
`R CMD check` completes with its existing `.codegraph` NOTE only, and the
incremental pkgdown site rebuild includes the new control-template reference.

## Current structure

- CodeGraph is initialized and synchronized for the repository's R, Python, and
  workflow source. The current index is used for symbol discovery and call-path
  review; Markdown, Rmd, Rd, data, and TOML still require exact text search.
- Package-maintainer documentation is indexed by `docs/codex/README.md`.
  Historical prompts, plans, audits, and long status logs are archived under
  `docs/codex/archive/`.
- Communication and regulatory activity implementations are separated into
  `R/analysis_communication.R` and `R/analysis_regulatory.R`.
- Repository benchmarks live together under `benchmarks/`. Development smoke
  scripts live under `scripts/smoke/` and are not installed with the package.

## Latest cleanup

The newest entries appear first. Older entries remain as point-in-time evidence;
historical validation counts and removed APIs do not describe the current
release gate.

- Added explicit local usage/timing tracking, its SQLite schema and privacy
  policy, a pkgdown article, a workflow-registry coverage assertion, and an
  overhead benchmark. Added the PBMC30k preparation and three-arm AutoZyme
  runners. The benchmark audit also found and fixed exact feature matching in
  program scoring and tightened relaxed AutoZyme gating so source trust is
  never disabled.

- Direct-reference audits found and repaired stale clustering reuse after count
  or metadata mutation, selection of an unrelated pre-existing SNN graph,
  hard-coded CCA/RPCA reduction provenance, split-count layer loss, unconditional
  HGNChelper admission, and missing blocked-HVG evidence. Six focused clustering
  regressions pass 34 assertions, and the no-batch Seurat stage projection
  matches its direct pipeline. Enrichment now passes its universe, adjustment,
  q-value, size, and exponent controls to the versioned upstream; rejects
  ambiguous numeric formulas and silent duplicate-rank collapse; preserves
  serious backend warnings; and no longer claims an fgsea patch intersection
  for clusterProfiler 4.20/enrichit. Direct ORA/GSEA plus storage-provenance
  pilots pass 10 assertions without warnings. The complete DE/enrichment file
  passes 92 assertions with three deliberately preserved upstream mapping/
  q-value diagnostics and one unavailable-COSG skip. The focused clustering
  regressions pass 34 assertions, preprocessing passes 77 assertions with ten
  optional-package skips, and the source build plus structural `R CMD check`
  finish with `Status: OK`. The pkgdown rebuild was attempted twice: the normal
  destination is not touchable by the current `pi` process because existing
  files are owned by `duansq`. A temporary-destination build successfully
  rendered the home/reference pages and the changed analysis-results,
  annotation, bulk-transcriptomics, and clustering articles; the all-article
  pass later stopped at the unrelated AutoZyme benchmark fixture. The checked-in
  `site/dev` tree was therefore not rewritten in this session.

- Added `shennong.dev/backend-conformance/v1`, installed JSON contracts, a
  static formals-to-parameter audit, a frozen legacy-method inventory, and a
  strict `backend-conformance` workflow. Six executable pilots compare
  Shennong with direct Seurat, cluster, edgeR, and clusterProfiler/enrichit
  calls while checking input immutability and wrapper-owned projections. The
  strict focused run passes 363 assertions with
  `FAIL 0 | WARN 0 | SKIP 0`; related preprocessing,
  metrics, bulk, and CI tests pass 483 assertions with 13 expected missing-
  optional-backend skips. The quick source build and structural `R CMD check`
  complete with `Status: OK`. A full-suite attempt in this reduced local
  dependency profile reaches the pre-existing Harmony-dependent clustering
  tests but cannot complete because Harmony is not installed; the dedicated
  conformance profile is complete.
  Pilot status is
  deliberately below admission: fixture hashes, fresh worker processes,
  real/OOD fixtures, evidence records, and extended version/platform matrices
  remain future C2--C4 work.

- `sn_run_cluster()` now expands conditional parameter grids for explicitly
  vectorized scalar controls and stores a versioned run/embedding/preprocessing
  manifest. The default multi-method path now creates UMAP only; t-SNE is
  explicit. Resolution variants recluster the existing graph and reuse the
  same UMAP. A real 40-cell smoke run covered 8 runs from two methods, two HVG
  counts, and two resolutions, yielding four unique embeddings/UMAPs and no
  t-SNE. `sn_compare_integrations()` now groups by preprocessing identity,
  exports each embedding once against its matching unintegrated baseline, uses
  one shared stratified cell set, and maps scores back to run-level rows. The
  focused two-group fake-backend test passed 36 assertions together with the
  legacy result-contract tests (38 assertions after embedding-level ranking was
  added). A real official PBMC1k/3k/4k run covered 8,274 cells and 19,721 genes
  with three methods, 1,000/2,000 HVGs, and resolutions 0.4/0.6: the object
  contains 12 run rows, six unique embeddings/UMAPs, zero t-SNE reductions, and
  unchanged raw count column sums. scib-metrics ran once per embedding in each
  preprocessing group; total scores for 1,000/2,000 HVGs were Harmony
  0.7615/0.7544, Coralysis 0.7135/0.7254, and unintegrated PCA 0.5785/0.5750.
  Both resolutions inherit the same embedding score, while within-group ranks
  are correctly 1/2/3. The saved qs2 result passes `sn_validate_result()`.
  The first full-suite pass reached 3,192 assertions before exposing an
  unrelated Seurat 5.4 `FeaturePlot()` signature change; the minimal assay
  compatibility repair then passed all 104 assertions in the original failing
  visualization file. The clean rerun of the original full-suite command passed
  3,237 assertions with three documented optional/local-fixture skips and only
  the existing BBKNN command-log warning. The final 0.3.0 source tarball built successfully;
  `_R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual --no-tests` completed
  with only the known `.codegraph` NOTE after an added namespace NOTE was fixed;
  and the complete pkgdown site rebuilt successfully under `site/dev` after two
  transient CRAN/Bioconductor DNS timeouts.

- Managed Python acceleration now uses the same pinned `zerostwo/autozyme`
  revision as the R integration. Static call-path review found real patch
  intersections only for BBKNN/Scanpy UMAP, scArches and stLearn Scanpy
  preprocessing, trajectory Scanpy preprocessing and dynamical scVelo,
  cell2location training, and CellPhoneDB statistical analysis. Those runners
  pin validated upstream versions, apply a second fail-closed runtime gate,
  activate guarded patches before their hot calls, and record activation
  evidence in backend manifests. Existing Shennong-managed pixi manifests gain
  the pinned dependency on their next preparation; custom workspaces are left
  unchanged. Backends without a matching patch remain explicitly unaccelerated.
  Focused Python tests passed 47 assertions; the full local suite passed 3,192
  assertions with the existing BBKNN warning and four optional/local-fixture
  skips. The source package built, `R CMD check --no-manual` completed with only
  the known `.codegraph` NOTE, and pkgdown rebuilt successfully under `site/dev`.
- `sn_compare_integrations()` now exports one sparse normalized matrix and the
  native reductions registered by a multi-method clustering run to a managed
  scib-metrics/JAX environment. It stores aggregate, per-metric, and ranking
  tables as a standard `integration_benchmark` result, records backend versions
  and JAX devices, excludes graph-only methods from the embedding comparison,
  and flags supervised reuse of the evaluation label. A real validation joined
  the official 10x PBMC 1k, 3k, and 4k matrices into 8,274 cells and 19,721
  features, retained unchanged raw count column sums, and completed
  unintegrated, Harmony, and Coralysis branches plus scib-metrics 0.6.0 over all
  8,274 cells and 2,000 shared features. The resulting total scores were
  Harmony 0.7544, Coralysis 0.7254, and unintegrated PCA 0.5750; these labels are
  coarse canonical-marker scores for runtime validation, not a curated gold
  standard. The saved 71.3 MiB qs2 object retains all three cluster columns and
  native/UMAP reductions, and its stored benchmark passes `sn_validate_result()`.
  The focused test file passes 27 assertions. On the local NVIDIA RTX 2000 Ada
  host, a clean managed GPU environment resolved the CUDA build of JAX and
  reported `jax.default_backend() == "gpu"` with device `cuda:0`; the temporary
  validation environment was removed afterward. The final full local suite
  passed 3,141 assertions with four optional/local-fixture skips and the
  existing BBKNN command-log warning. The 0.3.0 source tarball built
  successfully; `_R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual
  --no-tests` completed with only the known `.codegraph` NOTE; and the complete
  pkgdown site rebuilt successfully under `site/dev`.
- `sn_run_cluster()` now accepts an ordered vector of RNA analysis methods and
  reuses normalization, HVG selection, scaling, and PCA across them. The
  returned object keeps method-specific graphs, `<method>_clusters`,
  `umap.<method>`, and `tsne.<method>` results plus a discovery manifest under
  `object@misc$integration_comparison`. Per-method `integration_control` values
  are supported, including Harmony `theta` and Coralysis `icp_args`/`pca_args`.
  Real 48-cell `unintegrated + harmony` and 36-cell `unintegrated + coralysis`
  smoke runs completed with all expected components and unchanged raw counts.
  A separate `return_cluster = TRUE` smoke run confirmed the multi-method data
  frame contract.
  The focused clustering test file passed 515 assertions with only the existing
  BBKNN command-log warning. The full local suite passed 3,112 assertions with
  the same warning and four optional/local-fixture skips. The 0.3.0 source
  tarball built successfully; `_R_CHECK_FORCE_SUGGESTS_=false R CMD check
  --no-manual --no-tests` completed with only the known `.codegraph` NOTE; and
  the complete pkgdown site rebuilt successfully in 57.9 seconds after one
  transient CRAN DNS retry.
- HILCA's 183,547-cell Coralysis run exposed that Seurat v5 returned the
  normalized BPCells layer as `RenameDims`, which `Coralysis::PrepareData()`
  rejects. The adapter now subsets the selected integration features before
  sparse coercion. On the exact HILCA counts workflow, the 3,000-by-183,547
  input became a 76,692,818-nonzero `dgCMatrix` occupying about 0.875 GiB,
  while the complete 18,149-feature layer stayed BPCells-backed. The focused
  clustering suite passed 501 assertions; the full local suite passed 3,098
  assertions with four optional/local-fixture skips and the existing BBKNN
  command-log warning. The complete pkgdown site rebuilt successfully in 140.4
  seconds. A subsequent real rerun passed `PrepareData()` and built the
  Coralysis training set, but upstream `threads = 0` forked about 20 R workers
  during cluster-seed computation and exhausted 91 GiB RAM plus 8 GiB swap.
  Shennong and the HILCA script now default explicitly to one Coralysis worker.
  A fresh real-data rerun passed `PrepareData()`, built the training set, and
  entered cluster-seed computation with `Parallelism disabled, because threads
  = 1`; across 39 minutes it retained one R process at roughly 12--14 GiB RSS
  with 68--70 GiB host memory available and no new kernel OOM event. The
  rerun session ended before a completed Coralysis artifact was written, so
  this evidence closes the original type error and worker-fan-out OOM but is
  not an end-to-end output acceptance result.
- RNA integration is now layer-consistent across the Python backends:
  scVI/scANVI/totalVI and scPoli export the requested `assay`/`layer`, and the
  provenance records `source_layer`. `sn_run_cluster()` also supports a real
  scPoli latent backend and a BBKNN backend whose imported connectivity graph
  is used directly for clustering and graph-derived UMAP. Complete expression
  and protein inputs remain sparse; scPoli only densifies one bounded encoder
  minibatch, and low-dimensional PCA/latent/UMAP outputs are dense by design.
  Real pixi smoke tests completed with scArches/scPoli 0.6.1 (one training epoch,
  20-by-3 latent), BBKNN 1.6.0 (20-by-20 graph and 20-by-2 UMAP), and a sparse
  totalVI registration preserving both RNA and protein CSR matrices. Focused
  tests cover a deliberately distinct `decontaminated_counts` layer and prevent
  regression to hard-coded
  `counts`. AutoZyme Seurat scopes now select and log the actual operation
  (`runpca`, `scaledata`, and so on); unpatched operations no longer announce
  unrelated join/merge patches. The final local suite passed 3,094 assertions
  with four optional/local-fixture skips and one existing warning. The 0.3.0
  source tarball built successfully; `_R_CHECK_FORCE_SUGGESTS_=false R CMD
  check --no-manual` completed with only the known `.codegraph` hidden-directory
  NOTE; and the complete pkgdown site rebuilt under `site/dev` in 145.7 seconds.

- AutoZyme integration now pins `zerostwo/autozyme` commit
  `3fe148271a21db90a6cc47ed91d82ed834de6bc4`. UCell, LISI, the refreshed
  scDblFinder patch, and SoupX are direct fork providers; the old Shennong
  scDblFinder and SoupX source/kernels were removed. The Shennong ambient
  path calls standalone `decontX::decontX()` and requests the
  `decontx_standalone` hook. The local AutoZyme build currently registers that
  hook, while the pinned public fork `main` revision does not yet contain it;
  `sn_check_autozyme()` therefore remains the source of truth for whether the
  decontX fast path is available in a given installation.
- `sn_remove_ambient_contamination(method = "auto")` now detects one
  ADT/protein/CITE assay and routes it to `decontX::decontPro()`; explicit
  `method = "decontpro"` and `assay =` are available for reproducible CITE-seq
  runs. The direct decontX matrix path was smoke-tested with integer output,
  zero-count handling, and scoped local AutoZyme activation. decontPro is
  intentionally not described as accelerated until a combined fork revision
  registers a validated patch for it.

- Historical note, superseded by the strict 2026-08-21 admission policy:
  CellChat and call-safe NicheNetR communication workflows temporarily passed
  `strict = FALSE` to their scoped AutoZyme manager. They now require exact
  admitted upstream versions. The scDblFinder refresh branch is merged
  into `main` with its five-target sparse/IRLBA patch and updated
  source fingerprint; the already-merged Seurat branch has no remaining delta.
- Historical benchmark note, superseded by the current enrichment call graph:
  clusterProfiler 4.20 and fgsea 1.38 were once scoped independently inside
  `sn_enrich()`. The Shennong-bundled clusterProfiler patch reports one active
  `get_GO_data()` target and the official fgsea patch reports three active
  targets, but current GSEA uses enrichit and no longer intersects fgsea.
  Focused enrichment/acceleration/hook tests passed 407 assertions
  without warnings. The measured enrichment file completed in 161.8 seconds
  versus 224.5 seconds before activation in this same checkout; a final
  post-refactor rerun passed its 91 assertions in 136.8 seconds. Repeated cached GSEA
  calls fell from about 15 seconds on the first call to about 1 second. ORA
  remains dominated by the upstream enrichit statistical core, so this is a
  measured local improvement rather than a universal speed claim. The source
  tarball built successfully; the pre-merge isolated `R CMD check` passed
  2,619 assertions, and the post-merge check passed 2,629 assertions, with 15
  Seurat 5.5.1 SCTransform validation warnings and 14
  optional/repository-only skips. Check completed with zero errors/warnings and
  one known local `.codegraph` NOTE. The complete pkgdown site rebuilt in
  206.4 seconds.
- Ambient and doublet clustering are now explicit backend choices. decontX and
  scDblFinder default to their native automatic clustering, while callers can
  select `cluster_backend = "shennong"`; SoupX retains Shennong clustering.
  Seurat AutoZyme scopes now tolerate version-label drift and cover neighbor and
  cluster calls in addition to the earlier stages. A Seurat 5.5.1 smoke run
  completed with all scopes restored, and a paired 1,000-feature/400-cell run
  returned identical clusters with 6.213 s upstream versus 5.483 s accelerated
  in that single directional sample. This is compatibility evidence, not a
  generalized performance claim.
- The SoupX, scDblFinder, and Seurat merge/JoinLayers accelerator branches are
  integrated into `main`. A reproducible formal
  benchmark uses the validated 2,000-cell, 32,738-feature Kotliarov PBMC
  CITE-seq fixture, three fresh workers per condition, alternating order,
  operation elapsed time, and whole-worker peak RSS. Median speedups were
  3.08x for `merge()`, 5.28x for `JoinLayers()`, and 2.26x for
  `sn_find_doublets()`, with exact canonical outputs and inactive patches after
  all 18 runs. Peak RSS changed by +8.2%, -20.0%, and -8.1%, respectively;
  the merge memory regression is retained in the shipped JSON/CSV evidence and
  pkgdown article rather than generalized away. The merged full suite passed
  2,999 assertions with four optional/local-fixture skips, the complete pkgdown
  site rendered the new article, and the source tarball plus structural
  `R CMD check` completed with `Status: OK`.
- Automatic AutoZyme scopes now emit an INFO log after successful activation,
  naming the patches enabled for that workflow call while preserving the
  distinction between patch activation and an internal guarded fast-path hit.
  The focused acceleration tests passed 193 assertions, the workflow-hook tests
  passed 111 assertions, and the full local suite passed 2,980 assertions with
  four optional/local-fixture skips. The pkgdown site rebuilt successfully.
- The exact-scoped AutoZyme `seurat_merge` and `seurat_joinlayers` sources are
  now bundled under `inst/autozyme/patches/` with their scope, manifest,
  changelog, and benchmark evidence. Shennong admits only the exact SHA-256
  fingerprints and SeuratObject 5.4.0/5.4.0.9001, registers only
  `merge.Assay5` and `JoinLayers.Assay5`, and leaves broader Seurat/v3/SCT
  dispatch upstream. The method contracts passed 33 merge assertions, 18
  JoinLayers assertions, and six combined-chain assertions; the full local
  suite passed 2,978 assertions with four unrelated optional/local-fixture
  skips. A clean source build/install proved Shennong-owned registration and
  rollback against an AutoZyme registry without either entry, and a BPCells
  JoinLayers call matched upstream exactly. Structural `R CMD check` completed
  without errors (its two warnings are the expected missing `inst/doc` outputs
  from a deliberate `--no-build-vignettes` tarball), and the complete pkgdown
  site rebuilt successfully.
- The validated scDblFinder 1.27.6 patch is now bundled under
  `inst/autozyme/patches/scdblfinder/` and integrated into
  `sn_find_doublets()`. Shennong verifies the installed source SHA-256, scopes
  activation around compatible default `dgCMatrix` calls up to 33,000 cells,
  and restores the prior AutoZyme state. Focused acceleration tests passed 166
  assertions. A clean temporary installation verified the installed patch
  fingerprint, real registration of all four namespace targets, activation and
  rollback, and a 3,005-cell real SCE workflow run in 25.073 seconds with the
  patch inactive afterward.
- Eligible CellChat, NicheNetR, clusterProfiler, fgsea, Seurat, SoupX,
  tradeSeq, and WGCNA patches now activate lazily only inside compatible
  Shennong workflow calls, with the pre-call state restored after success or
  error. The automatic
  path requires the pinned AutoZyme build and installed, exactly validated
  upstream versions; absent or drifted dependencies and approximate patches are
  skipped safely. Option/environment opt-outs block automatic scopes without
  deactivating manually active patches, and explicit helpers ignore them.
  BPCells-backed Seurat layers bypass the Seurat fast patch to avoid
  `dgCMatrix` materialization; a manually active Seurat patch is suspended for
  that call and then restored. AutoZyme's `future.globals.maxSize` load-time
  mutation is restored before analysis. Result-producing scopes retain patches
  active for compatible calls in provenance. CellChat, tradeSeq, and other
  backend contracts remain separate BPCells-compatibility boundaries and may
  require controlled materialization or aggregation. The final local suite
  passed 2,847 assertions with four optional/local-fixture skips; the source
  package built as `Shennong_0.3.0.tar.gz`; packaged tests passed 2,439
  assertions with 14 repository-only/optional skips; and `R CMD check`
  completed without errors or warnings. Its sole NOTE is caused by the live
  CodeGraph daemon socket during local source staging; `.codegraph` is already
  excluded and a clean checkout has no socket. The complete pkgdown site also
  rebuilt successfully under `site/dev`.

- The SoupX ambient-correction path now scopes the exact AutoZyme `soupx`
  patch around `adjustCounts()` and restores the prior patch state on exit.
  Shennong now ships the patch R source and compiled native kernels, validates
  their SHA-256, and registers them through official AutoZyme at runtime when
  that package does not provide `soupx` itself. The installed patch directory
  also carries its exact `SCOPE.md` and a finalized benchmark table containing
  only equivalence-passing PBMC 10k/20k measurements; failed and baseline-less
  runs are excluded.
  Shennong requests the attested fractional correction and then applies the
  same stochastic integer rounding as `SoupX::adjustCounts(roundToInt = TRUE)`;
  a seeded `scToy` comparison returned an identical `dgCMatrix`. Local strict
  compatibility accepts the exact recorded Shennong-bundled patch SHA-256.
  A clean source installation found the installed bundled patch, compiled and
  loaded the Shennong native library, registered against an official AutoZyme
  build with no `soupx` entry, returned an identical result, and restored the
  patch to inactive. Focused acceleration and clustering suites passed 146 and
  460 assertions; the fresh-process full suite passed 2,875 assertions with
  four optional/local-fixture skips, and pkgdown rebuilt successfully.
- `sn_set_layer_backend()` now provides a bidirectional Seurat layer-storage
  boundary: selected layers can be staged and rebound to external BPCells
  directories or materialized directly as `dgCMatrix` without a dense
  intermediate. Multi-layer round-trip, preflight, source-retention, automatic
  `uint32_t`, and unsafe-conversion behavior are covered by focused tests. The
  standard pre-push path passed 2,613 tests with four optional/local-fixture
  skips, built the source tarball, completed structural `R CMD check` with
  status OK, and validated the pkgdown reference index; the incremental
  pkgdown site rebuild also completed successfully.
- Issue #7 tracks BPCells compatibility for Seurat initialization and
  sample-aware doublet detection. The implementation preserves BPCells
  `IterableMatrix` counts in `sn_initialize_seurat_object()`, converts BPCells
  chunks directly to `dgCMatrix` without a dense intermediate, and runs
  `scDblFinder()` on independently materialized `group_by` samples. A real
  BPCells/scDblFinder smoke test resolved 240/240 cells across two samples
  while leaving the returned counts layer as `RenameDims`. Focused tests passed
  467 assertions; the full suite passed 2,597 assertions with four expected
  optional-dependency/local-fixture skips; the source package built; and
  `_R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual` completed with no
  error or warning and only the existing `.codegraph` hidden-directory NOTE.
  The pkgdown site rebuilt successfully under `site/dev`.
- R-CMD-check run `30214842915` exposed one real package-metadata warning:
  `tests/testthat/test-figure-engine.R` conditionally loads `DOSE`, but the
  package had not declared it in `Suggests`. `DESCRIPTION` now declares the
  optional dependency; no hard installation dependency was introduced. The
  focused figure-engine suite passed 35 tests with one local-fixture skip, and
  the rebuilt source package passed 2,168 tests with 14 source-package skips,
  no warning, and no failure. Its sole local check NOTE is the live CodeGraph
  socket, which is absent from a clean GitHub checkout. Follow-up
  R-CMD-check `30215609091`, test-coverage `30215609075`, and pkgdown
  `30215609067` all passed for revision
  `e58f7fbfea1ea02cebc07efddcabb66d6e8bb48a`.
- The five-repository implementation milestone now has explicit source,
  commit, publication, image, deployment, and end-to-end states in
  `docs/codex/Ecosystem.md` and the machine-readable compatibility lock
  `docs/codex/ecosystem-lock.json`. The loop is governed by
  `shennong.dev/data-bundle/v1`,
  `shennong.dev/analysis-result-bundle/v1`, and
  `shennong.dev/runtime-r-toolchain/v1`, with OS PAT/project authority and DB
  exact/digest readback. All five implementation revisions are pushed, the
  three service images are published and inspected, and no service is recorded
  as deployed.
- Shennong and ShennongData package gates pass. OS revision
  `e3f421bcba70e85a82596634178f24f3270df621`, Runtime revision
  `1649f1da43b25d79e8d820bd1ef093cf49a77114`, and DB revision
  `24445c3b08bb708bcdc2d1cf2eccba1735815633` have green remote CI and Docker
  publication. Together they implement PAT/RBAC, exact toolchain admission,
  staged immutable input, Result Bundle validation, authoritative platform
  provenance, exact upload/readback, SHA-256, service-only no-redirect reads,
  and zero-byte Artifact support. Immutable OCI digests and run IDs are
  recorded in the ecosystem lock.
- The modality verdict remains conservative. scRNA-seq is the closest because
  Shennong has mature analysis and DB has a real PBMC3K headless fixture, but
  it has not crossed all five repositories. Bulk support means bulk RNA only;
  spatial and CITE-seq mean prepared-object/contract support; scATAC-seq has
  descriptor/transport contracts but no production analysis or five-service
  evidence. Project resource discovery is a bounded snapshot; stable cursor
  pagination is a follow-up.
- GitHub pkgdown run `30208397482` restored an exact dependency cache containing
  `stringfish` 0.19.0, then upgraded `RcppParallel` to 6.0.0. The cached
  `stringfish.so` still referenced the legacy TBB ABI and failed while pkgdown
  read `README.md`; this was a native dependency-cache defect, not a Shennong
  source failure. The pkgdown workflow now starts cache epoch 2, load-tests
  `stringfish` after dependency resolution, and rebuilds it from source only
  when the active `RcppParallel`/TBB ABI rejects the cached binary.
- Shennong now owns the local candidate-result handoff boundary
  `shennong.dev/analysis-result-bundle/v1`. `sn_build_result_bundle()` carries
  one validated canonical analysis result, immutable input
  identifier/revision/SHA-256 references, package/execution provenance, and candidate artifact
  roles/digests; `sn_export_result_bundle()` writes credential-free JSON and
  returns its SHA-256 digest. This does not imply OS authorization, byte
  upload, DB promotion, or an end-to-end ecosystem loop.
- `sn_list_methods(available = TRUE/FALSE)` now uses the caller value outside
  dplyr's data mask. Named Cox phenotypes in Scissor continue to retain their
  aligned sample identifiers, and the regression test now treats those names
  as part of the public alignment evidence rather than ignoring them.
- The original five-repository audit remains embedded as historical evidence
  in `docs/codex/Ecosystem.md`. Its pre-implementation gap table is explicitly
  marked historical; the current interface ledger above it is authoritative.
- The final 2026-07-26 full local package suite on the pre-existing dirty
  development tree completed with `PASS 2576`, `SKIP 4`, `WARN 0`, and
  `FAIL 0` in 542.9 seconds. The skips are explicit optional/local-fixture
  boundaries: `lisi`, `ROGUE`, `scmap`, and the unavailable local public-data
  figure fixture. The previously failing Scissor expectation now preserves
  and asserts the aligned sample identifiers, and the three Result Bundle
  exports are accounted for by the runtime-coverage governance test.
- The final clean staged source build and tarball audit exclude `.codegraph`
  and other repository-only state. The staged source package passed
  `_R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual`, including its
  complete packaged test suite, with `Status: OK` in 597.5 seconds. The
  check-safe AutoZyme metadata test now reads the installed package
  DESCRIPTION, and the corrected result-audit Rd markup passes `checkRd`.
  A subsequent `scripts/build-pkgdown.R --full` rebuilt the complete installed
  site, including the Result Bundle reference and discovery/retrieval example,
  in 266.8 seconds.
- Data distribution has moved out of Shennong. The package no longer exports
  `sn_load_data()`, `sn_list_datasets()`, `sn_download_zenodo()`, or
  `sn_upload_zenodo()`, and no longer ships the `pbmc_small` analysis fixtures.
  Real public-data fixtures for local validation and pkgdown rendering live
  under the ignored `SHENNONG_REAL_DATA_DIR` boundary. The four small runtime
  reference assets required by annotation, species mapping, gene filtering,
  and signature lookup remain package data. The current source, namespace, and
  data tree reflect that ownership change; the final full-suite,
  source-package, `R CMD check`, and clean non-real-data site verdicts are
  recorded above. A full real-data site rebuild remains a separate local-data
  gate.
- The local real-public-data matrix is prepared as four logical bundles and five
  checksummed artifacts: Kotliarov PBMC CITE-seq (RNA/ADT, 2,000 cells, 20
  samples); GSE72056 melanoma (6,000 features, 1,500 cells, 19 tumors) paired
  with TCGA-SKCM (3,000 finite features, 155 patients, 75 overall-survival
  events); Hermann spermatogenesis (2,000 features, 1,000 cells, spliced and
  unspliced layers); and fresh-frozen plus FFPE/CytAssist Visium lymph node (400
  real spots per section). The current validator accepts all five artifacts.
  Acquisition code, source/license metadata, deterministic subset rules, and
  coverage declarations are tracked; raw and derived data are not.
- `scripts/build-pkgdown.R --real` validates the complete local matrix before
  evaluating articles with real tables and figures, without downloading data.
  Article-level smoke rendering has already produced local real-output HTML and
  plot assets. `--full --real` remains the clean integrated site gate, and
  `--extended` is reserved for installed optional backends.
- The runtime-coverage runner traces declared analysis and visualization
  functions while their mapped articles render, blocks network calls, and
  writes ignored CSV/JSON evidence. The currently saved report is an
  intermediate diagnostic snapshot, not a release verdict: it exposed article
  failures and unobserved core rows that are being resolved before the final
  clean run. No final runtime-coverage count is recorded here until that rerun
  completes.
- Analytical result boundaries are now audited and normalized on schema
  `1.0.0`. New `sn_audit_results()` and `sn_upgrade_results()` APIs report and
  migrate canonical, legacy, and invalid analysis results; registered runtime
  payloads report as `artifact`, while every other populated top-level
  `object@misc` entry reports as `unregistered` and is never rewritten.
  Annotation (including label transfer), program scoring, trajectory,
  registered legacy collections, and standalone bulk workflows expose a
  canonical `tables$primary`. The durable coverage and exception inventory is
  recorded in `docs/codex/ResultContractAudit.md`.
  Unsupported or future schemas are never rewritten, critical contract fields
  use exact lookup, and canonical generic QC writes materialize their registered
  compatibility views before storage validation. Provenance types are checked
  before a result can be reported as valid or unified.
- Scissor now has a direct cell-first `sn_run_scissor()` entry point plus state,
  cell, sample, correlation, and reliability plots. It validates exact bulk
  sample alignment, retains all cell coefficients and model evidence, records
  the bulk sample as the inferential unit, and shares its backend implementation
  with the compatible `sn_prioritize_states(method = "scissor")` path. A missed
  upstream alpha-search cutoff is now retained as an explicit exploratory
  warning and diagnostic rather than silently accepted. The requested assay
  layer is now used consistently through preprocessing and backend execution;
  unavailable layers fail explicitly rather than falling back to counts.
- Bulk survival now retains adjusted Cox estimates, model failures, model
  performance, Kaplan-Meier and cumulative-hazard curves, group assignments,
  log-rank tests, risk tables, proportional-hazards tests, and scaled
  Schoenfeld-residual evidence. All six survival plot views build from retained
  tables without refitting models.
- Optional AutoZyme support is pinned to revision
  `718541d9489596c7c1d75f52e9b3a8b2a429d1f9` and is explicit, scoped, and
  transactional. Strict mode requires both that exact AutoZyme source and an
  exact validated upstream version; approximate patches require separate
  consent, and active patches are captured in analysis-result provenance. On
  the current maintainer host, the real-data benchmark recorded CellChat at
  1.264s versus 0.126s (10.0317x), with identical interaction output, and WGCNA
  at 7.781s versus 0.073s (106.5890x), with identical modules and eigengene
  difference `1.71e-07` plus trait-correlation difference `8.92e-08` within the
  declared `1e-05` tolerance. Patches were inactive after each scoped run;
  these timings are machine-specific measurements rather than general
  guarantees.
- Real execution exposed and fixed four cross-version/data-contract defects:
  CellChat now safely encodes numeric group labels and decodes its output;
  Seurat 5 label transfer requests the table-returning `TransferData()` form;
  `sn_plot_feature(assay = "ADT")` uses the selected assay for shared limits;
  and result-aware plots use exact lookup so backup fields cannot partially
  match `tables$primary`. Focused regressions cover each boundary; the final
  integrated package and site gates remain pending as described above.
- The post-`v0.2.0` development line now unifies single-cell and standalone
  bulk differential expression under `sn_find_de()`: Seurat objects select the
  single-cell path, while matrices, lists, and `SummarizedExperiment` inputs
  select bulk DE. `sn_find_bulk_de()` remains a compatibility wrapper over the
  same internal bulk implementation.
- RegVelo is registered and implemented as
  `sn_run_velocity(method = "regvelo")` in the managed trajectory pixi
  environment. The adapter accepts a named target-by-regulator matrix, an edge
  table, or a CSV prior GRN; standardizes velocity, latent time, confidence,
  embeddings, model artifacts, and manifest metadata into the existing
  velocity result contract; and documents downstream CellRank use. A real pixi
  solve succeeds with RegVelo 0.4.2, PyTorch 2.5.1 CPU, SciPy 1.15.2, and
  scvi-tools 1.1.6.post1.
- Shennong now ships a read-only stdio MCP server plus a validated
  `use-shennong-mcp` Agent Skill. Agents can list methods, inspect method
  status, read exact exported-function help, retrieve workflow guides, and
  inspect package entry points without arbitrary R evaluation or analysis-file
  mutation.
- Routine CI is split from release-quality validation. Push-time pkgdown builds
  are lazy and skip examples, while scheduled/manual full builds remain clean
  and exhaustive. R CMD check no longer reruns tests already covered by the
  coverage job unless a manual full-test input is selected. Path filters and
  concurrency cancellation avoid obsolete documentation-only or superseded
  runs. Locally, the complete pkgdown build took 147.8 seconds and the
  immediately following incremental build took 29.4 seconds, a roughly 80%
  reduction.
- The accumulated roadmap, workflow, compatibility, documentation, and
  publication-figure changes were published as GitHub Release `v0.2.0` from
  commit `5f75d4a`. The tag is non-draft and non-prerelease; `main` now advances
  to `0.2.0.9000`, and `NEWS.md` has a fresh `Unreleased` boundary so subsequent
  features cannot silently enter the released source. The release gate passes
  `FAIL 0 | WARN 0 | SKIP 6 | PASS 1906`, builds
  `Shennong_0.2.0.tar.gz`, completes structural
  `_R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual` with `Status: OK`,
  and validates the complete pkgdown reference index. Remote release checks
  pass in `R-CMD-check` run `29405957672`, coverage run `29405957618`, and
  pkgdown deployment run `29405957610`.
- The first remote `0.2.0` pkgdown run reached the evaluated annotation article
  but exposed a clean-runner dependency omission: UCell was not installed for
  the default program-scoring example. The website workflow now declares
  `bioc::UCell`, and the corrected clean-runner pkgdown build and deployment
  passes in GitHub Actions run `29400422339`.
- The corrected run passed the UCell-dependent article and then exposed an
  unconditional Zenodo download in `data-io-projects.Rmd`; Zenodo returned HTTP
  504. That article now requires the explicit
  `SHENNONG_RUN_NETWORK_VIGNETTES=true` opt-in so website builds are not coupled
  to third-party availability.
- The same clean-runner gate also exposed an R 4.6 development-runner failure
  where the multimodal roadmap test's four-formal `sn_run_cluster()` namespace
  mock leaked into the later clustering test file, plus two
  optional-backend checks that ran before dependency-independent input
  validation. The roadmap mock is now expression-scoped so the real export is
  restored immediately; its focused suite passes 35 assertions and a post-test
  check confirms all 18 public wrapper formals remain intact. Clustering also
  uses a compact public compatibility
  wrapper and a single-argument private implementation; its allowlisted tail
  resolver
  preserves named and positional calls, defaults, and explicit `NULL` values
  with explicit validation. Scissor/Symphony
  validate required inputs before checking their optional packages. The full
  clustering module passes `FAIL 0 | WARN 0 | SKIP 0 | PASS 442`; the combined
  clustering, abundance, and annotation regression previously passed
  `FAIL 0 | WARN 0 | SKIP 1 | PASS 505`, and the data-I/O article renders with
  network examples disabled.
- Full source-package checking exposed that excluding only CodeGraph contents
  could still leave an empty top-level `.codegraph` directory in the tarball.
  The root directory now has an explicit build-ignore rule, the active index
  socket was stopped before packaging, tarball audit confirms the directory is
  absent, and the rebuilt source completes full `R CMD check --no-manual` with
  `Status: OK`.
- The final pre-merge `testthat::test_local(stop_on_failure = TRUE)` run passes
  with `FAIL 0 | WARN 0 | SKIP 6 | PASS 1906`. Skips are limited to unavailable
  optional lisi, ROGUE, scmap, and zen4R dependencies.
- Pre-merge validation exposed one stale registry expectation and one benign
  `enrichit` qvalue fallback warning on a tiny deterministic ORA fixture. The
  registry test now matches the implemented Slingshot state, and the exact
  backend fallback is muffled without hiding other enrichment warnings;
  focused regression tests pass 116 assertions without warnings.
- The analysis/publication roadmap was completed on
  `feat/analysis-publication-roadmap`, merged into `main` as `1f85f80`, and the
  merged source passed `FAIL 0 | WARN 0 | SKIP 6 | PASS 1906` plus full
  `_R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual` with `Status: OK`.
  The merged topic branch was then pruned locally.
- `README.Rmd` now exposes a module-by-module one-command software matrix and
  regenerates `README.md`. The matrix covers 33 public workflow entry points
  and all 64 registry methods, and explicitly distinguishes direct R,
  Shennong-managed pixi/CLI, and external runner/result adapters from current
  local dependency availability.
- The final roadmap adapter pass adds the explicit multimodal entry point,
  keeps backend-specific annotation/trajectory/CellRank functions internal,
  adds direct optional Monocle 3 inference, standardizes external Palantir and
  scCODA/pertpy results, and leaves no method-registry entry marked
  unimplemented. Focused adapter tests pass with
  `FAIL 0 | WARN 0 | SKIP 0 | PASS 34`, including a real Monocle 3 principal
  graph and pseudotime run. The combined abundance/trajectory/adapter quick
  suite passes 97 assertions, its structural `R CMD check` reports
  `Status: OK`, and the complete pkgdown site rebuild includes the multimodal,
  abundance, and trajectory documentation.
- Milestone D now has six generic figure profiles, automatic specifications for
  500 through 5,000,000 simulated points, publication preflight QA,
  PDF/SVG/TIFF/PNG export, reproducible bundles, result-aware DE/enrichment/GSEA
  figures, core QC/clustering/integration diagnostics, and spec migration for
  the established core and bulk plot families. Figure-engine tests pass with
  `FAIL 0 | WARN 0 | SKIP 0 | PASS 34`; the quick structural check passes with
  `Status: OK`, and the full pkgdown site includes the publication-figure
  article and all new reference pages.
- Milestone C now includes standalone bulk QC, fixed/mixed-design differential
  expression through edgeR, DESeq2, limma, and dream, pathway scores, WGCNA,
  Cox survival, clinical associations, six result-aware plot families, and a
  matrix/list/`SummarizedExperiment` input contract. The focused suite passes
  with `FAIL 0 | WARN 0 | SKIP 0 | PASS 39`, including real DESeq2, dream, and
  WGCNA backend runs. The quick pre-push source build and structural
  `R CMD check` pass with `Status: OK`, and the complete pkgdown site includes
  the new article and all thirteen bulk reference pages.
- Milestone A now has a shipped method registry, availability diagnostics, and
  a common analysis-result schema with generic store/get/list/delete/validate
  APIs. Existing registered result collections are upgraded compatibly rather
  than copied.
- Milestone B1 now has the unified annotation API, marker/reference consensus,
  SingleR/CellTypist/Seurat/Symphony/scmap/scANVI adapters, confidence and
  runner-up calibration, hierarchical labels, a versioned Cell Ontology
  mapping, low-confidence review, and three result-aware diagnostic plots.
- Milestone B2 now has UCell/AUCell/GSVA/ssGSEA/mean program scoring, signature
  coverage diagnostics, cell metadata storage, sample-aware program tests, and
  result-aware activity/heatmap plots.
- Milestone B3 now has Slingshot topology/pseudotime/lineage probabilities,
  terminal-state and curve storage, optional tradeSeq dynamic/branch tests,
  fitted trends and convergence diagnostics, plus six result-aware plots.
- Milestone B4 now has unified Propeller/permutation/Milo abundance tests,
  sample-level evidence and contributions, Augur-style sample-held-out state
  separability, real Scissor bulk integration, RareQ topology discovery plus
  phenotype association, and two result-aware plots.
- Milestone B5 now has one comparable communication schema across LIANA,
  CellChat, CellPhoneDB, NicheNet, and MultiNicheNet; cross-method consensus and
  concordance; sample-level LR evidence and condition comparisons; retained
  ligand-target evidence; and seven result-aware plot modes.
- Milestone B6 now has multi-restart local NMF discovery, cNMF/Hotspot
  adapters, a real optional GENIE3 backend, explicit pySCENIC/SCENIC/GRNBoost2
  adapters, regulon activity and group-specificity summaries, and six
  result-aware program/GRN plot modes.
- Milestone B7 now has a managed scVelo/CellRank pixi environment, real
  spliced/unspliced velocity inference, projected vectors and transition
  evidence, GPCCA terminal states/fate probabilities/lineage-driver import,
  and two result-aware dynamics plots.
- Milestones C1/C2 now have one spatial dispatcher plus explicit feature,
  domain, neighborhood, deconvolution, mapping, integration, and communication
  APIs. The local path includes Moran's I with permutation evidence,
  memory-bounded KNN graphs, neighborhood enrichment/co-occurrence, and
  distance-constrained communication; nnSVG/BANKSY and existing Python
  backends remain discoverable optional paths.
- The CNV/malignancy/metabolism milestone now wraps inferCNVpy and CopyKAT in
  one stored result, exports chromosome evidence from the pixi backend, derives
  reference-calibrated malignancy/subclone/sample diagnostics, and adds a
  curated sample-aware metabolism workflow plus optional heavy-backend adapters.

- Removed internal helpers that had no callers and removed the now-unused
  `data.tree` and `later` dependencies.
- Consolidated three overlapping Codex design notes into `Governance.md`.
- Hardened the pre-push build path so temporary tarballs/check directories are
  cleaned and source-package contents are audited for local CodeGraph or plot
  artifacts.
- Preserved all pre-existing 2026-07-02 SCTransform `block_genes` source,
  documentation, skill, and test changes.
- Preserved the merged grouped-dotplot coordinate fix and the
  minimal-dependency vector fallback for feature plots without `ggrastr`.
- Aligned the pending data-server integration with the installed ShennongData
  0.2 client contract.
- Declared `msigdbr` in the pkgdown workflow's explicit website dependency
  list so the feature-annotation vignette can execute on a clean CI runner.
- Marked the CITE-seq WNN template as non-executing because the article does
  not construct or ship the required multimodal `pbmc_cite` input object.

## Validation

- The complete local suite passes with
  `FAIL 0 | WARN 0 | SKIP 6 | PASS 1960` in 535.7 seconds. Skips are limited to
  unavailable optional lisi, ROGUE, scmap, and zen4R packages. The focused MCP
  suite passes 19 assertions, including source-tree and installed-package help
  rendering.
- `R CMD build .` succeeds, and the resulting source package passes
  `_R_CHECK_FORCE_SUGGESTS_=false R CMD check --no-manual --no-tests` with
  `Status: OK`; the test phase is intentionally omitted there because the full
  testthat suite was run separately. The complete and incremental pkgdown
  builds both finish successfully.
- The RegVelo Python entry point parses successfully, its registry and manifest
  contracts pass focused tests, and the managed pixi dependency graph resolves
  on CPU. End-to-end biological model training is not part of this local
  validation pass.
- Spatial workflow tests pass with
  `FAIL 0 | WARN 0 | SKIP 0 | PASS 46`, covering spatial autocorrelation,
  adapters, domains, metadata storage, graph/enrichment/co-occurrence evidence,
  communication distance filtering, integration, the dispatcher, registry,
  and eight rendered plot types.
- Velocity/fate R contract tests pass with
  `FAIL 0 | WARN 0 | SKIP 0 | PASS 20`. A real CPU pixi smoke run using
  scVelo 0.3.4 and CellRank 2.3.2 completed end to end on 80 synthetic cells,
  returning 80 velocity vectors, 2,594 transition edges, and 80 CellRank fate
  probabilities.
- Program-discovery and GRN focused tests pass with
  `FAIL 0 | WARN 0 | SKIP 0 | PASS 79`, covering local NMF, external program
  adapters, GRN edge standardization, derived and supplied regulon activity,
  group specificity, method discovery, metadata storage, and all plots.
- CNV/metabolism tests pass with `FAIL 0 | WARN 0 | SKIP 0 | PASS 56`, covering
  the actual inferCNVpy import adapter, CopyKAT-style predictions, reference
  calibration, all CNV plots, curated and external metabolism backends,
  sample-level contrasts, registry discovery, and all metabolism plots.
- `scripts/check-prepush.R --filter=cnv-metabolism --quick` passes source build,
  structural `R CMD check` with `Status: OK`, and reference-index validation;
  the complete pkgdown site also rebuilds with the new article and five new
  reference pages.
- Communication tests pass with `FAIL 0 | WARN 0 | SKIP 0 | PASS 48` against
  real NicheNet and MultiNicheNet backends,
  synthetic CellChat/LIANA standardization, CellPhoneDB file parsing,
  sample-level contrasts, consensus/concordance, unified result storage, and
  all communication plots.
- Abundance and state-priority tests pass with
  `FAIL 0 | WARN 0 | SKIP 0 | PASS 34`. The suite runs real speckle Propeller,
  miloR, RareQ, and Scissor backends plus sample-label permutation and
  sample-held-out separability on deterministic synthetic data.
- The quick pre-push source build and structural check pass with `Status: OK`,
  and the complete pkgdown site rebuilt with the abundance/state-priority
  article and four new reference pages.
- Trajectory tests pass with `FAIL 0 | WARN 0 | SKIP 0 | PASS 29`, including
  deterministic branching and linear Slingshot paths, real tradeSeq GAM fits,
  result-contract retrieval, convergence/trend tables, and rendered ggplot
  grobs. The quick pre-push path passed source build and structural check; the
  initial Rd markup note was corrected before the final validation pass.
- The complete pkgdown site rebuilt successfully with the trajectory article
  and all seven new trajectory/dynamic-gene reference pages.
- Program-scoring tests currently pass with
  `FAIL 0 | WARN 0 | SKIP 0 | PASS 50`, including all five scoring backends,
  sparse mean scoring, coverage/drop behavior, sample-level inference, and
  rendered ggplot grobs.
- `scripts/check-prepush.R --filter=program-scoring --quick` passes the source
  build, structural R CMD check, and reference-index validation with
  `Status: OK`; the complete pkgdown site also rebuilds with the program
  scoring/comparison documentation and four new reference pages.
- Annotation workflow tests currently pass with
  `FAIL 0 | WARN 0 | SKIP 1 | PASS 39`; scmap is the single skip because the
  optional package is not installed locally. SingleR and Symphony CPU paths run
  against deterministic synthetic query/reference objects.
- `scripts/check-prepush.R --filter=annotation-workflow --quick` passes the
  annotation tests, source build, tarball audit, structural R CMD check, and
  pkgdown reference-index validation with `Status: OK`.
- The complete pkgdown site was rebuilt from an installed copy of the current
  source; the expanded annotation article and all eight new reference pages
  rendered successfully.
- The Milestone A result/registry plus interpretation compatibility tests pass
  with `FAIL 0 | WARN 0 | SKIP 0 | PASS 142`.
- pkgdown was rebuilt from an installed copy of the current source; the new
  analysis-method/result-contract article and all existing references rendered
  successfully.
- `scripts/check-prepush.R --filter=analysis-registry-result --quick` passed
  its targeted tests, source build, tarball audit, structural
  `_R_CHECK_FORCE_SUGGESTS_=false R CMD check`, and pkgdown reference-index
  validation with `Status: OK`.
- `devtools::document()` completed and updated regulatory Rd source pointers
  after the module split.
- Focused tests for signatures, communication/regulatory activity,
  interpretation, visualization, and shipped Codex skills passed with
  `FAIL 0 | WARN 0 | SKIP 0 | PASS 308`.
- The full local suite passed with
  `FAIL 0 | WARN 1 | SKIP 5 | PASS 1429`. The warning is the known tiny-input
  enrichment `qvalue` fallback; skips are for unavailable optional packages.
- The temporary-source pre-push build completed, and its tarball-content audit
  confirmed that `.codegraph`, `Rplots.pdf`, and the development smoke script
  were absent. The structural `R CMD check` completed with `Status: OK` after
  the CI compatibility fixes.
- The complete pkgdown site rebuilt successfully with the local pkgdown template
  override.
- CodeGraph re-synced successfully at 65 indexed files, 1,091 nodes, and 3,015
  edges; the new regulatory module is discoverable and the removed helper names
  have no exact source occurrences.
- PR #4 passed `R-CMD-check`, `test-coverage`, and both Codecov statuses before
  merging into `main`.
- Post-merge pkgdown deployments exposed a missing website-only `msigdbr`
  installation and an unseeded CITE-seq example after a transient jsDelivr SSL
  failure; both reproducibility fixes are being validated on GitHub Actions.

## Deferred local data

Ignored `dev/outputs/` and benchmark input/run caches contain several gigabytes
of untracked data and scripts. They were not deleted automatically; see
`Roadmap.md`.
