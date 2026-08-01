# Shennong Maintainer Status

Last updated: 2026-08-01

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
  that package does not provide `soupx` itself.
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
