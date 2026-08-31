All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

# Shennong (development version)

### Breaking changes

- Removed legacy `.qs` serialization, its exported rio adapters, and its
  GitHub auto-installer. Use `.qs2` with `qs2` throughout project IO and
  benchmark scripts. Requests for `.qs` fail with a migration message rather
  than falling through to rio. Existing files require conversion in a
  compatible older environment; renaming their extension is insufficient.

- Removed `method = "consensus"` and Shennong's homegrown consensus algorithm
  from `sn_run_annotation()`. The entry point is now purely reference-based:
  the default method is `singleR`, and `celltypist`, `seurat`, `symphony`,
  `scmap`, `scanvi`, and `popv` remain available. The exported helpers
  `sn_annotation_consensus()` and `sn_annotation_confidence()`, the bundled
  marker-scoring path, margin calibration, second-best labels, and
  supporting/conflicting marker fields are gone; `sn_plot_annotation_markers()`
  is removed because marker evidence no longer exists. Wrapper parameters
  `marker_database`, `consensus_reference_method`, `low_confidence_threshold`,
  and `margin_threshold` are dropped. Cluster summaries now report the modal
  predicted label per group with an `agreement_share`; cells without a finite
  backend score or with rejected labels (`unassigned`/`unknown`) are flagged
  `low_confidence`. Hierarchy levels and Cell Ontology mapping are retained.
- Replaced AutoZyme acceleration with the companion **ShennongOpt** package.
  Public controls are now `sn_check_acceleration()`, `sn_enable_acceleration()`,
  `sn_disable_acceleration()`, and `sn_with_acceleration()`; automatic scopes
  are disabled via `options(shennong.acceleration = FALSE)` or
  `SHENNONG_ACCELERATION_DISABLED=true`, and thread budgets via
  `shennong.opt.threads` / `SHENNONG_OPT_THREADS`. Covered patch families are
  Seurat RunPCA/ScaleData, scran, decontX, scDblFinder, Coralysis, UCell, LISI,
  and Rogue. Hot paths without a ShennongOpt counterpart (CellChat, NicheNetR,
  SoupX, clusterProfiler caches, tradeSeq, WGCNA, Seurat merge/JoinLayers) now
  execute plain upstream code until new patches land. The vendored Python
  autozyme integration (pixi dependencies, runner hooks, `_shared` startup
  hook, and the autozyme benchmark artifact/article) was removed; managed pixi
  environments run plain upstream Python.

### Added

- `sn_add_qc_metrics()` now defaults to `suffix = NULL`: selecting
  `decontaminated_counts` (including its dot-separated split layers)
  automatically writes `percent.mt_corrected`,
  `percent.ribo_corrected`, and `percent.hb_corrected`. Other layers
  retain the original column names. Explicit suffixes, including `""` to
  overwrite original columns, take precedence.

- `sn_add_qc_metrics()` refreshes mitochondrial, ribosomal, and hemoglobin
  percentages independently of initialization, with explicit assay/layer
  selection and optional output-column suffixes for before/after correction
  comparisons. Totals come from the current count layer rather than stale
  `nCount_*` metadata; sparse/BPCells storage is preserved. Initialization
  now reuses this helper with the same signatures and hemoglobin patterns.

### Fixed

- Repaired GitHub coverage and strict backend conformance: optional Python
  probes now honor the configured Shennong runtime instead of a maintainer's
  home directory. Python distribution versions are checked by the managed
  interpreter, and strict CI provisions pinned PopV/Scrublet environments
  before executing their candidate/oracle comparisons. Missing environments
  remain failures in strict mode and skips in ordinary package coverage.
  Shared Scrublet smoke tests use the same probe, and workflow assertions use
  base R without an undeclared YAML-parser dependency.

- Fixed the invalid scDesign3 deprecation-help link that caused GitHub
  R-CMD-check to fail on warnings; qualified annotation helper bindings and
  wrapped long example lines to remove the accompanying check notes.
  The local pre-push script now rejects check warnings even when R exits zero,
  matching the GitHub warning policy.

- Fixed a crash when printing Seurat command records whose parameters contain
  nested lists (`Error ... type 'list' cannot be handled by 'cat'`), first
  observed with `object@commands$sn_remove_ambient_contamination`. Shennong now
  stores logged commands as an internal `sn_seurat_command` subclass of
  `SeuratCommand` whose `show()` method renders nested list parameters as
  compact deparsed one-liners instead of failing. The stored `@params`
  structure is unchanged, so programmatic access such as
  `@params$requested$method` continues to work.

- The ShennongOpt bridge now expands umbrella patch keys (`seurat`, `scran`)
  to every registered granular patch sharing their prefix
  (`seurat_runpca`, `seurat_scaledata`, `seurat_findneighbors`,
  `scran_computeSumFactors`, ...), keeping legacy wrapper keys working after
  ShennongOpt moved to per-operation patch names.

### Changed

- Integrated the new ShennongOpt accelerations. `sn_run_cluster()`'s
  variable-feature selection, nearest-neighbor graph, and clustering stages now
  run through the guarded `seurat_findvariablefeatures`, `seurat_findneighbors`,
  and `seurat_findclusters` patches whenever ShennongOpt is installed; the
  umbrella `seurat` scope expands to every registered granular patch
  automatically. The `seuratobject_merge` and `seuratobject_joinlayers` patches
  are reachable through explicit `sn_enable_acceleration()` requests via the new
  `seurat_merge` / `seurat_joinlayers` scope keys and are never activated
  automatically. Parity evidence on a 240-cell synthetic run: identical
  variable features and clusters, sign-aligned PCA difference ~3e-4.

- Deprecated the 13 environment-specific `sn_call_*()` pixi aliases
  (`sn_call_scvi()`, `sn_call_scanvi()`, `sn_call_mmochi()`,
  `sn_call_scarches()`, `sn_call_scpoli()`, `sn_call_infercnvpy()`,
  `sn_call_trajectory()`, `sn_call_cellphonedb()`, `sn_call_cell2location()`,
  `sn_call_tangram()`, `sn_call_squidpy()`, `sn_call_spatialdata()`, and
  `sn_call_stlearn()`). They now emit a deprecation warning and forward to
  `sn_call_pixi_environment("<environment>", ...)`, which remains the single
  supported runtime primitive for direct managed-Python commands. The aliases
  will be removed in a future major release; object-level `sn_run_*()`
  workflows are unaffected.

- Renamed the remaining noun-first getters and off-family actions behind
  `.Deprecated()` forwarding shims (AUDIT-13): `sn_enrich()` →
  `sn_run_enrichment()`, `sn_deconvolve_bulk()` → `sn_run_bulk_deconvolution()`,
  `sn_metabolic_signatures()` → `sn_get_metabolic_signatures()`,
  `sn_method_status()` → `sn_get_method_status()`, `sn_figure_spec()` →
  `sn_get_figure_spec()`, `sn_integration_control_template()` →
  `sn_get_integration_control_template()`, `sn_pixi_paths()` →
  `sn_get_pixi_paths()`, `sn_pixi_config_path()` → `sn_get_pixi_config_path()`,
  `sn_mcp_server_config()` → `sn_get_mcp_server_config()`, and
  `sn_mcp_server()` → `sn_run_mcp_server()`. The old names keep working with a
  deprecation warning. Vignettes, pkgdown references, and the shipped Codex
  skills now teach only the new names.
- Demoted `sn_simulate_scdesign3()` to a deprecated shim;
  `sn_simulate(method = "scdesign3", ...)` is the supported entry point.
- Backend-conformance contracts now accept `unsupported`, `fallback`, and
  `waived` parameter roles, and `parameter_map` is reserved for the dots
  formal; the stricter gate reclassified three real catch-all parameter
  sweeps in the popv and bulk-DE contracts (AUDIT-14).

- Reorganized modules without behavior change (AUDIT-08/09): the Python-backed
  runners moved from `package_tools.R` to their owning domains
  (`sn_run_scarches`/`sn_run_scpoli` to clustering, `sn_run_infercnvpy` to CNV,
  `sn_run_cellphonedb` to communication, `sn_run_cell2location`/`tangram`/
  `squidpy`/`spatialdata`/`stlearn` to spatial), `sn_store_enrichment()` and the
  stored-result retrieval family moved next to the unified result contract, the
  shared `.sn_get_misc_result`/`.sn_validate_seurat_object` helpers moved to
  `utils.R`, `sn_simulate*` moved to a new `simulation.R`, and
  `sn_detect_accelerator()` moved to `acceleration.R`. `analysis_clustering.R`
  shrank from ~7,000 lines by extracting label transfer
  (`analysis_label_transfer.R`), simulation, and CellTypist modules. Twenty-two
  inline Seurat guards now route through the shared validator.

### Added

- The method registry is now load-bearing (AUDIT-06): a two-directional
  parity test asserts every `inst/methods/*.yml` entry is accepted by its
  owning function's dispatch surface and every dispatch value is registered or
  explicitly waived. `sn_run_velocity()` and `sn_run_fate()` now declare their
  choice sets in the `method` formal like every other workflow. The hardcoded
  integration, normalization, doublet, and pseudobulk-DE choice sets are pinned
  by a fail-closed snapshot.
- Added `sn_delete_artifact()`: explicit deletion path for registered workflow
  artifact collections (`clustering_stage_cache`, `integration_comparison`,
  `label_transfer`, ...), by member or whole container, failing closed on
  unknown types (AUDIT-12).
- Stored-result writers `sn_store_milo()`, `sn_store_deconvolution()`, and
  `sn_store_regulatory_activity()` now record a provenance block with package
  versions, timestamp, and an optional `random_seed`; `sn_run_milo()` gained a
  top-level `seed=` that flows into the stored result (AUDIT-12).

- Unified the analysis interface behind `object=` aliases (AUDIT-04): every
  x-first exported analysis function now also accepts the Seurat object or
  stored result as `object=`, covering all `sn_calculate_*` metrics,
  `sn_assess_integration()`, `sn_run_milo()`, `sn_enrich()`,
  `sn_deconvolve_bulk()`, `sn_run_celltypist()`, `sn_filter_genes()`,
  `sn_filter_cells()`, and the result-aware `sn_plot_*` family. Supplying
  both `x`/first argument and `object=` fails fast. `sn_run_grn()`,
  `sn_score_programs()`, and `sn_discover_programs()` accept the standard
  `store_name=` storage argument alongside their legacy `name=` formals
  (which keep working).
- Added top-level `seed=` and `verbose=` controls to workflow entry points
  that previously buried them in control bags (AUDIT-05):
  `sn_run_cluster()`, `sn_run_trajectory()`, `sn_run_velocity()`,
  `sn_run_fate()`, `sn_find_spatial_features()`, and
  `sn_find_spatial_domains()`. Precedence is `seed > control$seed > task
  default`, and the resolved seed is stamped into stored-result provenance.
- `sn_plot_de()` and `sn_plot_enrichment()` now follow the object-or-result
  pattern: pass a Seurat object plus `de_name=`/`enrichment_name=` to plot a
  stored DE or enrichment result directly (AUDIT-07).

- Added `method = "scrublet"` to `sn_find_doublets()`: Scrublet doublet
  detection through scanpy's native `sc.pp.scrublet()` wrapper, managed by a new
  `scrublet` pixi environment
  (`sn_prepare_pixi_environment("scrublet", install_environment = TRUE)`).
  The backend requires raw counts (non-count layers fail fast), ignores
  scDblFinder-specific arguments with a warning, writes `scrublet.class` and
  `scrublet.score` metadata columns, and is covered by a pilot
  backend-conformance contract (`doublets::scrublet`) whose C2 evidence shows
  exact score/call parity against a direct scanpy oracle on a committed pbmc3k
  fixture.

- Added `method = "popv"` to `sn_run_annotation()`: multi-algorithm majority-vote
  consensus annotation through the YosefLab PopV Python package, managed by a
  new `popv` pixi environment
  (`sn_prepare_pixi_environment("popv", install_environment = TRUE)`).
  The backend requires raw counts for the query and an annotated reference,
  exposes upstream knobs (`methods`, `hvg`, `n_samples_per_label`,
  `prediction_mode`, batch keys, seed) through
  `backend_control = list(popv = ...)`, stores per-cell consensus labels plus
  agreement-normalized scores, and is covered by an admitted backend-conformance
  contract (`annotation::popv`) with a direct-upstream oracle, committed pbmc3k
  fixtures, and C1/C2 evidence.

### Fixed

- Resolved the 2026-08-24 comprehensive-audit P0 findings: AUCell program and
  metabolism scoring now share one control-aware implementation with a direct
  upstream differential test; legacy spatial deconvolution/mapping aliases
  emit deprecation warnings in favor of `sn_run_cell2location()` and
  `sn_run_tangram()`; and WGCNA execution uses an explicit local namespace
  wrapper instead of attaching and detaching the package (AUDIT-01 to
  AUDIT-03).
- Hardened stored-result discovery and mutation: registered workflow artifacts,
  including integration-comparison state, can be discovered with
  `sn_list_results(include_artifacts = TRUE)`; artifact-reserved type names are
  rejected by `sn_store_result()`; and deleting the final generic result now
  removes empty containers (AUDIT-12). The pkgdown index now includes
  `sn_call_trajectory()` (AUDIT-11).
- Updated the ShennongOpt bridge for the companion package's standardized
  `sn_*` management API, while retaining runtime compatibility with an older
  installed build during migration. `sn_remove_ambient_contamination()` now
  stores its public call plus a schema-versioned parameter record in
  `object@commands`: requested and resolved method/assay/clustering choices,
  wrapper defaults, supplied and effective method-specific arguments, compact input-source provenance,
  backend function, and backend package version are available without copying
  a raw count matrix into the command log.
- Fixed `sn_run_annotation(method = "scmap")` with current scmap releases: the
  adapter now detects whether `scmapCluster()` returned a cells-by-references
  or references-by-cells label matrix, degrades rejected assignments to an
  explicit low-confidence `"unassigned"` label with zero score instead of
  crashing downstream consensus, and warns when upstream returns no usable
  labels.
- Fixed the coverage workflow dependency set so the executable MSigDB
  enrichment conformance tests have `msigdbr` available on GitHub Actions.
- Fixed `sn_install_shennong()` so an explicit local source is installed
  without remote version probes, while `channel = "auto"` falls back to the
  current Shennong source tree when both CRAN and GitHub version checks are
  unavailable.
- Fixed `sn_run_regulatory_activity(method = "progeny")`: the PROGENy model is
  now reshaped from its wide gene-by-pathway form into the long
  source/target/weight network that `decoupleR` expects, instead of failing
  with "subscript out of bounds". `progeny_top` still selects the top targets
  per pathway through the upstream model.
- Fixed `sn_run_cell_communication(method = "liana")` with an unset
  `resource`: the default LIANA Consensus resource is no longer shadowed by an
  explicit `NULL`, which previously failed with "argument is of length zero".
- Fixed `sn_run_metabolism(method = "scmetabolism")` on Seurat v5 objects:
  scoring now reads the bundled scMetabolism GMT gene sets directly and scores
  them with Shennong's own gene-set machinery, so the upstream helper no
  longer touches removed Assay slots. `scoring_method` gains `"aucell"` (the
  upstream default) for this backend and for `method = "geneset"`.

### Added

- Added opt-in SQLite workflow observability with
  `sn_enable_usage_tracking()`, `sn_disable_usage_tracking()`,
  `sn_with_usage_tracking()`, `sn_check_usage_tracking()`,
  `sn_list_usage_runs()`, and `sn_summarize_usage()`. Every exported
  non-control function family can record elapsed time, warning/error status,
  nested parentage, explicit
  development/production/test/benchmark mode, workflow invocation number, and
  a privacy-sanitized parameter fingerprint. Tracking is disabled at load,
  excludes objects, identifiers, paths, credentials, prompts, and free text,
  and never uploads without a separately consented flush. `sn_time_call()`
  provides the same timing display for an ad-hoc expression.
- Usage observability now covers every safely wrappable public export, including
  plot/get/list/store, IO, validation, installation, project and backend-call
  APIs: 257 of the current 267 function exports, with only ten usage-control
  functions excluded to prevent recursion. Because instrumentation replaces
  namespace bindings when enabled, references saved or imported beforehand
  remain outside interception and are reported as such by the status API.
  `sn_create_usage_store()`, `sn_confirm_usage_consent()`, and
  `sn_flush_usage_tracking()` add a managed DBI destination backed by a local
  SQLite outbox. Remote connections require explicit versioned research
  consent. Stable tamper-checked receipts isolate differently consented
  sessions in one outbox, consent categories gate every optional remote field,
  and unique session/run IDs make completed rows retry-safe. Scientific calls
  never wait for the network. Allowlisted method selectors passed through
  variables are resolved without evaluating objects, paths, text, or arbitrary
  call expressions.
- UCell scoring now uses a serial BiocParallel fallback when its exact
  AutoZyme patch is not active, avoiding worker-socket requirements in
  restricted runtimes. An active admitted patch retains the validated
  `BPPARAM = NULL` fast-path envelope, and an explicit user `BPPARAM` is
  preserved.
- Added a generated all-public-API formal/selector inventory, an explicit
  `sn_run_cluster()` method matrix covering all three normalization methods and
  eleven integration methods, and a 24-cell `sn_enrich()` dispatch matrix for
  ORA/GSEA across GO, KEGG, MSigDB and grouped/ungrouped inputs. Pairwise axes,
  high-risk cases, unsupported grouped GSEA, and current pilot status are
  explicit.
  These are static completeness and admission matrices, not evidence that every
  listed case has executed or passed a backend comparison.
- Reorganized pkgdown around a new research-workflow map with four real-data
  scientific narratives, module-specific visual checkpoints, and five article
  groups. The Bootstrap 5 theme now has responsive workflow cards, accessible
  light/dark palettes, keyboard focus, reduced-motion and print behavior, and
  no external runtime fonts, stylesheets, analytics, or CDN scripts.
  The obsolete tidytemplate Suggests/Remotes/website installation dependency
  was removed from package metadata and pkgdown CI.
- Added `scripts/clean-generated.R` and a dated repository redundancy audit.
  The cleanup command is dry-run by default and refuses tracked files,
  symlinks, scientific fixtures, research outputs, Git state, user history,
  project settings, and CodeGraph data.
- Added a fresh-process, three-arm AutoZyme workflow benchmark and a
  deterministic preparation script for a 30,000-cell, three-capture public
  PBMC fixture. The runner separates direct-upstream agreement, Shennong
  wrapper overhead, acceleration parity, elapsed time, peak worker RSS, patch
  eligibility, scoped activation, and the still-unknown internal fast-path-hit
  state. On the maintainer host with the pinned AutoZyme revision, the
  30,000-cell median speedups were 9.11x for LISI, 33.26x for UCell, 4.49x for
  Seurat normalization, 3.04x for Assay5 merge, 1.69x for JoinLayers, and 1.75x
  for scDblFinder, with all declared output comparators passing. scDblFinder
  traded about 290 MiB of additional peak worker RSS for that time reduction.
- `sn_enrich()` now exposes the upstream ORA/GSEA controls
  `p_adjust_method`, `qvalue_cutoff`, `universe`, `min_gs_size`,
  `max_gs_size`, and `gsea_exponent`. Stored enrichment results retain these
  effective parameters plus the clusterProfiler, enrichit, annotation, and
  MSigDB package versions used by the call.
- Added direct backend-conformance pilots for the no-batch Seurat clustering
  pipeline and for MSigDB-style ORA/GSEA against
  `clusterProfiler::enricher()` and `clusterProfiler::GSEA()`, including an
  explicit ORA universe and an enrichit gene-set-overlap boundary case.
- Parameter-grid `sn_run_cluster()` calls can now persist and resume completed
  combinations with `checkpoint_dir`, `resume`, and `checkpoint_compress`.
  Checkpoints are atomically published after each run, retain only the latest
  complete state, and are matched by a package/object/grid/argument signature.
  Each comparison result now records workflow and integration wall time, peak R
  heap memory, and, on Linux pixi backends, maximum backend process-tree RSS;
  `sn_compare_integrations()` joins these fields into its scIB tables.
- Added `sn_integration_control_template()` with complete executable templates
  for every supported integration backend, including shared pixi/runtime,
  accelerator, model/training, graph, CITE-seq, and method-specific controls.
- `sn_run_cluster()` now expands explicitly vectorized scalar controls such as
  `nfeatures`, `npcs`, `resolution`, `cluster_algorithm`, `rare_feature_n`, and
  Harmony `theta` into a conditional Cartesian grid. Natural vector inputs such
  as `dims`, `hvg_features`, regression covariates, and blocked genes remain
  intact. The versioned comparison manifest records run, embedding, and
  preprocessing identities; resolution variants reuse their graph and UMAP.
- `sn_compare_integrations()` now evaluates each unique native embedding once
  within its matching preprocessing group and unintegrated baseline, then maps
  scores back to every run sharing that embedding. Results include run,
  embedding, preprocessing, baseline, within-group rank, and overall rank.

- Added exact-scoped, single-core AutoZyme acceleration for
  `Coralysis::RunParallelDivisiveICP()`. Coralysis integration activates the
  patch only around the backend call, restores prior AutoZyme state on exit,
  and remains pinned to the validated Coralysis 0.99.10 and AutoZyme revision.
- Managed Python backends now pin the same `zerostwo/autozyme` revision and
  activate only patches that match executed hotspots: Scanpy for BBKNN,
  scArches, stLearn, and velocity preprocessing; scVelo for dynamical recovery;
  cell2location training; and CellPhoneDB statistical analysis. Activation is
  fail-closed outside the exact validated upstream versions and unsafe
  cell2location GPU/minibatch modes, falls back to the original implementation,
  honors
  `AUTOZYME_DISABLED`/`AUTOZYME_DISABLE`, and is recorded in backend manifests.
- `sn_run_cluster()` now accepts multiple RNA integration methods in one call,
  including an explicit `"unintegrated"` PCA baseline. Shared preprocessing is
  reused while each method retains its own low-dimensional reduction, neighbor
  graphs, cluster column, UMAP, and t-SNE. The result map is stored under
  `object@misc$integration_comparison`, and `integration_control` can be keyed
  by method for backend-specific parameters.
- `sn_run_cluster()` now accepts `integration_method = "scpoli"` and
  `"bbknn"`. scPoli trains a real scPoli model in the managed scArches pixi
  family and imports its latent representation. BBKNN runs in a dedicated
  managed pixi environment, imports the batch-balanced connectivity graph, and
  uses that graph directly for Seurat clustering and UMAP.
- Added `sn_compare_integrations()` to benchmark the native reductions retained
  by a multi-method `sn_run_cluster()` result with scib-metrics. Every method
  uses the same cells, features, batch labels, and biological labels; normalized
  expression remains sparse; UMAP/t-SNE coordinates are never metric inputs;
  graph-only methods and supervised label reuse are reported explicitly. The
  summary, per-metric values, ranking, backend versions, and JAX devices are
  stored as a discoverable `integration_benchmark` analysis result. Managed CPU
  and CUDA-specific JAX builds are selected explicitly, and a selected GPU
  environment that does not report a JAX GPU device is retained as a result
  warning.

### Fixed

- `sn_run_celltypist()` now exports Seurat layers as sparse MatrixMarket with
  exact gene/cell sidecars instead of materializing a dense genes-by-cells CSV.
  Both orientations are tested, executable discovery follows the documented
  option/`PATH` order, and existing path inputs no longer require Seurat. The
  unified annotation adapter now uses `counts` instead of silently
  double-normalizing the default `data` layer. XLSX-only output, CellTypist's
  automatic over-clustering, and one-column small-cell majority-vote results
  are handled according to the actual upstream output shape. User-supplied
  output directories are never mistaken for temporary directories and removed,
  and custom output prefixes are passed to the CLI as one quoted argument.
  Explicit normalized/scaled Seurat layers now fail closed instead of producing
  a known double-normalized result, and XLSX import dependencies are checked
  before the external process starts. When CellTypist supplies its probability
  matrix, Shennong stores the selected-label probability as confidence; missing
  confidence is conservatively scored as unresolved instead of fabricating 1.0.
- Real-data documentation now keeps core rendering dependency- and
  download-safe: built-in abundance permutation replaces an unconditional
  Propeller requirement, the feature-class example uses a frozen local
  resource, Louvain replaces an auto-installing Leiden call, and
  HGNChelper/Slingshot/tradeSeq examples are explicitly extended. Runtime
  coverage traces both the source namespace and attached exports and Git
  provenance commands accept the repository as an invocation-local safe
  directory. The final core audit observes all 94 declared functions across 15
  articles with no failures or network attempts.

- Automatic AutoZyme admission now contains only owned call shapes that are
  currently defensible. The automatic set is `cellchat`, `clusterprofiler`,
  `lisi`, `nichenetr`, `scdblfinder`, `seurat`, `seurat_merge`, `soupx`, and
  `ucell`; each still applies only to its documented operation and input
  envelope. Coralysis, standalone decontX, broad Seurat targets, JoinLayers,
  tradeSeq and WGCNA remain explicit-only until their guards and Shennong
  defaults pass new contracts. RPCA/CCA layer integration and default
  JoinLayers are forcibly suspended from a manually active broad patch;
  all-default scDblFinder is strict while non-default calls are forcibly
  upstream. Standard-count `sn_run_cluster()` no longer passes a redundant
  `layer` through NormalizeData dots, so the validated normalization fast path
  can execute, and owned label-transfer/simulation merges now scope the exact
  Assay5 merge patch.
- `sn_calculate_composition()` and `sn_calculate_roe()` now use a shared base-R
  contingency kernel for their large metadata count/totals stage while
  preserving factor, NA, `min_cells`, multi-group and tibble output semantics.
  A one-million-row profile showed the base count core faster than both dplyr
  and data.table; grouped QC and small result-table transformations remain in
  their clearer existing form because profiling did not justify replacement.
- Relaxed AutoZyme checks now relax only the tested upstream package version;
  every automatic or explicit activation still requires the pinned AutoZyme
  revision or an exact trusted Shennong-vendored patch source. An unverified
  same-version AutoZyme installation can no longer be silently admitted by a
  workflow using `strict = FALSE`.
- Removed `fgsea` from the automatic workflow-default inventory because no
  current Shennong call reaches its namespace: clusterProfiler 4.20 GSEA uses
  enrichit. The fgsea patch remains available for explicit management of
  direct fgsea calls.
- Program scoring now prioritizes exact feature names before case/version
  normalization. Signatures containing both an exact symbol and a
  version-suffixed feature no longer lose one member through a duplicated
  normalized lookup key, restoring direct UCell parity on the real PBMC
  benchmark.
- Standard single-`counts`-layer Seurat workflows no longer capture and restore
  redundant full-matrix temporary layer snapshots. On the 30,000-cell
  NormalizeData benchmark this reduced the unaccelerated wrapper boundary from
  3.45 to 3.01 seconds and cut its whole-worker peak RSS from about 3.25 GiB to
  2.27 GiB while preserving direct output parity.
- `sn_enrich()` now delegates `pvalueCutoff` and related significance rules to
  the versioned clusterProfiler/enrichit backend instead of forcing permissive
  upstream cutoffs and then retaining terms by raw p-value alone. GSEA rejects
  duplicate or non-finite ranked identifiers unless an explicit collapse policy
  is selected, numeric formula RHS values require an explicit `analysis`, and
  scientifically relevant mapping/q-value/invalid-p warnings remain visible.
- `sn_enrich()` no longer advertises the AutoZyme fgsea patch around
  clusterProfiler 4.20 GSEA calls, because that release executes enrichit and
  does not enter the fgsea namespace. The clusterProfiler GSON-cache scope and
  its acceleration provenance remain attached to stored results.
- `sn_run_cluster()` no longer requires `HGNChelper` for ordinary Seurat
  clustering or bundled blocked-gene signatures; the dependency is checked
  only when custom blocked gene symbols require validation.
- `sn_run_cluster()` reuse and persistent checkpoint signatures now include a
  blockwise SHA-256 of the selected expression layer plus relevant batch, HVG,
  regression, and integration-label metadata. Calls with the same names and
  dimensions but changed values no longer reuse stale normalization, PCA,
  graphs, or checkpoints.
- Clustering now selects the SNN graph created or overwritten by the current
  `FindNeighbors()` call instead of falling back to an unrelated pre-existing
  `_snn` graph. Seurat CCA/RPCA also honor a custom
  `integration_control$new.reduction` through downstream graph construction and
  provenance, while restoring the original split count layers exactly.
- Log-normalized clustering now records the actually removed blocked HVGs, and
  resolution changes no longer invalidate HVG/PCA/graph reuse when rare-feature
  selection is disabled.
- `sn_normalize_data(method = "scran")` and scran-backed
  `sn_run_cluster()` workflows now accept BPCells-backed count layers. The
  selected layer is materialized directly as a sparse `dgCMatrix` at the scran
  compatibility boundary, while the returned Seurat object's original BPCells
  counts layer remains on disk.
- Multi-method `sn_run_cluster()` calls now honor the existing
  `run_tsne = FALSE` default. UMAP remains the default projection and t-SNE is
  created only when explicitly requested with `run_tsne = TRUE`.
- `sn_plot_feature()` now selects its requested assay through the local Seurat
  object's default assay instead of forwarding the removed `assay` argument to
  `Seurat::FeaturePlot()`, restoring compatibility with Seurat 5.4.

- Coralysis integration now accepts BPCells-backed Seurat v5 assays by
  materializing only the selected integration features as a sparse
  `dgCMatrix` before `Coralysis::PrepareData()`. The complete expression layer
  is never converted to a dense matrix. Coralysis now defaults to one worker
  instead of `threads = 0` (all available workers), preventing forked R workers
  from multiplying the in-memory sparse matrix on large datasets; callers can
  still opt into more workers through `integration_control$icp_args$threads`.
- scVI, scANVI, and totalVI integration now export the requested `layer`
  instead of always reading `counts`. The selected source layer is recorded in
  backend configuration and `object@misc$integration`, so comparisons using
  inputs such as `decontaminated_counts` are auditable. scPoli follows the same
  explicit layer contract, while BBKNN derives its PCA from that selected layer.
- Python expression inputs remain SciPy sparse matrices. scPoli latent export
  densifies only one bounded neural-network minibatch at a time, totalVI keeps
  its protein matrix sparse, and the generic scArches fallback uses sparse
  scaling and truncated PCA rather than materializing a full expression matrix.
- Automatic Seurat acceleration now scopes and logs the actual operation, such
  as `runpca`, instead of reporting unrelated active `seurat_joinlayers` and
  `seurat_merge` patches for every wrapped call.

# Version 0.3.0

Released 2026-08-01.

### Added

- Added a reproducible real-data AutoZyme benchmark and pkgdown article for
  `seurat_merge`, `seurat_joinlayers`, and `scdblfinder`. The benchmark runs
  baseline and accelerated conditions in fresh R processes against the
  validated 2,000-cell Kotliarov PBMC fixture, records operation time and
  whole-worker peak RSS, requires exact outputs and patch rollback, and ships
  both raw JSON evidence and a compact CSV summary.
- Added `sn_set_layer_backend()` as the bidirectional Seurat layer-storage
  interface. Selected layers can be written and rebound to BPCells or safely
  materialized as in-memory `dgCMatrix` objects without changing assay/layer
  names. Count-like integer layers use BPCells `uint32_t` storage automatically,
  while `sn_convert_bpcells()` remains available as a compatible one-way
  wrapper.
- Added `sn_build_result_bundle()`, `sn_validate_result_bundle()`, and
  `sn_export_result_bundle()` as the package-owned,
  `shennong.dev/analysis-result-bundle/v1` JSON handoff boundary. Bundles carry
  one validated canonical result, immutable input
  identifier/revision/SHA-256 references, package and execution provenance,
  and candidate artifact roles/digests without service calls or credentials.
- Added `sn_audit_results()` and `sn_upgrade_results()` to inspect and migrate
  stored analytical results to the canonical `1.0.0` result contract. The
  audit distinguishes valid, safely upgradeable legacy, invalid, and
  intentionally out-of-contract runtime artifacts, and reports unknown
  top-level `object@misc` payloads as `unregistered`, without mutating objects.
- Added the direct `sn_run_scissor()` workflow and `sn_plot_scissor()` views for
  bulk-phenotype-guided cell selection. Unified results retain all-cell
  coefficients, state and sample summaries, model/correlation evidence, and
  optional bootstrap reliability output while the existing
  `sn_prioritize_states(method = "scissor")` entry point remains compatible.
- Expanded `sn_run_survival()` into a complete per-feature survival workflow
  with adjusted Cox models, Kaplan-Meier curves, log-rank tests, risk tables,
  cumulative hazards, model-performance statistics, proportional-hazards
  tests, and scaled Schoenfeld residuals. `sn_plot_survival()` now renders
  forest, Kaplan-Meier, risk-table, cumulative-hazard, scaled-residual, and
  proportional-hazards-test views.
- Added guarded AutoZyme acceleration support through
  `sn_check_autozyme()`, `sn_enable_autozyme()`, `sn_disable_autozyme()`, and
  `sn_with_autozyme()`. The strict automatic set now covers CellChat, NicheNetR,
  clusterProfiler, fgsea, Seurat, tradeSeq, and WGCNA. Shennong checks these
  patches lazily only when a compatible workflow needs them, then restores the
  pre-call patch state after success or error. Automatic activation normally
  requires the pinned AutoZyme build, an installed upstream package, and an
  exactly validated upstream version. Seurat scopes deliberately admit
  version-label drift behind runtime structure guards, while scDblFinder uses
  exact target-function fingerprints. Missing dependencies are skipped safely
  and approximate patches are never automatic. Result-producing scopes record the
  patches that were active for the compatible call in provenance before the
  automatic scope is restored.
- Added a local-only real-public-data harness under `scripts/real-data/`. Its
  four logical bundles cover Kotliarov PBMC CITE-seq, a GSE72056/TCGA-SKCM
  single-cell-to-bulk melanoma bridge, Hermann spermatogenesis spliced and
  unspliced counts, and two 10x Visium lymph-node sections. Source metadata,
  deterministic preparation, validation, and function-to-article coverage are
  tracked; raw and prepared data remain ignored and are never packaged.
- Added runtime article coverage and AutoZyme benchmark runners. Runtime
  coverage records which declared analysis and visualization functions were
  actually observed while real-data articles rendered, keeps optional extended
  backends explicit, and rejects network access. The benchmark compares
  baseline and scoped accelerated CellChat/WGCNA execution and checks output
  equivalence plus patch restoration.
- Added RegVelo as a managed `sn_run_velocity(method = "regvelo")` backend.
  Shennong accepts regulator-target edge tables, named target-by-regulator
  matrices, or CSV priors; the unified velocity result retains RegVelo latent
  time, projected vectors, transition evidence, trained-model metadata, and
  the H5AD artifact used by CellRank.
- Added a read-only stdio MCP server with `sn_mcp_server()` and
  `sn_mcp_server_config()`. Agents can discover methods, inspect exact installed
  function help, and read bundled workflow guides without arbitrary R execution
  or file mutation. A matching `use-shennong-mcp` Agent Skill is shipped and
  installable through `sn_install_codex_skill()`.

### Changed

- AutoZyme integration now pins the `zerostwo/autozyme` fork and scopes its
  `scDblFinder`, UCell, LISI, and SoupX accelerators around the corresponding
  package calls. The former Shennong-local scDblFinder and SoupX patches and
  SoupX C++ kernels were removed so the fork is the single provider for those
  implementations.
- `sn_remove_ambient_contamination()` now calls the standalone
  `decontX::decontX()` API directly instead of `celda::decontX()`. Its
  `method = "auto"` default detects a single ADT, protein, or CITE assay and
  routes it to `decontX::decontPro()`; use `method = "decontpro"` and
  `assay =` to make the CITE-seq choice explicit. The standalone decontX
  AutoZyme hook is scoped only around `decontX::decontX()`.
- `sn_remove_ambient_contamination(method = "decontx")` now defaults to
  decontX's native clustering instead of unconditionally running a complete
  `sn_run_cluster()` workflow first. The new `cluster_backend` argument selects
  `"native"` or `"shennong"`; explicit `cluster` labels still take precedence.
  SoupX keeps its method-specific `"shennong"` default because SoupX has no
  native clustering implementation.
- `sn_find_doublets()` now exposes the same `cluster_backend` choice. Its
  default remains scDblFinder's native automatic clustering; selecting
  `"shennong"` explicitly obtains assignments from `sn_run_cluster()`.
- Automatic Seurat AutoZyme scopes now admit registered non-approximate
  patches across upstream version labels and guard execution by the patch's
  runtime input/structure checks. `sn_run_cluster()` now scopes its neighbor
  and cluster stages as well as normalization, HVG, scaling, and PCA calls, so
  every currently supported Seurat fast path can participate. BPCells-backed
  objects still suppress the broad Seurat patch. The scDblFinder patch likewise
  uses exact function fingerprints rather than the package version string and
  fails closed when any targeted implementation changes.
- Enrichment AutoZyme scopes now activate clusterProfiler and fgsea
  independently, so failure of one patch cannot roll back the other. Shennong
  bundles a clusterProfiler 4.20/enrichit-compatible exact cache for repeated
  `get_GO_data()` requests, while ranked GSEA calls use AutoZyme's three-target
  fgsea patch with relaxed version-label gating.
- CellChat and call-safe NicheNetR communication workflows now request relaxed
  AutoZyme version gating. Their runtime input guards and upstream fallbacks
  remain active, but installed package version labels no longer prevent the
  compatible patch from starting.
- Successful automatic AutoZyme scopes now emit an INFO log naming the patches
  enabled for the current workflow call. The message is printed only after
  activation and state verification; it does not claim that every guarded
  internal fast path will accept the current input.
- Shennong now vendors the exact-scoped AutoZyme `seurat_merge` and
  `seurat_joinlayers` patches before their official AutoZyme release. They
  replace only `SeuratObject::merge.Assay5` and
  `SeuratObject::JoinLayers.Assay5`, verify bundled source fingerprints before
  registration, preserve captured-upstream fallback, and can be enabled with
  `sn_enable_autozyme(c("seurat_merge", "seurat_joinlayers"))`.
- `sn_find_doublets()` now lazily scopes Shennong's bundled exact AutoZyme
  `scdblfinder` patch around the validated default call.
  The fast path is body/formals guarded, limited to exact `dgCMatrix`
  inputs with 1--33,000 cells and positive finite library sizes, and restores
  the pre-call AutoZyme state. Non-default, BPCells-grouped, oversized, drifted,
  or otherwise unsupported calls continue through captured upstream code.
  The vendored patch now uses a non-materialized sparse transpose operator for
  exact generic IRLBA PCA and sparse-only artificial-doublet Poisson resampling;
  expanded inputs above 50,000 columns retain the upstream PCA path.
- `sn_remove_ambient_contamination(method = "soupx")` now lazily activates the
  exact-scoped AutoZyme `soupx` patch for `SoupX::adjustCounts()`. The patch and
  its native C++ kernels are now bundled by Shennong and registered through
  AutoZyme at runtime, so the official AutoZyme package no longer needs to ship
  a SoupX patch. The installed patch also includes its validated scope and the
  finalized, equivalence-passing PBMC 10k/20k benchmark summary. Shennong
  retains SoupX's stochastic integer-rounding behavior, restores the pre-call
  patch state, and falls back to upstream SoupX when AutoZyme is absent or
  incompatible.
- AutoZyme automatic activation is now lazy and scoped to compatible Shennong
  workflow calls rather than package loading or persistent process mutation. Set
  `options(shennong.autozyme = FALSE)`, `AUTOZYME_DISABLED=true`, or
  `AUTOZYME_DISABLE=true` to prevent automatic scopes. These opt-outs do not
  deactivate manually active patches and do not affect explicit
  `sn_enable_autozyme()` or `sn_with_autozyme()` calls. BPCells-backed Seurat
  layers never enter the Seurat fast patch because that path can materialize an
  on-disk matrix as `dgCMatrix`; an already active Seurat patch is suspended for
  the BPCells-backed call and restored afterward. Eligible NicheNetR use remains
  limited to its validated call-safe dense-prior path; analytical errors are
  never retried as unaccelerated work. Automatic scopes also restore the
  caller's `future.globals.maxSize` option after AutoZyme is loaded. This Seurat
  guard does not make every Shennong backend BPCells-native: backends such as
  CellChat and tradeSeq may still require a controlled sparse materialization or
  an aggregated input that fits memory.
- `sn_initialize_seurat_object()` now accepts all BPCells `IterableMatrix`
  subclasses and preserves their on-disk counts backend. `sn_find_doublets()`
  now handles BPCells-backed Seurat layers by materializing grouped samples
  independently; `ncores = 1` bounds count-matrix materialization to one sample
  at a time, while ungrouped BPCells input fails early with memory-safe
  guidance. BPCells-to-`dgCMatrix` conversion no longer passes through a dense
  matrix.
- `sn_list_methods(available = TRUE/FALSE)` now filters against the caller's
  requested value instead of resolving the argument through dplyr's data mask.
- Named Scissor Cox phenotypes explicitly preserve bulk-sample identifiers
  after alignment; the regression contract now verifies both reordered values
  and retained names.
- Analytical results now use `schema_version = "1.0.0"` and a canonical
  data-frame `tables$primary` across registered and generic result stores.
  Existing table aliases and specialized getters remain available as
  compatibility views, while non-analytical caches and backend manifests remain
  explicitly classified as artifacts.
- Reference label transfer now stores cell-level predictions as a canonical
  `annotation` result while retaining the existing compact
  `object@misc$label_transfer` manifest as a registered compatibility artifact.
- `scripts/build-pkgdown.R --real` now validates the complete local public-data
  matrix before evaluating article chunks with real tables and plots. It never
  downloads data; `--data-root` selects an existing local matrix, `--extended`
  opts into available external backends, and `--full --real` is the clean
  real-output release gate.
- `sn_find_de()` is now the common differential-expression entry point. It
  retains all Seurat marker, contrast, and pseudobulk behavior while
  automatically dispatching matrix, list, and `SummarizedExperiment` inputs to
  the standalone bulk engine. `sn_find_bulk_de()` remains a compatible wrapper.

### Removed

- Removed Shennong's data-distribution layer: `sn_load_data()`,
  `sn_list_datasets()`, `sn_download_zenodo()`, and `sn_upload_zenodo()` are no
  longer exported, and the bundled `pbmc_small` / `pbmc_small_raw` analysis
  datasets are no longer shipped. Dataset discovery, download, caching, and
  publication now belong to `ShennongData`; Shennong consumes materialized
  matrices, paths, and Seurat objects. The runtime reference assets used for
  annotation, species mapping, gene filtering, and signatures remain bundled.

### Fixed

- Declared `DOSE` in `Suggests` because the figure-engine tests conditionally
  load its S4 classes. This keeps the backend optional while satisfying
  `R CMD check`'s static test-dependency audit.
- Managed scVelo preprocessing now explicitly normalizes fractional splicing
  estimates and log1p-transforms the expression matrix before Seurat-flavor
  HVG selection. The choices are retained in velocity parameters and backend
  manifests, and CellRank computes its eigendecomposition before automatic
  macrostate selection when no state count is supplied.
- `sn_annotate_de_features()` now synchronizes the legacy `table` alias with
  `tables$primary`, preventing retrieval of stale differential-expression
  evidence after feature annotation.
- Result migration now rejects unsupported or future schema versions without
  rewriting them, uses exact field lookup so misspelled aliases cannot pass
  validation, and materializes the registered QC compatibility views when a
  canonical generic result is stored.
- `sn_run_scissor()` now uses the requested assay layer consistently for
  variable-feature selection, PCA, graph construction, and the Scissor backend;
  requesting an unavailable layer fails explicitly instead of silently falling
  back to counts.
- Stored spatial communication now retrieves the registered
  `cell_communication` result type instead of looking under an unreachable
  `communication` type.
- CellChat execution now encodes backend-unsafe group names such as the numeric
  label `"0"` and decodes sender/receiver labels in returned tables, preserving
  the user's original biological labels.
- Seurat label transfer now requests the table-returning `TransferData()` form.
  This avoids Seurat 5 returning a query object that was then incorrectly
  coerced to a data frame.
- `sn_plot_feature(assay = ...)` now uses the requested assay when calculating
  shared feature limits, so ADT and other non-default-assay features render
  without changing `DefaultAssay()`.
- Result-aware plotting uses exact nested-field lookup, preventing similarly
  named fields such as `tables_backup` from being partially matched as the
  canonical `tables$primary` result.
- Enrichment plotting now accepts the native S4 `enrichResult`, `gseaResult`,
  and `compareClusterResult` objects returned by clusterProfiler/DOSE, so
  standalone real-data enrichment results can be plotted without manual
  coercion.
- The BayesPrism adapter now restores the R session temporary directory if a
  backend removes it, preventing successful deconvolution from breaking later
  knitr graphics and file output in the same session.
- Cross-method communication concordance now counts finite rank pairs before
  computing Spearman correlation, and LIANA/NATMI `prod_weight` and
  `edge_specificity` output is recognized as quantitative evidence.
  LIANA/CellChat consensus therefore uses both backend rankings when available
  and retains an explicit `NA` diagnostic when overlap is not rankable instead
  of failing after both backends complete successfully.

### Performance

- On the current maintainer host, the reproducible real-data benchmark recorded
  CellChat at 1.264s baseline versus 0.126s accelerated (10.0317x), with
  identical 16-interaction output, and WGCNA at 7.781s versus 0.073s
  (106.5890x), with identical modules, maximum eigengene difference `1.71e-07`,
  and maximum trait-correlation difference `8.92e-08` (tolerance `1e-05`).
  These are measured host-specific results, not package-level performance
  guarantees; the ignored JSON artifact retains the exact inputs, timings,
  versions, source revision, and equivalence checks.

- Daily R package checks no longer repeat the full test suite already run by
  coverage; release audits can opt back in through `workflow_dispatch`.
- pkgdown push deployments preserve the previous site, rebuild only changed
  pages, and skip example execution. A scheduled or manually requested full
  build remains the clean release-quality gate, while
  `scripts/build-pkgdown.R` provides the same incremental/full split locally.
- Documentation-only changes avoid unnecessary coverage and package-check
  runs, and superseded runs are cancelled per branch.

# Version 0.2.0

Released 2026-07-15.

### Fixed

- The multimodal roadmap test now scopes its `sn_run_cluster()` mock to one
  expression, preventing the four-formal test double from leaking into later
  installed-package clustering tests.
- `sn_run_cluster()` now dispatches its long-standing tail controls through a
  compact, allowlisted compatibility wrapper. Named calls, positional tail
  calls, explicit `NULL` values, and defaults are preserved in one named
  argument bundle, avoiding long-formal dispatch on the R 4.6 development
  runner.
- Scissor and Symphony validate required user inputs before checking optional
  backend installations, so dependency-independent contract errors remain
  testable on minimal installations.
- The pkgdown deployment workflow now installs UCell before evaluating the
  program-scoring article, so the documented default `method = "ucell"`
  backend is available in clean GitHub Actions runners.
- The data-I/O article no longer performs unconditional Zenodo downloads during
  pkgdown deployment; its network-backed examples now require the explicit
  `SHENNONG_RUN_NETWORK_VIGNETTES=true` opt-in.
- `sn_plot_dot()` now renders named feature lists as Seurat's free-width marker
  facets without adding an incompatible fixed coordinate ratio; ordinary
  feature vectors retain the existing fixed-coordinate layout.
- `sn_plot_feature(raster = TRUE)` now falls back to vector point layers when
  the optional `ggrastr` package is unavailable, instead of asking Seurat to
  perform ordered rasterization and failing during minimal-dependency checks.
- `sn_load_data(backend = "api")` now uses the current ShennongData 0.2
  resource-client contract (`sn_connect()`, `sn_load_data()`, `sn_assay()`,
  and `collect()`) rather than the retired schema endpoint interface.
- `sn_run_cluster(normalization_method = "sctransform")` now applies
  `block_genes` to SCTransform-selected HVGs before PCA, preventing default
  ribosomal, mitochondrial, heat-shock, immunoglobulin, TCR, and pseudogene
  signatures from dominating SCT-based clustering unless explicitly forced via
  `hvg_features`.
- `sn_run_cluster(block_genes = ...)` now resolves each entry independently, so
  bundled signature queries such as `cellCycle.G2M`, `cellCycle.G1S`, `ribo`,
  `mito`, `heatshock`, and `pseudogenes` can be mixed with custom gene symbols
  and are removed from the final stored HVG set.
- `sn_write()` now creates missing parent directories before dispatching both
  rio and custom writers, so nested `.qs`, `.h5ad`, `.h5`, and BPCells outputs
  no longer fail only because the containing directory does not exist.
- `sn_write()` now auto-installs missing optional writer dependencies by
  default for custom formats such as `.qs2`, `.h5ad`, `.h5`, and BPCells.
  `.qs` output now attempts to install `qs` from the GitHub remote
  `qsbase/qs`, while `.qs2` remains the recommended new serialization format.
- `sn_run_cluster()` now uses `batch` as the only integration metadata
  argument; the older `batch_by` alias has been removed from this entry point.
- GitHub Actions now install `leidenbase` anywhere evaluated clustering
  examples or tests can request `cluster_algorithm = "leiden"`.
- `sn_run_cluster(cluster_algorithm = "leiden")` now checks for `leidenbase`
  before calling Seurat and auto-installs it through `sn_install_dependencies()`
  by default. Set `auto_install = FALSE` to keep the previous fail-fast
  behavior.
- `sn_run_cluster()` now skips redundant Seurat PCA work for integration
  backends that do not consume a Seurat PCA reduction, including Coralysis,
  scVI/scANVI, totalVI, MMoCHi, and CITE-seq protein-only Coralysis/MMoCHi
  paths.
- Native Coralysis clustering now stores the trained
  `SingleCellExperiment` under `object@misc$coralysis` by default, so the
  returned object can be used directly as a label-transfer reference. Set
  `integration_control = list(store_sce = FALSE)` only for clustering-only
  runs where the Coralysis reference object is not needed.
- `sn_run_cluster()` now calls configurable Seurat steps such as
  `FindClusters()` and `RunUMAP()` through symbolic object calls rather than
  `do.call(object = object, ...)`, preventing `object@commands` entries from
  serializing the full Seurat object into command-history call strings.
- Grouped HVG selection in `sn_run_cluster()` now skips cells with missing or
  empty `hvg_group_by` labels and runs temporary per-group HVG calls on the
  selected assay only, avoiding spurious ADT assay removal messages on
  multimodal objects.
- SCTransform workflows in `sn_run_cluster()` and `sn_normalize_data()` now
  temporarily raise `future.globals.maxSize` from the detected system memory and
  object size before calling Seurat, avoiding the default 500 MiB future export
  limit on large objects while restoring the caller's option afterwards.
- Package source builds now exclude local benchmark outputs and OmnipathR test
  log directories, keeping `R CMD build` tarballs small and free of generated
  validation artifacts.
- Package source builds and git status now ignore local `.codegraph/`,
  `.agents/`, and `.codex/` agent artifacts, keeping analysis state out of
  source builds and commits.
- Signature catalog add/update/delete helpers now edit the packaged Shennong
  signature tree directly, so custom signature maintenance no longer requires
  the optional upstream `SignatuR` package at runtime.
- The `.qs` writer parent-directory regression test now mocks optional writer
  dependency detection completely, avoiding live GitHub installation attempts
  during local tests.
- Local CIBERSORTx execution now uses structured `system2()` calls instead of a
  pasted shell command, and dry-run/stored command artifacts redact CIBERSORTx
  account email and token values.
- `scripts/check-prepush.R` now reports per-step timings, defaults local checks
  to `_R_CHECK_FORCE_SUGGESTS_=false`, skips duplicate test execution inside
  `R CMD check` after a successful full `test_local()` pass, and adds a
  `--quick` edit-loop mode.
- Stored-result writes through Shennong's shared `object@misc` helper now pass
  through a central collection registry and schema validator. Malformed stored
  DE/enrichment/interpretation/deconvolution/Milo/communication/regulatory/QC
  entries now fail early with an actionable schema error, and `sn_list_results()`
  also reports stored QC assessments.
- Enrichment wrappers now muffle known benign `clusterProfiler`/`fgsea` warnings
  produced by tiny deterministic examples, keeping local test and package-check
  output warning-clean while preserving hard errors, including the current
  `enrichit` qvalue fallback when a tiny result cannot estimate q values.

### Added

- Added the explicit `sn_run_multimodal()` CITE-seq entry point while keeping
  annotation, trajectory, and fate backend adapters internal to their unified
  workflow APIs.
- Added direct optional Monocle 3 trajectory inference plus standardized
  Palantir runner/result and scCODA/pertpy runner/result adapters. All shipped
  method-registry entries are now implemented, and scCODA retains biological
  samples as the inferential unit.
- Added a publication figure engine with generic screen/column/page/slide
  profiles, automatic size/point/raster/layout/pagination specifications,
  structured figure QA, deterministic PDF/SVG/TIFF/PNG export, and reproducible
  figure bundles containing source data, specs, sessions, manifests, and
  checksums.
- Added result-aware DE/GSEA/enrichment figures plus QC threshold, doublet,
  ambient correction, HVG, elbow, cluster tree, resolution sweep, integration,
  and reference-projection diagnostics. Core dimensional, feature, dot,
  heatmap, violin, box, bar, composition, Milo, and bulk plots now carry figure
  specifications without changing their native plot classes.
- Added a standalone bulk transcriptomics mainline for matrix/list/
  `SummarizedExperiment` inputs: sample QC, design-aware edgeR/DESeq2/limma/
  dream differential expression, pathway scoring, WGCNA module-trait analysis,
  Cox survival models, clinical associations, and result-aware plots.
- Added `sn_run_spatial()` as a dispatcher plus explicit spatial feature,
  domain, neighborhood, deconvolution, mapping, integration, and communication
  entry points. Local Moran's I and memory-bounded KNN workflows retain spatial
  graphs, permutation evidence, co-occurrence, coordinates, and diagnostics;
  nnSVG/BANKSY run when installed and heavyweight alternatives use explicit
  result adapters.
- Added spatial coordinate, feature, domain, SVG, neighborhood,
  deconvolution, and distance-aware communication plots that preserve tissue
  aspect ratio.
- Added `sn_run_velocity()` and `sn_plot_velocity()` with a managed scVelo
  pixi backend for spliced/unspliced preprocessing, projected velocity vectors,
  velocity pseudotime/confidence, transition edges, and retained H5AD evidence.
- Added `sn_run_fate()` and `sn_plot_fate()` with a managed CellRank GPCCA
  backend for terminal-state discovery, fate probabilities, optional lineage
  drivers, metadata storage, and explicit terminal-state controls.
- Added `sn_discover_programs()` and `sn_plot_discovered_programs()` for
  multi-restart NMF with reconstruction/stability diagnostics plus explicit
  cNMF and Hotspot result adapters. Discovered gene weights and per-cell
  activities use the shared result contract and can be stratified by metadata.
- Added `sn_run_grn()` and `sn_plot_regulon()` for real GENIE3 inference and
  explicit pySCENIC, legacy SCENIC, and GRNBoost2 adapters. Results standardize
  regulatory edges, regulons, per-cell activity, and transparent group
  specificity without hiding external motif databases or Python runtimes.
- Added `sn_run_cnv()` and `sn_plot_cnv()` as the unified inferCNVpy/CopyKAT
  workflow. Stored results now include reference-calibrated malignancy scores,
  malignant calls, subclones, sample summaries, chromosome-level CNV, optional
  CNV UMAP coordinates, and CNV-expression associations.
- Added `sn_metabolic_signatures()`, `sn_run_metabolism()`, and
  `sn_plot_metabolism()` for curated UCell/GSVA/ssGSEA/mean pathway activity,
  sample-level differential metabolism, scMetabolism, and standardized
  scFEA/Compass result adapters.
- Expanded `sn_run_cell_communication()` into a multi-backend communication
  workflow for LIANA, CellChat, CellPhoneDB, NicheNet, and MultiNicheNet. All
  backends now map to a shared ligand-receptor schema with method concordance,
  consensus ranks, sample-level expression evidence, condition contrasts,
  ligand-target links, retained backend artifacts, and a reserved spatial
  distance field. Added bubble, heatmap, network, chord, river,
  ligand-target, and differential-communication plots.
- Added `sn_test_abundance()` as the stable differential-abundance entry point
  for sample-level Propeller, transparent sample-label permutation, and
  neighborhood-level Milo. Results include standardized effects, adjusted
  significance, completed sample proportions, sample contributions, design
  data, backend evidence, and permutation nulls where applicable.
- Added `sn_prioritize_states()` for sample-aware Augur-style held-out
  separability, explicit bulk-input Scissor selection, and RareQ topology
  discovery followed by sample-level phenotype association. State rankings,
  cell scores, uncertainty, null distributions, and sample contributions use
  the shared result contract; `sn_plot_abundance()` and
  `sn_plot_state_priority()` render stored results.
- Added `sn_run_trajectory()` with Slingshot lineage inference and optional
  tradeSeq dynamic-gene, branch-pattern, differential-end, convergence, and
  fitted-trend outputs. Per-lineage pseudotime/probability, principal curves,
  terminal states, topology, diagnostics, and provenance use the unified
  analysis-result contract, while primary lineage/pseudotime are also added to
  Seurat metadata.
- Added result-aware trajectory, pseudotime, lineage-probability,
  dynamic-heatmap, gene-trend, and branch-comparison plots through
  `sn_plot_trajectory()`, `sn_plot_pseudotime()`,
  `sn_plot_lineage_probability()`, `sn_plot_dynamic_heatmap()`,
  `sn_plot_gene_trend()`, and `sn_plot_branch_comparison()`.
- Added `sn_score_programs()` with UCell (default per-cell), AUCell, GSVA,
  ssGSEA, and sparse-aware mean-expression backends. It records signature
  feature coverage, stores long-form scores in the unified result contract,
  and adds cell-level scores to Seurat metadata without silently converting
  large sparse matrices for GSVA.
- Added `sn_test_programs()` for condition comparisons that aggregate to the
  sample/patient level before inference when `sample_by` is supplied, plus
  `sn_plot_program_activity()` and `sn_plot_program_heatmap()` for stored
  score results.
- Added `sn_run_annotation()` as the stable annotation entry point with
  marker-only consensus plus optional SingleR, CellTypist, Seurat transfer,
  Symphony, scmap, and scANVI backends. It stores cell- and cluster-level
  hierarchical labels, calibrated confidence/margin, runner-up labels, marker
  support/conflicts, reference coverage, raw backend predictions, diagnostics,
  and provenance in the unified result contract.
- Added `sn_annotation_consensus()`, `sn_annotation_confidence()`,
  `sn_map_cell_ontology()`, and `sn_review_annotation()` for transparent
  evidence aggregation, a versioned bundled Cell Ontology snapshot, and
  explicit low-confidence review without allowing LLM output to overwrite
  computational labels.
- Added result-aware annotation confidence, marker-evidence, and confusion
  plots through `sn_plot_annotation_confidence()`,
  `sn_plot_annotation_markers()`, and `sn_plot_annotation_confusion()`.
- Added a shipped method registry with `sn_list_methods()` and
  `sn_method_status()`. It records runtime, optional dependency, installation
  action, input/output contract, CPU/GPU expectations, citations, and whether
  each current or planned backend is actually implemented and available.
- Added the versioned generic analysis-result contract and
  `sn_store_result()`, `sn_get_result()`, `sn_delete_result()`, and
  `sn_validate_result()`. Existing registered DE, enrichment, Milo,
  communication, deconvolution, regulatory, QC, and interpretation results are
  upgraded on write/read without breaking their specialized getters.
- `sn_list_results()` now accepts an optional `type` filter and includes new
  generic result types stored under `object@misc$analysis_results`.
- `sn_load_data()` can open lazy Shennong Data Server resources with
  `backend = "api"`, select assay/layer views through `api_args`, and
  materialize explicitly with `lazy = FALSE`.
- `sn_convert_bpcells()` now converts selected Seurat assay layers to
  BPCells-backed matrix directories and rebinds those layers in the returned
  object, helping large count or normalized-expression layers stay on disk.
- `sn_annotate_de_features()` now flags stored marker/DE genes that encode
  transcription factors, cell-surface or plasma-membrane proteins, cytokines,
  and chemokines. It can annotate direct DE tables or store annotated tables
  back under `object@misc$de_results` for retrieval with `sn_get_de_result()`.
- `sn_prepare_label_transfer_reference()` now creates compact reference
  objects for `sn_transfer_labels()`. For native Coralysis, it keeps only the
  trained Coralysis models, PCA model, feature names, and selected labels while
  dropping reference assays, reductions, and stored joint probabilities. For
  Seurat, scANVI, and scArches workflows, it returns a slim Seurat reference
  with selected assay layers and labels.
- `sn_run_cluster()` now accepts `umap_control`, a named list of
  `Seurat::RunUMAP()` arguments such as `n.neighbors`, `min.dist`, `spread`,
  `metric`, `seed.use`, and `reduction.name`, so users can tune embedding
  geometry without rerunning Coralysis, neighbor graph construction, or
  clustering.
- `sn_run_cluster(modality = "cite_seq")` now runs Seurat CITE-seq weighted
  nearest-neighbor clustering from paired RNA and ADT assays, including ADT CLR
  normalization, ADT PCA, `weighted.nn` / `wsnn` graph construction, clustering,
  and `wnn.umap` embedding.
- `sn_run_cluster(modality = "cite_seq", multimodal_method = ...)` now exposes
  a unified CITE-seq backend selector. In addition to Seurat WNN, users can run
  native Coralysis on the ADT protein assay, scvi-tools totalVI on paired RNA
  and ADT counts, or MMoCHi ADT landmark registration through a managed pixi
  backend.
- `sn_run_cluster(modality = "cite_seq", multimodal_method = "mmochi")` now
  supports single-sample CITE-seq runs with `batch = NULL` by passing a constant
  internal batch key to the MMoCHi backend instead of requiring a user-supplied
  batch column.
- `sn_run_cluster()` now records reusable stage signatures for normalization,
  cell-cycle scoring, HVG/rare-feature selection, PCA, integration, neighbor
  graph construction, clustering, and UMAP. Re-running on its own output
  reuses matching stages by default, so changing only `resolution` starts at
  clustering, changing `integration_method` starts at integration, and changing
  HVG controls starts at feature selection. Use `reuse = FALSE` or
  `rerun_from = "hvg"` / `"integration"` / another stage to force recompute.
- `sn_run_cluster()` now exposes key `Seurat::FindClusters()` controls,
  including `cluster_algorithm = "leiden"` / `"louvain"` / `"slm"`, custom
  cluster metadata names, random seed, start/iteration counts, singleton
  handling, and Leiden method/objective options.
- `sn_calculate_variance_explained()` now ranks metadata variables such as
  platform, study, tissue, and sample by weighted embedding variance explained,
  with single-variable and partial multi-variable modes for batch-effect
  diagnostics.
- `sn_calculate_roe()` now computes observed-over-expected enrichment for
  categorical composition tables from Seurat metadata or data frames, with long
  table output by default and optional matrix output for heatmaps.
- `sn_transfer_labels()` now wraps reference mapping with a query-first API
  for pipe-friendly workflows. The default Seurat anchor workflow is retained,
  and `method = "coralysis"` now projects queries onto Coralysis-trained
  references with native `Coralysis::ReferenceMapping()`. `method = "scanvi"`
  and `method = "scarches"` now provide semi-supervised scVI-family label
  transfer through the pixi-managed scverse backend.
- `sn_upload_zenodo()` now uploads reusable data files to Zenodo through
  `zen4R`, with a simple draft-first interface and an automatically uploaded
  Shennong manifest that records dataset version, package version, file sizes,
  md5 checksums, and sha256 checksums for reproducible reuse.
- `sn_download_zenodo()` now downloads reusable files from public Zenodo
  records without requiring a token, with optional token support for
  restricted records. `sn_load_data()` now uses this download layer and accepts
  multiple example datasets such as `dataset = c("pbmc1k", "pbmc3k")`; filtered
  datasets are returned as one merged Seurat object and raw datasets as a named
  list of sparse matrices.
- `sn_list_datasets()` now lists sample-level datasets available through the
  Shennong public Zenodo collection. `sn_load_data()` can load those samples
  from the `shennong_index.json` layout in Zenodo record `20044788`, downloading
  the study ZIP, extracting the requested sample's filtered/raw H5 or Cell
  Ranger metrics file, and optionally validating extracted files against
  `manifest.tsv`.
- `sn_simulate()` now provides a method-based simulation entry point.
  `method = "scdesign3"` wraps `scDesign3::scdesign3()` for Seurat or
  SingleCellExperiment inputs and can return simulated counts as a Seurat
  object, SingleCellExperiment, sparse matrix, or raw scDesign3 result.
  `sn_simulate_scdesign3()` remains as the backend-specific wrapper.
- `sn_plot_heatmap()` now draws focused heatmaps for user-selected genes, with
  cell-level and group-averaged modes, optional grouping/splitting, default
  rasterization, hidden cell names/ticks, 8 pt group labels, Paired group-bar
  colors, and automatic scaling of requested features when needed.
- `sn_run_cell_communication()` now wraps real cell-cell communication
  backends: CellChat, NicheNet (`nichenetr`), and LIANA. Results can be stored
  and retrieved with `sn_store_cell_communication()` and
  `sn_get_cell_communication_result()`.
- `sn_run_regulatory_activity()` now runs fast footprint-style activity
  inference with DoRothEA regulons or PROGENy pathway models through
  `decoupleR::run_ulm()`. Results can be stored and retrieved with
  `sn_store_regulatory_activity()` and `sn_get_regulatory_activity_result()`.
- `sn_run_cluster()` now accepts `integration_method` for batch workflows.
  In addition to the historical Harmony path, users can run Coralysis
  multi-level integration or Seurat layer integration with CCA/RPCA through
  the same clustering entry point.
- `sn_run_cluster()` now accepts `integration_method = "scvi"` and
  `"scanvi"`. These backends export selected count data to a pixi-managed
  scverse runtime under `~/.shennong/pixi/`, run the Python model, import the
  latent representation as a Seurat reduction, and then continue Shennong's
  neighbors/clustering/UMAP workflow. The scANVI path requires
  `integration_control = list(label_by = ...)`. Convenience wrappers
  `sn_run_scvi()` and `sn_run_scanvi()` expose the same workflows directly.
- New pixi helpers `sn_check_pixi()`, `sn_install_pixi()`,
  `sn_ensure_pixi()`, `sn_pixi_paths()`, `sn_list_pixi_environments()`,
  `sn_pixi_config_path()`, `sn_prepare_pixi_environment()`,
  `sn_call_pixi_environment()`, `sn_detect_accelerator()`, and
  `sn_configure_pixi_mirror()` expose the Python-runtime setup used by
  scVI/scANVI and future Python backends. Shennong now keeps pixi workspaces
  under `~/.shennong/pixi/`, renders package-bundled configs from
  `inst/pixi/`, can auto-install pixi when missing, selects CPU or CUDA pixi
  environments automatically, and can write China mirror configuration into
  the Shennong `PIXI_HOME`.
- Bundled pixi configs and command helpers are available for concrete Python
  method families including `scvi` (shared by scVI/scANVI), `scarches`
  (shared by scArches/scPoli), `infercnvpy`, `cellphonedb`, `cell2location`,
  `tangram`, `squidpy`, `spatialdata`, and `stlearn`. Use environment calls
  such as `sn_call_cell2location()` or analysis wrappers such as
  `sn_run_tangram()` to run commands inside those managed environments.
- `sn_run_infercnvpy(object = seurat_obj, ...)` now provides an object-level
  infercnvpy workflow: it exports the selected Seurat assay/layer with gene
  positions, runs infercnvpy in the managed pixi environment, and imports CNV
  metadata and optional CNV reductions back into the Seurat object.
- Packaged Python runner scripts now live with their pixi family configs under
  `inst/pixi/<family>/scripts/`. Object-level Seurat workflows are available
  for scArches/scPoli, CellPhoneDB, cell2location, Tangram, Squidpy,
  SpatialData, and stLearn through their corresponding `sn_run_*()` wrappers;
  method-specific Python settings can be supplied with `method_control`.
- `R CMD check` namespace diagnostics are tighter: optional `qs`/`qs2`
  serialization packages are now declared, and previous `ave`, `tail`,
  `target`, and `mor` code-analysis notes have been resolved.
- `sn_run_cluster()` now accepts `hvg_features`, a user-supplied feature list
  that is validated against the object and merged with internally selected
  HVGs and rare-aware features before scaling/PCA. This lets users force rare
  population marker genes into the clustering feature set when global HVG
  selection misses them.
- `sn_sweep_cluster_resolution()` now provides a formal resolution-sweep
  interface for empirically comparing candidate cluster counts across Seurat
  resolutions with metrics such as silhouette width, graph connectivity,
  cluster purity, clustering agreement, and optional ROGUE summaries.
- `sn_list_dependencies()` now reports the package's required and recommended
  R package surface with install status and expected source, and
  `sn_install_dependencies()` can install missing CRAN, Bioconductor, and
  GitHub dependencies in one step.
- `sn_list_10x_paths()` now scans a root directory for 10x Genomics outputs
  and can return `outs/` directories, filtered matrix paths, raw matrix paths,
  H5 files, or `metrics_summary.csv` paths. The default now returns `outs/`
  paths so the result can be passed directly to `sn_initialize_seurat_object()`.
  Returned vectors are now named with inferred sample identifiers.
- `sn_interpret_annotation()` now accepts `label_candidates` so sorted or
  enriched datasets can constrain annotation toward expected cell-type spaces
  such as `ILC1` / `ILC2` / `ILC3` instead of relying only on free-text
  background notes.
- `sn_interpret_annotation()` now supports `annotation_mode = "agentic"` for
  a two-stage workflow: broad lineage/state annotation followed by focused
  refinement on ambiguous or lineage-sensitive clusters. Annotation evidence
  now also exposes a `canonical_marker_snapshot` table so prompts can compare
  lineage-defining markers across clusters without relying only on top-ranked
  DE hits. The `ellmer` path now also uses native structured output for the
  final annotation table and can run a tool-assisted focused-comparison step
  before the refinement pass.
- Annotation heuristics and prompt priors are now more robust for ILC-rich
  datasets: blood ILC workflows now bias `KIT+ ILCP-like` over premature
  mature `ILC3` calls when the type-3 program is incomplete, mixed `T/NK` and
  `NK/ILC3` transitional states are called out more explicitly, and dominant
  hemoglobin programs can now surface as `erythroid contamination`.
- `sn_interpret_annotation()` and `sn_prepare_annotation_evidence()` now
  support `marker_selection = "specific"` and
  `enrichment_selection = "specific"` so annotation evidence can prefer
  cluster-restricted marker genes and pathway terms over generic top-ranked
  features.
- Annotation evidence now adds concise canonical lineage heuristic hints from
  known marker programs so the LLM can use deterministic guardrails such as
  `ILC2-like`, `KIT+ ILC-like`, `T-cell-like`, or `B-cell-like` when those
  programs are clearly supported.
- `sn_assess_qc()` now summarizes overall and per-sample QC status, reports
  current QC risk signals such as failed-QC fractions, doublet rates, and
  decontamination zero-count rates, and can compare a filtered object against a
  pre-filter reference to quantify low-quality-cell removal, doublet removal,
  and clean-cell retention. Reports can be stored under
  `object@misc$qc_assessments`.
- `sn_plot_composition()` now provides a composition-focused bar plot helper
  for grouped proportions, counts, QC pass/fail summaries, and similar
  categorical tables.
- `sn_compare_composition()` now compares sample-level composition between two
  groups and reports mean proportions, differences, log2 fold changes, and
  optional Wilcoxon/FDR statistics per category.
- `sn_run_milo()` now provides a Shennong wrapper around miloR for
  neighborhood-level differential abundance testing between two sample groups
  from a Seurat embedding.
- `sn_list_palettes()` and `sn_get_palette()` now expose the package palette
  registry directly, and `sn_list_palettes()` now includes preview-oriented
  display output plus the `OkabeIto` palette from `ggokabeito`.

### Changed

- Standardized public metadata-selector arguments on the current API names and
  removed old compatibility aliases such as `group`, `group_col`,
  `sample_col`, `label_col`, `labels_key`, `annotation_col`, `condition_col`,
  `cluster_col`, `reference_key`, `cell_type_key`, `cell_type_col`,
  `cell_state_col`, `groupby`, and `cnv_score_groupby`. Backend-local config
  keys may still use backend names when required by the external Python tools.
- `sn_run_cluster()` now uses a simpler rare-feature interface. The supported
  automatic rare feature methods are `gini` and `local_markers`; less common
  `local_hvg` and `ciara` modes were removed from the clustering wrapper.
  Advanced thresholds are consolidated into `rare_feature_control =
  list(group_max_fraction = ..., group_max_cells = ..., gene_max_fraction =
  ..., min_cells = ...)`; the old scalar threshold arguments have been
  removed.
- `sn_plot_feature()` now silently replaces Seurat's default expression color
  scale when a Shennong palette is requested, avoiding the noisy duplicate
  colour-scale message.
- `sn_plot_feature()` now exposes additional Seurat 5.5 `FeaturePlot()`
  arguments including `assay`, `dims`, `cells`, `alpha`, `stroke_size`,
  `min_cutoff`, and `raster_dpi`. With `ggrastr` available, `raster = TRUE`
  rasterizes a regular ggplot point layer so `pt_size` behaves like
  `raster = FALSE`, fixing overly large rasterized feature points.
- `sn_plot_feature(raster = TRUE)` now falls back to the vector point layer
  when optional package `ggrastr` is unavailable instead of failing during
  Seurat's ordered rasterization.
- `sn_plot_dim()` now exposes `label_halo` so users can disable the white label
  halo/background, and `label = TRUE, repel = TRUE` now keeps a repel-aware
  label layer instead of replacing it with fixed-position shadow text.
- `sn_plot_dot()` now uses black colorbar frame/tick styling and suppresses the
  duplicate colour-scale replacement message when applying Shennong palettes.
- `sn_list_dependencies()` and `sn_install_dependencies()` now classify
  `anndataR` and `tidytemplate` as GitHub-hosted optional dependencies and
  `Nebulosa` as a Bioconductor dependency instead of routing them through CRAN.
  Legacy `.qs` support now remains opportunistic when `qs` is already
  installed, while new one-step dependency installation uses `qs2` and avoids
  attempting the archived, R 4.6-incompatible `qs` package.
- `sn_install_dependencies()` now installs required dependencies for
  GitHub-hosted optional packages by default and stops with the package names
  that remain missing after installer warnings, making partial installation
  failures easier to diagnose.
- Reworked the pkgdown article set around a PBMC3k tutorial path, adding
  explicit data/project and visualization articles and rewriting workflow
  articles to explain why each Shennong function is used before showing the
  code. Heavy or credentialed chunks now remain opt-in through
  `SHENNONG_RUN_VIGNETTES=true` so local website builds stay fast.
- `sn_run_cluster(normalization_method = "sctransform", batch = ...)` now runs
  SCTransform followed by Harmony integration instead of rejecting
  SCTransform-based integration workflows.

- `sn_run_cluster()` now applies `rare_feature_n` per selected
  `rare_feature_method` before de-duplicating the combined rare-aware feature
  set, matching the documented contract. The stored
  `object@misc$rare_feature_selection` record now keeps both the requested
  `rare_feature_n`, the resolved `rare_feature_control`, and the realized
  `selected_rare_feature_n`.
- `sn_calculate_rogue()` now avoids materializing the full matrix before
  optional subsampling, skips redundant entropy work when grouped ROGUE scores
  are requested, and returns tidy per-cluster or per-sample-per-cluster tables
  when grouping metadata are supplied.
- `sn_list_palettes()` now renders palette names and swatches without overlap
  in plot mode and includes built-in viridis-family palettes in the shared
  palette registry. `sn_plot_dim()`, `sn_plot_feature()`, and the other
  `sn_plot_*()` wrappers now reconcile `aspect_ratio` with `panel_widths` /
  `panel_heights` automatically instead of erroring when fixed panel sizes are
  requested, and Seurat reduction plots now apply axis hiding more reliably
  while `sn_plot_dim(label = TRUE)` adds a white halo behind labels for better
  legibility. `sn_plot_feature()` now also supports
  `mode = "density"` for Nebulosa-style embedding density maps with a
  galaxy-like default theme and shared colorbar collection across multi-feature
  plots. `sn_find_doublets()` now records skipped cells as
  `unresolved` instead of `NA` in the stored class column and orders doublet
  classes as `singlet`, `doublet`, then other levels.
- `sn_initialize_seurat_object()` now accepts the named character vectors
  returned by `sn_list_10x_paths()` and imports all detected 10x samples in one
  call, returning a named list of Seurat objects. Annotation-aware
  `sn_filter_genes()` warnings now report example unmatched feature names, and
  `sn_plot_*()` legends use safer non-negative spacing so legend text does not
  overlap plotted panels.
- `sn_initialize_seurat_object()` now recognizes typical 10x Genomics `outs/`
  directories, reads the filtered matrix automatically, and stores discovered
  source metadata such as `raw_feature_bc_matrix` paths and
  `metrics_summary.csv` contents in `Seurat::Misc(object, "input_source")`.
  `sn_remove_ambient_contamination()` now reuses that stored raw path
  automatically when the selected method can use background droplets and the
  caller leaves `raw = NULL`. When `sample_name = NULL` and the input path is a
  named 10x path, the inferred sample identifier is now written into
  `meta.data$sample` automatically.
- `sn_enrich()` now reuses in-session caches for repeated MSigDB term-table
  loads and repeated SYMBOL-to-ENTREZ conversions, which reduces repeated
  overhead during enrichment-heavy test and analysis sessions.
- `sn_interpret_annotation()` now supports cluster-level functional evidence
  through `enrichment_name`, can incorporate cluster QC summaries into the
  prompt, requests structured annotation JSON by default, stores the parsed
  cluster annotation table, can map normalized cell-type labels plus
  confidence/risk fields back onto Seurat metadata, and now defaults to an
  `ellmer`-backed provider path rather than Shennong-managed local provider
  config. It also now resolves `de_name` automatically when omitted,
  preferring a stored `default` marker result, then a single available DE
  result, and otherwise the most recent marker result. Annotation prompts now
  include the full cluster_by evidence table instead of truncating at eight rows,
  request one record per cluster, use more conservative evidence-grounded
  label selection, can inject candidate-label priors for sorted datasets, and
  can attach cluster-neighborhood geometry from reductions such as UMAP.
  Prompt assembly is now more explicitly markdown-structured, which keeps the
  system/task/evidence sections easier to iterate on as prompt templates.
  `sn_interpret_de()`, `sn_interpret_enrichment()`, `sn_write_results()`,
  `sn_write_figure_legend()`, and `sn_write_presentation_summary()` now share
  the same step-wise progress logging and elapsed-time reporting surface.
- LLM-provider integration is now centered on `ellmer`. The old
  `sn_configure_llm_provider()`, `sn_list_llm_providers()`,
  `sn_get_llm_provider()`, `sn_make_openai_provider()`, and
  `sn_make_sub2api_provider()` compatibility shims have been removed, along
  with the legacy `~/.shennong` provider/history workflow. The supported entry
  point is now `sn_make_ellmer_provider()`, which can explicitly forward
  `reasoning_effort` to compatible GPT-5 chat-completions endpoints. Default
  environment-variable discovery is now limited to `OPENAI_*`; the temporary
  `SUB2API_*` compatibility path has been removed.
- Annotation metadata write-back is now lean by default. High-level
  interpretation writes only the core fields needed for visualization and
  grouping (`label`, `broad_label`, `confidence`, `status`, `risk_flags`)
  into Seurat metadata, while detailed supporting markers/functions/notes stay
  in the stored interpretation result table under `object@misc`.
- `sn_calculate_composition()` now supports multi-column `group_by` values,
  can return proportions, counts, or both through `measure`, preserves factor
  columns from the source metadata, filters returned composition categories by
  `min_cells`, and can sort a single grouping column by a chosen category level
  such as WT proportion.
- `sn_initialize_project()` is now the single project bootstrap entry point. It
  now also writes a repository `.gitignore` plus a generated project `.Rproj`
  file into initialized analysis repositories.
- The built-in `Paired` discrete palette is now overridden by Shennong so its
  brightest yellow swatch uses `#ECD577` instead of Brewer's `#FFFF99`.
- Visualization helpers now share a common discrete-palette resolver. Named
  palettes such as `\"Paired\"` expand automatically when more categories are
  present than the base palette length, and key plotting helpers now expose
  `panel_widths` / `panel_heights` plus consistent axis-label handling.
- `sn_plot_dim()` and `sn_plot_feature()` now default to hiding coordinate
  axes, choose point sizes automatically when `pt_size = NULL`, and keep a
  shared point-size heuristic for small versus large datasets. `sn_plot_dot()`
  now defaults to a warmer-high / cooler-low color direction, keeps the
  Z-score legend ahead of the percent legend, supports `Min`/`Max` legend
  labels, uses hollow black-edged percent legend dots, accepts
  `legend_position`, and uses thicker black colorbar ticks.
- Continuous-color helpers now use the same palette registry through
  `sn_get_palette(..., palette_type = "continuous")`, and expression-oriented
  plotting functions now apply continuous palettes through the shared internal
  resolver instead of separate ad hoc `scale_*_distiller()` logic.
- `sn_plot_barplot()` now supports automatic summary bars for repeated
  observations, optional SD/SE error bars, and optional jittered raw points,
  making it suitable for sample-level effect summaries as well as simple
  identity bars.
- `sn_compare_composition()` now adds a `change` factor with levels
  `Increase` and `Decrease` derived from the sign of `log2_fc`.
- `sn_find_doublets()` now skips zero-count and low-feature cells before
  running `scDblFinder()` on corrected layers, records corrected-layer results
  with `_corrected` suffixes, and works with the zero-count flags produced by
  ambient-RNA correction.
- `sn_filter_cells()` now validates its `method` argument, keeps constant-value
  QC groups when MAD collapses to zero, and checks plotting dependencies
  explicitly before rendering diagnostics. `sn_filter_genes()` now validates
  `min_cells` and keeps its threshold summary stable when the requested
  threshold exceeds the number of cells in the object.
- `sn_run_cluster()` now defaults `hvg_group_by` to the `batch` column when
  batch integration is requested and the user does not explicitly supply a
  separate HVG grouping variable.
- `sn_standardize_gene_symbols()` now also accepts character vectors of gene
  symbols or gene IDs and returns the standardized vector directly, while
  preserving the existing matrix and Seurat-object behavior.
- `sn_run_cluster(normalization_method = "scran", batch = ...)` now runs
  scran normalization before the selected batch-integration backend instead of
  rejecting batch workflows.
- The repository now ships `scripts/check-prepush.R` so maintainers can run
  documentation, targeted tests, the full test suite, `R CMD build`, and
  `R CMD check --no-manual` in one local pre-push command.
- Internal DE result storage now reuses the shared misc-result helper instead
  of maintaining a second collection-specific implementation.
- `sn_check_version()` and `sn_install_shennong()` now use the unified
  `source` / `ref` arguments for GitHub and local source paths; the old
  `github_repo`, `github_ref`, and `local_path` aliases have been removed.

### Fixed

- Fixed `sn_list_10x_paths()` so detected samples are returned in deterministic
  sample-name order instead of depending on platform-specific filesystem or
  `find` traversal order.
- Fixed `sn_standardize_gene_symbols()` so unresolved or ambiguous
  `HGNChelper` suggestions no longer propagate `NA` row names or drop
  otherwise valid original symbols. Truly missing or empty feature names are
  still removed before duplicate symbols are aggregated.
- Fixed IO edge cases where `sn_read(row_names = "column")` failed for column
  names, detected 10x spatial directories were not dispatched to a custom
  reader, `sn_write()` failed for existing `SingleCellExperiment` h5ad exports,
  and `qs2` serialization called `qs2::qs_save()` with the wrong argument name.
- Fixed `sn_run_celltypist()` path inputs so precomputed CellTypist inputs return
  a prediction table instead of trying to write metadata onto a character path.
- Fixed `sn_remove_ambient_contamination(method = "soupx")` so Seurat returns
  now add `nCount_<assay>_corrected` and `nFeature_<assay>_corrected` metadata
  from the SoupX-corrected layer, matching the decontX writeback behavior.
- Fixed the pkgdown reference index by adding the exported
  `sn_sweep_cluster_resolution()` topic.
- Fixed a malformed hidden R chunk in the clustering vignette that prevented
  pkgdown from rendering articles.
- Fixed local and CI test helpers so Seurat fixtures used in `de_enrich` and
  `utils` tests no longer depend on implicit species inference or non-returned
  normalization calls.
- Fixed pkgdown reference indexing so new helpers such as `sn_assess_qc()` and
  `sn_list_10x_paths()` are included in the generated site configuration.

# Version 0.1.2

Released 2026-03-25.

### Added

- Bundled human and mouse GENCODE gene-annotation data was added for gene-level
  filtering workflows. The snapshot now also stores genome-context fields such
  as sequence name, source, coordinates, strand, gene status/source, and
  annotation level in addition to gene identifiers and gene types.
- Small built-in PBMC example assets were added as `pbmc_small` and
  `pbmc_small_raw`, sampled from the packaged `pbmc1k` / `pbmc3k` references
  for check-safe examples and README workflows.
- Integration and cluster-diagnostics metrics were expanded with
  `sn_calculate_silhouette()`, `sn_calculate_graph_connectivity()`,
  `sn_calculate_pcr_batch()`, `sn_calculate_clustering_agreement()`,
  `sn_calculate_isolated_label_score()`, `sn_calculate_cluster_entropy()`,
  `sn_calculate_cluster_purity()`, `sn_identify_challenging_groups()`, and the
  aggregate `sn_assess_integration()` wrapper.
- Rare-cell-aware clustering support was added through
  `sn_detect_rare_cells()` and new `sn_run_cluster()` parameters that can
  append rare-aware features such as Gini-selected genes, local HVGs, local
  markers, or optional CIARA-derived features before PCA and Harmony.
- Bulk RNA-seq deconvolution support was added through
  `sn_deconvolve_bulk()`, `sn_store_deconvolution()`, and
  `sn_get_deconvolution_result()`, covering local BayesPrism runs plus
  local CIBERSORTx container workflows and result import.
- pkgdown documentation is now reorganized around workflow stages rather than
  source files alone. New end-to-end articles cover preprocessing and QC,
  clustering and integration, metrics and diagnostics, annotation and pathways,
  composition analysis, and interpretation/reporting.
- `sn_find_de()` now supports pseudobulk differential expression with
  `limma` in addition to `DESeq2` and `edgeR`.
- `sn_find_de()` now supports marker discovery with `COSGR` when the optional
  GitHub package is installed.
- `sn_initialize_project()` now scaffolds user analysis repositories from the
  shipped `inst/codex/project-template/` assets, creating a governed
  project layout with `AGENTS.md`, `memory/`, `docs/standards/`, `skills/`,
  `config/`, `data/`, `scripts/`, `notebooks/`, `runs/`, and `results/`.
- `sn_initialize_project()` is now a convenience wrapper over the packaged
  project-template initializer.
- `sn_get_codex_skill_path()` now exposes packaged Codex asset paths for the
  Codex root, package-usage skills, project template, and project-template
  skills.
- `sn_install_codex_skill()` now installs package-usage skills, project
  governance skills, or both from the packaged Codex asset layout.
- Signature catalog helpers were added:
  `sn_list_signatures()`, `sn_add_signature()`, `sn_update_signature()`, and
  `sn_delete_signature()`.
- Stored-result discovery and retrieval helpers were added:
  `sn_list_results()`, `sn_get_de_result()`, `sn_get_enrichment_result()`, and
  `sn_get_interpretation_result()`.
- The interpretation layer now supports user-supplied background context and
  dual output styles for either model-facing prompt bundles or human-readable
  summaries.

### Changed

- Internal sparse-matrix workflows were optimized for better runtime and lower
  peak memory use. Pseudobulk DE aggregation now groups columns without
  materializing dense matrices, split Seurat assay layers are combined through
  sparse triplet assembly instead of repeated indexed writes, exact kNN
  fallbacks now use blockwise distance evaluation instead of constructing a
  full cell-by-cell distance matrix, and gene-symbol standardization now reuses
  bundled annotation data plus grouped row aggregation instead of external CSV
  reads and per-column duplicate collapsing.
- Developer-facing informational notifications are now routed through internal
  package logging helpers instead of a mix of ad hoc `logger`, `message()`,
  and `cli` progress calls. Progress wording is now more consistent across
  initialization, clustering, enrichment, and preprocessing workflows, while
  real `warning()` and `stop()` conditions retain their existing semantics.
- `sn_enrich()` now uses a single `x` input with automatic dispatch across
  gene vectors, ranked named vectors, data frames, and Seurat-stored DE
  results. Its `gene_clusters` formula now drives both grouped ORA
  (`gene ~ cluster`) and ranked GSEA (`gene ~ log2fc`) workflows, and
  multi-database requests such as `database = c("H", "GOBP", "KEGG",
  "C2:CP:REACTOME")` are supported in one call.
- `sn_enrich()` now aligns MSigDB arguments with `msigdbr` by supporting
  `collection` / `subcollection`, while still accepting collection strings in
  `database` such as `H` or `C2:CP:REACTOME`.
- `sn_enrich()` now filters returned enrichment tables by raw p-value rather
  than adjusted p-value and writes one result file per requested database with
  stable `prefix` / `outdir` naming.
- Dataset documentation is now being consolidated into a single `R/data.R`
  source, and user-facing examples are shifting from simulated matrices or
  network-backed PBMC downloads to the built-in small PBMC package data.
- `sn_find_de()` now uses a single `method` argument instead of separate
  `test_use` and `pseudobulk_method` parameters.
- `sn_find_de()` now uses `layer` consistently for Seurat v5 workflows and no
  longer exposes the legacy `slot` argument.
- `sn_filter_genes()` now supports annotation-aware filtering through
  `gene_class` and exact `gene_type` values backed by the bundled GENCODE
  snapshot for human and mouse.
- Signature registry maintenance now delegates add/update/delete operations to
  the upstream `SignatuR` package API instead of maintaining a parallel custom
  tree-editing implementation inside Shennong.
- Harmony-backed integration now targets the `immunogenomics/harmony`
  `harmony2` developer branch in package metadata and CI instead of the CRAN
  release line.
- Integration metrics now prefer stored Seurat neighbor graphs when available
  and otherwise fall back to Annoy-based approximate kNN or exact distance
  search, keeping the default assessment path fast enough for routine use.
- `sn_get_signatures()` now reads from a package-owned signature snapshot built
  from the full `SignatuR` tree during development, so runtime signature
  retrieval is stable, tree-structured, and no longer depends on the installed
  `SignatuR` package version.
- Signature build assets now center on `data/shennong_signature_catalog.rda`,
  which is rebuilt directly from the upstream `SignatuR` dataset during
  development, replacing the opaque `R/sysdata.rda` storage used by the earlier
  snapshot prototype.
- `sn_enrich()` now stores enrichment results in
  `object@misc$enrichment_results[[store_name]]` when a Seurat object is
  supplied, aligning enrichment with the existing stored DE workflow.
- pkgdown articles, shipped Codex skill references, and `NEWS.md` are now
  treated as required deliverables for any user-facing workflow change.

### Fixed

- Fixed `sn_install_shennong(channel = "github")` so it no longer tries to
  resolve runtime-optional GitHub `Suggests` by default during installation.
- Removed the optional rare-cell backends `FiRE`, `CellSIUS`, and `EDGE` from
  Shennong's supported dependency surface and `sn_detect_rare_cells()`
  interface.
- Fixed `sn_plot_dot()` theme handling so the optional `catplot` theme no
  longer tries to impose an additional aspect ratio on top of `coord_fixed()`.
- Fixed Rd example line-width failures in `R CMD check` and declared the
  runtime `data.tree` dependency explicitly.

# Version 0.1.1

Released 2026-03-19.

### Added

- Automatic human/mouse species inference through `sn_get_species()` using
  `hom_genes` and mitochondrial naming patterns.
- A layer-aware pkgdown article covering inferred species, non-default count
  layers, and stored differential-expression metadata.
- Additional regression tests for species inference, DE metadata, and
  BPCells-backed Seurat layers.

### Changed

- `sn_initialize_seurat_object()` now attempts species inference before deciding
  whether to compute species-specific QC metrics.
- `sn_get_signatures()` now uses package-local fallback signatures for the core
  categories needed by Shennong workflows when `SignatuR` is unavailable.
- Stored DE results now include schema and provenance metadata such as package
  version, timestamp, assay/layer context, and threshold settings.
- `sn_plot_*()` helpers now treat `catplot` as an optional enhancement rather
  than a mandatory dependency.

### Fixed

- Fixed `sn_remove_ambient_contamination()` for BPCells-backed Seurat layers by
  materializing BPCells matrices before passing them to `decontX`.
- Fixed optional dependency handling so missing `SignatuR` no longer breaks core
  initialization and clustering paths.
- Fixed SoupX validation order so missing `raw` input is reported before package
  availability issues.
- Fixed the `sn_deconvolve_bulk()` example so the `cibersortx` dry-run path no
  longer fails package examples for missing credentials.
- Fixed GitHub Actions dependency installation so check and coverage jobs keep
  a minimal, solvable dependency set instead of trying to install unsupported
  GitHub-only optional backends during lockfile generation.

# Version 0.1.0

Released 2026-03-18.

### Added

- `sn_load_data()` as the primary example-data loader.
- Consolidated clustering into `sn_run_cluster()`.
- A unified ambient contamination interface for SoupX and decontX.
- Expanded test coverage across composition, data loading, clustering,
  utilities, visualization, and ambient contamination.
- Missing help pages for exported functions, CI scaffolding, a richer README,
  and updated package metadata for current R syntax requirements.

### Changed

- Removed the old `sn_load_pbmc()` wrapper in favor of the primary
  `sn_load_data()` entry point.
- Refreshed generated documentation and package metadata as part of the
  modernization effort.

### Fixed

- Addressed multiple modernization issues in the package build, test, and
  documentation pipeline.
