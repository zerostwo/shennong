# Shennong Modernization Decisions

Last updated: 2026-05-03

## 2026-05-03

- Python single-cell methods should enter Shennong through the existing
  `sn_run_cluster()` integration-method contract rather than through a parallel
  public clustering API. scVI and scANVI therefore behave like other batch
  integration backends from the user's perspective and still return a regular
  Seurat object.
- Shennong should not make the package itself depend on Python or reticulate
  for scVI/scANVI. The R side writes simple MatrixMarket/CSV interchange files,
  pixi owns the Python environment under `~/.shennong/pixi/<family>/`, and the
  Python backend writes explicit artifacts that can be inspected or reused.
- scVI/scANVI outputs should be imported as reductions, not as a replacement
  assay. The learned latent representation is the integration product used for
  neighbors, clusters, and UMAP, while counts and normalized assay data remain
  under Seurat's existing assay/layer structure.
- Python environments should be generated under the user-level Shennong home,
  not inside the analysis project and not inside the installed R package. The
  package vendors backend scripts beside reusable pixi config templates under
  `inst/pixi/<family>/scripts/`; materialized pixi manifests, lockfiles,
  environments, runtime inputs, and outputs belong under `~/.shennong/`.
- The pixi binary is a user-level tool, but pixi runtime state should be
  isolated when Shennong launches Python backends. `sn_run_cluster()` therefore
  ensures pixi is available, then runs with a Shennong-level `PIXI_HOME` so
  mirror configuration and caches are discoverable and do not overwrite global
  user settings.
- GPU support should be opt-out and CPU fallback should remain reliable. The
  scVI/scANVI manifest defines separate CPU and CUDA environments; automatic
  selection uses CUDA only when `nvidia-smi` reports an NVIDIA GPU, while other
  systems run the CPU environment unless the caller explicitly overrides it.
- Pixi configs should represent installable environment families, not every
  user-facing method alias. `scanvi` therefore resolves to the shared `scvi`
  environment, and `scpoli` resolves to the shared `scarches` environment.
  Spatial analysis is represented by concrete tools such as `cell2location`,
  `tangram`, `squidpy`, `spatialdata`, and `stlearn`, rather than a generic
  `spatial` environment.
- Public `sn_run_*` Python method wrappers should prefer object-level
  contracts whenever Shennong can define a stable Seurat input/output schema.
  Command-mode wrappers remain available as a fallback, but user-facing
  analysis APIs such as `sn_run_infercnvpy(object = ...)` should handle export,
  pixi execution, and import back into Seurat.
- Runner scripts belong with their pixi family config under
  `inst/pixi/<family>/scripts/`, not in a separate global `inst/python/`
  directory. This keeps each installable Python family self-contained.

## 2026-03-30

- Shennong should use `ellmer`'s native structured-output surface when the
  task contract is inherently tabular. Cluster annotation now prefers
  `chat_structured()` with an explicit schema and only falls back to text JSON
  parsing when the provider does not expose structured data directly.
- Shennong should not keep temporary proxy-specific environment-variable
  branches once the backend is standardized. The LLM setup path now recognizes
  only `OPENAI_*`, leaving provider-specific gateway handling to the explicit
  `base_url` passed into `sn_make_ellmer_provider()`.
- Prompt templates should stay in files and be assembled into explicit
  markdown sections rather than one long concatenated blob. This keeps the
  interpretation layer closer to `ellmer`'s prompt-design guidance and makes
  future prompt iteration easier to review.
- Tool calling should be applied narrowly as a workflow primitive, not as a
  replacement for annotation logic. The current use is a focused
  cluster-evidence lookup step in agentic annotation before the final
  structured refinement pass.

## 2026-03-29

- When Shennong initializes a Seurat object from a recognized 10x Genomics
  `outs/` directory, it should retain that provenance on the object instead of
  treating the input as an anonymous count matrix. Storing the discovered
  filtered/raw paths and metrics summary under `Seurat::Misc(object,
  "input_source")` keeps ambient-RNA correction and downstream reporting
  ergonomic without introducing new user-facing mandatory arguments.
- The public path-discovery helper should match the dominant workflow instead
  of exposing a broad but vague directory scanner. Shennong now uses a
  10x-specific `sn_list_10x_paths()` entry point that defaults to returning
  `outs/` directories, with opt-in selection of filtered/raw/H5/metrics assets
  for downstream automation.
- Repeated enrichment support lookups should be cached within an R session.
  MSigDB term tables and SYMBOL-to-ENTREZ mappings are deterministic for a
  given input set and package state, so caching them reduces avoidable latency
  in test suites and interactive enrichment workflows without changing results.
- Annotation-oriented LLM workflows should be structured, not free-form. When
  Shennong asks a provider to annotate clusters, it should request a stable
  JSON schema, normalize returned cell-type names into one user-selected style,
  and, when possible, map the parsed cluster-level labels back onto each cell's
  metadata for downstream filtering and plotting.
- Interpretation helpers should follow the same stored-result ergonomics as
  enrichment helpers. When annotation workflows need a stored marker table and
  the caller omits `de_name`, Shennong should prefer the stored `default`
  result, then a single available DE result, and otherwise the most recent
  marker result rather than failing on a missing required argument.
- Cluster annotation prompts should optimize first for evidence completeness,
  then for brevity. Shennong now allows annotation prompts to include the full
  cluster evidence table instead of a generic eight-row preview so LLM-based
  annotation cannot silently drop later clusters in larger Seurat objects.
- Shennong should not maintain its own persistent LLM provider registry or
  request-history layer when `ellmer` already owns transport, credentials, and
  provider semantics. High-level interpretation now uses `ellmer` as the only
  supported backend and no longer depends on `~/.shennong` provider files.
- Annotation improvements should stay at the level of comparative evidence and
  canonical lineage guardrails rather than a hand-maintained subtype score
  table. The current guardrails now explicitly cover `KIT+ ILCP-like` versus
  mature `ILC3` in blood, mixed `T/NK` states, `NK/ILC3` transitional states,
  and `erythroid contamination`, but they do so through marker-program priors
  and context-sensitive prompt guidance rather than a closed manual label map.
- Sorted or enriched single-cell datasets need explicit annotation priors. For
  cases such as blood ILC-sorted data, Shennong now exposes
  `label_candidates` so prompts can constrain the plausible label space
  without forcing a hard whitelist.
- Annotation evidence should favor specificity over raw rank when possible.
  Shennong now allows marker and enrichment summaries to prefer
  cluster-restricted signals and to attach reduction-neighborhood summaries as
  supporting evidence, while keeping geometry subordinate to marker evidence.
- Annotation logic should encode an evidence-comparison workflow, not a
  hand-maintained subtype scorecard. Shennong now keeps only broad lineage
  guardrails and moves accuracy improvements into a staged `agentic`
  annotation workflow that compares clusters against each other, reuses
  canonical-marker snapshots, and asks the model to refine ambiguous clusters
  in a second pass.
- Seurat metadata should stay compact and plot-oriented. Cluster-level
  annotation evidence such as supporting markers/functions, notes, and
  recommended checks now remains in the stored interpretation table by
  default, while metadata only receives the core fields needed for grouping
  and visualization.
- 10x path-discovery helpers should return workflow-ready sample labels instead
  of anonymous path vectors. Shennong now infers sample names from the
  directory structure and carries those names into Seurat initialization when
  `sample_name = NULL`, reducing repetitive manual sample annotation in common
  multi-sample 10x workflows.
- Embedding-style plots should optimize first for readability on small
  datasets. `sn_plot_dim()` and `sn_plot_feature()` now hide coordinate axes by
  default and use shared automatic point sizes unless the caller overrides
  `pt_size`, while `sn_plot_dot()` uses a more opinionated legend layout and
  colorbar styling tuned for publication-oriented defaults.

## 2026-03-25

- Default collaboration for this repository should stop at local commits unless
  the user explicitly asks to push to GitHub. Do not push after a successful
  fix or refactor by default; wait for an explicit "push"/"submit to GitHub"
  instruction.
- `NEWS.md` updates belong under an unreleased section while work is still in
  progress. Once a version is marked as released, treat that section as frozen
  and record subsequent changes in the current unreleased section instead of
  appending to the released version.
- Pre-push validation should be encoded in one repository-local command rather
  than left to memory. Shennong now ships `scripts/check-prepush.R` so the
  standard maintainer path can run documentation, tests, build, and check in a
  predictable order before code reaches GitHub Actions.
- QC reporting should be recomputable from merged Seurat metadata rather than
  depend on fragile per-object history alone. `sn_assess_qc()` therefore uses
  sample-level metadata aggregation by default, can optionally compare against a
  pre-filter reference object for before/after judgments, and stores the
  resulting report as a convenience snapshot rather than as the sole source of
  truth.
- QC filtering should fail safe on degenerate thresholds instead of classifying
  an entire group as outliers. When MAD collapses to zero, Shennong should
  treat values exactly on the derived bound as passing, and it should validate
  the supported outlier method plus plotting dependencies at the public entry
  point rather than relying on downstream failures.
- When Shennong has to fall back to exact neighbor search, it should optimize
  for bounded memory rather than for the shortest implementation. The exact
  kNN path now computes distances blockwise so larger embeddings do not require
  an `n x n` dense distance matrix just to recover the top-k neighbors.
- Sparse grouped aggregation should stay sparse all the way through DE and
  preprocessing helpers whenever the underlying algorithm allows it. Internal
  grouped count aggregation and duplicate-feature collapsing now use matrix
  multiplication / grouped row sums instead of dense transposes or per-column
  `apply()` loops.
- Split Seurat assay layers should be merged by sparse assembly, not by
  repeated subassignment into a preallocated sparse matrix. This avoids
  unnecessary reallocations when the package needs a temporary combined layer
  for downstream workflows.
- Repeated nearest-neighbor plumbing should live in one shared helper layer.
  Annoy neighbor extraction and self-neighbor cleanup are now centralized so
  clustering and metrics use the same implementation and future tuning applies
  consistently across both workflows.
- Grouped HVG selection should reduce peak memory without changing ranking
  semantics. Shennong now iterates through one metadata group at a time rather
  than materializing all split Seurat objects up front; the ranking logic
  remains inline because it is specific to this one feature-selection path.
- Developer-facing progress notifications should use one lightweight internal
  logging surface. Shennong now routes informational messages through internal
  glue-aware wrappers so progress wording and formatting stay consistent,
  without changing `warning()` or `stop()` semantics for real runtime
  conditions.
- Integration assessment should be a first-class workflow, not an ad hoc
  collection of standalone helpers. Shennong now exposes `sn_assess_integration()`
  as the aggregate entry point while still exporting lower-level metric
  functions for users who want one score at a time.
- The default performance strategy for graph-based metrics is now:
  reuse an existing Seurat neighbor graph when available, otherwise fall back
  to Annoy-based approximate kNN through Seurat, and only then use an exact
  distance-based neighbor search. This gives a practical speed/accuracy balance
  without introducing a new heavy dependency just for nearest-neighbor search.
- Rare groups and poorly separated groups should be reported explicitly instead
  of being inferred indirectly from a single global score. The metrics layer
  therefore now includes `sn_identify_challenging_groups()` with group size,
  neighbor purity, connectivity, and silhouette-derived diagnostics.
- Integration assessment should include explicit small-population and
  within-cluster mixing diagnostics rather than relying only on global LISI or
  PCR summaries. The default workflow now includes isolated-label preservation,
  cluster purity against biological labels, and cluster entropy against batch
  labels so rare or weakly separated groups remain visible in routine reports.
- Rare-cell support should enter the workflow at two points: as a diagnostic
  layer (`sn_detect_rare_cells()`) and as a feature-selection enhancement for
  clustering/integration (`sn_run_cluster(rare_feature_method = ...)`). This is
  more controllable than replacing the main HVG logic with a single opaque
  rare-cell algorithm.
- Bulk RNA-seq deconvolution should be exposed as a first-class Shennong
  workflow rather than as scattered export scripts. The package now uses a
  single `sn_deconvolve_bulk()` entry point with backend-specific behavior:
  BayesPrism runs locally through its R API, while CIBERSORTx is handled as a
  local container workflow with explicit credentials and reproducible commands
  instead of a manual website-upload step.
- User-facing documentation should be organized by analysis stage rather than
  by source-file boundaries. The primary taxonomy is now:
  preprocessing/QC; clustering/integration; diagnostics/benchmarking;
  annotation/markers/pathways; composition/comparative analysis; and
  interpretation/reporting. pkgdown reference groups and articles should follow
  that workflow-oriented structure.
- Because Shennong is still explicitly experimental, current refactors do not
  need to preserve backward compatibility with earlier 0.x interfaces. Legacy
  wrappers, deprecated aliases, and "backwards-compatible" transitional
  arguments should be removed rather than carried forward.
- `sn_enrich()` should prefer one primary input contract over parallel
  `x`/`object`/`gene_clusters` entry points. `x` is now the only primary input;
  stored-result enrichment on Seurat objects should flow through `x` plus
  `source_de_name`, while direct enrichment of vectors/data frames returns raw
  results and can be stored later with `sn_store_enrichment()`.
- Enrichment intent should be inferred from the formula RHS whenever practical:
  categorical/grouping RHS values mean grouped ORA, while numeric RHS values
  mean ranked GSEA. This keeps `sn_enrich()` aligned with the user's table
  shape instead of forcing separate argument vocabularies.
- Multi-database enrichment is part of the public workflow, not a caller-side
  loop. `sn_enrich()` should therefore accept a vectorized `database` argument
  spanning GO, KEGG, and MSigDB collection strings in one call, with
  per-database storage/output names generated deterministically.
- Upstream `clusterProfiler` p-value cutoffs should not define the final user
  contract for Shennong enrichment. Shennong now runs enrichment with permissive
  upstream cutoffs and applies its own raw p-value filtering afterward so the
  exported `pvalue_cutoff` argument has the intended semantics.
- Dataset documentation should be centralized in one `R/data.R` source rather
  than split across multiple small `data_*.R` files. Small built-in PBMC
  example assets should be packaged from sampled `pbmc1k`/`pbmc3k` data so
  README and examples can avoid network downloads.

## 2026-03-24

- Gene-annotation filtering should be driven by bundled package data built from
  local GENCODE GTF sources rather than ad hoc runtime imports. `sn_filter_genes()`
  now accepts a coarse `gene_class` selector and exact `gene_type` values, with
  species-aware matching against either gene symbols or GENCODE gene IDs.
  The packaged annotation snapshot should retain enough genome context to stay
  useful outside simple biotype filtering, so it now keeps seqname/source,
  coordinates, strand, gene status/source, level, and version-stripped IDs.
- Shennong should not maintain a second independent implementation of
  signature-registry editing when `SignatuR` already provides the authoritative
  mutation API. Signature add/update/delete helpers now bridge into a transient
  `SignatuR` object, call the upstream functions, and then serialize back to
  Shennong's packaged snapshot.
- Harmony integration should track the upstream `immunogenomics/harmony`
  `harmony2` developer branch explicitly in package metadata and CI, rather
  than implicitly relying on the CRAN release line. The runtime namespace
  remains `harmony`, so the clustering implementation can stay on
  `harmony::RunHarmony()` while installation sources are pinned to the intended
  development branch.

## 2026-03-21

- The package root must remain a standard R package repository. Analysis-project governance directories such as `memory/`, `runs/`, `results/`, and raw/processed project data directories belong only in initialized user projects, not in the package root.
- Initialized-project governance is now shipped under `inst/codex/project-template/`. Package-usage skills are now shipped under `inst/codex/package-skills/`. Maintainership docs remain under `docs/codex/`.
- The shipped project template should expose one canonical initializer:
  `sn_initialize_project()`.
- Empty directories inside the shipped project template should be created from an explicit manifest (`inst/codex/project-template/directories.txt`), not hidden `.gitkeep` placeholders. This keeps the scaffold maintainable while avoiding hidden-file notes in `R CMD check`.
- Coverage expansion should prioritize durable, stateful behavior over superficial line hits. For this pass, the preferred targets were Seurat layer/state helpers, stored-result contracts in `object@misc`, packaged project scaffolding, IO dispatch paths, and external-tool integration seams such as CellTypist.
- Signature retrieval should not depend on an optional upstream package at runtime. Shennong now ships a meaningful `shennong_signature_catalog` dataset under `data/`, rebuilt directly from the upstream `SignatuR` dataset during development, and user-facing docs should describe bundled signature paths rather than `SignatuR` installation requirements.

## 2026-03-18

- The package will be modernized in small, validated milestones rather than a repo-wide rewrite.
- The current exported API in `NAMESPACE`, plus documented datasets, is the default compatibility boundary until tests or documentation justify a change.
- `sn_quick_cluster()` is currently treated as a documented-but-not-exported inconsistency that requires explicit handling in a later milestone rather than an immediate behavior change.
- `docs/codex/` is the canonical project memory location for this effort; these files remain outside package builds through the existing `^docs$` rule in `.Rbuildignore`.
- Pre-existing working tree changes in `.Rprofile`, `NAMESPACE`, and removed legacy files are considered user state and will not be reverted.
- The first technical modernization slice after scaffolding will prioritize additional tests and metadata alignment because they reduce risk for every later refactor.
- For `sn_calculate_composition()`, the documented behavior was treated as the compatibility source of truth: inconsistent `additional_cols` values now trigger a warning while still using the first value per group.
- The repository-level `docs` ignore rule was narrowed so `docs/codex/` can be versioned while generated pkgdown output under `docs/` remains ignored.
- For the current task, `sn_load_data()` will become the primary public example-data loader. `sn_load_pbmc()` should be retained only as a compatibility bridge if that can be done without re-introducing duplicate implementation.
- `sn_quick_cluster()` will no longer remain as a separate clustering implementation; its non-integration workflow should be absorbed by `sn_run_cluster()`.
- `sn_load_pbmc()` was retained as a deprecated wrapper around `sn_load_data()` to avoid an unnecessary hard break while still moving the public examples toward the generalized loader.
- `sn_run_cluster()` now owns both single-dataset clustering and Harmony integration. The integration path was simplified to a shared HVG/PCA workflow before Harmony instead of the previous split/join-plus-`SelectIntegrationFeatures5()` path because the simpler flow validated reliably on real PBMC data.
- `sn_remove_ambient_contamination()` now uses a common API for SoupX and decontX: shared input coercion, shared optional cluster injection, SoupX-specific `raw` requirement, and Seurat-object round-tripping only when the caller actually supplied a Seurat object.
- `DESCRIPTION` was updated to reflect the packages touched in this task, and the minimum R version was raised to `R (>= 4.1.0)` because the package already uses base pipe syntax.
- The package README is now maintained through `README.Rmd` and rendered to `README.md`; examples are shown but not executed during README rendering to keep the build stable.
- The clustering vignette is now guarded by the `SHENNONG_RUN_VIGNETTES` environment variable so package checks remain deterministic while preserving the documented workflow.
- `sn_read()` and `sn_write()` were moved off `rio:::` internals. Custom format handling now uses explicit dispatch helpers around public `rio::import()` and `rio::export()` calls.
- `sn_add_data_from_anndata()` now uses a CSV-specific fallback reader for AnnData-exported metadata and embedding tables because `rio::import()` does not preserve row names reliably for one-column CSV metadata.
- `sn_run_celltypist()` now resolves its executable path from an explicit argument first, then `getOption("shennong.celltypist_path")`, then `Sys.which("celltypist")`, removing the hidden dependency on a deleted global path object.
- Seurat workflows that conceptually analyze a count matrix should expose `assay` and `layer` when practical instead of silently assuming `RNA/counts`. This was applied to clustering, doublet detection, gene filtering, scran normalization, CellTypist export, and ROGUE scoring.
- For merged Seurat v5 objects that store split layers such as `counts.sample1` and `counts.sample2`, Shennong now combines matching layers internally before downstream analysis rather than rejecting the object or only using the first split layer.
- `sn_run_cluster()` preserves the original `counts` layer even when the analysis is driven by a different layer. The selected layer is copied into a temporary `counts` layer only for the duration of the workflow and restored before returning the Seurat object.
- GitHub-only optional backends should not be added to the baseline CI dependency solve unless the corresponding examples, tests, or pkgdown builds actually execute them. Some upstream repos are missing, unsolved on current Bioconductor stacks, or depend on unavailable packages, so the safer policy is to keep them in `Suggests` and install them only in targeted environments.
- Explicit `methods::slot` and `slot<-` imports are required for the Seurat command-logging pattern used throughout the package; relying on implicit availability leaves a persistent `R CMD check` note.
- decontX correction should not invent synthetic counts to rescue cells that round to zero. The safer compatibility policy is to preserve the original counts for those cells by default, emit a warning with the affected cell count, and offer an explicit `remove_zero_count_cells` switch for callers who prefer to discard them.
- When ambient-contamination correction returns a Seurat object, the corrected per-cell totals should be surfaced explicitly in metadata as assay-scoped `nCount_*_corrected` and `nFeature_*_corrected` columns so downstream QC decisions remain inspectable.
- The old `pipeline` argument on `sn_run_cluster()` was too implementation-specific and ambiguous. Public APIs should instead describe the normalization strategy directly with `normalization_method`, using the same vocabulary across `sn_run_cluster()` and `sn_normalize_data()`.
- `sn_initialize_seurat_object()` should not infer study/sample metadata from `project`; explicit `sample_name` and `study` parameters are clearer and keep metadata creation opt-in.
- Cell-cycle scoring is treated as an opportunistic clustering enhancement, not a hard requirement. If the selected assay lacks sufficient marker overlap for the requested species, the workflow should warn and continue instead of failing inside Seurat.
- Shennong should expose package-maintenance helpers for users directly in the API. `sn_check_version()` reports whether the installed package is current relative to CRAN or GitHub, and `sn_install_shennong()` provides an explicit installation entry point for either channel.
- pkgdown output must not reuse `docs/`, because this repository uses `docs/codex/` for persistent modernization memory. The generated site now writes to `site/`, which is ignored in git and excluded from package builds.
- GitHub Actions should avoid generic installation of all Suggests when that pulls in GitHub-only or non-standard packages such as `BPCells`. CI now installs hard dependencies plus an explicit, validated subset of optional packages instead of relying on `pak` to solve every Suggests entry.
- Automatically installing missing GitHub packages at runtime is an unsafe side effect for an analysis package. `check_installed_github()` now fails with a clear `remotes::install_github()` instruction instead of modifying the user library implicitly.
- Seurat command logging should be centralized in one internal helper, but the recorded command names must remain the public `sn_*` function names rather than the helper name so reproducibility metadata stays meaningful.
- Batch-aware integration should not hard-code a single HVG strategy. `sn_run_cluster()` now uses `hvg_group_by` instead of an abstract `hvg_method` enum: `NULL` means global HVGs, while any metadata column name triggers per-group HVG discovery followed by ranked merging.
- pkgdown articles should demonstrate real outputs, not code-only skeletons. The PBMC workflow article is now evaluated during pkgdown builds and uses `pbmc1k`/`pbmc3k` to show actual clustering, integration, composition, LISI, marker, and enrichment results while remaining check-safe outside pkgdown.
- The repository should not vendor the `SignatuR` R6 `Node` object as package data. It breaks lazydata installation and duplicates upstream state. Signature lookups now load the dataset from the installed `SignatuR` package at runtime instead.
- End-user Codex materials should ship inside the installed package under `inst/codex/project-template/` and `inst/codex/package-skills/`, while repository-only planning and modernization memory remain in `docs/` and `AGENTS.md`.
- Bundled Codex-skill helpers are part of the public API and therefore must also follow the package naming convention. The installed-package entry points are `sn_get_codex_skill_path()` and `sn_install_codex_skill()`.
- User-facing differential-expression workflows should converge on `sn_find_de()` instead of ad hoc Seurat calls. Stored DE results now live in `object@misc$de_results`, which makes downstream marker visualization through `sn_plot_dot(features = "top_markers")` reproducible.
- Enrichment support should include both ORA and ranked-list GSEA so marker interpretation can operate on either thresholded gene sets or full ranked tables.
- Coverage reporting should be treated as a first-class maintenance signal: the repository now exposes a Codecov badge in the README and validates `covr::package_coverage()` through a dedicated GitHub Actions workflow.
- pkgdown deployment should be validated from a clean `site/` rebuild, because stale generated files can otherwise mask removed topics or old aliases in the published site.
- Species handling should be permissive for the common human/mouse case. When `species` is omitted, Shennong now attempts to infer it from feature names using `hom_genes` and mitochondrial naming patterns; only genuinely ambiguous inputs still require an explicit species.
- Optional packages that enhance core workflows should not block the default path. `SignatuR` and `catplot` are now treated as optional enhancements rather than mandatory dependencies for initialization, clustering, and plotting.
- Stored DE results should include lightweight provenance metadata (`schema_version`, package version, timestamp, thresholds, and assay/layer context) so downstream plotting and interpretation have a stable object contract.
- Seurat v5-facing DE APIs should expose `layer`, not `slot`. Shennong now maps the requested analysis layer onto whatever Seurat backend layer name the selected method requires internally, while keeping the public interface layer-based.
- Differential-expression backends should be selected with a single `method` argument scoped by `analysis`, rather than separate parameters for Seurat tests and pseudobulk engines. This keeps the public contract smaller and avoids redundant state.
- The experimental DE API can change incompatibly before the first stable release. The package version is therefore advanced to `0.2.0` to reflect the intentional public API change to `sn_find_de()`.
- pkgdown deployment in CI should not assume the package is already installed in the runner library. The workflow now performs an explicit `R CMD INSTALL .` before `pkgdown::deploy_to_branch()`.
- The first LLM integration should be a thin interpretation layer, not a statistical-analysis layer. Shennong now treats DE, enrichment, and clustering outputs as upstream evidence and builds prompts or narrative summaries on top of them.
- LLM features must remain provider-agnostic at the package boundary. The public API returns prompt bundles by default and only calls a model when the user supplies a provider function explicitly.
- Interpretation artifacts should be stored alongside other analysis products inside `object@misc`, using dedicated collections such as `enrichment_results` and `interpretation_results` instead of inventing a parallel state system.
- The `R/` directory should now be organized by durable workflow domains rather than by a mix of historical implementation slices. Preprocessing-related public APIs and helpers are co-located in `preprocessing.R`; analysis code is split into `analysis_clustering.R`, `analysis_de.R`, `analysis_enrichment.R`, and `analysis_metrics.R`; example data, IO, signatures, and package-management helpers each have dedicated modules.
- User-facing analysis features are not considered complete until their usage is discoverable. Any new or changed public workflow must update both pkgdown articles and the shipped Codex skill references in the same change set, especially when stored-result schemas or retrieval helpers change.
- pkgdown is part of the delivery surface, not just a derived artifact. After changing user-facing functionality, the site must be rebuilt locally to verify that the rendered reference and articles reflect the new behavior before the task is considered complete.
- Release notes are also part of the delivery contract. Any user-facing feature change must update `NEWS.md` in the same change set as the code, tests, pkgdown articles, and shipped skill references.
- Signature retrieval should be runtime-stable. Shennong now ships its own internal signature snapshot generated directly from `SignatuR` during development, and `sn_get_signatures()` reads that bundled data instead of the user's installed `SignatuR` package.
- End-user projects need a package-level way to bootstrap a governed analysis
  repository from packaged assets. Shennong exposes
  `sn_initialize_project()` as that entry point.
- Package examples must be check-safe even when they exercise credentialed integrations. The `sn_deconvolve_bulk(method = "cibersortx", cibersortx_dry_run = TRUE)` example now passes explicit placeholder credentials so dry-run validation stays local and deterministic.
- CI should install only the validated optional packages that are actually exercised in that workflow. Generic check and coverage jobs should not attempt to solve rare-cell or deconvolution backends that are only runtime-optional and are guarded by `check_installed()` / `check_installed_github()`.
- `sn_install_shennong(channel = "github")` should default to a conservative install path (`dependencies = FALSE`, `upgrade = "never"`) so end users can install the package even when some runtime-optional GitHub backends are temporarily unavailable upstream.
- The optional rare-cell backends `FiRE`, `CellSIUS`, and `EDGE` should be removed entirely rather than kept as nominally supported integrations. Their upstream installation paths add maintenance risk without being required for the default Shennong workflow, while `gini`, `gapclust`, `sca`, and `challenging_groups` still cover the rare-cell API surface.
- Composition summaries should preserve the source metadata's factor levels by default. Users often rely on those levels for plotting and reporting, so `sn_calculate_composition()` should not silently coerce grouped categorical columns to character unless a workflow truly requires it.
- Composition plotting is a common enough package workflow to justify a dedicated helper. `sn_plot_composition()` now covers stacked/dodged bar charts for proportions and counts instead of requiring each vignette and user script to rebuild the same `ggplot2` scaffolding.
- Composition significance should be computed on sample-level proportions, not pooled cell counts. Shennong now treats `sn_compare_composition()` as the correct API for group-wise composition fold changes and tests, avoiding pseudo-replication from per-cell tests.
- Neighborhood differential abundance should be exposed as a separate workflow from global composition. Shennong now wraps miloR through `sn_run_milo()` for local abundance shifts, while `sn_compare_composition()` remains the sample-level global proportion interface.
- Discrete color handling should be centralized rather than letting each plotting helper interpret palette names on its own. Shennong now resolves named palettes through a shared helper, with `Paired` as the default discrete palette and automatic interpolation when category counts exceed the base palette size.
- Palette management should be usable both interactively and programmatically. Printing helpers such as `show_all_palettes()` remain for discovery, but scripts should use `sn_list_palettes()` and `sn_get_palette()` as the stable API surface.
- Discrete and continuous palette resolution should follow one registry and one direction convention. Shennong now treats `sn_get_palette()` as the common palette API and keeps plot helpers thin wrappers around shared palette resolvers instead of mixing manual/discrete and distiller/continuous code paths.
- `sn_plot_barplot()` should cover the common statistical bar-chart use case directly instead of forcing users to drop down to raw `ggplot2` for every sample-level summary. Automatic replicate summarization plus optional error bars and raw-point overlays are now part of that helper's intended scope.

## 2026-03-28

- Project scaffolding should present one clear public entry point:
  `sn_initialize_project()`.
- Initialized analysis repositories should feel ready for immediate GitHub and RStudio use. The scaffold now includes a generated `<project>.Rproj` file plus a repository `.gitignore` alongside the existing governance template.
- Dependency discovery and installation should be package-driven rather than
  buried in `DESCRIPTION` or scattered across workflow docs. Shennong now
  exposes `sn_list_dependencies()` for inspection and
  `sn_install_dependencies()` for one-click installation of missing CRAN,
  Bioconductor, and GitHub dependencies.

- Default workflow remains local commit only; for installer ergonomics, `sn_install_shennong()` should prefer unified `source` / `ref` arguments and keep legacy aliases only for compatibility.

- Composition filtering and composition comparison need separate thresholds: `sn_calculate_composition()` should filter returned categories by count, while `sn_compare_composition()` should keep `min_cells` as a sample-level replicate filter and expose a simple ordered direction label derived from `log2_fc`.
- Any new exported function must update `_pkgdown.yml` in the same change set. Pkgdown reference-index omissions are now a known recurring failure mode, so local pre-push validation should include `pkgdown::build_reference_index()`.
- IO helpers should honor their documented contracts across supported custom formats. Named `row_names`, detected 10x spatial directories, existing `SingleCellExperiment` h5ad exports, and current `qs2` save signatures are maintained as regression-tested behavior.
- `sn_run_celltypist()` has two return modes by input type: Seurat inputs are updated with prediction metadata, while path inputs return the parsed prediction table because there is no Seurat object available for metadata writeback.
- Pkgdown articles should use PBMC3k as the main teaching dataset and explain the reasoning behind each Shennong call before showing code. Heavy analysis, external providers, and credentialed/containerized backends should be opt-in with `SHENNONG_RUN_VIGNETTES=true` so routine website builds remain fast.
- Visualization is a first-class Shennong feature and should have a dedicated article that demonstrates `sn_plot_*()` wrappers, palette APIs, stored-marker dot plots, and composition plots rather than leaving plotting examples scattered only across analysis articles.
- `sn_run_cluster()` should allow user feature priors through `hvg_features` rather than forcing all PCA features to come from internal HVG discovery. User-supplied features are validated against the object, recorded in `object@misc$hvg_selection`, and merged after blocked-gene filtering so intentional rare-population markers are not silently removed.
- SCTransform is a valid normalization path for Harmony integration in `sn_run_cluster()`; only scran remains single-dataset-only for now.
- Gene-symbol standardization should never produce `NA` feature names. When `HGNChelper` cannot provide a clean unambiguous replacement but the original symbol is present, preserve the original symbol and only drop truly missing or empty names.
- Rare-aware feature discovery in `sn_run_cluster()` should stay small and understandable. The clustering wrapper now keeps `gini` for expression-sparsity discovery and `local_markers` for rare-group marker discovery, while removing less common `local_hvg` and PCA-dependent `ciara` modes from this path.
- Manual and automatic feature priors should meet at one PCA feature set. `hvg_features` remains the explicit user prior, while `rare_feature_method` supplies automatic additions; both are merged with internal HVGs before scaling and PCA.
- Reference label projection is common enough to deserve a Shennong wrapper. `sn_transfer_labels()` keeps the standard Seurat anchor workflow, uses the query object as the first argument for pipe-friendly calls, and records metadata/provenance in one reproducible call.
- Coralysis reference mapping belongs in `sn_transfer_labels()` rather than a separate public function. The same query-first API now selects `method = "coralysis"` and expects a Coralysis-trained reference stored in `reference@misc$coralysis`, preserving the object-priority style while using the real Coralysis mapping backend.
- Simulation should be method-based because scDesign3 is not the only possible backend. `sn_simulate(method = "scdesign3")` is the public entry point, while `sn_simulate_scdesign3()` remains the backend-specific implementation wrapper.
- Label readability in `sn_plot_dim()` should be controllable. The white label halo/background is now optional through `label_halo`, and repel labels should remain repel-aware instead of being converted to fixed-position shadow text.
- Focused marker heatmaps are a first-class visualization helper. `sn_plot_heatmap()` is intentionally gene-list driven and validates/scales the selected features instead of trying to infer markers automatically. It supports both cell-level inspection and group-averaged summaries because users need both raw heterogeneity and compact marker-panel comparisons.
- Rasterized feature plots should preserve the ordinary ggplot point-size contract. `sn_plot_feature(raster = TRUE)` therefore prefers `ggrastr` post-rasterization when available instead of Seurat's `scattermore` branch, whose `pointsize` units do not match regular `geom_point(size = ...)` closely enough for tiny values such as `pt_size = 0.01`.
- Communication inference must delegate to real backends rather than improvised ligand-receptor joins. `sn_run_cell_communication()` wraps CellChat, NicheNet, and LIANA, and the NicheNet path requires explicit prior resources so provenance remains clear.
- Fast regulatory inference should use footprint-style resources rather than SCENIC-style GRN reconstruction. `sn_run_regulatory_activity()` uses decoupleR ULM with DoRothEA for TF activity and PROGENy for pathway activity, with optional group averaging for routine Seurat reports.
- Batch integration in `sn_run_cluster()` should be explicit rather than inferred as Harmony forever. `integration_method` now selects Harmony, Coralysis, Seurat CCA, or Seurat RPCA. Coralysis is restricted to log-normalized workflows because its public API expects a SingleCellExperiment with `logcounts`, while SCTransform integration remains Harmony-only for now.
