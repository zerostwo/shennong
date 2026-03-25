# Shennong Modernization Decisions

Last updated: 2026-03-25

## 2026-03-25

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
- `sn_initialize_codex_project()` is the package entry point that materializes the shipped project template into a user project. `sn_initialize_project()` remains as a convenience wrapper.
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
- The GitHub Actions dependency-resolution failure should be fixed at the package metadata level rather than by weakening CI. GitHub-only optional dependencies are now declared in `DESCRIPTION` under `Remotes`, and this was validated locally with `pak::lockfile_create()`.
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
- End-user projects need a package-level way to bootstrap a governed analysis repository from packaged assets. Shennong now exposes `sn_initialize_codex_project()` for the full project template and retains `sn_initialize_project()` as a wrapper for convenience.
