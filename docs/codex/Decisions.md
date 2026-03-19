# Shennong Modernization Decisions

Last updated: 2026-03-18

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
- End-user Codex skill materials should ship inside the installed package under `inst/codex/skills/shennong/`, while repository-only planning and modernization memory remain in `docs/` and `AGENTS.md`.
- Bundled Codex-skill helpers are part of the public API and therefore must also follow the package naming convention. The installed-package entry points are `sn_get_codex_skill_path()` and `sn_install_codex_skill()`.
- User-facing differential-expression workflows should converge on `sn_find_de()` instead of ad hoc Seurat calls. Stored DE results now live in `object@misc$de_results`, which makes downstream marker visualization through `sn_plot_dot(features = "top_markers")` reproducible.
- Enrichment support should include both ORA and ranked-list GSEA so marker interpretation can operate on either thresholded gene sets or full ranked tables.
- Coverage reporting should be treated as a first-class maintenance signal: the repository now exposes a Codecov badge in the README and validates `covr::package_coverage()` through a dedicated GitHub Actions workflow.
- pkgdown deployment should be validated from a clean `site/` rebuild, because stale generated files can otherwise mask removed topics or old aliases in the published site.
- Species handling should be permissive for the common human/mouse case. When `species` is omitted, Shennong now attempts to infer it from feature names using `hom_genes` and mitochondrial naming patterns; only genuinely ambiguous inputs still require an explicit species.
- Optional packages that enhance core workflows should not block the default path. `SignatuR` and `catplot` are now treated as optional enhancements rather than mandatory dependencies for initialization, clustering, and plotting.
- Stored DE results should include lightweight provenance metadata (`schema_version`, package version, timestamp, thresholds, and assay/layer context) so downstream plotting and interpretation have a stable object contract.
