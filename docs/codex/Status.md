# Shennong Modernization Status

Last updated: 2026-03-19
Current milestone: DE API consolidation and CI deployment hardening

## Snapshot

- Repository audit completed.
- Codex memory files and repository guidance are in place and now trackable in git.
- The task-specific API consolidation around `sn_load_data()`, `sn_remove_ambient_contamination()`, and `sn_run_cluster()` is now implemented and validated.
- Seurat analysis helpers that operate on count-like layers now accept explicit `assay`/`layer` inputs where that behavior matters operationally.
- Public Seurat preprocessing and clustering APIs now expose explicit normalization-method controls instead of overloading the old `pipeline` wording.
- decontX-based ambient-RNA removal now preserves original counts for zeroed cells by default, can optionally drop them, and writes corrected per-cell count/feature metadata.
- The package now exposes explicit version-checking and install helpers for CRAN-or-GitHub release management.
- pkgdown now builds to `site/` instead of `docs/`, so the documentation site no longer collides with `docs/codex/` project memory.
- Exported functions now have generated help coverage through roxygen-managed `.Rd` files.
- The package now has a generated `README.md`, guarded vignette execution, repository-local CI workflows, and a passing `R CMD check --no-manual`.
- The local test surface now covers composition, data loading, clustering, IO/visualization helpers, preprocessing, and utility helpers.
- A coverage workflow and README coverage badge are now in place, and local `covr::package_coverage()` succeeds.
- pkgdown has been revalidated with a clean rebuild to `site/`, including the rendered PBMC article and bundled Codex-skill reference pages.
- Species detection now falls back to automatic inference from feature names using `hom_genes` plus mitochondrial naming patterns, so core workflows no longer require an explicit `species` argument in common human/mouse inputs.
- pkgdown now includes a second evaluated article covering layer-aware workflows and stored DE metadata.
- The differential-expression API has now been consolidated around `analysis`, `layer`, and `method`, removing the old mixed `slot` / `test_use` / `pseudobulk_method` vocabulary.
- `sn_find_de()` now supports `limma` pseudobulk analysis and optional `COSGR` marker discovery.
- The pkgdown workflow now installs the package before `deploy_to_branch()`, addressing the Actions failure where the deployment worktree could not `library(Shennong)`.
- A first LLM-ready interpretation layer now exists on top of the existing analysis stack, including enrichment storage, evidence preparation, prompt construction, and optional provider-backed writing helpers.
- The `R/` source tree has now been reorganized around stable domains: preprocessing, clustering, DE, enrichment, metrics, example data, IO, package tools, signatures, visualization, utilities, and interpretation.
- Stored analysis outputs are now treated as a first-class user interface: Seurat objects can expose DE, enrichment, and interpretation artifacts through explicit listing and retrieval helpers instead of requiring direct access to `object@misc`.

## Concise Audit

### Current package architecture

- `R/` contains 13 broad modules, with major concentration in `visualization.R`, `integrate.R`, `quality_control.R`, and `io.R`.
- `NAMESPACE` defines the main exported API surface.
- `man/` is present but incomplete relative to exports.
- `tests/testthat/` contains a single file covering `sn_calculate_composition()`.
- `vignettes/clustering.Rmd` exists and exercises a remote-data workflow.
- `_pkgdown.yml` exists, but no CI workflow directory is present in the checked-out repository.

### Anti-patterns and inconsistencies

- Domain-oriented source layout still mixes public APIs and internal helpers inside large files.
- The working tree is already dirty before modernization work begins.

### Risky areas

- IO helpers and adapter registration in `R/io.R`
- Seurat-heavy analysis helpers across preprocessing, QC, integration, and plotting
- Hidden optional dependencies not clearly separated from core package requirements

### Missing tests

- Enrichment, differential-expression, and CellTypist workflows still have lighter coverage than clustering, IO, and preprocessing.
- Heavy optional-package paths are mostly validated through package checks and targeted skips rather than exhaustive matrix tests.

### Documentation gaps

- Dataset-level narrative docs and long-form articles can still be expanded.
- Some runtime warnings from upstream packages remain visible in tests, even though validation passes.

### Likely public API surface

- Exported functions from `NAMESPACE`
- Documented datasets: `marker_genes` and `hom_genes`

### Proposed target architecture

- Keep domain-oriented source files initially, but make public-versus-internal boundaries explicit
- Normalize `DESCRIPTION`, roxygen imports, and `NAMESPACE`
- Add lightweight tests before refactoring high-risk code paths
- Make documentation, pkgdown, and vignettes safe for standard package workflows

## Validation Log

- Passed: `Rscript -e 'testthat::test_local(stop_on_failure = TRUE)'`
- Passed: `Rscript -e 'testthat::test_local(filter = "composition", stop_on_failure = TRUE)'`
- Passed: `Rscript -e 'testthat::test_local(stop_on_failure = TRUE)'` after the composition change
- Passed: `Rscript -e 'devtools::document()'`
- Passed: `Rscript -e 'testthat::test_local(filter = "data_loading|clustering", stop_on_failure = TRUE)'`
- Passed: `Rscript -e 'testthat::test_local(stop_on_failure = TRUE)'`
- Passed: `Rscript inst/scripts/smoke_cluster_pbmc.R`
- Passed: `R CMD build .`
- Passed: `Rscript -e 'rmarkdown::render("README.Rmd", output_format = "github_document", quiet = TRUE)'`
- Passed: `R CMD check --no-manual Shennong_0.1.0.tar.gz`
- Passed: `Rscript -e 'testthat::test_local(filter = "clustering|preprocessing", stop_on_failure = TRUE)'`
- Passed: `Rscript -e 'testthat::test_local(stop_on_failure = TRUE)'` after `assay`/`layer` support
- Passed: `Rscript -e 'if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak", repos = "https://cloud.r-project.org"); pak::lockfile_create(".", lockfile = tempfile("pkg-", fileext = ".lock"), upgrade = TRUE)'`
- Passed: `R CMD build .` after `methods::slot` namespace cleanup
- Passed: `R CMD check --no-manual Shennong_0.1.0.tar.gz` with `Status: OK`
- Passed: `Rscript -e 'testthat::test_local(filter = "clustering|preprocessing", stop_on_failure = TRUE)'` after the decontX and normalization API refinement
- Passed: `Rscript -e 'testthat::test_local(stop_on_failure = TRUE)'` after the decontX and normalization API refinement
- Passed: `R CMD build .` after the decontX and normalization API refinement
- Passed: `R CMD check --no-manual Shennong_0.1.0.tar.gz` with `Status: OK` after the decontX and normalization API refinement
- Passed: `Rscript -e 'devtools::document()'` after adding version/install helpers and pkgdown metadata updates
- Passed: `Rscript -e 'testthat::test_local(filter = "versioning|preprocessing|clustering|utils", stop_on_failure = TRUE)'`
- Passed: `Rscript -e 'rmarkdown::render("README.Rmd", output_format = "github_document", quiet = TRUE)'`
- Passed: `Rscript -e 'pkgdown::build_site(new_process = FALSE, install = FALSE)'`
- Passed: `Rscript -e 'testthat::test_local(stop_on_failure = TRUE)'`
- Passed: `R CMD build .`
- Passed: `R CMD check --no-manual Shennong_0.1.0.tar.gz` with `Status: OK` after the CI/pkgdown hardening
- Passed: `Rscript -e 'testthat::test_local(filter = "io_visualization|codex_skill|signatures", stop_on_failure = TRUE)'`
- Passed: `Rscript -e 'cov <- covr::package_coverage(); cat(covr::percent_coverage(cov), "\n")'` with `53.49386`
- Passed: `Rscript -e 'Sys.setenv(SHENNONG_RUN_VIGNETTES = "true"); pkgdown::build_site(new_process = FALSE, install = TRUE)'`
- Passed: `R CMD build .` after coverage and pkgdown updates
- Passed: `R CMD check --no-manual Shennong_0.1.0.tar.gz` with `Status: OK` after coverage and pkgdown updates
- Passed: `Rscript -e 'testthat::test_local(filter = "utils|preprocessing|clustering|de_enrich", stop_on_failure = TRUE)'`
- Passed: `Rscript -e 'testthat::test_local(stop_on_failure = TRUE)'` with `PASS 112`
- Passed: `R CMD build .` after species inference and article updates
- Passed: `R CMD check --no-manual Shennong_0.1.0.tar.gz` with `Status: OK` after species inference and article updates
- Passed: `Rscript -e 'devtools::document()'` after the DE API consolidation
- Passed: `Rscript -e 'testthat::test_local(filter = "de_enrich", stop_on_failure = TRUE)'`
- Passed: `Rscript -e 'testthat::test_local(stop_on_failure = TRUE)'` with `PASS 118`
- Passed: `Rscript -e 'Sys.setenv(SHENNONG_RUN_VIGNETTES = "true"); pkgdown::build_site(new_process = FALSE, install = TRUE)'`
- Passed: `Rscript -e 'testthat::test_local(filter = "interpretation", stop_on_failure = TRUE)'`
- Passed: `Rscript -e 'testthat::test_local(stop_on_failure = TRUE)'` with `PASS 139`
- Passed: `R CMD build .` after adding the interpretation layer

## Completed In This Iteration

- Replaced the placeholder `AGENTS.md` with repository-specific development guidance.
- Added the codex memory files under `docs/codex/`.
- Narrowed `.gitignore` so `docs/codex/` is not hidden behind the repo-wide `docs` ignore rule.
- Updated `sn_calculate_composition()` so it now warns when `additional_cols` are not constant within a `group_by` group, matching the existing documentation contract.
- Added regression tests for inconsistent `additional_cols` behavior and the `min_cells` no-results error path.
- Added `sn_load_data()` as the generalized example-data loader and retained `sn_load_pbmc()` as a compatibility wrapper.
- Removed the separate `sn_quick_cluster()` implementation and folded single-dataset clustering into `sn_run_cluster()`.
- Reworked `sn_remove_ambient_contamination()` into a unified SoupX/decontX interface with shared input and cluster handling.
- Added new unit tests for data loading, clustering, and ambient contamination, plus a PBMC smoke script covering `pbmc1k` and `pbmc3k`.
- Corrected `sn_initialize_seurat_object()` so plain base matrices are accepted as count inputs.
- Updated `DESCRIPTION`, `NAMESPACE`, `.Rd` files, and the clustering vignette to match the new APIs.
- Added roxygen coverage for previously undocumented exported helpers, then regenerated `man/` and `NAMESPACE`.
- Rebuilt `README` from `README.Rmd`, expanded the user-facing quick start, and removed misleading coverage-badge drift.
- Added `CONTRIBUTING.md`, a pkgdown GitHub Actions workflow, and Conventional Commit guidance in `AGENTS.md`.
- Reworked `sn_read()`/`sn_write()` away from `rio:::` internals onto public `rio` APIs plus explicit custom-format dispatch.
- Fixed multiple stale roxygen mismatches and example failures so package documentation now survives `R CMD check --no-manual`.
- Added internal helpers for validating, reading, combining, and temporarily activating Seurat assay layers, including split-layer support on merged Seurat v5 objects.
- Updated `sn_run_cluster()`, `sn_find_doublets()`, `sn_filter_genes()`, `sn_normalize_data()`, `sn_run_celltypist()`, and `sn_calculate_rogue()` to accept explicit `assay`/`layer` parameters.
- Added regression coverage proving non-default layers can drive clustering, doublet detection, gene filtering, and normalization without overwriting the original `counts` layer.
- Declared GitHub remotes for `BPCells`, `catplot`, `lisi`, `ROGUE`, and `SignatuR`, then locally reproduced a passing `pak::lockfile_create()` to match the failing Actions dependency-resolution step.
- Imported `methods::slot` and `slot<-` explicitly so `R CMD check` no longer reports namespace notes for command logging.
- Refined `sn_remove_ambient_contamination()` so decontX-generated zero-count cells no longer receive synthetic counts. Affected cells are restored from the original counts by default, can optionally be removed, and now produce corrected `nCount_*` / `nFeature_*` metadata columns.
- Updated `sn_initialize_seurat_object()` to accept optional `sample_name` and `study` metadata at object creation time.
- Replaced the overloaded `pipeline` argument in `sn_run_cluster()` with `normalization_method`, while retaining a deprecated compatibility alias. The same normalization vocabulary is now used by `sn_normalize_data()`.
- Hardened cell-cycle scoring so clustering skips that step with a warning when the requested species markers do not overlap the selected assay features, avoiding a low-level `sample.int()` failure.
- Added `sn_check_version()` and `sn_install_shennong()` so users can check whether their installation is current and install either the CRAN release or GitHub development version once both channels exist.
- Changed pkgdown output to `site/`, expanded the reference index to cover all user-facing documented topics, and locally validated a successful site build.
- Adjusted GitHub Actions dependency installation so CI no longer tries to solve or unpack the GitHub-only `BPCells` Suggests dependency through `pak` during the generic dependency bootstrap step.
- Replaced repeated Seurat command-logging snippets with a shared internal helper while preserving the public `sn_*` command names recorded in `object@commands`.
- Removed automatic GitHub package installation side effects from `check_installed_github()`. Missing GitHub-only optional packages now produce an explicit install instruction instead of mutating the user library during analysis.
- Added `hvg_method` to `sn_run_cluster()` so integration workflows can now choose between global HVG selection and per-batch HVG selection followed by ranked merging.
- Reduced package-side test warnings to zero in the clustering, preprocessing, and IO/visualization slices by removing avoidable coercion, metadata-overwrite, assay-default, and upstream signature warnings.
- Replaced `hvg_method` with `hvg_group_by` on `sn_run_cluster()`. The clustering API now treats group-aware HVG selection as a direct metadata-column choice rather than an abstract method enum.
- Added a PBMC workflow smoke script and a pkgdown article that run real `pbmc1k`/`pbmc3k` analyses covering single-sample clustering, Harmony integration, composition summaries, LISI scoring, marker discovery, and GO enrichment with concrete outputs.
- Removed the vendored `SignatuR` package dataset from `data/` because its `data.tree::Node` / R6 structure is not compatible with package lazydata installation. `sn_get_signatures()` now loads the optional dataset directly from the installed `SignatuR` package when needed.
- Added a distributable end-user Codex skill under `inst/codex/skills/shennong/`, including package overview, function map, object model, workflow recipes, and packaged helper scripts. Repository-only memory remains in `docs/codex/`.
- Added `sn_get_codex_skill_path()` and `sn_install_codex_skill()` so installed-package users can discover and copy the bundled skill into their local agent skill directory.
- Expanded the analysis layer around `sn_find_de()`, `sn_plot_dot()`, and `sn_enrich()`:
  `sn_find_de()` now supports marker discovery, direct contrasts, and pseudobulk DE with stored results in `object@misc$de_results`; `sn_plot_dot()` can reuse stored top markers; `sn_enrich()` now supports both ORA and GSEA.
- The PBMC vignette now compares pre- and post-integration embeddings, uses Shennong plotting helpers, and shows stored-marker dot plots plus GSEA-driven interpretation.
- Added targeted unit tests for bundled Codex skill installation, stored-marker/GSEA workflows, signature retrieval, and Seurat plotting helpers, bringing the local full-suite pass count to 105.
- Added a GitHub Actions coverage workflow and surfaced the Codecov badge in the README.
- Performed a clean pkgdown rebuild into `site/` to remove stale pages and confirm the article/reference site matches the current package sources.
- Added species inference through `sn_get_species()` for Seurat objects, matrices, and feature vectors, with `sn_initialize_seurat_object()` now attempting inference automatically before deciding whether to compute species-specific QC metrics.
- Made built-in signature retrieval resilient to missing `SignatuR` by providing package-local fallback signatures for the categories required by the core workflows.
- Expanded the stored DE schema with package/version/timestamp metadata and exposed it in a new evaluated pkgdown article, `layered-workflows`.
- Consolidated `sn_find_de()` onto a Seurat v5-style `layer` interface, removed the old `slot` / `test_use` / `pseudobulk_method` arguments, and introduced a single method selector covering Seurat, pseudobulk, and COSGR-backed marker workflows.
- Added `limma` as a supported pseudobulk backend and added targeted tests for both `limma` and optional `COSGR` marker discovery.
- Hardened the pkgdown GitHub Actions deployment job by installing the package before `pkgdown::deploy_to_branch()`.
- Added an interpretation-layer module that keeps LLM prompting separate from analysis. It introduces stored enrichment metadata, evidence bundles for annotation/DE/pathway/writing tasks, prompt builders, and optional provider-backed writing helpers.
- Consolidated the source layout so workflow-adjacent functions now live together: `preprocessing.R` owns initialization, normalization, QC, doublet detection, ambient correction, and species inference; analysis files are split into clustering, DE, enrichment, and metrics; package-maintenance helpers are grouped under `package_tools.R`.

## Remaining High-Priority Work

- Review the remaining stale worktree deletions and doc/export mismatches outside this task, especially removed `grn` and legacy upstream files.
- Decide whether to suppress or explicitly document known Seurat/HGNChelper runtime warnings in tests and smoke paths.
- Decide whether `sn_remove_ambient_contamination()` should also gain an explicit input `assay`/`layer` pair for Seurat objects, matching the newer downstream layer-aware APIs.
- Add deeper coverage for enrichment, DE, CellTypist, and other heavyweight optional integrations when stable fixtures are available.
- Revisit `scDblFinder`-originating warnings in tests; they are upstream warnings today, but some deprecation notices may require argument updates in a future compatibility pass.
