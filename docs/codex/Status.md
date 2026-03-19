# Shennong Modernization Status

Last updated: 2026-03-18
Current milestone: CI, pkgdown, and package-maintenance hardening

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
- Documented datasets: `marker_genes`, `hom_genes`, and `SignatuR`

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

## Remaining High-Priority Work

- Review the remaining stale worktree deletions and doc/export mismatches outside this task, especially removed `grn` and legacy upstream files.
- Decide whether to suppress or explicitly document known Seurat/HGNChelper runtime warnings in tests and smoke paths.
- Decide whether `sn_remove_ambient_contamination()` should also gain an explicit input `assay`/`layer` pair for Seurat objects, matching the newer downstream layer-aware APIs.
- Add deeper coverage for enrichment, DE, CellTypist, and other heavyweight optional integrations when stable fixtures are available.
- Revisit `scDblFinder`-originating warnings in tests; they are upstream warnings today, but some deprecation notices may require argument updates in a future compatibility pass.
