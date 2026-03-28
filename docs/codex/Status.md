# Shennong Modernization Status

Last updated: 2026-03-25

## 2026-03-25

- Added a repository-local pre-push helper at `scripts/check-prepush.R`. It
  now runs the maintainer-standard validation sequence in one place:
  `devtools::document()`, an optional targeted `testthat` pass, the full local
  test suite, `R CMD build`, and `R CMD check --no-manual`.
- Added `sn_assess_qc()` as a lightweight QC-reporting layer for Seurat
  objects. It now summarizes per-sample QC status, reports current risk signals
  such as failed-QC fractions, doublet rates, and decontamination zero-count
  rates, and can compare a filtered object against a reference object to
  quantify low-quality-cell removal, doublet removal, and clean-cell
  retention. Reports can be stored under `object@misc$qc_assessments`.
- Reviewed `sn_filter_genes()` / `sn_filter_cells()` and fixed several QC-edge
  cases. Cell filtering now validates the supported method, treats zero-MAD
  groups as stable rather than failing every cell, and checks plotting
  dependencies explicitly; gene-filter threshold plots now stay well-defined
  even when `min_cells` exceeds the number of cells.
- Reworked several internal sparse-matrix hot paths for better performance and
  lower peak memory use. Pseudobulk DE aggregation now uses sparse grouped
  matrix multiplication instead of `t(as.matrix(...))` + `rowsum()`, exact kNN
  fallbacks now compute neighbors blockwise instead of allocating a full
  distance matrix, split Seurat layer merges now assemble one sparse matrix
  from triplets rather than repeated indexed writes, and gene-symbol
  standardization now uses bundled annotation data plus grouped row
  aggregation.
- Continued the refactor by removing duplicated kNN cleanup / Annoy lookup
  logic across clustering and metrics. Shared helpers now own self-neighbor
  removal and Annoy neighbor extraction, and grouped HVG selection now
  processes one metadata level at a time instead of materializing
  `Seurat::SplitObject()` output for every group at once.
- Standardized package progress logging through internal wrappers instead of a
  mix of ad hoc `logger::log_*()`, `message()`, and `cli` progress calls.
  Informational notifications now flow through shared helpers with glue-based
  interpolation, while direct warnings/errors remain reserved for actual
  conditions that callers may need to handle.
- Added regression coverage for the new sparse aggregation helpers and the
  blockwise exact-kNN helper in `tests/testthat/test_utils.R`.
- Expanded `R/analysis_metrics.R` from a minimal LISI/ROGUE/composition module
  into a broader integration-assessment layer. New exported helpers now cover
  silhouette widths, graph connectivity, PCR batch scoring, clustering
  agreement, isolated-label preservation, cluster entropy/purity diagnostics,
  rare/difficult-group diagnostics, and an aggregate `sn_assess_integration()`
  wrapper that returns summary, per-cell, and per-group outputs.
- Added rare-cell-aware clustering support in `R/analysis_clustering.R`. The
  package now exposes `sn_detect_rare_cells()` plus optional
  `sn_run_cluster(rare_feature_method = ...)` feature augmentation using native
  Gini scoring, local rare-group HVGs/markers, and an optional CIARA backend.
- Added structured metrics tests that exercise graph reuse, Annoy fallback,
  exact-kNN fallback, PCR improvement against a baseline reduction, aggregate
  summary scoring, and rare-group flagging on synthetic integrated embeddings.
- Updated the README, PBMC clustering vignette, and shipped single-cell skill
  reference so the new metrics workflow is visible at the user level instead of
  remaining an undocumented helper layer.
- Added bulk deconvolution support in a new `R/analysis_deconvolution.R`
  module. Shennong now supports local BayesPrism runs and local CIBERSORTx
  container workflows, with stored results available under
  `object@misc$deconvolution_results`.
- Reorganized pkgdown around workflow stages instead of mixed implementation
  topics. Added new end-to-end articles for preprocessing/QC, metrics and
  diagnostics, annotation/pathway analysis, and composition/comparative
  summaries, while retitling the clustering and interpretation articles to
  match the new taxonomy.
- Refactored `sn_enrich()` around a single primary `x` input plus a formula
  contract for enrichment semantics. `gene_clusters = gene ~ cluster` now
  drives grouped ORA, while `gene_clusters = gene ~ log2fc` or a named numeric
  vector drives ranked GSEA without a separate `score_col`-only path.
- Expanded `sn_enrich()` database dispatch so one call can run against mixed
  requests such as `c("H", "GOBP", "KEGG", "C2:CP:REACTOME")`, while also
  exposing `collection` and `subcollection` arguments aligned with `msigdbr`.
- Fixed enrichment result persistence so `prefix` / `outdir` save one `.rds`
  per requested database with stable filenames, and post-filtered all returned
  results by raw p-value rather than the upstream adjusted-p cutoff behavior.
- Added targeted enrichment tests for formula auto-detection, Seurat-stored DE
  reuse, multi-database storage/output, helper validation, and raw p-value
  filtering. Local file-level coverage for `R/analysis_enrichment.R` now
  reaches `78.77095%`, up from the earlier ~44% baseline.
- Removed the legacy `sn_load_pbmc()` wrapper and started collapsing package
  data documentation into a single `R/data.R` source. New built-in example data
  objects `pbmc_small` and `pbmc_small_raw` are now built from sampled
  `pbmc1k`/`pbmc3k` assets through `data-raw/build_shennong_example_data.R`.
- README is being rewritten around a standard package-facing structure, using
  the new built-in PBMC sample data instead of governance text or network-backed
  examples.

## 2026-03-24

- Added a new `data-raw` build path for bundled human/mouse GENCODE gene
  annotations and wired that package data into `sn_filter_genes()` through
  `gene_class` and exact `gene_type` filtering. The packaged snapshot now keeps
  not only gene names/IDs/types but also genome context fields such as seqname,
  source, coordinates, strand, status/source tags, and annotation level.
- Replaced the custom signature-registry tree-editing logic with wrappers that
  materialize the packaged snapshot as a `SignatuR` object and then call the
  upstream `SignatuR::AddNode()`, `AddSignature()`, and `RemoveSignature()`
  APIs before serializing the result back to Shennong's packaged snapshot.
- Switched all Harmony installation paths from the CRAN package line to the
  upstream `immunogenomics/harmony@harmony2` developer branch. This includes
  `DESCRIPTION` remotes, GitHub Actions dependency setup, README install
  guidance, and the clustering vignette note so local, CI, and pkgdown runs all
  target the same integration backend.

## 2026-03-21

- Fixed the failing `R CMD check` example/doc issues in the interpretation layer by regenerating `man/` from the current roxygen source. `sn_list_results()` examples now normalize before requesting the `data` layer, and the `background` / `output_format` arguments are documented again for the interpretation-writing wrappers.
- Reworked the packaged project-template scaffold so empty analysis directories are created from an explicit `inst/codex/project-template/directories.txt` manifest instead of hidden `.gitkeep` files. This removes the package hidden-file NOTE while preserving initialized-project directory creation.
- Refactored the Codex architecture into a clean package-vs-project split. The package root remains an R package, shipped initialized-project governance now lives under `inst/codex/project-template/`, shipped package-usage skills now live under `inst/codex/package-skills/`, and `sn_initialize_codex_project()` now scaffolds user projects from the packaged template assets.
- Expanded the test suite around Seurat object state, stored-result contracts, IO dispatch, packaged project scaffolding, and metrics/clustering integrations. Local `covr::package_coverage()` now reaches 70.30%, up from the high-60s baseline during this task.
- Replaced the runtime `SignatuR` dependency with a bundled `shennong_signature_catalog` dataset stored under `data/` and rebuilt directly from the upstream `SignatuR` dataset during development. Signature retrieval is now package-stable, tree-structured, and no longer depends on opaque `sysdata` storage.
- Added signature-management interfaces for package maintenance: `sn_list_signatures()`, `sn_add_signature()`, `sn_update_signature()`, and `sn_delete_signature()`. They now operate on the packaged snapshot and `data-raw/build_shennong_signatures.R` rebuilds that snapshot directly from `SignatuR`.

## 2026-03-19

- Fixed `sn_plot_dot()` so the optional `catplot` theme does not add a second aspect-ratio constraint on top of `coord_fixed()`.
- Replaced runtime `SignatuR` lookups with a bundled Shennong signature snapshot plus build script so maintainers can regenerate package data reproducibly from the upstream source.
- Added `sn_initialize_project()` as a convenience wrapper around the packaged initialized-project template. When used with agent governance, the created project now uses `AGENTS.md`, `memory/`, `docs/standards/`, `skills/`, `runs/`, and `results/` rather than the older flat `docs/codex/` layout inside user projects.
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
- Passed: `Rscript -e 'testthat::test_local(filter = "utils|preprocessing|interpretation|io_visualization", stop_on_failure = TRUE)'`
- Passed: `Rscript -e 'testthat::test_local(filter = "data_io_contracts|composition|clustering", stop_on_failure = TRUE)'`
- Passed: `Rscript -e 'testthat::test_local(filter = "codex_skill", stop_on_failure = TRUE)'`
- Passed: `Rscript -e 'testthat::test_local(filter = "codex_skill", stop_on_failure = TRUE)'` after the project-template directory-manifest change
- Passed: `Rscript -e 'testthat::test_local(filter = "de_enrich", stop_on_failure = TRUE)'`
- Passed: `Rscript -e 'testthat::test_local(filter = "metrics", stop_on_failure = TRUE)'`
- Passed: `Rscript -e 'testthat::test_local(filter = "data_loading", stop_on_failure = TRUE)'`
- Passed: `Rscript -e 'testthat::test_local(filter = "data_io_contracts", stop_on_failure = TRUE)'`
- Passed: `Rscript -e 'testthat::test_local(stop_on_failure = TRUE)'` after the coverage expansion
- Passed: `Rscript -e 'cov <- covr::package_coverage(type = c("tests")); cat(covr::percent_coverage(cov), "\n")'` with `70.29831`
- Passed: `Rscript -e 'if (requireNamespace("devtools", quietly = TRUE)) devtools::document() else stop("devtools not installed")'`
- Passed: `R CMD build .`
- Passed: `R CMD check --no-manual -o /tmp/shennong-check-current Shennong_0.1.1.tar.gz`
- Passed: `Rscript -e 'testthat::test_local(stop_on_failure = TRUE)'` after the CI-fix/template-manifest change (`PASS 372`)
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
- Added a distributable Codex asset layer inside the package. Package-usage skills now live under `inst/codex/package-skills/`, while initialized-project governance and skills now live under `inst/codex/project-template/`. Repository-only memory remains in `docs/codex/`.
- Added `sn_get_codex_skill_path()` and `sn_install_codex_skill()` so installed-package users can discover and copy the bundled skill into their local agent skill directory.
- Tightened the GitHub Actions dependency bootstrap back to the minimal baseline needed for checks and coverage after `pak::lockfile_create()` failed on unsupported optional GitHub backends such as `FiRE` and `GapClust`.
- Verified locally that `pak::lockfile_create(".", lockfile = ".github/pkg.lock", upgrade = TRUE)` now succeeds again with the revised CI dependency set.
- Updated `sn_install_shennong(channel = "github")` so GitHub installs default to `dependencies = FALSE` and `upgrade = "never"`, preventing optional GitHub `Suggests` from breaking package installation.
- Added regression tests covering the default and overridden GitHub-install argument paths for `sn_install_shennong()`.
- Removed the `FiRE`, `CellSIUS`, and `EDGE` rare-cell backends from `sn_detect_rare_cells()`, package metadata, tests, and vignette references so Shennong no longer depends on those unstable optional integrations.
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
- Fixed the `sn_deconvolve_bulk()` example so the `cibersortx` dry-run example path supplies placeholder credentials during checks.
- Added the missing GitHub-hosted optional packages (`BayesPrism`, `CIARA`, `CellSIUS`, `EDGE`, `FiRE`, `GapClust`) to `DESCRIPTION`/CI installation paths so Actions can satisfy the same optional backends referenced in package code and examples.
- Qualified the remaining `setNames()` call under `stats::setNames()` so `R CMD check` no longer reports a global-function note in the rare-cell workflow.
- `sn_calculate_composition()` now preserves factor columns from source metadata and can order a single grouping column by a chosen category level, reducing the amount of manual post-processing needed before plotting.
- Added `sn_plot_composition()` as a higher-level plotting helper for stacked composition-style bar charts so users no longer need to rebuild the same `ggplot2::geom_col()` pattern for composition, QC status, or similar categorical summaries.
- Removed a duplicate DE-specific misc-storage helper and now route DE result persistence through the shared `.sn_store_misc_result()` path used by other stored-result workflows.
- Added `sn_compare_composition()` for the statistically correct replicate-aware composition workflow: compute sample-level proportions first, then estimate group differences, log2 fold changes, and optional Wilcoxon/FDR summaries per category.
- Added `sn_run_milo()` as a neighborhood differential abundance interface over miloR, using Seurat embeddings plus sample/group metadata and returning neighborhood-level DA statistics with optional annotation labels.
- Standardized discrete palette handling across plotting helpers so string palettes such as `Paired` work consistently for both fill and color scales, including automatic interpolation when the number of categories exceeds the base palette length. The main plot helpers also now share panel-size and axis-label options.
- Added user-facing palette discovery and retrieval helpers, `sn_list_palettes()` and `sn_get_palette()`, so scripts can reuse the same palette resolution logic as the plotting helpers instead of relying on printed output only.
- Extended palette normalization to continuous color scales as well, so `sn_plot_feature()` and `sn_plot_dot()` now consume the same palette registry and direction semantics as the discrete plotting helpers.
- Expanded `sn_plot_barplot()` from a thin `geom_col()` wrapper into a more generally useful grouped-summary plot helper with automatic repeated-observation summarization, optional SD/SE error bars, and optional jittered raw points.

## Remaining High-Priority Work

- Review the remaining stale worktree deletions and doc/export mismatches outside this task, especially removed `grn` and legacy upstream files.
- Decide whether to suppress or explicitly document known Seurat/HGNChelper runtime warnings in tests and smoke paths.
- Decide whether `sn_remove_ambient_contamination()` should also gain an explicit input `assay`/`layer` pair for Seurat objects, matching the newer downstream layer-aware APIs.
- Add deeper coverage for enrichment, DE, CellTypist, and other heavyweight optional integrations when stable fixtures are available.
- Revisit `scDblFinder`-originating warnings in tests; they are upstream warnings today, but some deprecation notices may require argument updates in a future compatibility pass.
