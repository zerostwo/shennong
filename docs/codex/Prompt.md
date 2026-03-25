# Shennong Modernization Prompt

Last updated: 2026-03-18
Status: Frozen for this effort unless explicitly revised in `docs/codex/Decisions.md`

## Objective

Modernize the `Shennong` R package incrementally and safely toward current R package best practices while preserving public behavior unless a change is explicitly justified, documented, and validated inside this repository.

## Current Requested Scope

- Generalize example-data loading by introducing `sn_load_data()` and moving the current PBMC loader behavior behind it.
- Extend `sn_remove_ambient_contamination()` to support both SoupX and decontX with one unified API surface.
- Consolidate clustering into `sn_run_cluster()` by folding in the single-dataset workflow now handled by `sn_quick_cluster()`.
- Validate `sn_run_cluster()` against `pbmc1k` and `pbmc3k` example data for both single-dataset clustering and integration paths.

## Non-Negotiable Constraints

- Do not do a blind rewrite.
- Work in small, reviewable milestones.
- Preserve exported behavior where practical.
- Do not introduce unnecessary dependencies.
- Do not remove legacy behavior without tests or documentation supporting the change.
- Do not claim correctness without running repository-local validation commands.
- Update `docs/codex/Status.md` and `docs/codex/Decisions.md` after each meaningful change.

## Audit Snapshot

### Current package architecture

- The package is organized into large domain files under `R/`: IO, preprocessing, QC, integration, visualization, metrics, composition, enrichment, GRN, data, and utilities.
- Public API is primarily determined by `NAMESPACE`, with additional documented datasets in `R/data.R`.
- Documentation is roxygen-driven, with generated output in `man/`.
- Tests are minimal: one `testthat` file focused on `sn_calculate_composition()`.
- One vignette exists: `vignettes/clustering.Rmd`.
- `_pkgdown.yml` is present, but no CI workflow files are currently checked into the repository.

### Anti-patterns and inconsistencies

- `DESCRIPTION` declares only a small subset of packages actually used in `R/`.
- Several exported functions are missing `.Rd` files, while `sn_quick_cluster()` is documented but not exported.
- `R/io.R` depends heavily on non-exported `rio:::` internals and dynamically assigns importer functions into `globalenv()`.
- Some functions perform runtime package installation helpers or assume optional packages without clear metadata separation.
- The vignette performs remote downloads and loads heavyweight packages unconditionally, which is not check-safe.
- The repository currently has unrelated working tree changes in `.Rprofile`, `NAMESPACE`, and removed legacy files; those changes must not be overwritten blindly.

### Risky areas

- `R/io.R`: format guessing, import/export hooks, global assignment, and private `rio` internals.
- Seurat-dependent workflows in `R/preprocessing.R`, `R/quality_control.R`, `R/integrate.R`, and `R/visualization.R`.
- Hidden dependency usage across the package, including `Seurat`, `SeuratObject`, `Matrix`, `cli`, `curl`, `clusterProfiler`, `SoupX`, `scran`, `scDblFinder`, and others not currently declared in `DESCRIPTION`.

### Missing tests

- No tests currently cover IO helpers, preprocessing, QC, clustering, plotting wrappers, enrichment, or dataset loading.
- Current tests do not verify edge cases such as empty filtered results, metadata joins, or optional-column behavior beyond one function.

### Documentation gaps

- Missing manual pages for several exported functions, including `sn_read()`, `sn_write()`, `sn_enrich()`, `sn_set_path()`, plotting helpers, and multiple rio adapter exports.
- Examples and vignette chunks are not consistently safe for package checks.
- README badges reference CI that is not present in the checked-out repository state.

### Likely public API surface

- Data IO: `sn_read()`, `sn_write()`, `sn_load_data()`, `sn_add_data_from_anndata()`
- Preprocessing and QC: `sn_initialize_seurat_object()`, `sn_score_cell_cycle()`, `sn_normalize_data()`, `sn_standardize_gene_symbols()`, `sn_filter_genes()`, `sn_filter_cells()`, `sn_find_doublets()`, `sn_remove_ambient_contamination()`
- Clustering and annotation: `sn_run_cluster()`, `sn_run_celltypist()`, `sn_run_pyscenic()`
- Metrics and summaries: `sn_calculate_composition()`, `sn_calculate_lisi()`, `sn_calculate_rogue()`
- Visualization: `sn_plot_dim()`, `sn_plot_violin()`, `sn_plot_dot()`, `sn_plot_feature()`, `sn_plot_barplot()`, `sn_plot_boxplot()`, `show_all_palettes()`
- Signatures and utilities: `sn_get_signatures()`, `sn_get_species()`, `sn_check_file()`, `sn_set_path()`, `sn_enrich()`

## Target Architecture

- Keep the current domain-oriented layout initially, but make boundaries clearer between exported APIs and internal helpers.
- Normalize package metadata so `DESCRIPTION`, roxygen imports, and `NAMESPACE` match actual runtime usage.
- Treat IO adapters as internal implementation details unless there is a documented reason to keep them exported.
- Strengthen tests around critical behavior before refactoring high-risk internals.
- Keep vignettes, examples, and pkgdown configuration compatible with non-interactive package checks.
- The current package is still experimental, so simplifying the public surface
  is preferred over retaining deprecated wrappers or transitional aliases.

## Validation Policy

- Prefer repository-local commands already aligned with the package toolchain:
  - `Rscript -e 'testthat::test_local(stop_on_failure = TRUE)'`
  - `Rscript -e 'testthat::test_local(filter = "<target>", stop_on_failure = TRUE)'`
  - `Rscript -e 'if (requireNamespace("devtools", quietly = TRUE)) devtools::document() else stop("devtools not installed")'`
  - `R CMD build .`
  - `R CMD check --no-manual Shennong_*.tar.gz`
- Run the narrowest relevant validation after each meaningful change.
- Run broader package validation when metadata, exports, or documentation generation changes.
