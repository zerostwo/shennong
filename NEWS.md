All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

# Unreleased

### Added

- `sn_find_de()` now supports pseudobulk differential expression with
  `limma` in addition to `DESeq2` and `edgeR`.
- `sn_find_de()` now supports marker discovery with `COSGR` when the optional
  GitHub package is installed.
- Stored-result discovery and retrieval helpers were added:
  `sn_list_results()`, `sn_get_de_result()`, `sn_get_enrichment_result()`, and
  `sn_get_interpretation_result()`.
- The interpretation layer now supports user-supplied background context and
  dual output styles for either model-facing prompt bundles or human-readable
  summaries.

### Changed

- `sn_find_de()` now uses a single `method` argument instead of separate
  `test_use` and `pseudobulk_method` parameters.
- `sn_find_de()` now uses `layer` consistently for Seurat v5 workflows and no
  longer exposes the legacy `slot` argument.
- `sn_enrich()` now stores enrichment results in
  `object@misc$enrichment_results[[store_name]]` when a Seurat object is
  supplied, aligning enrichment with the existing stored DE workflow.
- pkgdown articles, shipped Codex skill references, and `NEWS.md` are now
  treated as required deliverables for any user-facing workflow change.
- The pkgdown deployment workflow now installs the current package before
  publishing the site, so reference generation can load `Shennong` inside the
  temporary deployment worktree.

### Fixed

- Fixed pkgdown GitHub Actions deployment when `deploy_to_branch()` runs in a
  clean runner without a preinstalled copy of `Shennong`.

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

- Retained `sn_load_pbmc()` as a deprecated compatibility wrapper around
  `sn_load_data()`.
- Refreshed generated documentation and package metadata as part of the
  modernization effort.

### Fixed

- Addressed multiple modernization issues in the package build, test, and
  documentation pipeline.
