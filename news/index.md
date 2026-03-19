# Changelog

## Version 0.1.1

Released 2026-03-19.

#### Added

- Automatic human/mouse species inference through
  [`sn_get_species()`](https://songqi.org/shennong/reference/sn_get_species.md)
  using `hom_genes` and mitochondrial naming patterns.
- A layer-aware pkgdown article covering inferred species, non-default
  count layers, and stored differential-expression metadata.
- Additional regression tests for species inference, DE metadata, and
  BPCells-backed Seurat layers.

#### Changed

- [`sn_initialize_seurat_object()`](https://songqi.org/shennong/reference/sn_initialize_seurat_object.md)
  now attempts species inference before deciding whether to compute
  species-specific QC metrics.
- [`sn_get_signatures()`](https://songqi.org/shennong/reference/sn_get_signatures.md)
  now uses package-local fallback signatures for the core categories
  needed by Shennong workflows when `SignatuR` is unavailable.
- Stored DE results now include schema and provenance metadata such as
  package version, timestamp, assay/layer context, and threshold
  settings.
- `sn_plot_*()` helpers now treat `catplot` as an optional enhancement
  rather than a mandatory dependency.

#### Fixed

- Fixed
  [`sn_remove_ambient_contamination()`](https://songqi.org/shennong/reference/sn_remove_ambient_contamination.md)
  for BPCells-backed Seurat layers by materializing BPCells matrices
  before passing them to `decontX`.
- Fixed optional dependency handling so missing `SignatuR` no longer
  breaks core initialization and clustering paths.
- Fixed SoupX validation order so missing `raw` input is reported before
  package availability issues.

## Version 0.1.0

Released 2026-03-18.

#### Added

- [`sn_load_data()`](https://songqi.org/shennong/reference/sn_load_data.md)
  as the primary example-data loader.
- Consolidated clustering into
  [`sn_run_cluster()`](https://songqi.org/shennong/reference/sn_run_cluster.md).
- A unified ambient contamination interface for SoupX and decontX.
- Expanded test coverage across composition, data loading, clustering,
  utilities, visualization, and ambient contamination.
- Missing help pages for exported functions, CI scaffolding, a richer
  README, and updated package metadata for current R syntax
  requirements.

#### Changed

- Retained
  [`sn_load_pbmc()`](https://songqi.org/shennong/reference/sn_load_pbmc.md)
  as a deprecated compatibility wrapper around
  [`sn_load_data()`](https://songqi.org/shennong/reference/sn_load_data.md).
- Refreshed generated documentation and package metadata as part of the
  modernization effort.

#### Fixed

- Addressed multiple modernization issues in the package build, test,
  and documentation pipeline.
