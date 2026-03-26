All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
- `sn_initialize_codex_project()` now scaffolds user analysis repositories from
  the shipped `inst/codex/project-template/` assets, creating a governed
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
- The pkgdown deployment workflow now installs the current package before
  publishing the site, so reference generation can load `Shennong` inside the
  temporary deployment worktree.

### Fixed

- Fixed pkgdown GitHub Actions deployment when `deploy_to_branch()` runs in a
  clean runner without a preinstalled copy of `Shennong`.
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
- Fixed GitHub Actions dependency installation for GitHub-hosted optional rare
  cell and deconvolution backends used during checks and coverage runs.

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
