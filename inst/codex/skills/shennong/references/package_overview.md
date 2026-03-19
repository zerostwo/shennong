# Package Overview

## What Shennong Is For

Shennong is an R package for single-cell analysis workflows built around
Seurat objects. It provides a coherent end-user API for:

- loading packaged PBMC example datasets
- initializing Seurat objects
- quality control and ambient RNA handling
- normalization, clustering, and Harmony integration
- annotation and marker discovery
- composition summaries and integration metrics
- visualization, enrichment, and GSEA

## Intended User

This skill is aimed at analysts who already have `Shennong` installed in an R
environment and want a stable, package-level interface instead of stitching
together many raw Seurat and Bioconductor calls.

## API Shape

The public API uses strict `sn_verb_noun` naming. Common families are:

- `sn_load_*`: packaged example data and installation helpers
- `sn_initialize_*`, `sn_normalize_*`, `sn_filter_*`: preprocessing and QC
- `sn_run_*`: clustering, annotation, and integration workflows
- `sn_find_*`: differential expression and doublet detection
- `sn_calculate_*`: composition and metric summaries
- `sn_plot_*`: plotting wrappers for common single-cell outputs
- `sn_enrich()`: enrichment and GSEA

## Package Boundaries

- Shennong’s main in-memory object is a Seurat object.
- Count-like matrix workflows increasingly expose explicit `assay` and `layer`
  parameters.
- Differential expression results can be stored back into the Seurat object in
  `object@misc$de_results`.

## Typical Session

1. `library(Shennong)`
2. `object <- sn_load_data("pbmc3k")`
3. `object <- sn_filter_cells(object, ...)`
4. `object <- sn_filter_genes(object, ...)`
5. `object <- sn_run_cluster(object, ...)`
6. `object <- sn_find_de(object, analysis = "markers", ...)`
7. `sn_plot_dim(object, ...)`
8. `sn_plot_dot(object, features = "top_markers", ...)`
9. `sn_enrich(...)`
