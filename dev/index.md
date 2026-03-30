# Shennong

`Shennong` is an experimental R package for single-cell transcriptomics
workflows built around Seurat objects. It focuses on practical analysis
steps that are commonly repeated across projects: preprocessing,
clustering, batch integration, differential expression, enrichment,
visualization, and interpretation-ready result storage.

## Installation

Install the current development version from GitHub:

``` r
install.packages("remotes")
remotes::install_github("zerostwo/shennong")
```

List Shennong’s required and recommended R packages, then install the
missing ones in one step:

``` r
deps <- sn_list_dependencies()
deps

sn_install_dependencies(scope = "required")
```

## What Shennong Covers

- Seurat object initialization and QC-aware preprocessing
- Gene and cell filtering, including bundled human/mouse GENCODE gene
  classes
- Doublet detection and ambient RNA correction
- Single-dataset clustering and Harmony-based integration with
  [`sn_run_cluster()`](https://songqi.org/shennong/dev/reference/sn_run_cluster.md)
- Multi-metric integration assessment with
  [`sn_assess_integration()`](https://songqi.org/shennong/dev/reference/sn_assess_integration.md)
- Marker detection, pseudobulk differential expression, and enrichment
  analysis
- Bulk deconvolution workflows with BayesPrism and local CIBERSORTx
  containers
- Stored-result workflows for downstream interpretation and reporting

## Built-In Example Data

The package ships a small built-in PBMC example derived from the
`pbmc1k` and `pbmc3k` assets:

- `pbmc_small`: a Seurat object with sample metadata
- `pbmc_small_raw`: a matching raw count matrix with extra droplets

For larger example datasets,
[`sn_load_data()`](https://songqi.org/shennong/dev/reference/sn_load_data.md)
can still download the full PBMC references on demand.

## Quick Start

``` r
library(Shennong)

data("pbmc_small", package = "Shennong")

pbmc <- sn_filter_cells(
  pbmc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  plot = FALSE
)
pbmc <- sn_filter_genes(pbmc, min_cells = 3, plot = FALSE)
pbmc <- sn_run_cluster(
  pbmc,
  normalization_method = "seurat",
  resolution = 0.6
)

sn_plot_dim(pbmc, group_by = "seurat_clusters", label = TRUE)
```

## Integration Example

The same built-in PBMC example can be used for batch-aware integration:

``` r
pbmc_integrated <- sn_run_cluster(
  pbmc,
  batch = "sample",
  normalization_method = "seurat",
  resolution = 0.6
)

sn_plot_dim(pbmc_integrated, group_by = "sample")
sn_plot_dim(pbmc_integrated, group_by = "seurat_clusters", label = TRUE)
```

You can summarize integration quality and surface rare or difficult
groups directly from the integrated object:

``` r
metrics <- sn_assess_integration(
  pbmc_integrated,
  batch = "sample",
  cluster = "seurat_clusters",
  reduction = "harmony",
  baseline_reduction = "pca"
)

metrics$summary
metrics$per_group$isolated_label_score
metrics$per_group$cluster_entropy
metrics$per_group$cluster_purity
metrics$per_group$challenging_groups
```

## Differential Expression And Enrichment

``` r
pbmc <- sn_find_de(
  pbmc,
  analysis = "markers",
  group_by = "seurat_clusters",
  layer = "data",
  store_name = "cluster_markers",
  return_object = TRUE,
  verbose = FALSE
)

pbmc <- sn_enrich(
  x = pbmc,
  source_de_name = "cluster_markers",
  gene_clusters = gene ~ cluster,
  database = c("GOBP", "H"),
  species = "human",
  store_name = "cluster_pathways",
  pvalue_cutoff = 0.05
)
```

## Documentation

Longer workflow articles are available in the package site and
vignettes, including:

- clustering and integration
- layer-aware preprocessing
- stored-result interpretation workflows

## Status

`Shennong` is still experimental. The package currently prioritizes a
clean and consistent workflow surface over backward compatibility across
early versions.
