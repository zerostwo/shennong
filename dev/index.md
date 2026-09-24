# Shennong

Shennong is an R package for single-cell, multimodal, spatial, and bulk
transcriptomics. It brings preprocessing, clustering, differential
analysis, pathway analysis, and plotting into a common workflow, and
keeps analytical results with their inputs and parameters on a Seurat
object.

**New here?** Follow [Get
started](https://zerostwo.github.io/shennong/dev/articles/get-started.html)
for a complete example using bundled data. No data download or Python
setup is needed. The package is in development, and documented API
changes may be breaking.

## Install

``` r

install.packages(c("remotes", "Seurat"))
remotes::install_github("zerostwo/shennong")
library(Shennong)
```

Install optional backends only for the workflows you need. Use
[`sn_list_methods()`](https://zerostwo.github.io/shennong/dev/reference/sn_list_methods.md)
and
[`sn_get_method_status()`](https://zerostwo.github.io/shennong/dev/reference/sn_get_method_status.md)
to inspect availability; [Choose a
backend](https://zerostwo.github.io/shennong/dev/articles/method-catalog.html)
lists methods and their execution requirements.

## A first analysis

This small example uses PBMC counts bundled with SeuratObject. Its
clusters illustrate the API and are not cell-type annotations or
condition-level tests.

``` r

data("pbmc_small", package = "SeuratObject")
pbmc <- SeuratObject::CreateSeuratObject(
  SeuratObject::GetAssayData(pbmc_small, assay = "RNA", layer = "counts")
)
pbmc <- sn_run_cluster(
  pbmc, integration_method = "unintegrated",
  nfeatures = 100, npcs = 10, dims = 1:10, resolution = 1.2, seed = 717,
  verbose = FALSE
)
sn_plot_dim(pbmc, reduction = "umap", group_by = "seurat_clusters")

pbmc <- sn_find_de(
  pbmc, analysis = "markers", group_by = "seurat_clusters",
  result_id = "markers", verbose = FALSE
)
sn_get_de_result(pbmc, "markers", direction = "up", top_n = 3)
```

`top_n = 3` selects up to three genes **per group**. Use
`top_scope = "all"` for three genes overall. Keep the complete result
for later use:

``` r

sn_list_results(pbmc, type = "de")
result <- sn_get_result(pbmc, type = "de", result_id = "markers")
head(result$tables$primary)
saveRDS(pbmc, "pbmc-analysed.rds")
```

Assign an object-returning workflow back to `pbmc` to retain its result.
Use `return_object = FALSE` on DE, enrichment, scoring, or Milo to
receive a unified result directly. See [Parameters and
results](https://zerostwo.github.io/shennong/dev/articles/parameters-and-results.html)
for return values, result IDs, filters, and safe reruns.

## Find the right guide

| I want to… | Guide |
|----|----|
| Learn the core workflow with runnable examples | [Get started](https://zerostwo.github.io/shennong/dev/articles/get-started.html) |
| Understand parameters, grouping, and stored results | [Parameters and results](https://zerostwo.github.io/shennong/dev/articles/parameters-and-results.html) |
| Import data, filter cells, and normalize | [Data IO](https://zerostwo.github.io/shennong/dev/articles/data-io-projects.html) · [Preprocessing](https://zerostwo.github.io/shennong/dev/articles/preprocessing-qc.html) |
| Cluster cells or integrate batches | [Clustering and integration](https://zerostwo.github.io/shennong/dev/articles/clustering.html) |
| Find markers, compare conditions, and enrich pathways | [Differential expression and enrichment](https://zerostwo.github.io/shennong/dev/articles/differential-expression.html) |
| Score signatures or annotate cells | [Annotation and pathways](https://zerostwo.github.io/shennong/dev/articles/annotation-pathways.html) |
| Compare composition or neighborhood abundance | [Composition and Milo](https://zerostwo.github.io/shennong/dev/articles/composition-analysis.html) |
| Plot expression and analytical results | [Visualization](https://zerostwo.github.io/shennong/dev/articles/visualization.html) |
| Work with spatial, bulk, or cell dynamics data | [All workflows](https://zerostwo.github.io/shennong/dev/articles/index.html) |
| Look up a specific function or backend | [Function reference](https://zerostwo.github.io/shennong/dev/reference/index.html) · [Backend catalog](https://zerostwo.github.io/shennong/dev/articles/method-catalog.html) |

The introductory examples run on bundled data. The larger real-data
workflows state their required files and optional dependencies; they do
not download study data during a normal website build. Use your own data
with the documented input structure, or follow the [research workflow
map](https://zerostwo.github.io/shennong/dev/articles/research-workflow-map.html)
for the public-data narratives. Dataset discovery and materialization
live in [ShennongData](https://github.com/zerostwo/shennong-data).

See the [release
notes](https://zerostwo.github.io/shennong/dev/news/index.html) for
breaking changes. Report reproducible problems in the [issue
tracker](https://github.com/zerostwo/shennong/issues).
