# Calculate ROGUE score for Seurat Object

This function calculates ROGUE score based on Seurat object.

## Usage

``` r
sn_calculate_rogue(
  x,
  cluster = NULL,
  sample = NULL,
  span = 0.9,
  assay = "RNA",
  layer = "counts"
)
```

## Arguments

- x:

  A Seurat object.

- cluster:

  Column name in metadata specifying cluster labels.

- sample:

  Column name in metadata specifying sample labels.

- span:

  The span parameter for rogue estimation.

- assay:

  Assay used for ROGUE calculation. Defaults to `"RNA"`.

- layer:

  Layer used as the input count matrix. Defaults to `"counts"`.

## Value

A data.frame of ROGUE score per cluster/sample or per cluster.

## Examples

``` r
if (FALSE) { # \dontrun{
data("pbmc_small", package = "Shennong")
pbmc <- sn_run_cluster(pbmc_small, normalization_method = "seurat", verbose = FALSE)
rogue_tbl <- sn_calculate_rogue(pbmc, cluster = "seurat_clusters")
head(rogue_tbl)
} # }
```
