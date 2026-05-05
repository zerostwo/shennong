# Calculate ROGUE score for Seurat Object

This function calculates ROGUE score based on Seurat object.

## Usage

``` r
sn_calculate_rogue(
  x,
  cluster_by = NULL,
  sample_by = NULL,
  span = 0.9,
  assay = "RNA",
  layer = "counts",
  cells = NULL,
  max_cells = 3000,
  stratify_by = NULL,
  seed = 717,
  min_cells = 10,
  min_genes = 10,
  cluster = NULL,
  sample = NULL
)
```

## Arguments

- x:

  A Seurat object.

- cluster_by:

  Column name in metadata specifying cluster_by labels.

- sample_by:

  Column name in metadata specifying sample labels.

- span:

  The span parameter for rogue estimation.

- assay:

  Assay used for ROGUE calculation. Defaults to `"RNA"`.

- layer:

  Layer used as the input count matrix. Defaults to `"counts"`.

- cells:

  Optional character vector of cell names to include.

- max_cells:

  Optional integer cap used to subsample cells before running ROGUE.
  Defaults to `3000`.

- stratify_by:

  Optional metadata column used to preserve representation during
  subsampling. Defaults to `sample` when supplied, otherwise `cluster`.

- seed:

  Random seed used when `max_cells` triggers subsampling.

- min_cells:

  Minimum cells retained by the upstream `ROGUE::matr.filter()` step.

- min_genes:

  Minimum detected genes retained by the upstream `ROGUE::matr.filter()`
  step.

- cluster:

  Deprecated alias for `cluster_by`.

- sample:

  Deprecated alias for `sample_by`.

## Value

When neither `cluster` nor `sample` is supplied, returns a single
numeric ROGUE score for the selected matrix. When `cluster` is supplied,
returns a data frame with per-cluster ROGUE scores. When both `cluster`
and `sample` are supplied, returns a tidy data frame with one row per
sample-cluster combination.

## Examples

``` r
if (FALSE) { # \dontrun{
data("pbmc_small", package = "Shennong")
pbmc <- sn_run_cluster(pbmc_small, normalization_method = "seurat", verbose = FALSE)
rogue_tbl <- sn_calculate_rogue(pbmc, cluster_by = "seurat_clusters")
head(rogue_tbl)
} # }
```
