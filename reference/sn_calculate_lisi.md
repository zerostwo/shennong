# Calculate LISI scores from a Seurat embedding

This function calculates the Local Inverse Simpson's Index (LISI) for
one or more metadata labels from a Seurat reduction. It is commonly used
to assess batch mixing or label separation after integration.

## Usage

``` r
sn_calculate_lisi(
  x,
  reduction = "pca",
  label = "sample",
  dims = NULL,
  cells = NULL,
  max_cells = NULL,
  stratify_by = label[[1]],
  seed = 717
)
```

## Arguments

- x:

  A Seurat object.

- reduction:

  Reduction name used to extract embeddings. Defaults to `"pca"`.

- label:

  Character vector of metadata column names passed to
  [`lisi::compute_lisi()`](https://rdrr.io/pkg/lisi/man/compute_lisi.html).

- dims:

  Optional integer vector of embedding dimensions to retain.

- cells:

  Optional character vector of cell names to score.

- max_cells:

  Optional integer cap used to subsample cells before running LISI. When
  `NULL`, use all selected cells.

- stratify_by:

  Optional metadata column used to preserve representation during
  subsampling. Defaults to the first requested `label`.

- seed:

  Random seed used when `max_cells` triggers subsampling.

## Value

A data frame with one row per retained cell. The first column is
`cell_id`; each requested label contributes one LISI score column.

## Examples

``` r
if (FALSE) { # \dontrun{
data("pbmc_small", package = "Shennong")
pbmc <- sn_run_cluster(
  pbmc_small,
  batch = "sample",
  species = "human",
  verbose = FALSE
)
lisi_tbl <- sn_calculate_lisi(
  pbmc,
  reduction = "harmony",
  label = "sample"
)
head(lisi_tbl)
} # }
```
