# Rank metadata variables by embedding variance explained

Quantifies how much variation in a dimensional reduction is explained by
one or more metadata variables. This is useful for identifying whether
`platform`, `study`, `tissue`, `sample`, or another covariate is the
dominant driver of residual batch structure.

## Usage

``` r
sn_calculate_variance_explained(
  x,
  variables,
  reduction = .sn_default_metric_reduction(x),
  dims = NULL,
  method = c("single", "partial"),
  cells = NULL,
  max_cells = NULL,
  stratify_by = NULL,
  seed = 717,
  return_dim_data = FALSE
)
```

## Arguments

- x:

  A Seurat object.

- variables:

  Character vector of metadata columns to rank.

- reduction:

  Reduction name used for the score. Defaults to `"harmony"` when
  present, otherwise `"pca"`.

- dims:

  Optional integer vector of embedding dimensions to retain.

- method:

  One of `"single"` or `"partial"`. Defaults to `"single"`.

- cells:

  Optional character vector of cell names to include.

- max_cells:

  Optional integer cap used to subsample cells before running the
  metric.

- stratify_by:

  Optional metadata column used to preserve representation during
  subsampling.

- seed:

  Random seed used when `max_cells` triggers subsampling.

- return_dim_data:

  Logical; if `TRUE`, return a list containing the summary table and
  per-dimension results.

## Value

A data frame ranked by `variance_explained`. When
`return_dim_data = TRUE`, a list with `summary` and `dim_data` is
returned.

## Details

Two modes are available. `method = "single"` fits one model per variable
and reports the weighted R-squared across the selected dimensions.
`method = "partial"` fits all variables together and reports the
incremental variance explained by each variable after the others.
Partial estimates are helpful but can be ambiguous when variables are
nested or confounded, such as one platform per study.

## Examples

``` r
if (FALSE) { # \dontrun{
variance_tbl <- sn_calculate_variance_explained(
  seu,
  variables = c("platform", "study", "tissue", "sample"),
  reduction = "pca",
  dims = 1:30
)
variance_tbl
} # }
```
