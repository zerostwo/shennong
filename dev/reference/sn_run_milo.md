# Run neighborhood differential abundance testing with miloR

This wrapper runs the standard miloR neighborhood differential abundance
workflow on a Seurat object using one embedding and one sample-level
group contrast. Cells are grouped into neighborhoods, counted per
sample, and then tested for differential abundance between the two
sample groups.

## Usage

``` r
sn_run_milo(
  x,
  sample_col,
  group_col,
  contrast = NULL,
  reduction = "pca",
  dims = NULL,
  cells = NULL,
  max_cells = NULL,
  stratify_by = sample_col,
  k = 20,
  d = NULL,
  prop = 0.1,
  refined = TRUE,
  refinement_scheme = "reduced_dim",
  covariates = NULL,
  annotation_col = NULL,
  fdr_weighting = c("k-distance", "neighbour-distance", "max", "graph-overlap", "none"),
  min_mean = 0,
  norm_method = c("TMM", "RLE", "logMS"),
  store_name = NULL,
  return_object = FALSE,
  return_intermediate = FALSE,
  verbose = TRUE
)
```

## Arguments

- x:

  A Seurat object.

- sample_col:

  Metadata column defining biological samples.

- group_col:

  Metadata column defining the sample-level comparison group.

- contrast:

  Optional character vector of length 2 giving the comparison as
  `c(case, control)`. When omitted, `group_col` must contain exactly two
  levels in the selected cells.

- reduction:

  Reduction name used to build neighborhoods. Defaults to `"pca"`.

- dims:

  Optional integer vector of embedding dimensions to retain.

- cells:

  Optional character vector of cells to include.

- max_cells:

  Optional integer cap used to subsample cells before neighborhood
  construction.

- stratify_by:

  Optional metadata column used to preserve representation during
  subsampling. Defaults to `sample_col`.

- k:

  Number of neighbors for graph and neighborhood construction.

- d:

  Number of embedding dimensions passed to miloR. Defaults to the
  selected embedding dimensionality.

- prop:

  Proportion of cells sampled as neighborhood indices.

- refined:

  Logical; if `TRUE`, use miloR's refined neighborhood sampling.

- refinement_scheme:

  Refinement scheme passed to `miloR::makeNhoods()`.

- covariates:

  Optional sample-level covariates added to the DA design formula
  alongside `group_col`.

- annotation_col:

  Optional cell-level metadata column used to annotate neighborhoods
  with `miloR::annotateNhoods()`.

- fdr_weighting:

  FDR weighting strategy passed to `miloR::testNhoods()`.

- min_mean:

  Minimum mean count threshold passed to `miloR::testNhoods()`.

- norm_method:

  Normalization method passed to `miloR::testNhoods()`.

- store_name:

  Optional name used under `object@misc$milo_results`. When supplied,
  the milo result is stored on the Seurat object.

- return_object:

  Logical; when `TRUE` and `store_name` is supplied, return the updated
  Seurat object.

- return_intermediate:

  Logical; if `TRUE`, return a list with the DA table, design data, and
  milo object.

- verbose:

  Logical; if `TRUE`, emit progress logs.

## Value

By default, a data frame of neighborhood-level DA statistics. When
`return_intermediate = TRUE`, a list with `table`, `design_df`, and
`milo` is returned.

## Examples

``` r
if (FALSE) { # \dontrun{
da_tbl <- sn_run_milo(
  seu,
  sample_col = "sample",
  group_col = "condition",
  contrast = c("treated", "control"),
  reduction = "pca"
)
} # }
```
