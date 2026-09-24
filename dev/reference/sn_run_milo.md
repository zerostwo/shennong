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
  sample_by = NULL,
  group_by = NULL,
  contrast = NULL,
  reduction = "pca",
  dims = NULL,
  cells = NULL,
  max_cells = NULL,
  stratify_by = NULL,
  k = 20,
  d = NULL,
  prop = 0.1,
  refined = TRUE,
  refinement_scheme = "reduced_dim",
  covariates = NULL,
  annotation_by = NULL,
  fdr_weighting = c("k-distance", "neighbour-distance", "max", "graph-overlap", "none"),
  min_mean = 0,
  norm_method = c("TMM", "RLE", "logMS"),
  result_id = NULL,
  return_object = TRUE,
  keep_model = FALSE,
  verbose = TRUE,
  seed = 717,
  overwrite = FALSE,
  object = NULL
)
```

## Arguments

- x:

  A Seurat object.

- sample_by:

  Metadata column defining biological samples.

- group_by:

  Metadata column defining the sample-level comparison group.

- contrast:

  Optional character vector of length 2 giving the comparison as
  `c(case, control)`. When omitted, `group_by` must contain exactly two
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
  subsampling. Defaults to `sample_by`.

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
  alongside `group_by`.

- annotation_by:

  Optional cell-level metadata column used to annotate neighborhoods
  with `miloR::annotateNhoods()`.

- fdr_weighting:

  FDR weighting strategy passed to `miloR::testNhoods()`.

- min_mean:

  Minimum mean count threshold passed to `miloR::testNhoods()`.

- norm_method:

  Normalization method passed to `miloR::testNhoods()`.

- result_id:

  Optional stable identifier. A unique `milo` ID is allocated when
  omitted.

- return_object:

  Return the updated Seurat object (default), or the unified analysis
  result when `FALSE`.

- keep_model:

  Retain the fitted Milo object in `models$milo`.

- verbose:

  Logical; if `TRUE`, emit progress logs.

- seed:

  Random seed used for cell and neighborhood sampling, applied locally
  without changing the caller's RNG state. `NULL` uses the caller's
  state; an effective numeric seed is recorded in provenance.

- overwrite:

  Explicitly replace an existing result; requires an explicit ID.

- object:

  Alias for `x`; supply only one of `x` and `object`.

## Value

A Seurat object or unified result. Retrieve neighborhood statistics with
[`sn_get_milo_result()`](https://zerostwo.github.io/shennong/dev/reference/sn_get_milo_result.md);
design data are in `tables$design`.

## Examples

``` r
if (FALSE) { # \dontrun{
seu <- sn_run_milo(
  seu,
  sample_by = "sample",
  group_by = "condition",
  contrast = c("treated", "control"),
  reduction = "pca"
)
} # }
```
