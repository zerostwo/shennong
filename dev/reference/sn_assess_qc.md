# Assess overall QC status and before/after filtering outcomes

This function summarizes the QC status of a Seurat object, optionally
compares it against a reference object captured before filtering, and
produces per-sample and overall quality scores. It is designed to answer
whether low-quality cells and called doublets were removed while
retaining clean cells, and whether ambient-RNA correction introduced
problematic zero-count cells.

## Usage

``` r
sn_assess_qc(
  object,
  reference = NULL,
  sample_by = NULL,
  store_name = "default",
  return_object = FALSE,
  verbose = TRUE,
  sample_col = NULL
)
```

## Arguments

- object:

  A Seurat object to assess.

- reference:

  Optional Seurat object representing the pre-filter state. When
  supplied, the function compares the retained cells in `object` against
  the QC flags and doublet calls recorded in `reference`.

- sample_by:

  Optional metadata column defining samples. When `NULL`, the function
  uses `sample` or `orig.ident` when available and otherwise treats the
  object as one sample.

- store_name:

  Name used when storing the assessment under
  `object@misc$qc_assessments`.

- return_object:

  Logical; when `TRUE`, store the assessment in the Seurat object and
  return the updated object.

- verbose:

  Logical; when `TRUE`, print a concise QC summary.

- sample_col:

  Deprecated alias for `sample_by`.

## Value

A list with `overall`, `by_sample`, `comparison`, and `messages` when
`return_object = FALSE`; otherwise the updated Seurat object with the
stored assessment.

## Examples

``` r
if (requireNamespace("Seurat", quietly = TRUE)) {
  data("pbmc_small", package = "Shennong")
  qc_report <- sn_assess_qc(pbmc_small, verbose = FALSE)
  qc_report$overall
}
#>   n_samples n_cells qc_score qc_label retention_fraction
#> 1         2     200 78.06004     good                 NA
#>   low_quality_removed_fraction doublet_removed_fraction
#> 1                           NA                       NA
```
