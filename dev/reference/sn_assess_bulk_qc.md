# Assess bulk transcriptomics sample quality

Computes library size, detected features, expression distributions,
sample PCA, sample correlations, and robust multivariate outlier flags.
The analysis-scale expression of the selected variable features is
retained in `tables$expression` for auditable sample-to-sample scatter
plots.

## Usage

``` r
sn_assess_bulk_qc(
  object,
  metadata = NULL,
  assay = NULL,
  top_variable = 2000L,
  outlier_z = 3.5,
  result_id = "bulk_qc"
)
```

## Arguments

- object:

  A feature-by-sample matrix, a `SummarizedExperiment`, or a list
  containing `counts`/`expression` and optional `metadata`.

- metadata:

  Optional sample metadata with rows matching sample names.

- assay:

  Assay name for `SummarizedExperiment` input.

- top_variable:

  Number of variable features used for PCA/correlation.

- outlier_z:

  Robust z-score threshold used for sample flags.

- result_id:

  Stable identifier for the returned bulk-QC result.

## Value

A validated Shennong bulk-QC result.
