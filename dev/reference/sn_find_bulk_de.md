# Find differential expression in bulk transcriptomics data

Find differential expression in bulk transcriptomics data

## Usage

``` r
sn_find_bulk_de(
  object,
  metadata = NULL,
  design = ~condition,
  contrast,
  method = c("auto", "edger", "deseq2", "limma", "dream"),
  assay = NULL,
  store_name = "bulk_de",
  backend_control = list()
)
```

## Arguments

- object:

  Bulk input accepted by
  [`sn_assess_bulk_qc()`](https://songqi.org/shennong/dev/reference/sn_assess_bulk_qc.md).

- metadata:

  Optional sample metadata.

- design:

  A fixed- or mixed-effects model formula.

- contrast:

  Character triple: variable, numerator, denominator.

- method:

  One of `auto`, `edger`, `deseq2`, `limma`, or `dream`.

- assay:

  Assay name for `SummarizedExperiment` input.

- store_name:

  Result name.

- backend_control:

  Backend options or a custom `runner`/precomputed `result`.

## Value

A validated Shennong bulk-DE result.
