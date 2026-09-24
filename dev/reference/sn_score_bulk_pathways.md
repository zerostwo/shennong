# Score pathways in bulk expression samples

Score pathways in bulk expression samples

## Usage

``` r
sn_score_bulk_pathways(
  object,
  signatures,
  method = c("mean", "gsva", "ssgsea"),
  metadata = NULL,
  assay = NULL,
  min_genes = 2L,
  result_id = "bulk_pathways",
  backend_control = list()
)
```

## Arguments

- object:

  Bulk input accepted by
  [`sn_assess_bulk_qc()`](https://zerostwo.github.io/shennong/dev/reference/sn_assess_bulk_qc.md).

- signatures:

  Named list or two-column data frame of gene sets.

- method:

  Scoring method: mean, GSVA, or ssGSEA.

- metadata:

  Optional sample metadata.

- assay:

  Assay name for `SummarizedExperiment` input.

- min_genes:

  Minimum matched genes per pathway.

- result_id:

  Stable identifier for the returned pathway-score result.

- backend_control:

  Backend-specific controls.

## Value

A validated bulk pathway-score result.
