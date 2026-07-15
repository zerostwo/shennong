# Run Cox proportional-hazards models for bulk features

Run Cox proportional-hazards models for bulk features

## Usage

``` r
sn_run_survival(
  object,
  time,
  event,
  features,
  covariates = NULL,
  metadata = NULL,
  assay = NULL,
  store_name = "bulk_survival"
)
```

## Arguments

- object:

  Bulk input accepted by
  [`sn_assess_bulk_qc()`](https://songqi.org/shennong/dev/reference/sn_assess_bulk_qc.md).

- time, event:

  Metadata columns containing follow-up time and event status.

- features:

  Expression features or numeric metadata columns.

- covariates:

  Optional adjustment variables.

- metadata:

  Optional sample metadata.

- assay:

  Assay name for `SummarizedExperiment` input.

- store_name:

  Result name.

## Value

A validated survival result with one adjusted Cox model per feature.
