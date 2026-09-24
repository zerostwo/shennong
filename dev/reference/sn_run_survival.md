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
  result_id = "bulk_survival",
  group_method = c("median", "quantile", "fixed", "none"),
  group_quantile = 0.5,
  group_cutpoint = NULL,
  group_labels = c("Low", "High"),
  ties = c("efron", "breslow", "exact"),
  risk_times = NULL
)
```

## Arguments

- object:

  Bulk input accepted by
  [`sn_assess_bulk_qc()`](https://zerostwo.github.io/shennong/dev/reference/sn_assess_bulk_qc.md).

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

- result_id:

  Stable identifier for the returned survival result.

- group_method:

  Feature grouping used for Kaplan-Meier analysis.

- group_quantile:

  Quantile used when `group_method = "quantile"`.

- group_cutpoint:

  Fixed scalar or feature-named cutpoints.

- group_labels:

  Labels ordered as lower/equal and higher than the cutpoint.

- ties:

  Cox partial-likelihood tie method.

- risk_times:

  Optional non-negative times shown in the risk table. The default uses
  at most eight deterministic pretty breaks per feature.

## Value

A validated survival result with one adjusted Cox model per feature.
