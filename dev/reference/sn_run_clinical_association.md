# Associate bulk features with clinical variables

Numeric clinical variables use linear models; two-level categorical
variables use Welch tests and multi-level variables use ANOVA.

## Usage

``` r
sn_run_clinical_association(
  object,
  features,
  clinical_vars,
  covariates = NULL,
  metadata = NULL,
  assay = NULL,
  result_id = "bulk_clinical"
)
```

## Arguments

- object:

  Bulk input accepted by
  [`sn_assess_bulk_qc()`](https://zerostwo.github.io/shennong/dev/reference/sn_assess_bulk_qc.md).

- features:

  Expression features or numeric metadata columns.

- clinical_vars:

  Clinical metadata variables to test.

- covariates:

  Optional adjustment variables for linear models.

- metadata:

  Optional sample metadata.

- assay:

  Assay name for `SummarizedExperiment` input.

- result_id:

  Stable identifier for the returned clinical-association result.

## Value

A validated clinical-association result.
