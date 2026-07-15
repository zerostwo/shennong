# Run a bulk transcriptomics workflow

Run a bulk transcriptomics workflow

## Usage

``` r
sn_run_bulk(
  object,
  workflow = c("qc", "de", "pathway", "network", "survival"),
  ...
)
```

## Arguments

- object:

  Bulk input accepted by
  [`sn_assess_bulk_qc()`](https://songqi.org/shennong/dev/reference/sn_assess_bulk_qc.md).

- workflow:

  Workflow to run.

- ...:

  Arguments passed to the selected workflow.

## Value

A validated Shennong analysis result.
