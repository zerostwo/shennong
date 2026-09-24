# Plot QC assessment summaries

Plot QC assessment summaries

## Usage

``` r
sn_plot_qc(
  x,
  metric = c("qc_score", "n_cells", "retention_fraction"),
  object = NULL
)
```

## Arguments

- x:

  A result from
  [`sn_assess_qc()`](https://zerostwo.github.io/shennong/dev/reference/sn_assess_qc.md)
  or its `by_sample` table.

- metric:

  Metric shown on the y-axis.

- object:

  Alias for `x`; supply only one of `x` and `object`.

## Value

A `ggplot` object with figure metadata.
