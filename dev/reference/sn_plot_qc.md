# Plot QC assessment summaries

Plot QC assessment summaries

## Usage

``` r
sn_plot_qc(x, metric = c("qc_score", "n_cells", "retention_fraction"))
```

## Arguments

- x:

  A result from
  [`sn_assess_qc()`](https://songqi.org/shennong/dev/reference/sn_assess_qc.md)
  or its `by_sample` table.

- metric:

  Metric shown on the y-axis.

## Value

A `ggplot` object with figure metadata.
