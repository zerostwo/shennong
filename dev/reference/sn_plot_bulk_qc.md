# Plot bulk sample quality metrics

Plot bulk sample quality metrics

## Usage

``` r
sn_plot_bulk_qc(
  x,
  metric = c("library_size", "detected_features", "mean_correlation")
)
```

## Arguments

- x:

  A bulk-QC result.

- metric:

  Sample metric to plot.

## Value

A `ggplot` object.
