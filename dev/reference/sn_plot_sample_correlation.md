# Plot bulk sample correlation

Plot bulk sample correlation

## Usage

``` r
sn_plot_sample_correlation(
  x,
  view = c("heatmap", "scatter"),
  sample_x = NULL,
  sample_y = NULL,
  features = NULL,
  method = c("spearman", "pearson", "kendall"),
  transform = c("none", "log1p"),
  add_fit = TRUE,
  label = FALSE,
  object = NULL
)
```

## Arguments

- x:

  A bulk-QC result.

- view:

  Correlation heatmap or a feature-level scatter plot comparing two
  samples.

- sample_x, sample_y:

  Sample names used by `view = "scatter"`. When both are omitted, the
  first two samples are used.

- features:

  Optional feature subset for the scatter view.

- method:

  Correlation method shown by the scatter view.

- transform:

  Optional axis transformation for the scatter view.

- add_fit:

  Add a least-squares trend line to the scatter view.

- label:

  Label features in the scatter view.

- object:

  Alias for `x`; supply only one of `x` and `object`.

## Value

A `ggplot` correlation heatmap or sample-to-sample scatter plot.
