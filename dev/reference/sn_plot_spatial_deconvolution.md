# Plot spatial deconvolution proportions

Plot spatial deconvolution proportions

## Usage

``` r
sn_plot_spatial_deconvolution(
  x,
  spatial_cols = c("spatial_x", "spatial_y"),
  location_col = "cell",
  cell_type_col = "cell_type",
  proportion_col = "proportion"
)
```

## Arguments

- x:

  A data frame with location, cell type, and proportion columns.

- spatial_cols:

  Coordinate column names.

- location_col, cell_type_col, proportion_col:

  Column names.

## Value

A faceted `ggplot` object.
