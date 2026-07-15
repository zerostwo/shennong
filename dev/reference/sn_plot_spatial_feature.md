# Plot expression in spatial coordinates

Plot expression in spatial coordinates

## Usage

``` r
sn_plot_spatial_feature(
  object,
  features,
  spatial_cols = NULL,
  assay = NULL,
  layer = "data",
  point_size = 1.5
)
```

## Arguments

- object:

  A Seurat object.

- features:

  Features to plot.

- spatial_cols:

  Coordinate metadata columns.

- assay, layer:

  Expression assay and layer.

- point_size:

  Point size.

## Value

A faceted `ggplot` object.
