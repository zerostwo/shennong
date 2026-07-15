# Plot spatial coordinates colored by metadata

Plot spatial coordinates colored by metadata

## Usage

``` r
sn_plot_spatial(object, group_by = NULL, spatial_cols = NULL, point_size = 1.5)
```

## Arguments

- object:

  A Seurat object.

- group_by:

  Optional metadata column.

- spatial_cols:

  Coordinate metadata columns.

- point_size:

  Point size.

## Value

A `ggplot` object preserving coordinate aspect ratio.
