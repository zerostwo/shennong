# Plot RNA velocity vectors

Plot RNA velocity vectors

## Usage

``` r
sn_plot_velocity(
  x,
  name = NULL,
  color_by = "pseudotime",
  arrow_scale = 1,
  point_size = 0.6
)
```

## Arguments

- x:

  A Seurat object or velocity result.

- name:

  Stored result name.

- color_by:

  Cell-table field used to color points.

- arrow_scale:

  Multiplicative arrow-length scale.

- point_size:

  Point size.

## Value

A `ggplot` object.
