# Plot RNA velocity vectors

Plot RNA velocity vectors

## Usage

``` r
sn_plot_velocity(
  x,
  result_id = NULL,
  color_by = "pseudotime",
  arrow_scale = 1,
  point_size = 0.6,
  object = NULL
)
```

## Arguments

- x:

  A Seurat object or velocity result.

- result_id:

  Stored result result_id.

- color_by:

  Cell-table field used to color points.

- arrow_scale:

  Multiplicative arrow-length scale.

- point_size:

  Point size.

- object:

  Alias for `x`; supply only one of `x` and `object`.

## Value

A `ggplot` object.
