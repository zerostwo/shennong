# Plot spatial deconvolution proportions or abundances

Values are displayed as supplied, without normalization or a forced
zero-to-one range, so abundance estimates remain distinct from
proportions.

## Usage

``` r
sn_plot_spatial_deconvolution(
  x,
  spatial_cols = c("spatial_x", "spatial_y"),
  location_col = "cell",
  cell_type_col = "cell_type",
  proportion_col = "proportion",
  split_by = NULL,
  point_size = 1.5,
  pt_size = NULL,
  cols = NULL,
  palette = "viridis",
  alpha = 1,
  na_value = "grey80",
  limits = NULL,
  order = TRUE,
  show_legend = TRUE,
  show_axis = FALSE,
  show_border = TRUE,
  title = NULL,
  ncol = NULL,
  flip_y = TRUE,
  raster = FALSE,
  raster_dpi = 512,
  style = "classic"
)
```

## Arguments

- x:

  A data frame with location, cell type, and value columns.

- spatial_cols:

  Coordinate column names in x, y order.

- location_col, cell_type_col, proportion_col:

  Column names.

- split_by:

  Optional table column identifying sections.

- point_size:

  Point size in mm; retained for positional compatibility.

- pt_size:

  Alias matching
  [`sn_plot_dim()`](https://zerostwo.github.io/shennong/dev/reference/sn_plot_dim.md).
  Supply only one size argument.

- cols:

  Explicit colors; named vectors match categorical labels. For feature
  maps these are gradient colors.

- palette:

  Shennong palette name or color vector.

- alpha:

  Point opacity between zero and one.

- na_value:

  Color for missing group or feature values.

- limits:

  Optional common color limits. Out-of-range colors are clamped;
  original values remain in the plot data and exported source table.

- order:

  Draw higher values last and missing values first.

- show_legend, show_axis, show_border:

  Show the legend, coordinate axes or panel border.

- title:

  Optional plot title.

- ncol:

  Number of facet columns.

- flip_y:

  Reverse the y axis for image coordinates. Set FALSE for Cartesian
  coordinates.

- raster:

  Rasterize only point layers with ggrastr; labels stay vector.

- raster_dpi:

  Raster resolution (scalar or vector, whose maximum is used).

- style:

  Presentation preset; see
  [`sn_list_plot_styles()`](https://zerostwo.github.io/shennong/dev/reference/sn_list_plot_styles.md).
  Coordinates, shared limits and source values are unchanged.

## Value

A faceted `ggplot` with a common continuous scale.
