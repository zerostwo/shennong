# Plot spatial coordinates colored by metadata

Spatial maps use the same palette, compact typography and label
conventions as
[`sn_plot_dim()`](https://zerostwo.github.io/shennong/dev/reference/sn_plot_dim.md),
while retaining equal coordinate units on both axes.

## Usage

``` r
sn_plot_spatial(
  object,
  group_by = NULL,
  spatial_cols = NULL,
  point_size = 1.5,
  cells = NULL,
  split_by = NULL,
  images = NULL,
  cols = NULL,
  palette = "Paired",
  pt_size = NULL,
  alpha = 1,
  na_value = "grey80",
  label = FALSE,
  label_size = 8 * 0.36,
  label_color = "black",
  label_box = FALSE,
  label_halo = TRUE,
  repel = FALSE,
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

- object:

  A Seurat object.

- group_by:

  Optional categorical metadata column. NULL shows all locations.

- spatial_cols:

  Two numeric metadata columns in x, y order. Recognized metadata
  coordinates take precedence over native images unless `images` is
  supplied. Coordinates are never rescaled or registered across
  sections.

- point_size:

  Point size in mm; retained for positional compatibility.

- cells:

  Optional cell names to plot, in the requested order.

- split_by:

  Metadata column identifying sections or panels. Required to separate
  sections when using metadata coordinates; no sample column is guessed.

- images:

  Optional Seurat image/FOV names. Native images are used when metadata
  coordinates are absent. Each image is plotted separately using
  full-resolution Visium coordinates or FOV centroids; histology and
  cell segmentation outlines are not drawn. Overlapping image cell sets
  are rejected.

- cols:

  Explicit colors; named vectors match categorical labels. For feature
  maps these are gradient colors.

- palette:

  Shennong palette name or color vector.

- pt_size:

  Alias matching
  [`sn_plot_dim()`](https://zerostwo.github.io/shennong/dev/reference/sn_plot_dim.md).
  Supply only one size argument.

- alpha:

  Point opacity between zero and one.

- na_value:

  Color for missing group or feature values.

- label:

  Show group labels anchored at observed locations near medians.

- label_size, label_color, label_box, label_halo, repel:

  Label controls as in
  [`sn_plot_dim()`](https://zerostwo.github.io/shennong/dev/reference/sn_plot_dim.md).
  Repulsion requires ggrepel; halos use shadowtext when available.

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

An editable `ggplot` with fixed coordinate aspect and figure export
metadata.
