# Create a dimensionality reduction plot for categorical data

This function creates a dimensionality reduction plot for categorical
data using Seurat and ggplot2. It allows for the selection of the
reduction method, grouping, and splitting variables, as well as the
visualization of labels, rasterization, and color palette. The
sn_plot_dim() function is intended to be used as a wrapper around
Seurat's DimPlot() function.

## Usage

``` r
sn_plot_dim(
  object,
  dims = c(1, 2),
  cells = NULL,
  cols = NULL,
  pt_size = NULL,
  reduction = NULL,
  group_by = NULL,
  split_by = NULL,
  shape_by = NULL,
  order = NULL,
  shuffle = FALSE,
  seed = 717,
  label = FALSE,
  label_size = 8 * 0.36,
  label_color = "black",
  label_box = FALSE,
  label_halo = TRUE,
  repel = FALSE,
  cells_highlight = NULL,
  cols_highlight = "#DE2D26",
  sizes_highlight = 1,
  na_value = "grey50",
  ncol = NULL,
  combine = TRUE,
  raster = TRUE,
  raster_dpi = c(512, 512),
  show_legend = TRUE,
  show_axis = FALSE,
  show_border = TRUE,
  title = NULL,
  palette = "Paired",
  aspect_ratio = 1,
  panel_widths = NULL,
  panel_heights = NULL,
  style = c("classic", "nebula", "glass"),
  style_control = list(),
  camera = NULL,
  interactive = FALSE,
  ...
)
```

## Arguments

- object:

  A Seurat object containing categorical data.

- dims:

  The dimensions to plot. Default is c(1, 2).

- cells:

  The cells to plot. Default is NULL.

- cols:

  The columns to plot. Default is NULL.

- pt_size:

  The size of the points on the plot. When `NULL`, Shennong chooses a
  value automatically based on the number of cells.

- reduction:

  The dimensionality reduction method. Default is NULL.

- group_by:

  The variable to group data by. Default is NULL.

- split_by:

  The variable to split data by. Default is NULL.

- shape_by:

  The variable to shape data by. Default is NULL.

- order:

  The order to plot the data in. Default is NULL.

- shuffle:

  Logical value indicating whether to shuffle the data before plotting.
  Default is FALSE.

- seed:

  The random seed to use for shuffling the data. Default is 1.

- label:

  Logical value indicating whether to show labels on the plot. Default
  is FALSE.

- label_size:

  The size of the labels on the plot. Default is 8 \* 0.36.

- label_color:

  The color of the labels on the plot. Default is "black".

- label_box:

  Logical value indicating whether to show a box around the labels on
  the plot. Default is FALSE.

- label_halo:

  Logical value indicating whether Shennong should add a white
  halo/background behind text labels. Set to `FALSE` to keep Seurat's
  native label layer unchanged.

- repel:

  Logical value indicating whether to use point repulsion to avoid
  overlapping labels. Default is TRUE.

- cells_highlight:

  The cells to highlight on the plot. Default is NULL.

- cols_highlight:

  The columns to highlight on the plot. Default is NULL.

- sizes_highlight:

  The sizes to highlight on the plot. Default is NULL.

- na_value:

  The value to use for missing data. Default is "grey50".

- ncol:

  The number of columns to use for the plot. Default is NULL.

- combine:

  Logical value indicating whether to combine the plots into a single
  plot. Default is TRUE.

- raster:

  Logical value indicating whether to use rasterization for improved
  performance. Default is TRUE.

- raster_dpi:

  The DPI to use for rasterization. Default is c(512, 512).

- show_legend:

  Logical value indicating whether to show the legend on the plot.
  Default is TRUE.

- show_axis:

  Logical value indicating whether to show the axis on the plot. Default
  is FALSE.

- show_border:

  Logical value indicating whether to show the panel and axis borders on
  the plot. Default is TRUE.

- title:

  The title for the plot. Default is NULL.

- palette:

  The color palette to use for the plot. Default is "Paired".

- aspect_ratio:

  Optional panel aspect ratio. Defaults to `1`. When used together with
  `panel_widths` or `panel_heights`, Shennong derives the missing panel
  dimension automatically.

- panel_widths, panel_heights:

  Optional panel size arguments forwarded to
  [`catplot::theme_cat()`](https://rdrr.io/pkg/catplot/man/theme_cat.html)
  when available.

- style:

  Rendering style: `"classic"` preserves the existing 2D plot;
  `"nebula"` and `"glass"` use planar density contours with the default
  `dims = c(1, 2)`, or real three-dimensional embeddings with
  `dims = 1:3`. Static 3D output is a ggplot with the complete
  point/surface scene rasterized at 600 dpi by default; labels and
  legends remain vector elements. Set `raster_dpi = 600` explicitly when
  exporting. Requires misc3d and htmlwidgets; static export additionally
  uses chromote, png, and Chrome/Chromium.

- style_control:

  Named 3D controls: `surface_alpha` (glass 0.16, nebula 0.035),
  `point_alpha` (0.85 for dimension plots), `glow` (glass 0.3, nebula
  0.85), `surface_mass` (0.95), `bandwidth` (0.75), `grid_size` (48,
  integer 16–64), `background` (glass `"#03030C"`, nebula `"#101322"`),
  and `auto_rotate` (FALSE, browser only). Surfaces are binned Gaussian
  KDE isosurfaces of group coordinates, not expression contours or
  biological boundaries. Groups with fewer than five cells or
  rank-deficient coordinates retain points without a surface.

- camera:

  Camera list or downloaded JSON accepted by
  [`sn_get_plot_camera()`](https://zerostwo.github.io/shennong/dev/reference/sn_get_plot_camera.md).

- interactive:

  If TRUE, return a local WebGL htmlwidget with drag rotation,
  shift-drag pan, scroll zoom, and camera export controls. Requires
  htmlwidgets; currently one panel only. FALSE returns a ggplot
  compatible with `ggsave()`. Browser and PDF use the same WebGL shader;
  GPU/CPU antialiasing and vector text layout can differ. 3D styles have
  a square, borderless panel, white labels and numeric feature legends.
  Seurat-only shape, highlight, repel, blend, density-mode and extra
  `...` options are rejected.

- ...:

  Additional parameters to be passed to the DimPlot() function in
  Seurat.

## Value

A ggplot2/patchwork object, a list when `combine = FALSE`, or an
htmlwidget when `interactive = TRUE` for a 3D style.

## Examples

``` r
if (FALSE) { # \dontrun{
pbmc <- qs2::qs_read(file.path(
  Sys.getenv("SHENNONG_REAL_DATA_DIR"), "single-cell", "kotliarov_pbmc.qs2"
))
pbmc <- sn_run_cluster(pbmc, normalization_method = "seurat", verbose = FALSE)
sn_plot_dim(
  object = pbmc,
  reduction = "umap",
  group_by = "seurat_clusters",
  palette = "Set1"
)
} # }
```
