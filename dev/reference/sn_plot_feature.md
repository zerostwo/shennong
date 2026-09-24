# Plot feature expression in reduced dimensions

This function plots feature expression in reduced dimensions using the
FeaturePlot function from the Seurat package.

## Usage

``` r
sn_plot_feature(
  object,
  features,
  reduction = NULL,
  assay = NULL,
  dims = c(1, 2),
  cells = NULL,
  label = FALSE,
  split_by = NULL,
  label_size = 8 * 0.36,
  pt_size = NULL,
  alpha = 1,
  stroke_size = NULL,
  layer = "data",
  min_cutoff = NA,
  max_cutoff = NA,
  shape_by = NULL,
  blend = FALSE,
  blend_threshold = 0.5,
  ncol = NULL,
  coord_fixed = FALSE,
  by_col = TRUE,
  mode = c("expression", "density"),
  density_method = c("wkde", "ks"),
  density_adjust = 1,
  density_style = c("galaxy", "plain"),
  raster = TRUE,
  raster_dpi = c(512, 512),
  seed = 717,
  title = NULL,
  legend_title = NULL,
  show_legend = TRUE,
  show_axis = FALSE,
  show_border = TRUE,
  palette = "YlOrRd",
  direction = 1,
  legend_labels = c("text", "numeric"),
  keep_scale = c("all", "feature", "none"),
  collect_legend = TRUE,
  aspect_ratio = 1,
  panel_widths = NULL,
  panel_heights = NULL,
  x_label = NULL,
  y_label = NULL,
  style = c("classic", "nebula", "glass"),
  style_control = list(),
  camera = NULL,
  interactive = FALSE,
  group_by = NULL,
  ...
)
```

## Arguments

- object:

  A Seurat object containing the data to plot.

- features:

  A character vector of feature names to plot.

- reduction:

  A character string specifying the dimensionality reduction to use
  (e.g., "PCA", "UMAP", "tSNE"). Defaults to NULL.

- assay:

  Optional assay passed to Seurat's `FeaturePlot()` and used for
  density-mode expression retrieval.

- dims:

  Two dimensions to plot. Defaults to `c(1, 2)`.

- cells:

  Optional cells to include.

- label:

  A character vector specifying the labels to use for each cell group.
  Defaults to label.

- split_by:

  A character vector specifying the cell groups to split the plot by.
  Defaults to NULL.

- label_size:

  A numeric value specifying the size of the labels. Defaults to 8 \*
  0.36.

- pt_size:

  A numeric value specifying the size of the points. When `NULL`,
  Shennong chooses a value automatically based on the number of cells.

- alpha:

  Point alpha passed to Seurat's `FeaturePlot()`.

- stroke_size:

  Point stroke size passed to Seurat's `FeaturePlot()`.

- layer:

  A character string specifying which layer in the Seurat object to use
  (for example, `"data"` or `"scale.data"`). Defaults to `"data"`.

- min_cutoff:

  A numeric value specifying the minimum expression cutoff. Defaults to
  `NA`.

- max_cutoff:

  A numeric value specifying the maximum expression cutoff. Defaults to
  NA.

- shape_by, blend, blend_threshold, ncol, coord_fixed, by_col,
  raster_dpi:

  Snake-case wrappers for the corresponding Seurat `FeaturePlot()`
  arguments.

- mode:

  One of `"expression"` or `"density"`. Density mode computes a
  Nebulosa-style weighted feature density over the selected embedding
  and renders it with a galaxy-like theme by default.

- density_method:

  Density estimator used when `mode = "density"`. One of `"wkde"` or
  `"ks"`. Defaults to `"wkde"`.

- density_adjust:

  Bandwidth adjustment forwarded to the density estimator when
  `mode = "density"`. Larger values smooth more.

- density_style:

  One of `"galaxy"` or `"plain"`. Defaults to `"galaxy"`.

- raster:

  A logical value specifying whether to use raster graphics. Defaults to
  TRUE. When `ggrastr` is installed, Shennong rasterizes the regular
  ggplot point layer after plotting so `pt_size` keeps the same behavior
  as `raster = FALSE`; otherwise it falls back to Seurat's native raster
  backend.

- seed:

  An integer value specifying the random seed. Defaults to 717.

- title:

  A character string specifying the plot title. Defaults to NULL.

- legend_title:

  A character string specifying the legend title. Defaults to NULL.

- show_legend:

  A logical value specifying whether to show the legend. Defaults to
  TRUE.

- show_axis:

  A logical value specifying whether to show the plot axis. Defaults to
  FALSE.

- show_border:

  A logical value specifying whether to show the plot border. Defaults
  to TRUE.

- palette:

  A character string specifying the color palette to use. Defaults to
  "YlOrRd".

- direction:

  A numeric value specifying the direction of the color palette.
  Defaults to 1.

- legend_labels:

  One of `"text"` to show `"Min"` / `"Max"` at the colorbar ends or
  `"numeric"` to show numeric break labels. Defaults to `"text"`.

- keep_scale:

  Passed to Seurat's `FeaturePlot(keep.scale = ...)`. Defaults to
  `"all"` so multi-feature plots share one comparable color scale and
  can collect a single legend.

- collect_legend:

  Logical; when `TRUE`, collect a shared legend for patchwork
  multi-feature plots. Defaults to `TRUE`.

- aspect_ratio:

  Optional panel aspect ratio. Defaults to `1`. When used together with
  `panel_widths` or `panel_heights`, Shennong derives the missing panel
  dimension automatically.

- panel_widths, panel_heights:

  Optional panel size arguments forwarded to
  [`catplot::theme_cat()`](https://rdrr.io/pkg/catplot/man/theme_cat.html)
  when available.

- x_label, y_label:

  Optional axis labels.

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

- group_by:

  For 3D styles, metadata column defining geometry groups and labels;
  defaults to active identities. Expression colors still represent the
  requested feature and assay layer, not the geometry group.

- ...:

  Additional parameters to pass to FeaturePlot.

## Value

A ggplot2/patchwork object, or an htmlwidget when a 3D style uses
`interactive = TRUE`.

## Examples

``` r
if (FALSE) { # \dontrun{
sn_plot_feature(x = mySeuratObject, features = c("CD3D", "CD8A", "CD4"), reduction = "UMAP")
} # }
```
