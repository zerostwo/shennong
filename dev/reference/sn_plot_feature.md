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
  slot = "data",
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

- slot:

  A character string specifying which slot in the Seurat object to use
  (e.g., "data", "scale.data", "integrated"). Defaults to "data".

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

- ...:

  Additional parameters to pass to FeaturePlot.

## Value

A ggplot2 object.

## Examples

``` r
if (FALSE) { # \dontrun{
sn_plot_feature(x = mySeuratObject, features = c("CD3D", "CD8A", "CD4"), reduction = "UMAP")
} # }
```
