# Plot a heatmap for selected genes or features

`sn_plot_heatmap()` draws a compact heatmap for user-selected genes. In
`mode = "cells"` it renders cells with Seurat's heatmap backend, while
`mode = "average"` first averages expression within the requested
metadata groups and displays a group-level heatmap.

## Usage

``` r
sn_plot_heatmap(
  object,
  features,
  group_by = "ident",
  split_by = NULL,
  assay = NULL,
  layer = "scale.data",
  slot = NULL,
  mode = c("cells", "average"),
  scale_data = TRUE,
  max_cells = NULL,
  seed = 717,
  palette = "RdBu",
  direction = -1,
  disp_min = -2.5,
  disp_max = NULL,
  col_min = -2.5,
  col_max = 2.5,
  average_fun = c("mean", "median"),
  average_scale = TRUE,
  group_bar = TRUE,
  label = TRUE,
  raster = TRUE,
  label_size = 8,
  label_angle = 45,
  angle = NULL,
  group_palette = "Paired",
  group_colors = NULL,
  show_cell_names = FALSE,
  show_ticks = FALSE,
  show_legend = TRUE,
  title = NULL,
  collect_legend = TRUE,
  aspect_ratio = NULL,
  panel_widths = NULL,
  panel_heights = NULL,
  ...
)
```

## Arguments

- object:

  A Seurat object.

- features:

  Character vector of genes or features to display.

- group_by:

  Metadata column used to order and label cells. Use `NULL`, `""`, or
  `"ident"` for current Seurat identities.

- split_by:

  Optional metadata column. When supplied, one heatmap is drawn for each
  split level and combined with patchwork.

- assay:

  Assay used for expression retrieval.

- layer:

  Data layer used for expression retrieval. Defaults to `"scale.data"`
  for cell-level heatmaps. In `mode = "average"`, use `"data"` to
  average normalized expression before optional feature-wise scaling.

- slot:

  Deprecated alias for `layer`.

- mode:

  One of `"cells"` or `"average"`. `"cells"` displays cells; `"average"`
  averages expression by `group_by` and optional `split_by`.

- scale_data:

  If `TRUE` and `layer = "scale.data"`, scale the requested features
  when they are missing from the scaled layer. Used only in
  `mode = "cells"`.

- max_cells:

  Optional maximum number of cells to show per cell-level heatmap panel.

- seed:

  Optional seed used when `max_cells` downsamples cells.

- palette:

  Continuous color palette name or vector.

- direction:

  Palette direction.

- disp_min, disp_max:

  Display limits passed to
  [`Seurat::DoHeatmap()`](https://satijalab.org/seurat/reference/DoHeatmap.html).

- col_min, col_max:

  Display limits for `mode = "average"`.

- average_fun:

  Summary used in `mode = "average"`.

- average_scale:

  Whether to feature-wise z-score average expression in
  `mode = "average"`.

- group_bar, label, raster:

  Passed to
  [`Seurat::DoHeatmap()`](https://satijalab.org/seurat/reference/DoHeatmap.html)
  for `mode = "cells"`. Rasterization is enabled by default.

- label_size, label_angle:

  Size in points and rotation angle for the cell-level group labels
  drawn above the heatmap.

- angle:

  Deprecated alias for `label_angle`.

- group_palette, group_colors:

  Palette or explicit colors for the cell-level group label bar.
  Defaults to `"Paired"`.

- show_cell_names, show_ticks:

  Whether to show cell-level x-axis labels and axis ticks in
  `mode = "cells"`. Both default to `FALSE`.

- show_legend:

  Whether to keep the heatmap legend.

- title:

  Optional plot title. When `split_by` is supplied, split levels are
  used as panel titles unless `title` is also supplied.

- collect_legend:

  Whether to collect patchwork legends for split plots.

- aspect_ratio, panel_widths, panel_heights:

  Optional sizing controls forwarded through the Shennong/catplot theme
  helper.

- ...:

  Additional arguments passed to
  [`Seurat::DoHeatmap()`](https://satijalab.org/seurat/reference/DoHeatmap.html).

## Value

A ggplot2 or patchwork object.

## Examples

``` r
if (FALSE) { # \dontrun{
sn_plot_heatmap(
  seurat_obj,
  features = c("CD3D", "MS4A1", "LYZ"),
  group_by = "seurat_clusters"
)
} # }
```
