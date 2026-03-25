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
  pt_size = 2,
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
  show_axis = TRUE,
  show_border = TRUE,
  title = NULL,
  palette = "Paired",
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

  The size of the points on the plot. Default is 1.

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
  is TRUE.

- show_border:

  Logical value indicating whether to show the panel and axis borders on
  the plot. Default is TRUE.

- title:

  The title for the plot. Default is NULL.

- palette:

  The color palette to use for the plot. Default is "Paired".

- ...:

  Additional parameters to be passed to the DimPlot() function in
  Seurat.

## Value

A ggplot2 object containing the dimensionality reduction plot.

## Examples

``` r
if (FALSE) { # \dontrun{
data("pbmc_small", package = "Shennong")
pbmc <- sn_run_cluster(pbmc, normalization_method = "seurat", verbose = FALSE)
sn_plot_dim(
  object = pbmc,
  reduction = "umap",
  group_by = "seurat_clusters",
  palette = "Set1"
)
} # }
```
