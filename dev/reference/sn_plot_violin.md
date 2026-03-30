# Plot a violin plot with categorical groups

This function plots a violin plot with categorical groups using the
VlnPlot function from the Seurat package.

## Usage

``` r
sn_plot_violin(
  object,
  features,
  pt_size = 0,
  sort = FALSE,
  group_by = NULL,
  split_by = NULL,
  show_legend = FALSE,
  angle_x = 0,
  palette = "Paired",
  aspect_ratio = 0.5,
  panel_widths = NULL,
  panel_heights = NULL,
  x_label = NULL,
  y_label = NULL,
  title = NULL
)
```

## Arguments

- object:

  A Seurat object containing the data to plot.

- features:

  A character vector of feature names to plot.

- pt_size:

  The size of the points to plot.

- sort:

  Whether to sort the features by their mean expression or not.

- group_by:

  A character vector specifying the grouping variable. Defaults to
  "ident".

- split_by:

  A character vector specifying the splitting variable.

- show_legend:

  Whether to show the legend or not. Defaults to FALSE.

- angle_x:

  The angle of the x-axis labels. Defaults to 0.

- palette:

  The color palette to use. Defaults to `"Paired"`.

- aspect_ratio:

  The aspect ratio of the plot. Defaults to 0.5. When used together with
  `panel_widths` or `panel_heights`, Shennong derives the missing panel
  dimension automatically.

- panel_widths, panel_heights:

  Optional panel size arguments forwarded to
  [`catplot::theme_cat()`](https://rdrr.io/pkg/catplot/man/theme_cat.html)
  when available.

- x_label, y_label, title:

  Optional plot labels.

## Value

A ggplot2 object.

## Examples

``` r
if (FALSE) { # \dontrun{
sn_plot_violin(object = mySeuratObject, features = c("CD3D", "CD8A", "CD4"))
} # }
```
