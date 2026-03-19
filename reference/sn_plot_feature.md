# Plot feature expression in reduced dimensions

This function plots feature expression in reduced dimensions using the
FeaturePlot function from the Seurat package.

## Usage

``` r
sn_plot_feature(
  object,
  features,
  reduction = NULL,
  label = label,
  split_by = NULL,
  label_size = 8 * 0.36,
  pt_size = 1,
  slot = "data",
  max_cutoff = NA,
  raster = TRUE,
  seed = 717,
  title = NULL,
  legend_title = NULL,
  show_legend = TRUE,
  show_axis = TRUE,
  show_border = TRUE,
  palette = "YlOrRd",
  direction = 1,
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

  A numeric value specifying the size of the points. Defaults to 1.

- slot:

  A character string specifying which slot in the Seurat object to use
  (e.g., "data", "scale.data", "integrated"). Defaults to "data".

- max_cutoff:

  A numeric value specifying the maximum expression cutoff. Defaults to
  NA.

- raster:

  A logical value specifying whether to use raster graphics. Defaults to
  TRUE.

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
  TRUE.

- show_border:

  A logical value specifying whether to show the plot border. Defaults
  to TRUE.

- palette:

  A character string specifying the color palette to use. Defaults to
  "YlOrRd".

- direction:

  A numeric value specifying the direction of the color palette.
  Defaults to 1.

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
