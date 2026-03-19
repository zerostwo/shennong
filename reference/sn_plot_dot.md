# Plot a dot plot with categorical groups

This function plots a dot plot with categorical groups using the DotPlot
function from the Seurat package.

## Usage

``` r
sn_plot_dot(
  x,
  assay = NULL,
  features,
  de_name = "default",
  n = 3,
  marker_groups = NULL,
  col_min = -2.5,
  col_max = 2.5,
  dot_min = 0,
  dot_scale = 4,
  idents = NULL,
  group_by = NULL,
  split_by = NULL,
  cluster_idents = FALSE,
  scale = TRUE,
  scale_by = "radius",
  scale_min = NA,
  scale_max = NA,
  palette = "RdBu"
)
```

## Arguments

- x:

  A Seurat object containing the data to plot.

- assay:

  The assay to plot. Defaults to the default assay.

- features:

  A character vector of feature names to plot, or `"top_markers"` to
  automatically use the top stored DE markers from
  `object@misc$de_results[[de_name]]`.

- de_name:

  Name of the stored DE result to use when `features = "top_markers"`.
  Defaults to `"default"`.

- n:

  Number of genes to select per group when `features = "top_markers"`.
  Defaults to `3`.

- marker_groups:

  Optional subset of DE result groups to include when
  `features = "top_markers"`.

- col_min:

  The minimum value for the color scale. Defaults to -2.5.

- col_max:

  The maximum value for the color scale. Defaults to 2.5.

- dot_min:

  The minimum value for the dot size scale. Defaults to 0.

- dot_scale:

  The size of the dots to plot.

- idents:

  A character vector of identities to plot. Defaults to all identities.

- group_by:

  A character vector specifying the grouping variable. Defaults to
  "ident".

- split_by:

  A character vector specifying the splitting variable.

- cluster_idents:

  Whether to cluster the identities or not. Defaults to FALSE.

- scale:

  Whether to scale the dot size or not. Defaults to TRUE.

- scale_by:

  The variable to scale the dot size by. Defaults to "radius".

- scale_min:

  The minimum value for the dot size scale. Defaults to NA.

- scale_max:

  The maximum value for the dot size scale. Defaults to NA.

- palette:

  The diverging palette used for the color scale.

## Value

A ggplot2 object.

## Examples

``` r
if (FALSE) { # \dontrun{
sn_plot_dot(x = mySeuratObject, features = c("CD3D", "CD8A", "CD4"))
} # }
```
