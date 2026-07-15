# Plot bulk sample PCA

Plot bulk sample PCA

## Usage

``` r
sn_plot_bulk_pca(x, metadata = NULL, color_by = NULL, pc_x = 1L, pc_y = 2L)
```

## Arguments

- x:

  A bulk-QC result.

- metadata:

  Optional sample metadata used for color labels.

- color_by:

  Optional metadata column.

- pc_x, pc_y:

  Principal components to display.

## Value

A `ggplot` object.
