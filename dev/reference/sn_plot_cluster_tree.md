# Plot cluster transitions across resolutions

Plot cluster transitions across resolutions

## Usage

``` r
sn_plot_cluster_tree(object, resolution_cols = NULL)
```

## Arguments

- object:

  A Seurat object with multiple resolution metadata columns, or a data
  frame containing those assignments.

- resolution_cols:

  Assignment columns ordered from low to high resolution.

## Value

A cluster-transition graph.
