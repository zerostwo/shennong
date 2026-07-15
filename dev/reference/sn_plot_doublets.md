# Plot doublet classifications in an embedding

Plot doublet classifications in an embedding

## Usage

``` r
sn_plot_doublets(object, class_col = NULL, score_col = NULL, reduction = NULL)
```

## Arguments

- object:

  A Seurat object after
  [`sn_find_doublets()`](https://songqi.org/shennong/dev/reference/sn_find_doublets.md).

- class_col, score_col:

  Metadata columns for doublet call and score.

- reduction:

  Reduction to plot.

## Value

A doublet embedding plot.
