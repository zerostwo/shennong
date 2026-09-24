# Plot sample-aware differential communication effects

Plot sample-aware differential communication effects

## Usage

``` r
sn_plot_communication_comparison(x, result_id = NULL, n = 30L, object = NULL)
```

## Arguments

- x:

  A Seurat object or unified communication result.

- result_id:

  Stored result result_id when `x` is a Seurat object.

- n:

  Maximum number of effects.

- object:

  Alias for `x`; supply only one of `x` and `object`.

## Value

A `ggplot` object.
