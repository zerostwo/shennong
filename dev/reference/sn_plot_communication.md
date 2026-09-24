# Plot standardized cell-cell communication results

Plot standardized cell-cell communication results

## Usage

``` r
sn_plot_communication(
  x,
  result_id = NULL,
  type = c("bubble", "heatmap", "network", "chord", "river"),
  n = 50L,
  table = "primary",
  object = NULL
)
```

## Arguments

- x:

  A Seurat object, unified communication result, or standardized table.

- result_id:

  Stored result result_id when `x` is a Seurat object.

- type:

  Plot type: bubble, heatmap, network, chord, or river.

- n:

  Maximum number of top-ranked interactions to display.

- table:

  Stored table to plot.

- object:

  Alias for `x`; supply only one of `x` and `object`.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) sn_plot_communication(object, "communication", type = "bubble") # \dontrun{}
```
