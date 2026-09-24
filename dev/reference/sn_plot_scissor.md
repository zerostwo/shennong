# Plot a unified Scissor result

Plot a unified Scissor result

## Usage

``` r
sn_plot_scissor(
  x,
  result_id = "scissor",
  type = c("states", "cells", "samples", "correlations", "reliability"),
  n = 5000L,
  object = NULL
)
```

## Arguments

- x:

  A Scissor state-priority result or a Seurat object containing one.

- result_id:

  Stored result result_id when `x` is a Seurat object.

- type:

  Plot state ranking, cell coefficients, sample contributions,
  cell-level correlation summaries, or optional reliability output.

- n:

  Maximum states or cells shown.

- object:

  Alias for `x`; supply only one of `x` and `object`.

## Value

A ggplot object.

## Examples

``` r
if (FALSE) { # \dontrun{
sn_plot_scissor(object, result_id = "scissor", type = "states")
sn_plot_scissor(object, result_id = "scissor", type = "correlations")
} # }
```
