# Plot metabolic pathway activity and differential results

Plot metabolic pathway activity and differential results

## Usage

``` r
sn_plot_metabolism(
  x,
  result_id = NULL,
  type = c("activity", "heatmap", "differential", "sample"),
  pathways = NULL,
  n = 30L,
  object = NULL
)
```

## Arguments

- x:

  A Seurat object or unified metabolism result.

- result_id:

  Stored result result_id when `x` is a Seurat object.

- type:

  Activity distribution, sample heatmap, differential effects, or
  sample-level pathway summary.

- pathways:

  Optional pathways to retain.

- n:

  Maximum pathways shown.

- object:

  Alias for `x`; supply only one of `x` and `object`.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) sn_plot_metabolism(object, "metabolism", type = "differential") # \dontrun{}
```
