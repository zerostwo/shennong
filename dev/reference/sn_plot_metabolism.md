# Plot metabolic pathway activity and differential results

Plot metabolic pathway activity and differential results

## Usage

``` r
sn_plot_metabolism(
  x,
  name = NULL,
  type = c("activity", "heatmap", "differential", "sample"),
  pathways = NULL,
  n = 30L
)
```

## Arguments

- x:

  A Seurat object or unified metabolism result.

- name:

  Stored result name when `x` is a Seurat object.

- type:

  Activity distribution, sample heatmap, differential effects, or
  sample-level pathway summary.

- pathways:

  Optional pathways to retain.

- n:

  Maximum pathways shown.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) sn_plot_metabolism(object, "metabolism", type = "differential") # \dontrun{}
```
