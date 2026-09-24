# Plot a gene regulatory network result

Plot a gene regulatory network result

## Usage

``` r
sn_plot_regulon(
  x,
  result_id = NULL,
  type = c("network", "activity", "specificity"),
  regulons = NULL,
  n = 50L,
  object = NULL
)
```

## Arguments

- x:

  A Seurat object or GRN result.

- result_id:

  Stored result result_id when `x` is a Seurat object.

- type:

  Network edges, regulon activity, or group specificity.

- regulons:

  Optional regulators/regulons to retain.

- n:

  Maximum network edges.

- object:

  Alias for `x`; supply only one of `x` and `object`.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) sn_plot_regulon(object, "grn_genie3", type = "specificity") # \dontrun{}
```
