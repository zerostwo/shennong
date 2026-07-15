# Plot a gene regulatory network result

Plot a gene regulatory network result

## Usage

``` r
sn_plot_regulon(
  x,
  name = NULL,
  type = c("network", "activity", "specificity"),
  regulons = NULL,
  n = 50L
)
```

## Arguments

- x:

  A Seurat object or GRN result.

- name:

  Stored result name when `x` is a Seurat object.

- type:

  Network edges, regulon activity, or group specificity.

- regulons:

  Optional regulators/regulons to retain.

- n:

  Maximum network edges.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) sn_plot_regulon(object, "grn_genie3", type = "specificity") # \dontrun{}
```
