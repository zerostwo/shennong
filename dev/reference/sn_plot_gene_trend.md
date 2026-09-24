# Plot fitted expression trends for selected genes

Plot fitted expression trends for selected genes

## Usage

``` r
sn_plot_gene_trend(x, result_id = NULL, features, object = NULL)
```

## Arguments

- x:

  A Seurat object or trajectory result.

- result_id:

  Stored trajectory result_id.

- features:

  Features to plot.

- object:

  Alias for `x`; supply only one of `x` and `object`.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) sn_plot_gene_trend(object, "development", c("MKI67", "GZMB")) # \dontrun{}
```
