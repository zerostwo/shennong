# Plot fitted expression trends for selected genes

Plot fitted expression trends for selected genes

## Usage

``` r
sn_plot_gene_trend(x, name = NULL, features)
```

## Arguments

- x:

  A Seurat object or trajectory result.

- name:

  Stored trajectory name.

- features:

  Features to plot.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) sn_plot_gene_trend(object, "development", c("MKI67", "GZMB")) # \dontrun{}
```
