# Plot WGCNA modules or trait associations

Plot WGCNA modules or trait associations

## Usage

``` r
sn_plot_wgcna(x, type = c("traits", "modules"), object = NULL)
```

## Arguments

- x:

  A bulk-network result.

- type:

  Module sizes or module-trait association heatmap.

- object:

  Alias for `x`; supply only one of `x` and `object`.

## Value

A `ggplot` object.
