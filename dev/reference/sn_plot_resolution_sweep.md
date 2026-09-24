# Plot a clustering resolution sweep

Plot a clustering resolution sweep

## Usage

``` r
sn_plot_resolution_sweep(x, metrics = NULL, object = NULL)
```

## Arguments

- x:

  Result from
  [`sn_sweep_cluster_resolution()`](https://zerostwo.github.io/shennong/dev/reference/sn_sweep_cluster_resolution.md)
  or its summary table.

- metrics:

  Numeric metric columns to display.

- object:

  Alias for `x`; supply only one of `x` and `object`.

## Value

A resolution-quality plot.
