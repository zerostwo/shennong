# Plot integration assessment metrics

Plot integration assessment metrics

## Usage

``` r
sn_plot_integration(x, aggregate = FALSE, object = NULL)
```

## Arguments

- x:

  Result from
  [`sn_assess_integration()`](https://zerostwo.github.io/shennong/dev/reference/sn_assess_integration.md)
  or its summary table.

- aggregate:

  Include aggregate rows.

- object:

  Alias for `x`; supply only one of `x` and `object`.

## Value

An integration score plot.
