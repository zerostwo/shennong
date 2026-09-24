# Deprecated direct scDesign3 wrapper

`sn_simulate_scdesign3()` is a deprecated compatibility entry point. Use
[`sn_simulate()`](https://zerostwo.github.io/shennong/dev/reference/sn_simulate.md)
with `method = "scdesign3"` instead.

## Usage

``` r
sn_simulate_scdesign3(object, ...)
```

## Arguments

- object:

  A Seurat or SingleCellExperiment object.

- ...:

  Additional arguments passed to the scDesign3 backend.

## Value

Simulated data in the backend's requested format.
