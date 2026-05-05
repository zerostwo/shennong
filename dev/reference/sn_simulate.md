# Simulate single-cell counts with scDesign3

`sn_simulate()` provides a method-based simulation entry point.
Currently `method = "scdesign3"` delegates to
[`sn_simulate_scdesign3()`](https://songqi.org/shennong/dev/reference/sn_simulate_scdesign3.md).

## Usage

``` r
sn_simulate(object, method = c("scdesign3"), ...)
```

## Arguments

- object:

  A Seurat or SingleCellExperiment object.

- method:

  Simulation backend. Currently supports `"scdesign3"`.

- ...:

  Additional arguments passed to the selected backend.

## Value

Simulated data in the backend's requested format.

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- sn_simulate(seurat_obj, method = "scdesign3", celltype = "cell_type")
} # }
```
