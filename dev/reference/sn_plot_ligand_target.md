# Plot NicheNet or MultiNicheNet ligand-target evidence

Plot NicheNet or MultiNicheNet ligand-target evidence

## Usage

``` r
sn_plot_ligand_target(x, result_id = NULL, n = 50L, object = NULL)
```

## Arguments

- x:

  A Seurat object or unified communication result.

- result_id:

  Stored result result_id when `x` is a Seurat object.

- n:

  Maximum number of ligand-target links.

- object:

  Alias for `x`; supply only one of `x` and `object`.

## Value

A `ggplot` object.
