# Plot ambient RNA correction totals

Plot ambient RNA correction totals

## Usage

``` r
sn_plot_ambient_correction(
  before,
  after = NULL,
  assay = NULL,
  layer_before = "counts",
  layer_after = "decontaminated_counts",
  n = 15L
)
```

## Arguments

- before:

  Raw/uncorrected feature-by-cell matrix, Seurat object, or list
  containing `before` and `after` matrices.

- after:

  Corrected matrix when not supplied in `before`.

- assay, layer_before, layer_after:

  Seurat assay/layers.

- n:

  Number of most changed features labeled.

## Value

A before/after feature-total plot.
