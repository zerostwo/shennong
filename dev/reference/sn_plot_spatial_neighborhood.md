# Plot spatial-neighborhood enrichment or co-occurrence

Plot spatial-neighborhood enrichment or co-occurrence

## Usage

``` r
sn_plot_spatial_neighborhood(
  x,
  result_id = NULL,
  type = c("enrichment", "cooccurrence"),
  object = NULL
)
```

## Arguments

- x:

  A Seurat object or spatial-neighborhood result.

- result_id:

  Stored result result_id.

- type:

  Enrichment or co-occurrence.

- object:

  Alias for `x`; supply only one of `x` and `object`.

## Value

A `ggplot` object.
