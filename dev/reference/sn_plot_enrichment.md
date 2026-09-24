# Plot enrichment results

Plot enrichment results

## Usage

``` r
sn_plot_enrichment(
  result,
  type = c("dot", "bar", "ridge", "network", "emap"),
  n = 20L,
  minimum_overlap = 0.1,
  object = NULL,
  result_id = "default"
)
```

## Arguments

- result:

  A Shennong enrichment result or enrichment table.

- type:

  Dot, bar, ridge, network, or enrichment-map view.

- n:

  Maximum pathways.

- minimum_overlap:

  Minimum Jaccard overlap for network edges.

- object:

  Optional Seurat object holding a stored enrichment result; supply
  either `result` or `object`, not both.

- result_id:

  Stored enrichment result name used when `object` is supplied.

## Value

A result-aware `ggplot` with source data and figure specification.
