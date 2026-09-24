# Plot fitted dynamic-gene trends as a heatmap

Plot fitted dynamic-gene trends as a heatmap

## Usage

``` r
sn_plot_dynamic_heatmap(
  x,
  result_id = NULL,
  features = NULL,
  lineage = NULL,
  scale_rows = TRUE,
  object = NULL
)
```

## Arguments

- x:

  A Seurat object or trajectory result.

- result_id:

  Stored trajectory result_id.

- features:

  Optional features to display.

- lineage:

  Optional lineage number or label.

- scale_rows:

  Standardize each feature within lineage.

- object:

  Alias for `x`; supply only one of `x` and `object`.

## Value

A `ggplot` heatmap.

## Examples

``` r
if (FALSE) sn_plot_dynamic_heatmap(object, "development") # \dontrun{}
```
