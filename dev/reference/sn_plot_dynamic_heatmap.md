# Plot fitted dynamic-gene trends as a heatmap

Plot fitted dynamic-gene trends as a heatmap

## Usage

``` r
sn_plot_dynamic_heatmap(
  x,
  name = NULL,
  features = NULL,
  lineage = NULL,
  scale_rows = TRUE
)
```

## Arguments

- x:

  A Seurat object or trajectory result.

- name:

  Stored trajectory name.

- features:

  Optional features to display.

- lineage:

  Optional lineage number or label.

- scale_rows:

  Standardize each feature within lineage.

## Value

A `ggplot` heatmap.

## Examples

``` r
if (FALSE) sn_plot_dynamic_heatmap(object, "development") # \dontrun{}
```
