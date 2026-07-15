# Plot lineage assignment probability

Plot lineage assignment probability

## Usage

``` r
sn_plot_lineage_probability(x, name = NULL, lineage, point_size = 0.7)
```

## Arguments

- x:

  A Seurat object or trajectory result.

- name:

  Stored trajectory name.

- lineage:

  Lineage to display.

- point_size:

  Point size.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) sn_plot_lineage_probability(object, "development", "Lineage1") # \dontrun{}
```
