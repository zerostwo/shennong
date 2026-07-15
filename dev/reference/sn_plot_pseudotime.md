# Plot pseudotime on the trajectory embedding

Plot pseudotime on the trajectory embedding

## Usage

``` r
sn_plot_pseudotime(x, name = NULL, lineage = NULL, point_size = 0.7)
```

## Arguments

- x:

  A Seurat object or trajectory result.

- name:

  Stored trajectory name.

- lineage:

  Optional lineage. Defaults to primary pseudotime.

- point_size:

  Point size.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) sn_plot_pseudotime(object, "development") # \dontrun{}
```
