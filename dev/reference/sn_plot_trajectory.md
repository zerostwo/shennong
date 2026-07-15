# Plot an inferred trajectory on its embedding

Plot an inferred trajectory on its embedding

## Usage

``` r
sn_plot_trajectory(x, name = NULL, color_by = "cluster", point_size = 0.7)
```

## Arguments

- x:

  A Seurat object or trajectory result.

- name:

  Stored trajectory name when `x` is a Seurat object.

- color_by:

  Cell-table column used for color.

- point_size:

  Point size.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) sn_plot_trajectory(object, "development") # \dontrun{}
```
