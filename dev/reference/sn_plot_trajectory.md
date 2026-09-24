# Plot an inferred trajectory on its embedding

Plot an inferred trajectory on its embedding

## Usage

``` r
sn_plot_trajectory(
  x,
  result_id = NULL,
  color_by = "cluster",
  point_size = 0.7,
  object = NULL
)
```

## Arguments

- x:

  A Seurat object or trajectory result.

- result_id:

  Stored trajectory result_id when `x` is a Seurat object.

- color_by:

  Cell-table column used for color.

- point_size:

  Point size.

- object:

  Alias for `x`; supply only one of `x` and `object`.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) sn_plot_trajectory(object, "development") # \dontrun{}
```
