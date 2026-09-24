# Plot pseudotime on the trajectory embedding

Plot pseudotime on the trajectory embedding

## Usage

``` r
sn_plot_pseudotime(
  x,
  result_id = NULL,
  lineage = NULL,
  point_size = 0.7,
  object = NULL
)
```

## Arguments

- x:

  A Seurat object or trajectory result.

- result_id:

  Stored trajectory result_id.

- lineage:

  Optional lineage. Defaults to primary pseudotime.

- point_size:

  Point size.

- object:

  Alias for `x`; supply only one of `x` and `object`.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) sn_plot_pseudotime(object, "development") # \dontrun{}
```
