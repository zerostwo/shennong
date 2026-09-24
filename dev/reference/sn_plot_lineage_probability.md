# Plot lineage assignment probability

Plot lineage assignment probability

## Usage

``` r
sn_plot_lineage_probability(
  x,
  result_id = NULL,
  lineage,
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

  Lineage to display.

- point_size:

  Point size.

- object:

  Alias for `x`; supply only one of `x` and `object`.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) sn_plot_lineage_probability(object, "development", "Lineage1") # \dontrun{}
```
