# Plot CellRank fate probabilities

Plot CellRank fate probabilities

## Usage

``` r
sn_plot_fate(
  x,
  result_id = NULL,
  states = NULL,
  point_size = 0.7,
  object = NULL
)
```

## Arguments

- x:

  A Seurat object or fate result.

- result_id:

  Stored result result_id.

- states:

  Optional terminal states to retain.

- point_size:

  Point size.

- object:

  Alias for `x`; supply only one of `x` and `object`.

## Value

A faceted `ggplot` object.
