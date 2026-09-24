# Plot cell-state priority scores

Plot cell-state priority scores

## Usage

``` r
sn_plot_state_priority(x, result_id = NULL, n = 30L, object = NULL)
```

## Arguments

- x:

  A Seurat object or state-priority result.

- result_id:

  Stored result result_id when `x` is a Seurat object.

- n:

  Maximum states to show.

- object:

  Alias for `x`; supply only one of `x` and `object`.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) sn_plot_state_priority(object, "priority") # \dontrun{}
```
