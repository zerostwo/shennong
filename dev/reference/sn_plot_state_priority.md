# Plot cell-state priority scores

Plot cell-state priority scores

## Usage

``` r
sn_plot_state_priority(x, name = NULL, n = 30L)
```

## Arguments

- x:

  A Seurat object or state-priority result.

- name:

  Stored result name when `x` is a Seurat object.

- n:

  Maximum states to show.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) sn_plot_state_priority(object, "priority") # \dontrun{}
```
