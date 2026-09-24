# Plot branch-specific dynamic-gene evidence

Plot branch-specific dynamic-gene evidence

## Usage

``` r
sn_plot_branch_comparison(
  x,
  result_id = NULL,
  test = c("pattern", "differential_end"),
  n = 20L,
  object = NULL
)
```

## Arguments

- x:

  A Seurat object or trajectory result.

- result_id:

  Stored trajectory result_id.

- test:

  Branch test to display.

- n:

  Maximum number of features.

- object:

  Alias for `x`; supply only one of `x` and `object`.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) sn_plot_branch_comparison(object, "development") # \dontrun{}
```
