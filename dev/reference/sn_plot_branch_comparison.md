# Plot branch-specific dynamic-gene evidence

Plot branch-specific dynamic-gene evidence

## Usage

``` r
sn_plot_branch_comparison(
  x,
  name = NULL,
  test = c("pattern", "differential_end"),
  n = 20L
)
```

## Arguments

- x:

  A Seurat object or trajectory result.

- name:

  Stored trajectory name.

- test:

  Branch test to display.

- n:

  Maximum number of features.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) sn_plot_branch_comparison(object, "development") # \dontrun{}
```
