# Plot differential-abundance effects

Plot differential-abundance effects

## Usage

``` r
sn_plot_abundance(x, name = NULL, n = 30L, adjusted_p_value = NULL)
```

## Arguments

- x:

  A Seurat object or differential-abundance result.

- name:

  Stored result name when `x` is a Seurat object.

- n:

  Maximum number of rows to display.

- adjusted_p_value:

  Optional adjusted-p-value cutoff.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) sn_plot_abundance(object, "abundance") # \dontrun{}
```
