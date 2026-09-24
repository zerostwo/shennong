# Plot differential-abundance effects

Plot differential-abundance effects

## Usage

``` r
sn_plot_abundance(
  x,
  result_id = NULL,
  n = 30L,
  adjusted_p_value = NULL,
  object = NULL
)
```

## Arguments

- x:

  A Seurat object or differential-abundance result.

- result_id:

  Stored result result_id when `x` is a Seurat object.

- n:

  Maximum number of rows to display.

- adjusted_p_value:

  Optional adjusted-p-value cutoff.

- object:

  Alias for `x`; supply only one of `x` and `object`.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) sn_plot_abundance(object, "abundance") # \dontrun{}
```
