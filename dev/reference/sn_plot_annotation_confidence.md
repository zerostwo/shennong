# Plot annotation confidence

Plot annotation confidence

## Usage

``` r
sn_plot_annotation_confidence(
  x,
  result_id = "annotation",
  level = c("cluster", "cell")
)
```

## Arguments

- x:

  A Seurat object or annotation result.

- result_id:

  Stored annotation name.

- level:

  Plot cluster- or cell-level confidence.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) sn_plot_annotation_confidence(object) # \dontrun{}
```
