# Plot annotation confidence

Plot annotation confidence

## Usage

``` r
sn_plot_annotation_confidence(
  x,
  store_name = "annotation",
  level = c("cluster", "cell")
)
```

## Arguments

- x:

  A Seurat object or annotation result.

- store_name:

  Stored annotation name.

- level:

  Plot cluster- or cell-level confidence.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) sn_plot_annotation_confidence(object) # \dontrun{}
```
