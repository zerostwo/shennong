# Plot CNV, malignancy, and subclone results

Plot CNV, malignancy, and subclone results

## Usage

``` r
sn_plot_cnv(
  x,
  result_id = NULL,
  type = c("heatmap", "umap", "score", "sample", "association"),
  n = 100L,
  object = NULL
)
```

## Arguments

- x:

  A Seurat object or unified CNV result.

- result_id:

  Stored result result_id when `x` is a Seurat object.

- type:

  Plot type: chromosome heatmap, CNV UMAP, score distribution, sample
  summary, or CNV-expression association.

- n:

  Maximum cells or features shown.

- object:

  Alias for `x`; supply only one of `x` and `object`.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) sn_plot_cnv(object, "cnv", type = "heatmap") # \dontrun{}
```
