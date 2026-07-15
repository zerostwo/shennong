# Plot cell-level QC thresholds

Plot cell-level QC thresholds

## Usage

``` r
sn_plot_qc_thresholds(
  x,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  thresholds = list(),
  sample_by = NULL,
  max_cells = 20000L
)
```

## Arguments

- x:

  A Seurat object or cell metadata data frame.

- features:

  QC metadata columns.

- thresholds:

  Named list with one/two numeric lower/upper thresholds.

- sample_by:

  Optional sample column used for color.

- max_cells:

  Maximum plotted cells.

## Value

A faceted QC distribution plot.
