# Plot miloR neighborhood differential-abundance results

Plot miloR neighborhood differential-abundance results

## Usage

``` r
sn_plot_milo(
  x,
  milo_name = "default",
  annotation_col = NULL,
  fdr_col = c("SpatialFDR", "FDR"),
  fdr_cutoff = 0.1,
  logfc_cutoff = NULL,
  palette = "Paired",
  title = NULL,
  x_label = "logFC",
  y_label = expression(-log[10]("FDR")),
  panel_widths = NULL,
  panel_heights = NULL
)
```

## Arguments

- x:

  A milo result data frame, or a Seurat object containing a stored milo
  result.

- milo_name:

  Name of the stored milo result when `x` is a Seurat object.

- annotation_col:

  Optional column used to color points.

- fdr_col:

  FDR column used on the y-axis. Defaults to `"SpatialFDR"`.

- fdr_cutoff:

  Horizontal significance threshold.

- logfc_cutoff:

  Optional vertical threshold for effect size.

- palette:

  Discrete palette used when `annotation_col` is supplied.

- title, x_label, y_label:

  Optional plot labels.

- panel_widths, panel_heights:

  Optional panel size arguments forwarded to
  [`catplot::theme_cat()`](https://rdrr.io/pkg/catplot/man/theme_cat.html)
  when available.

## Value

A ggplot object.

## Examples

``` r
milo_df <- data.frame(
  logFC = c(1.2, -0.8, 0.4),
  SpatialFDR = c(0.01, 0.2, 0.04),
  cell_type = c("T", "B", "Myeloid")
)
sn_plot_milo(milo_df, annotation_col = "cell_type")

```
