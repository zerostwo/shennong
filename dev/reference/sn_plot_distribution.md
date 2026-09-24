# Plot numeric distributions from Seurat metadata, expression, or tables

`sn_plot_distribution()` consolidates routine violin, box, histogram,
density, ridge, and QC-distribution views under one object-first
interface. When `sample_by` is supplied, observations are summarized to
biological samples before plotting.

## Usage

``` r
sn_plot_distribution(
  object,
  features,
  group_by = NULL,
  sample_by = NULL,
  assay = NULL,
  layer = "data",
  view = c("violin", "box", "histogram", "density", "ridge"),
  aggregate_fun = c("mean", "median", "sum"),
  thresholds = list(),
  bins = 30L,
  show_points = FALSE,
  point_alpha = 0.6,
  jitter_width = 0.15,
  palette = "Paired",
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  aspect_ratio = NULL,
  panel_widths = NULL,
  panel_heights = NULL
)
```

## Arguments

- object:

  A Seurat object or data frame.

- features:

  Numeric metadata/expression features or data-frame columns.

- group_by:

  Optional metadata/data-frame grouping column.

- sample_by:

  Optional biological-sample column used before plotting.

- assay, layer:

  Seurat expression source.

- view:

  One of `"violin"`, `"box"`, `"histogram"`, `"density"`, or `"ridge"`.

- aggregate_fun:

  Sample summary when `sample_by` is supplied.

- thresholds:

  Optional named list of lower/upper reference values for the requested
  features.

- bins:

  Histogram bins.

- show_points:

  Overlay observations on violin or box plots.

- point_alpha, jitter_width:

  Point transparency and horizontal jitter.

- palette:

  Named Shennong palette or explicit colors.

- title, x_label, y_label:

  Optional plot labels.

- aspect_ratio, panel_widths, panel_heights:

  Optional panel sizing controls.

## Value

A faceted `ggplot` object with source data and figure metadata.

## Examples

``` r
values <- data.frame(
  sample = rep(c("S1", "S2"), each = 3),
  condition = rep(c("control", "treated"), each = 3),
  score = c(1, 2, 3, 2, 4, 5)
)
sn_plot_distribution(values, "score", group_by = "condition", view = "box")
```
