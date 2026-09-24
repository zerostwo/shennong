# Plot a numeric association or sample-to-sample expression correlation

This is the canonical scatter/correlation entry point for tabular,
Seurat, and bulk expression inputs. Seurat observations can optionally
be aggregated to biological samples before correlation. For a
feature-by-sample matrix, `x` and `y` name the two sample columns and
each point is a feature.

## Usage

``` r
sn_plot_association(
  object,
  x,
  y,
  sample_by = NULL,
  group_by = NULL,
  assay = NULL,
  layer = "data",
  aggregate_fun = c("mean", "median", "sum"),
  method = c("spearman", "pearson", "kendall"),
  features = NULL,
  transform = c("none", "log1p"),
  add_fit = TRUE,
  label = FALSE,
  max_points = NULL,
  seed = 717,
  point_size = 1.5,
  point_alpha = 0.7,
  palette = "Paired",
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  aspect_ratio = 1,
  panel_widths = NULL,
  panel_heights = NULL
)
```

## Arguments

- object:

  A Seurat object, data frame, numeric feature-by-sample matrix,
  `SummarizedExperiment`, bulk input list, or bulk-QC result containing
  `tables$expression`.

- x, y:

  Numeric metadata/feature names, data-frame columns, or sample columns
  for matrix-like input.

- sample_by:

  Optional Seurat metadata or data-frame column used to aggregate
  cell/observation values to biological samples.

- group_by:

  Optional metadata/data-frame column mapped to color. With `sample_by`,
  it must be constant within each sample.

- assay, layer:

  Seurat expression source used when `x` or `y` is a feature.

- aggregate_fun:

  Sample summary: `"mean"`, `"median"`, or `"sum"`.

- method:

  Correlation method.

- features:

  Optional feature subset for matrix-like input.

- transform:

  Optional transformation applied to both numeric axes.

- add_fit:

  Add an ordinary least-squares trend line.

- label:

  Label observations/features.

- max_points:

  Optional deterministic maximum number of displayed points.

- seed:

  Seed used when `max_points` downsamples observations.

- point_size, point_alpha:

  Point size and transparency.

- palette:

  Named Shennong palette or explicit colors.

- title, x_label, y_label:

  Optional plot labels.

- aspect_ratio, panel_widths, panel_heights:

  Optional panel sizing controls.

## Value

A `ggplot` object with source data and figure metadata.

## Examples

``` r
expression <- matrix(
  c(3, 5, 8, 4, 6, 9), nrow = 3,
  dimnames = list(c("G1", "G2", "G3"), c("sample_a", "sample_b"))
)
sn_plot_association(expression, x = "sample_a", y = "sample_b")
```
