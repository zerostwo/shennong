# Plot sample-level composition summaries with string column names

Plot sample-level composition summaries with string column names

## Usage

``` r
sn_plot_sample_bar(
  data,
  x_by,
  fill_by,
  sample_by,
  panel_widths = NULL,
  panel_heights = NULL,
  ...
)
```

## Arguments

- data:

  A Seurat object or metadata/summary data frame.

- x_by:

  Column name defining the x axis.

- fill_by:

  Column name string for the composition category.

- sample_by:

  Column name identifying biological samples.

- panel_widths, panel_heights:

  Optional positive numeric panel dimensions in points (pt), applied
  directly through ggplot2 for every composition type. Scalars repeat
  across facets; vectors specify facet column widths or row heights.
  Style C accepts scalars only to preserve circular pies. If only one
  dimension and `aspect_ratio` are supplied, the other is derived; C
  derives it from its coordinate ratio even when `aspect_ratio` is
  omitted. Explicit width and height take precedence over
  `aspect_ratio`. These sizes exclude outer labels/margins; choose a
  sufficiently large export canvas.

- ...:

  Additional named arguments to
  [`sn_plot_composition()`](https://zerostwo.github.io/shennong/dev/reference/sn_plot_composition.md),
  including `sankey_layout`, `stratum_gap`, `show_stratum_boxes`,
  `stratum_label_size`, `show_pies`, `node_palette`, `axis_labels`,
  `facet_row_by`, `facet_col_by`, `y_by` for summarized weights, and
  `data_kind`.

## Value

A ggplot with Shennong figure metadata.
