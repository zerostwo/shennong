# Plot a Sankey diagram with a publication-style preset

Dedicated composition entry points use ordinary strings for column
names. They share data preparation, counting, and figure metadata with
[`sn_plot_composition()`](https://zerostwo.github.io/shennong/dev/reference/sn_plot_composition.md).
They do not duplicate statistical implementations. Text defaults to 8
pt. For style C, `show_pies = FALSE` hides bottom pies without changing
the ribbons; TRUE shows them (the default for C).

## Usage

``` r
sn_plot_sankey(
  data,
  flow_by,
  fill_by = NULL,
  style = c("A", "B", "C"),
  panel_widths = NULL,
  panel_heights = NULL,
  show_pies = NULL,
  ...
)
```

## Arguments

- data:

  A Seurat object or metadata/summary data frame.

- flow_by:

  Character vector of metadata columns in stage order.

- fill_by:

  Optional color-group column name string. Defaults to the first stage
  for A/C, and the second stage for B. C requires the first stage.

- style:

  `"A"` hierarchy, `"B"` bracketed annotation comparison, or `"C"`
  vertical source/target composition with pies. B/C require two stages.

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

- show_pies:

  Show source-composition pies below style C target nodes; NULL enables
  them only for C. All pies have equal size, with slices computed from
  the same retained weights as the flows within each facet.

- ...:

  Additional named arguments to
  [`sn_plot_composition()`](https://zerostwo.github.io/shennong/dev/reference/sn_plot_composition.md),
  including `sankey_layout`, `stratum_gap`, `show_stratum_boxes`,
  `stratum_label_size`, `show_pies`, `node_palette`, `axis_labels`,
  `facet_row_by`, `facet_col_by`, `y_by` for summarized weights, and
  `data_kind`.

## Value

A ggplot with Shennong figure metadata and source data attached.

## Examples

``` r
d <- data.frame(original = c("T", "T", "NK"), final = c("T", "NK", "NK"))
sn_plot_sankey(d, flow_by = c("original", "final"), style = "B")
```
