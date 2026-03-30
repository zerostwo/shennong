# Plot grouped composition-style bar charts

This helper is designed for stacked or dodged categorical bar plots such
as cell composition, QC pass/fail summaries, doublet-class proportions,
or other grouped category tables produced by Shennong.

## Usage

``` r
sn_plot_composition(
  data,
  x,
  y = NULL,
  fill,
  facet_row = NULL,
  facet_col = NULL,
  position = c("stack", "fill", "dodge"),
  order_by = NULL,
  order_value = NULL,
  order_desc = FALSE,
  palette = "Paired",
  angle_x = 45,
  show_legend = TRUE,
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  panel_widths = NULL,
  panel_heights = NULL
)
```

## Arguments

- data:

  A data frame.

- x:

  Column mapped to the x-axis.

- y:

  Column mapped to the y-axis. If omitted, `proportion` is used when
  present, otherwise `count`.

- fill:

  Column mapped to the fill aesthetic.

- facet_row, facet_col:

  Optional columns used for faceting.

- position:

  One of `"stack"`, `"fill"`, or `"dodge"`.

- order_by:

  Optional y-column used to reorder the x-axis. Defaults to the selected
  `y` column.

- order_value:

  Optional level of `fill` used when ordering x-axis levels.

- order_desc:

  Logical; if `TRUE`, order x-axis levels in descending order.

- palette:

  Optional named or unnamed vector passed to
  [`ggplot2::scale_fill_manual()`](https://ggplot2.tidyverse.org/reference/scale_manual.html).

- angle_x:

  Rotation angle for x-axis labels.

- show_legend:

  Logical; if `FALSE`, hide the legend.

- title:

  Optional plot title.

- x_label, y_label:

  Optional axis labels.

- panel_widths, panel_heights:

  Optional panel size arguments forwarded to
  [`catplot::theme_cat()`](https://rdrr.io/pkg/catplot/man/theme_cat.html)
  when available.

## Value

A ggplot object.

## Examples

``` r
plot_data <- data.frame(
  sample = c("A", "A", "B", "B"),
  cell_type = c("T", "B", "T", "B"),
  proportion = c(60, 40, 30, 70)
)
sn_plot_composition(
  plot_data,
  x = sample,
  fill = cell_type
)

```
