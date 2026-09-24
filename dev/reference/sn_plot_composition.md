# Plot metadata-driven composition summaries

`sn_plot_composition()` accepts a Seurat object, cell-level metadata, or
an already summarized composition table. It can draw stacked counts or
proportions, sample-level replicate summaries, donor/sample
distributions, and two- or multi-stage alluvial (Sankey-style) diagrams.
For Seurat input, all variables are read from `object[[]]` and
expression matrices are not accessed.

## Usage

``` r
sn_plot_composition(
  data,
  x = NULL,
  y = NULL,
  fill = NULL,
  type = "bar",
  data_kind = c("auto", "metadata", "summary"),
  sample_by = NULL,
  unit_by = NULL,
  flow_by = NULL,
  min_cells = 0,
  position = c("stack", "fill", "dodge"),
  summary_fun = c("mean", "median"),
  errorbar = c("se", "sd", "none"),
  show_points = TRUE,
  point_alpha = 0.7,
  jitter_width = 0.15,
  bins = 30,
  alluvium_alpha = 0.8,
  show_stratum_labels = TRUE,
  facet_row = NULL,
  facet_col = NULL,
  order_by = NULL,
  order_value = NULL,
  order_desc = FALSE,
  palette = "Paired",
  angle_x = 45,
  show_legend = TRUE,
  title = NULL,
  x_label = NULL,
  y_label = NULL,
  aspect_ratio = NULL,
  panel_widths = NULL,
  panel_heights = NULL,
  sankey_layout = c("parallel", "expanded"),
  stratum_gap = 0.02,
  show_stratum_boxes = FALSE,
  stratum_label_size = 8/ggplot2::.pt,
  x_by = NULL,
  y_by = NULL,
  fill_by = NULL,
  facet_row_by = NULL,
  facet_col_by = NULL,
  style = c("A", "B", "C"),
  show_pies = NULL,
  node_palette = "Set3",
  axis_labels = NULL
)
```

## Arguments

- data:

  A Seurat object or data frame. A data frame may contain cell-level
  metadata or a table with a `proportion` or `count` column.

- x:

  Metadata/table column mapped to the x-axis. It may be omitted for an
  alluvial plot when `flow_by` is supplied.

- y:

  Optional value column for summarized input. For metadata input it must
  be `proportion` or `count`; the default is `proportion` for filled
  bars and sample-level plots, and `count` otherwise.

- fill:

  Optional composition category or fill column. It is required for
  stacked and sample-level composition plots. For alluvial plots it
  defaults to the first `flow_by` column for A/C and the last column for
  B.

- type:

  Plot type: `"bar"`, `"sample_bar"`, `"sample_boxplot"`, `"histogram"`,
  or `"alluvial"`. The aliases `"stacked_bar"`, `"errorbar"`, and
  `"sankey"` are also accepted.

- data_kind:

  One of `"auto"`, `"metadata"`, or `"summary"`. In auto mode, a data
  frame containing `proportion`/`count` (or an explicit `y`) is treated
  as summarized; Seurat input is always metadata.

- sample_by:

  Metadata column defining biological samples for `sample_bar` and
  `sample_boxplot`.

- unit_by:

  Optional unique unit for `bar`, `histogram`, or `alluvial` plots, for
  example a donor ID repeated across cells. Defaults to `sample_by`
  outside sample-level plot types.

- flow_by:

  Character vector of two or more metadata columns forming the alluvial
  axes, for example `c("level1", "level2", "level3")`.

- min_cells:

  Minimum total number of cells required per biological sample for
  sample-level plots. For an ordinary metadata bar, it is the minimum
  count required for a returned x/fill combination.

- position:

  Bar position: `"stack"`, `"fill"`, or `"dodge"`.

- summary_fun:

  Center shown by `sample_bar`: `"mean"` or `"median"`.

- errorbar:

  Error interval for `sample_bar`: `"none"`, `"sd"`, or `"se"`.

- show_points:

  Logical; overlay sample-level observations on a summary bar or
  boxplot.

- point_alpha, jitter_width:

  Point transparency and horizontal jitter.

- bins:

  Number of bins used by `histogram`.

- alluvium_alpha:

  Transparency of alluvial flows (default 0.8 for A/B, 0.35 for C).

- show_stratum_labels:

  Logical; label alluvial strata.

- facet_row, facet_col:

  Optional metadata/table columns used for faceting.

- order_by:

  Optional value column used to reorder the x-axis.

- order_value:

  Optional `fill` level used when ordering the x-axis.

- order_desc:

  Logical; order x levels in descending order.

- palette:

  Named Shennong palette or an explicit color vector.

- angle_x:

  Rotation angle for x-axis labels.

- show_legend:

  Logical; show the legend.

- title, x_label, y_label:

  Optional plot labels.

- aspect_ratio:

  Optional panel aspect ratio.

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

- sankey_layout:

  Sankey spacing: `"parallel"` (default) gives all axes the same total
  gap budget (single-node axes have no gaps); `"expanded"` uses the same
  gap between nodes, so axes with more categories grow taller. Ribbon
  thickness always retains the original count/weight scale; expansion
  adds whitespace only. Style C scales each facet to a common horizontal
  span; compare proportions within facets, not absolute counts between
  facets.

- stratum_gap:

  Gap as a fraction of total flow weight (default 0.02 for A/C, 0.008
  for B).

- show_stratum_boxes:

  Draw outlined Sankey node boxes (default FALSE). Boxes are white for
  A/B; C retains colored bars and adds outlines. Labels sit outside the
  flows in either mode.

- stratum_label_size:

  Sankey label size in mm for compatibility with ggplot2 text layers.
  The default `8 / ggplot2::.pt` renders at 8 pt. Composition titles,
  axes, facet strips, and legends also default to 8 pt.

- x_by, y_by, fill_by:

  Character scalar column names for the x axis, summarized values, and
  color grouping. These are the preferred interface; string variables
  are supported. Do not combine with the corresponding legacy `x`, `y`,
  or `fill` argument.

- facet_row_by, facet_col_by:

  Character scalar column names for facets.

- style:

  Sankey preset: `"A"` separated annotation hierarchy, `"B"` two-axis
  horizontal annotation comparison with brackets and destination-colored
  flows, or `"C"` two-axis top-to-bottom composition with colored node
  bars and source-composition pies. B/C require two axes.

- show_pies:

  Show source-composition pies below style C target nodes; NULL enables
  them only for C. All pies have equal size, with slices computed from
  the same retained weights as the flows within each facet.

- node_palette:

  Palette for style C target nodes, separate from the source palette
  used for flows and pies.

- axis_labels:

  Optional character vector of axis titles for styles A/B. Style C
  labels source categories directly; use `title` for its heading.

## Value

A ggplot object with its plot-ready source table attached to the
Shennong figure specification.

## Details

Sample-level plots first calculate each biological sample's composition
and only then summarize across conditions. This keeps samples, rather
than cells, as the visual replicate and inserts zeroes for cell types
absent from an otherwise retained sample.

Sankey axes respect each input column's factor levels independently,
from top to bottom; character columns use first appearance order. Unused
levels are omitted. Missing paths are excluded. Sankey plots use a clean
theme and separated stage panels, and need no optional plotting backend.

## Examples

``` r
metadata <- data.frame(
  sample = rep(paste0("S", 1:4), each = 20),
  condition = rep(c("Control", "Treated"), each = 40),
  cell_type = rep(c("T", "B"), 40)
)
sn_plot_composition(metadata, x = sample, fill = cell_type)
sn_plot_composition(
  metadata, x = condition, fill = cell_type,
  type = "sample_boxplot", sample_by = "sample"
)
if (FALSE) { # \dontrun{
sn_plot_composition(
  seu, type = "alluvial",
  flow_by = c("cell_type_level1", "cell_type_level2", "cell_type_level3")
)
} # }
```
