# Plot a supplied quantitative matrix or bubble matrix

A presentation-only entry point for Ro/e, correlation, module activity,
pathway scores and other already computed matrices. Values are never
scaled, normalized or tested. Missing values remain missing, not zero.
Clustering is opt-in and uses Euclidean distance with complete linkage
on finite matrices.

## Usage

``` r
sn_plot_matrix(
  x,
  type = c("heatmap", "bubble"),
  size = NULL,
  style = "classic",
  palette = NULL,
  midpoint = NULL,
  limits = NULL,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  row_split = NULL,
  column_split = NULL,
  show_values = FALSE,
  digits = 2L,
  max_size = 4,
  value_label = "Value",
  size_label = "Size",
  na_color = "grey90",
  title = NULL,
  panel_widths = NULL,
  panel_heights = NULL,
  compare = NULL,
  row_annotation = NULL,
  column_annotation = NULL,
  row_summary = NULL,
  column_summary = NULL,
  summary_label = "Supplied summary"
)
```

## Arguments

- x:

  Numeric matrix with unique row and column names.

- type:

  `heatmap` or `bubble`.

- size:

  Numeric non-negative matrix for bubble area, with identical dimensions
  and dimnames to `x`. Zero has zero area. Required for bubbles.

- style:

  Presentation preset; see
  [`sn_list_plot_styles()`](https://zerostwo.github.io/shennong/dev/reference/sn_list_plot_styles.md).

- palette:

  Optional palette name or explicit colors. Defaults to a sequential
  palette, or a diverging palette when `midpoint` is supplied.

- midpoint:

  Optional numeric midpoint for a diverging color scale, for example
  zero for correlations or one for Ro/e. No transformation is applied.

- limits:

  Optional two-number increasing color limits. Out-of-range values are
  squished for display; source values remain attached unchanged.

- cluster_rows, cluster_cols:

  Whether to cluster rows or columns.

- row_split, column_split:

  Optional named vectors assigning rows/columns to facet blocks. Names
  must exactly match the corresponding matrix names.

- show_values:

  Whether to print the supplied values in tiles.

- digits:

  Number of significant digits in tile labels.

- max_size:

  Maximum bubble diameter in mm.

- value_label, size_label:

  Legend titles, including units where relevant.

- na_color:

  Fill used for missing values.

- title:

  Optional title.

- panel_widths, panel_heights:

  Optional panel sizes in points.

- compare:

  Second matrix with identical axes for `split_triangle`: `x` fills the
  upper triangle and `compare` the lower; diagonal omitted.

- row_annotation, column_annotation:

  Data frames with row names exactly matching the corresponding axis;
  columns are categorical annotation tracks.

- row_summary:

  Named numeric vector rendered beside `heatmap_bar`.

- column_summary:

  Named numeric vector drawn above `heatmap_density`. Values are
  supplied summaries, not inferred densities or expression means.

- summary_label:

  Axis label for the supplied summary track.

## Value

A ggplot or patchwork with unchanged input tables attached. Structural
styles are function-specific; see
[`sn_list_plot_styles()`](https://zerostwo.github.io/shennong/dev/reference/sn_list_plot_styles.md).
`dendrogram` clusters both axes using Euclidean distance and complete
linkage. Triangle styles require identical row/column names and forbid
independent clustering. `circular` places columns around a ring and rows
from inner to outer radius. Annotation tracks are categorical, with
their own legends.

## Examples

``` r
m <- matrix(c(0.5, 1.5, 2, 0.8), 2,
            dimnames = list(c("NK", "B"), c("Lung", "Skin")))
sn_plot_matrix(m, midpoint = 1, style = "atlas_module", value_label = "Ro/e")
```
