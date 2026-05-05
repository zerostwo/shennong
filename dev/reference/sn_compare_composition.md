# Compare sample-level composition between groups

Computes sample-level composition for a categorical variable and then
compares the resulting per-sample proportions between two groups. This
avoids treating individual cells as independent replicates and is
therefore the recommended way to estimate composition fold changes and
significance when replicate samples are available.

## Usage

``` r
sn_compare_composition(
  x,
  sample_by = NULL,
  group_by = NULL,
  variable,
  contrast,
  min_cells = 20,
  pseudocount = 0.5,
  test = c("wilcox", "none"),
  adjust_method = "BH",
  additional_cols = NULL,
  return_sample_data = FALSE,
  sample_col = NULL,
  group_col = NULL
)
```

## Arguments

- x:

  A Seurat object or a data frame containing cell-level metadata.

- sample_by:

  Column defining biological samples.

- group_by:

  Column defining the group or condition to compare between samples.

- variable:

  Column whose composition should be compared, for example cell type or
  cell-cycle phase.

- contrast:

  Character vector of length 2 giving the group levels as
  `c(case, control)`.

- min_cells:

  Minimum number of cells required per sample before that sample is
  retained for comparison. Defaults to `20`.

- pseudocount:

  Small value added to group means before computing `log2_fc`. Defaults
  to `0.5`.

- test:

  Statistical test to apply to sample-level proportions. One of
  `"wilcox"` or `"none"`. Defaults to `"wilcox"`.

- adjust_method:

  Multiple-testing correction method passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Defaults
  to `"BH"`.

- additional_cols:

  Optional sample-level columns to carry into the intermediate
  composition table before comparison.

- return_sample_data:

  Logical; if `TRUE`, return both the summary table and the completed
  sample-level composition table.

- sample_col:

  Deprecated alias for `sample_by`.

- group_col:

  Deprecated alias for `group_by`.

## Value

A data frame with one row per `variable` level. When
`return_sample_data = TRUE`, a list with `summary` and `sample_data` is
returned.

## Examples

``` r
if (FALSE) { # \dontrun{
comparison <- sn_compare_composition(
  seu,
  sample_by = "sample",
  group_by = "condition",
  variable = "cell_type",
  contrast = c("treated", "control")
)
} # }
```
