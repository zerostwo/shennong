# Calculate composition proportions

Calculates the proportion of different categories within groups from
metadata.

## Usage

``` r
sn_calculate_composition(
  x,
  group_by,
  variable,
  min_cells = 20,
  measure = c("proportion", "count", "both"),
  additional_cols = NULL,
  sort_by = NULL,
  sort_value = NULL,
  sort_desc = FALSE
)
```

## Arguments

- x:

  A Seurat object or a data frame.

- group_by:

  Column name or character vector of column names used to define groups.

- variable:

  Column name whose proportions should be calculated.

- min_cells:

  Minimum number of cells required for a returned `group_by + variable`
  category to be retained. Defaults to `20`.

- measure:

  One of `"proportion"`, `"count"`, or `"both"`. Defaults to
  `"proportion"`.

- additional_cols:

  Optional character vector of additional columns to carry into the
  output.

- sort_by:

  Optional summary column used to order the primary `group_by` column in
  the returned table. One of `NULL`, `"proportion"`, `"count"`, or
  `"group_total"`. Sorting is only applied when `group_by` has length 1.

- sort_value:

  Optional level of `variable` used for sorting when `sort_by` is
  supplied. Defaults to the first observed level.

- sort_desc:

  Logical; if `TRUE`, sort in descending order.

## Value

A data frame with per-group proportions.

## Details

This function takes either a Seurat object or a data frame (like cell
metadata) and computes the percentage composition of a given `variable`
(for example, cell type) for each category specified by `group_by` (for
example, sample or condition). When `group_by` contains multiple
columns, proportions are calculated within each unique combination of
those grouping columns.

## Examples

``` r
if (FALSE) { # \dontrun{
composition_df <- sn_calculate_composition(
  x = seu,
  group_by = c("sample_id", "cell_type"),
  variable = "Phase",
  min_cells = 10,
  measure = "both"
)
composition_df <- sn_calculate_composition(
  x = seu,
  group_by = "sample_id",
  variable = "cell_type",
  min_cells = 10
)
print(composition_df)
} # }
```
