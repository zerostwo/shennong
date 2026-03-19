# Calculate Composition Proportions

Calculates the proportion of different categories within groups from
metadata.

## Usage

``` r
sn_calculate_composition(
  x,
  group_by,
  variable,
  min_cells = 20,
  additional_cols = NULL
)
```

## Arguments

- x:

  A Seurat object or a data frame (e.g., metadata).

- group_by:

  A character string specifying the column name in the metadata to group
  by (e.g., "sample_id", "treatment").

- variable:

  A character string specifying the column name in the metadata for
  which to calculate proportions (e.g., "cell_type", "cluster_id").

- min_cells:

  An integer specifying the minimum number of cells required for a group
  (defined by `group_by`) to be included in the analysis. Defaults
  to 20. Set to 0 to include all groups.

- additional_cols:

  A character vector specifying the names of additional metadata columns
  to include in the output table. These columns should ideally have
  values that are constant within each `group_by` category. If values
  are not constant, the function will use the first value encountered
  for each group and issue a warning. Defaults to NULL.

## Value

A data frame with the following columns:

- group_by:

  The grouping variable categories.

- variable:

  The variable categories whose proportions are calculated.

- proportion:

  The calculated proportion (percentage) for each variable category
  within each group.

- ...:

  Additional columns specified by `additional_cols`.

## Details

This function takes either a Seurat object or a data frame (like cell
metadata) and computes the percentage composition of a given `variable`
(e.g., cell type) for each category specified by `group_by` (e.g.,
sample, condition). It allows filtering out groups with fewer than a
minimum number of cells and optionally adding other metadata columns to
the results.

## Examples

``` r
if (FALSE) { # \dontrun{
composition_df <- sn_calculate_composition(
  x = seu,
  group_by = "sample_id",
  variable = "cell_type",
  min_cells = 10
)
print(composition_df)
} # }
```
