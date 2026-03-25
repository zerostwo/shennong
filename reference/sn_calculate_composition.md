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
  additional_cols = NULL
)
```

## Arguments

- x:

  A Seurat object or a data frame.

- group_by:

  Column name used to define groups.

- variable:

  Column name whose proportions should be calculated.

- min_cells:

  Minimum number of cells required for a `group_by` category to be
  retained. Defaults to `20`.

- additional_cols:

  Optional character vector of additional columns to carry into the
  output.

## Value

A data frame with per-group proportions.

## Details

This function takes either a Seurat object or a data frame (like cell
metadata) and computes the percentage composition of a given `variable`
(for example, cell type) for each category specified by `group_by` (for
example, sample or condition).

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
