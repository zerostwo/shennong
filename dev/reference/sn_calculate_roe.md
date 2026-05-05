# Calculate observed-over-expected enrichment

Calculates observed-over-expected (RO/E) enrichment for a categorical
`variable` across metadata groups. RO/E is useful for asking whether a
cell type, cluster, state, or annotation is over-represented in a
condition relative to the marginal distribution of the full table.

## Usage

``` r
sn_calculate_roe(
  x,
  group_by,
  variable,
  pseudocount = 0,
  return_matrix = FALSE,
  matrix_value = c("roe", "log2_roe", "observed", "expected")
)
```

## Arguments

- x:

  A Seurat object or a data frame.

- group_by:

  Column name or character vector of column names used to define the
  rows/groups of the contingency table.

- variable:

  Column name defining the categories to test for enrichment.

- pseudocount:

  Numeric scalar added to both observed and expected counts before
  computing the ratio. Defaults to `0`.

- return_matrix:

  Logical; if `TRUE`, return a numeric matrix with `group_by` levels as
  rows and `variable` levels as columns.

- matrix_value:

  Value to place in the matrix when `return_matrix = TRUE`. One of
  `"roe"`, `"log2_roe"`, `"observed"`, or `"expected"`.

## Value

A data frame with `observed`, `expected`, totals, `roe`, and `log2_roe`
columns. When `return_matrix = TRUE`, a numeric matrix is returned.

## Examples

``` r
if (FALSE) { # \dontrun{
roe_tbl <- sn_calculate_roe(
  seu,
  group_by = "condition",
  variable = "cell_type"
)
roe_mat <- sn_calculate_roe(
  seu,
  group_by = "condition",
  variable = "cell_type",
  return_matrix = TRUE
)
} # }
```
