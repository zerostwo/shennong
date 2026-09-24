# Run phenotype-guided Scissor cell selection

`sn_run_scissor()` is the direct Scissor entry point. It returns a
cell-first unified Scissor result while sharing the same backend
implementation used by
[`sn_prioritize_states()`](https://zerostwo.github.io/shennong/dev/reference/sn_prioritize_states.md).
The result also carries state and sample summaries, correlation
summaries, model metadata, and optional reliability output.

## Usage

``` r
sn_run_scissor(
  object,
  bulk_expression,
  bulk_phenotype,
  state_by = NULL,
  sample_by = NULL,
  family = c("binomial", "gaussian", "cox"),
  phenotype = "bulk_phenotype",
  assay = NULL,
  layer = "data",
  result_id = "scissor",
  seed = 717L,
  backend_control = list(),
  return_object = TRUE
)
```

## Arguments

- object:

  A Seurat object.

- bulk_expression:

  Numeric gene-by-bulk-sample expression matrix with unique gene row
  names and sample column names.

- bulk_phenotype:

  A binary or numeric phenotype vector, or a two-column time/event
  matrix for `family = "cox"`. Named inputs are aligned strictly to
  `bulk_expression`; unnamed inputs retain positional compatibility.

- state_by:

  Optional cell-state metadata column. When omitted, the direct result
  remains cell-first and uses one descriptive `all_cells` state summary.

- sample_by:

  Optional biological-sample metadata column.

- family:

  Scissor response family.

- phenotype:

  Descriptive label stored with the result.

- assay, layer:

  Single-cell expression source.

- result_id:

  Stored result name under `scissor`.

- seed:

  Random seed recorded in provenance.

- backend_control:

  Direct Scissor controls, or a list containing a `scissor` sub-list.
  Set `reliability = TRUE` to run bootstrap reliability; set
  `retain_correlation_matrix = TRUE` to retain the long sample-cell
  correlation table subject to `max_correlation_rows`. Scissor's
  `cutoff` is an alpha-search stopping criterion rather than a
  guaranteed upper bound; Shennong warns and marks the result
  exploratory when no candidate alpha achieves it.

- return_object:

  Return the updated object instead of the result.

## Value

A Seurat object or unified Scissor state-priority result.

## Examples

``` r
if (FALSE) { # \dontrun{
result <- sn_run_scissor(
  object,
  bulk_expression = bulk_matrix,
  bulk_phenotype = bulk_group,
  state_by = "cell_type",
  sample_by = "sample",
  return_object = FALSE
)
result$tables$cells
result$tables$states
} # }
```
