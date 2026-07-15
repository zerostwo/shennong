# Prioritize phenotype-responsive or rare cell states

Prioritize phenotype-responsive or rare cell states

## Usage

``` r
sn_prioritize_states(
  object,
  method = c("augur", "scissor", "rareq"),
  phenotype,
  sample_by = NULL,
  state_by = NULL,
  contrast = NULL,
  assay = NULL,
  layer = "data",
  features = NULL,
  max_features = 500L,
  max_cells_per_state = 500L,
  folds = 3L,
  permutations = 100L,
  reduction = "pca",
  dims = 1:20,
  bulk_expression = NULL,
  bulk_phenotype = NULL,
  family = c("binomial", "gaussian", "cox"),
  store_name = "priority",
  seed = 717L,
  backend_control = list(),
  return_object = TRUE
)
```

## Arguments

- object:

  A Seurat object.

- method:

  State-priority backend.

- phenotype:

  Cell metadata phenotype for Augur/RareQ, or a descriptive label when
  Scissor receives explicit bulk inputs.

- sample_by:

  Optional biological sample column. Strongly recommended and required
  for RareQ phenotype association.

- state_by:

  Existing state/cell-type column for Augur and Scissor.

- contrast:

  Binary phenotype labels ordered as `c(case, control)`.

- assay, layer:

  Expression source for Augur.

- features, max_features:

  Features used for separability scoring.

- max_cells_per_state:

  Maximum balanced cells per state.

- folds:

  Cross-validation folds when samples are unavailable.

- permutations:

  Number of label permutations.

- reduction, dims:

  Reduction used to construct RareQ neighbors.

- bulk_expression:

  Gene-by-bulk-sample matrix required by Scissor.

- bulk_phenotype:

  Bulk phenotype vector or survival matrix for Scissor.

- family:

  Scissor phenotype family.

- store_name:

  Stored result name.

- seed:

  Random seed.

- backend_control:

  Backend-specific options.

- return_object:

  Return the updated object instead of the result.

## Value

A Seurat object or unified state-priority result.

## Examples

``` r
if (FALSE) { # \dontrun{
object <- sn_prioritize_states(
  object, phenotype = "condition", sample_by = "sample",
  state_by = "cell_type", contrast = c("treated", "control")
)
sn_get_result(object, "state_priority", "priority")
} # }
```
