# Test differential abundance across biological samples

Provides one stable entry point for sample-level Propeller or
permutation tests and neighborhood-level Milo testing. Cells are never
treated as independent replicates.

## Usage

``` r
sn_test_abundance(
  object,
  method = c("propeller", "milo", "sccoda", "permutation"),
  sample_by,
  condition_by,
  cell_type_by,
  design = NULL,
  contrast = NULL,
  store_name = "abundance",
  transform = c("logit", "asin"),
  permutations = 1000L,
  seed = 717L,
  backend_control = list(),
  return_object = TRUE
)
```

## Arguments

- object:

  A Seurat object.

- method:

  Abundance backend.

- sample_by:

  Biological sample column.

- condition_by:

  Sample-level condition column.

- cell_type_by:

  Cell-type/state column.

- design:

  Optional no-intercept formula or sample-level covariate names for
  Propeller; covariate names are forwarded to Milo.

- contrast:

  Two condition labels ordered as `c(case, control)`.

- store_name:

  Stored differential-abundance result name.

- transform:

  Propeller proportion transformation.

- permutations:

  Number of sample-label permutations.

- seed:

  Random seed for permutation testing.

- backend_control:

  Named backend argument lists.

- return_object:

  Return the updated object instead of the result.

## Value

A Seurat object or unified differential-abundance result.

## Examples

``` r
if (FALSE) { # \dontrun{
object <- sn_test_abundance(
  object, sample_by = "sample", condition_by = "condition",
  cell_type_by = "cell_type", contrast = c("treated", "control")
)
sn_get_result(object, "differential_abundance", "abundance")
} # }
```
