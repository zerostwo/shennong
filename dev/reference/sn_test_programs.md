# Test program activity between conditions

Aggregates cell-level scores to the sample level before inference when
`sample_by` is supplied. This prevents cells from being treated as
independent biological replicates.

## Usage

``` r
sn_test_programs(
  object,
  score_name,
  condition_by,
  sample_by = NULL,
  group_by = NULL,
  contrast = NULL,
  method = c("wilcox", "limma"),
  store_name = NULL,
  return_object = TRUE
)
```

## Arguments

- object:

  A Seurat object containing a stored program-scoring result.

- score_name:

  Stored scoring result name.

- condition_by:

  Condition metadata column.

- sample_by:

  Sample/patient metadata column. Strongly recommended for inference.

- group_by:

  Optional cell-type or state column used for stratified tests.

- contrast:

  Optional two condition levels; defaults to the first two.

- method:

  `"wilcox"` or `"limma"`.

- store_name:

  Result name.

- return_object:

  Return the object or unified result.

## Value

A Seurat object or unified program-comparison result.

## Examples

``` r
if (FALSE) { # \dontrun{
object <- sn_test_programs(
  object, "immune_programs", condition_by = "condition",
  sample_by = "patient", group_by = "cell_type"
)
} # }
```
