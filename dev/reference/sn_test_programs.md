# Test program activity between conditions

Aggregates cell-level scores to the sample level before inference when
`sample_by` is supplied. This prevents cells from being treated as
independent biological replicates. Complete pairs sharing the same
sample ID across conditions use a paired Wilcoxon test or a
sample-adjusted limma model. A mixture of paired and unpaired samples is
rejected. The result table records `paired`.

## Usage

``` r
sn_test_programs(
  object,
  source_result_id,
  condition_by,
  sample_by = NULL,
  group_by = NULL,
  contrast = NULL,
  method = c("wilcox", "limma"),
  result_id = NULL,
  return_object = TRUE
)
```

## Arguments

- object:

  A Seurat object containing a stored program-scoring result.

- source_result_id:

  Stored scoring result identifier.

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

- result_id:

  Stable identifier for the stored test result.

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
