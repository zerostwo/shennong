# Score gene programs in cells or aggregated samples

Provides one stable interface for sparse-aware UCell, AUCell, GSVA,
ssGSEA, and dependency-free mean-expression scoring. Per-cell scores are
added to Seurat metadata; aggregated scores remain in the stored result.
Metadata names are unique across programs and existing columns. The
stored `tables$metadata_columns` maps each program to its metadata
column.

## Usage

``` r
sn_score_programs(
  object,
  signatures,
  method = c("ucell", "aucell", "gsva", "ssgsea", "mean"),
  assay = NULL,
  layer = "data",
  group_by = NULL,
  aggregate = NULL,
  species = NULL,
  min_genes = 1L,
  backend_control = list(),
  seed = 717,
  return_object = TRUE,
  result_id = NULL,
  overwrite = FALSE
)
```

## Arguments

- object:

  A `Seurat` object.

- signatures:

  Named list, program/gene data frame, named gene vector, or bundled
  signature query vector.

- method:

  Scoring backend. UCell is the default for per-cell data.

- assay, layer:

  Expression source.

- group_by:

  Optional metadata column identifying groups to aggregate.

- aggregate:

  Required with `group_by`: `"expression"` averages expression before
  scoring, while `"scores"` averages per-cell scores after scoring.
  These differ for nonlinear scoring methods.

- species:

  Species required for bundled signature queries.

- min_genes:

  Minimum matched features required per signature.

- backend_control:

  Named backend-specific control list, keyed by method.

- seed:

  Random seed applied locally during scoring. Use `NULL` to use the
  caller's random state. Effective seeds are recorded in provenance.

- return_object:

  Return the updated object or the unified result.

- result_id:

  Stable identifier for the stored program-scoring result and its
  metadata prefix. When omitted, a unique method-derived ID is
  allocated.

- overwrite:

  Explicitly replace an existing result. Requires an explicit
  `result_id`; obsolete score columns owned by that result are removed.

## Value

A Seurat object or unified program-scoring result.

## Examples

``` r
if (FALSE) { # \dontrun{
object <- sn_score_programs(
  object,
  signatures = list(T_cell = c("CD3D", "CD3E")),
  method = "ucell",
  result_id = "immune_programs"
)
sn_get_result(object, "program_scoring", "immune_programs")
} # }
```
