# Score gene programs in cells or aggregated samples

Provides one stable interface for sparse-aware UCell, AUCell, GSVA,
ssGSEA, and dependency-free mean-expression scoring. Per-cell scores are
added to Seurat metadata; aggregated scores remain in the stored result.

## Usage

``` r
sn_score_programs(
  object,
  signatures,
  method = c("ucell", "aucell", "gsva", "ssgsea", "mean"),
  assay = NULL,
  layer = "data",
  name = NULL,
  group_by = NULL,
  species = NULL,
  min_genes = 1L,
  backend_control = list(),
  return_object = TRUE
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

- name:

  Stored-result and metadata prefix. A stable method-derived name is
  used when omitted.

- group_by:

  Optional metadata column used to average expression before GSVA/ssGSEA
  or other group-level scoring.

- species:

  Species required for bundled signature queries.

- min_genes:

  Minimum matched features required per signature.

- backend_control:

  Named backend-specific control list.

- return_object:

  Return the updated object or the unified result.

## Value

A Seurat object or unified program-scoring result.

## Examples

``` r
if (FALSE) { # \dontrun{
object <- sn_score_programs(
  object,
  signatures = list(T_cell = c("CD3D", "CD3E")),
  method = "ucell",
  name = "immune_programs"
)
sn_get_result(object, "program_scoring", "immune_programs")
} # }
```
