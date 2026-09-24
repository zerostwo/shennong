# Delete a registered workflow artifact

Registered artifact collections (for example `"clustering_stage_cache"`,
`"integration_comparison"`, or `"label_transfer"`) are removed either
member-by-member or as a whole container. Because these historical
collections live at top level in `object@misc`, their names alone cannot
prove Shennong ownership; every deletion therefore requires explicit
confirmation. Unknown artifact types fail closed; unified analysis
results must be deleted through
[`sn_delete_result`](https://zerostwo.github.io/shennong/dev/reference/sn_delete_result.md)
instead.

## Usage

``` r
sn_delete_artifact(object, artifact_type, artifact_id = NULL, confirm = FALSE)
```

## Arguments

- object:

  A Seurat object.

- artifact_type:

  Registered artifact type, i.e. the collection name or its
  `"*_artifact"` alias (for example `"integration"` or
  `"integration_artifact"`).

- artifact_id:

  Optional artifact identifier. When omitted, the entire artifact
  collection container is removed only when `confirm = TRUE`.

- confirm:

  Explicit confirmation required before deleting a member or an entire
  top-level artifact collection. This protects user-owned payloads that
  use a legacy collection name also registered by Shennong.

## Value

The modified Seurat object.

## Examples

``` r
if (FALSE) { # \dontrun{
obj <- sn_delete_artifact(
  obj, "integration_comparison", "pbmc_grid", confirm = TRUE
)
obj <- sn_delete_artifact(obj, "label_transfer", confirm = TRUE)
} # }
```
