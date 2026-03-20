# Retrieve a stored enrichment result from a Seurat object

Retrieve a stored enrichment result from a Seurat object

## Usage

``` r
sn_get_enrichment_result(
  object,
  enrichment_name = "default",
  top_n = NULL,
  groups = NULL,
  with_metadata = FALSE
)
```

## Arguments

- object:

  A `Seurat` object.

- enrichment_name:

  Name of the stored enrichment result.

- top_n:

  Optional number of top terms to keep.

- groups:

  Optional subset of cluster/group labels when the stored table includes
  a `Cluster` column.

- with_metadata:

  If `TRUE`, return the full stored result list instead of just the term
  table.

## Value

A tibble or stored-result list.

## Examples

``` r
if (FALSE) { # \dontrun{
terms <- sn_get_enrichment_result(seurat_obj, enrichment_name = "cluster_gsea", top_n = 10)
} # }
```
