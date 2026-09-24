# Retrieve an enrichment result table

Retrieve an enrichment result table

## Usage

``` r
sn_get_enrichment_result(
  object,
  result_id = NULL,
  top_n = NULL,
  groups = NULL,
  top_scope = c("group", "all"),
  p_adjusted_cutoff = NULL
)
```

## Arguments

- object:

  A `Seurat` object or unified enrichment result.

- result_id:

  Identifier of the stored enrichment result. May be omitted when
  exactly one enrichment result is stored.

- top_n:

  Optional number of top terms to keep.

- groups:

  Optional subset of cluster/group labels when the stored table includes
  a `Cluster` column.

- top_scope:

  Whether to retain the top terms within each group or overall.

- p_adjusted_cutoff:

  Optional maximum adjusted p-value, from zero to one.

## Value

A filtered tibble. Use `sn_get_result(object, "enrichment")` for the
complete result and its metadata.

## Examples

``` r
if (FALSE) { # \dontrun{
terms <- sn_get_enrichment_result(seurat_obj, result_id = "cluster_gsea", top_n = 10)
} # }
```
