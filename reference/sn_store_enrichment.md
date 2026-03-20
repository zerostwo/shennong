# Store an enrichment result on a Seurat object

This helper stores enrichment output inside
`object@misc$enrichment_results[[store_name]]` so interpretation and
writing helpers can reuse it later.

## Usage

``` r
sn_store_enrichment(
  object,
  result,
  store_name = "default",
  analysis = c("ora", "gsea"),
  database = "GOBP",
  species = NULL,
  source_de_name = NULL,
  gene_col = "gene",
  score_col = NULL,
  return_object = TRUE
)
```

## Arguments

- object:

  A `Seurat` object.

- result:

  An enrichment result object or data frame coercible with
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html).

- store_name:

  Name used under `object@misc$enrichment_results`.

- analysis:

  One of `"ora"` or `"gsea"`.

- database:

  Database used for enrichment, for example `"GOBP"`.

- species:

  Species label used in the enrichment run.

- source_de_name:

  Optional stored DE result name that produced the input ranked gene
  list or gene set.

- gene_col:

  Column containing gene symbols when the enrichment input came from a
  data frame.

- score_col:

  Column containing ranking scores for GSEA inputs.

- return_object:

  If `TRUE`, return the updated Seurat object.

## Value

A `Seurat` object or a stored-result list.

## Examples

``` r
if (requireNamespace("Seurat", quietly = TRUE)) {
  counts <- matrix(rpois(10 * 12, lambda = 1), nrow = 10, ncol = 12)
  rownames(counts) <- c("CD3D", "CD3E", "TRAC", "LTB", "MS4A1", "CD79A", "HLA-DRA", "LYZ", "ACTB", "MALAT1")
  colnames(counts) <- paste0("cell", 1:12)
  obj <- sn_initialize_seurat_object(counts, species = "human")
  enrich_tbl <- tibble::tibble(
    ID = c("GO:0001", "GO:0002"),
    Description = c("immune response", "lymphocyte activation"),
    NES = c(2.1, 1.7),
    p.adjust = c(0.01, 0.03)
  )
  obj <- sn_store_enrichment(obj, enrich_tbl, store_name = "demo_gsea")
  names(obj@misc$enrichment_results)
}
#> INFO [2026-03-20 02:56:34] Initializing Seurat object for project: Shennong
#> INFO [2026-03-20 02:56:34] Running QC metrics for human ...
#> Maps last updated on: Sat Nov 16 10:35:32 2024
#> Maps last updated on: Sat Nov 16 10:35:32 2024
#> INFO [2026-03-20 02:56:34] Seurat object initialization complete.
#> [1] "demo_gsea"
```
