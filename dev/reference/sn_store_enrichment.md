# Store an enrichment result on a Seurat object

This helper stores enrichment output inside the canonical Shennong
result registry so interpretation and writing helpers can reuse it
later.

## Usage

``` r
sn_store_enrichment(
  object,
  result,
  result_id = "default",
  analysis = c("ora", "gsea"),
  database = "GOBP",
  species = NULL,
  source_de_result_id = NULL,
  gene_col = "gene",
  score_col = NULL,
  parameters = list(),
  return_object = TRUE
)
```

## Arguments

- object:

  A `Seurat` object.

- result:

  An enrichment result object or data frame coercible with
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html).

- result_id:

  Stable identifier for the stored enrichment result.

- analysis:

  One of `"ora"` or `"gsea"`.

- database:

  Database used for enrichment, for example `"GOBP"`.

- species:

  Species label used in the enrichment run.

- source_de_result_id:

  Optional stored DE result name that produced the input ranked gene
  list or gene set.

- gene_col:

  Column containing gene symbols when the enrichment input came from a
  data frame.

- score_col:

  Column containing ranking scores for GSEA inputs.

- parameters:

  Named list of effective enrichment parameters retained for discovery
  and reproducibility.

- return_object:

  If `TRUE`, return the updated Seurat object.

## Value

A `Seurat` object or a stored-result list.

## Examples

``` r
if (requireNamespace("Seurat", quietly = TRUE)) {
  counts <- matrix(rpois(10 * 12, lambda = 1), nrow = 10, ncol = 12)
  rownames(counts) <- c(
    "CD3D", "CD3E", "TRAC", "LTB", "MS4A1",
    "CD79A", "HLA-DRA", "LYZ", "ACTB", "MALAT1"
  )
  colnames(counts) <- paste0("cell", 1:12)
  obj <- sn_initialize_seurat_object(counts, species = "human")
  enrich_tbl <- tibble::tibble(
    ID = c("GO:0001", "GO:0002"),
    Description = c("immune response", "lymphocyte activation"),
    NES = c(2.1, 1.7),
    p.adjust = c(0.01, 0.03)
  )
  obj <- sn_store_enrichment(obj, enrich_tbl, result_id = "demo_gsea")
  sn_list_results(obj, type = "enrichment")
  sn_get_result(obj, type = "enrichment", result_id = "demo_gsea")
}
#> INFO [2026-09-28 08:50:02] Initializing Seurat object for project: Shennong.
#> INFO [2026-09-28 08:50:02] Running QC metrics for human.
#> INFO [2026-09-28 08:50:02] Seurat object initialization complete.
#> $schema_version
#> [1] "2.0.0"
#> 
#> $method
#> [1] "ora"
#> 
#> $backend
#> [1] "clusterProfiler"
#> 
#> $models
#> $models$backend_results
#> $models$backend_results$GOBP
#> # A tibble: 2 × 4
#>   ID      Description             NES p.adjust
#>   <chr>   <chr>                 <dbl>    <dbl>
#> 1 GO:0001 immune response         2.1     0.01
#> 2 GO:0002 lymphocyte activation   1.7     0.03
#> 
#> 
#> 
#> $analysis
#> [1] "ora"
#> 
#> $database
#> [1] "GOBP"
#> 
#> $species
#> NULL
#> 
#> $source_de_result_id
#> NULL
#> 
#> $gene_col
#> [1] "gene"
#> 
#> $score_col
#> NULL
#> 
#> $parameters
#> list()
#> 
#> $provenance
#> $provenance$package_versions
#> $provenance$package_versions$Shennong
#> [1] "0.3.0.9000"
#> 
#> $provenance$package_versions$R
#> [1] "4.6.1"
#> 
#> 
#> $provenance$random_seed
#> [1] NA
#> 
#> $provenance$timestamp
#> [1] "2026-09-28 08:50:02 UTC"
#> 
#> $provenance$result_id
#> [1] "demo_gsea"
#> 
#> $provenance$analysis_type
#> [1] "enrichment"
#> 
#> 
#> $analysis_type
#> [1] "enrichment"
#> 
#> $result_id
#> [1] "demo_gsea"
#> 
#> $input
#> list()
#> 
#> $tables
#> $tables$primary
#> # A tibble: 2 × 4
#>   ID      Description             NES p.adjust
#>   <chr>   <chr>                 <dbl>    <dbl>
#> 1 GO:0001 immune response         2.1     0.01
#> 2 GO:0002 lymphocyte activation   1.7     0.03
#> 
#> 
#> $embeddings
#> list()
#> 
#> $graphs
#> list()
#> 
#> $diagnostics
#> list()
#> 
#> $warnings
#> character(0)
#> 
```
