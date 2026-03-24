# List stored Shennong analysis and interpretation results on a Seurat object

List stored Shennong analysis and interpretation results on a Seurat
object

## Usage

``` r
sn_list_results(object)
```

## Arguments

- object:

  A `Seurat` object.

## Value

A tibble describing stored DE, enrichment, and interpretation results.

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
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- sn_find_de(
    obj,
    analysis = "markers",
    group_by = NULL,
    layer = "data",
    min_pct = 0,
    logfc_threshold = 0,
    return_object = TRUE,
    verbose = FALSE
  )
  sn_list_results(obj)
}
#> INFO [2026-03-24 17:45:58] Initializing Seurat object for project: Shennong
#> INFO [2026-03-24 17:45:58] Running QC metrics for human ...
#> INFO [2026-03-24 17:45:58] Seurat object initialization complete.
#> Warning: No DE genes identified
#> Warning: The following tests were not performed: 
#> Warning: When testing Shennong versus all:
#>  Cells in one or both identity groups are not present in the data requested
#> # A tibble: 1 × 8
#>   collection type  name    analysis method created_at              n_rows source
#>   <chr>      <chr> <chr>   <chr>    <chr>  <chr>                    <int> <chr> 
#> 1 de_results de    default markers  wilcox 2026-03-24 17:45:59 UTC      0 NA    
```
