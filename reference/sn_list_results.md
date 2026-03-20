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
  rownames(counts) <- c("CD3D", "CD3E", "TRAC", "LTB", "MS4A1", "CD79A", "HLA-DRA", "LYZ", "ACTB", "MALAT1")
  colnames(counts) <- paste0("cell", 1:12)
  obj <- sn_initialize_seurat_object(counts, species = "human")
  obj <- sn_find_de(obj, analysis = "markers", group_by = NULL, layer = "data", min_pct = 0, logfc_threshold = 0, return_object = TRUE, verbose = FALSE)
  sn_list_results(obj)
}
#> INFO [2026-03-20 02:56:23] Initializing Seurat object for project: Shennong
#> INFO [2026-03-20 02:56:23] Running QC metrics for human ...
#> Maps last updated on: Sat Nov 16 10:35:32 2024
#> Maps last updated on: Sat Nov 16 10:35:32 2024
#> INFO [2026-03-20 02:56:23] Seurat object initialization complete.
#> Error in .sn_validate_seurat_assay_layer(object = object, assay = assay,     layer = source_layer): Layer 'data' was not found in assay 'RNA'.
```
