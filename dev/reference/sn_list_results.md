# List stored Shennong analysis and interpretation results on a Seurat object

List stored Shennong analysis and interpretation results on a Seurat
object

## Usage

``` r
sn_list_results(object, type = NULL, include_artifacts = FALSE)
```

## Arguments

- object:

  A `Seurat` object.

- type:

  Optional analysis type used to filter the result inventory.

- include_artifacts:

  Include registered workflow artifacts that do not implement the
  unified analysis-result contract.

## Value

A tibble describing registered Shennong stored-result collections,
including DE, enrichment, interpretation, deconvolution, Milo,
communication, regulatory activity, and QC assessment entries when
present. The canonical lookup key is reported in `result_id`.

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
```
