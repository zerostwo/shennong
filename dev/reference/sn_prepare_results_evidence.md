# Prepare manuscript-style results evidence

Prepare manuscript-style results evidence

## Usage

``` r
sn_prepare_results_evidence(
  object,
  cluster_de_name = NULL,
  contrast_de_name = NULL,
  enrichment_name = NULL,
  cluster_by = NULL,
  n_markers = 5,
  n_terms = 10,
  cluster_col = NULL
)
```

## Arguments

- object:

  A `Seurat` object.

- cluster_de_name:

  Optional stored cluster-marker result.

- contrast_de_name:

  Optional stored contrast or pseudobulk result.

- enrichment_name:

  Optional stored enrichment result.

- cluster_by:

  Metadata column containing cluster labels.

- n_markers:

  Number of marker genes to retain per cluster.

- n_terms:

  Number of enrichment terms to retain.

- cluster_col:

  Deprecated alias for `cluster_by`.

## Value

A structured list ready for prompt construction.

## Examples

``` r
if (requireNamespace("Seurat", quietly = TRUE)) {
  counts <- matrix(rpois(10 * 24, lambda = 1), nrow = 10, ncol = 24)
  rownames(counts) <- c(
    "CD3D", "CD3E", "TRAC", "LTB", "MS4A1",
    "CD79A", "HLA-DRA", "LYZ", "ACTB", "MALAT1"
  )
  colnames(counts) <- paste0("cell", 1:24)
  obj <- sn_initialize_seurat_object(counts, species = "human")
  obj$cell_type <- rep(c("Tcell", "Bcell"), each = 12)
  Seurat::Idents(obj) <- obj$cell_type
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- sn_find_de(obj, analysis = "markers", group_by = "cell_type",
    layer = "data", min_pct = 0, logfc_threshold = 0,
    store_name = "celltype_markers", return_object = TRUE, verbose = FALSE
  )
  obj <- sn_store_enrichment(
    obj,
    tibble::tibble(ID = "GO:0001", Description = "immune response", NES = 2, p.adjust = 0.01),
    store_name = "demo_gsea"
  )
  evidence <- sn_prepare_results_evidence(
    obj,
    cluster_de_name = "celltype_markers",
    enrichment_name = "demo_gsea",
    cluster_col = "cell_type"
  )
  names(evidence)
}
#> INFO [2026-05-05 20:17:34] Initializing Seurat object for project: Shennong.
#> INFO [2026-05-05 20:17:34] Running QC metrics for human.
#> INFO [2026-05-05 20:17:34] Seurat object initialization complete.
#> Warning: No DE genes identified
#> Warning: `cluster_col` is deprecated; use `cluster_by` instead.
#> [1] "task"               "dataset"            "cluster_summary"   
#> [4] "cluster_markers"    "enrichment_summary"
```
