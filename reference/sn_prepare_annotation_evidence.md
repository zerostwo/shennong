# Prepare cluster-annotation evidence from a Seurat object

Prepare cluster-annotation evidence from a Seurat object

## Usage

``` r
sn_prepare_annotation_evidence(
  object,
  de_name,
  cluster_col = "seurat_clusters",
  n_markers = 10
)
```

## Arguments

- object:

  A `Seurat` object.

- de_name:

  Name of a stored marker result in `object@misc$de_results`.

- cluster_col:

  Metadata column containing cluster labels.

- n_markers:

  Number of top markers to retain per cluster.

## Value

A structured list ready for prompt construction.

## Examples

``` r
if (requireNamespace("Seurat", quietly = TRUE)) {
  counts <- matrix(rpois(20 * 24, lambda = 1), nrow = 20, ncol = 24)
  rownames(counts) <- c(
    paste0("GENE", 1:14),
    "CD3D", "CD3E", "TRAC", "MS4A1", "CD79A", "HLA-DRA"
  )
  colnames(counts) <- paste0("cell", 1:24)
  counts[c("CD3D", "CD3E", "TRAC"), 1:12] <- counts[c("CD3D", "CD3E", "TRAC"), 1:12] + 20
  counts[c("MS4A1", "CD79A", "HLA-DRA"), 13:24] <- counts[c("MS4A1", "CD79A", "HLA-DRA"), 13:24] + 20
  obj <- sn_initialize_seurat_object(counts, species = "human")
  obj$cell_type <- rep(c("Tcell", "Bcell"), each = 12)
  Seurat::Idents(obj) <- obj$cell_type
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- sn_find_de(obj, analysis = "markers", group_by = "cell_type",
    layer = "data", min_pct = 0, logfc_threshold = 0,
    store_name = "celltype_markers", return_object = TRUE, verbose = FALSE
  )
  evidence <- sn_prepare_annotation_evidence(obj, de_name = "celltype_markers", cluster_col = "cell_type")
  names(evidence)
}
#> INFO [2026-03-20 02:56:28] Initializing Seurat object for project: Shennong
#> INFO [2026-03-20 02:56:28] Running QC metrics for human ...
#> Maps last updated on: Sat Nov 16 10:35:32 2024
#> Maps last updated on: Sat Nov 16 10:35:32 2024
#> INFO [2026-03-20 02:56:28] Seurat object initialization complete.
#> [1] "task"             "cluster_col"      "source_de_name"   "analysis_method" 
#> [5] "species"          "cluster_summary"  "top_marker_table" "caveats"         
```
