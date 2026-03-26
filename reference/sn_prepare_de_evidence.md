# Prepare differential-expression evidence from a stored DE result

Prepare differential-expression evidence from a stored DE result

## Usage

``` r
sn_prepare_de_evidence(object, de_name, n_genes = 15)
```

## Arguments

- object:

  A `Seurat` object.

- de_name:

  Name of a stored DE result in `object@misc$de_results`.

- n_genes:

  Number of top genes to include per direction or group.

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
  counts[c("CD3D", "CD3E", "TRAC"), 1:12] <-
    counts[c("CD3D", "CD3E", "TRAC"), 1:12] + 20
  counts[c("MS4A1", "CD79A", "HLA-DRA"), 13:24] <-
    counts[c("MS4A1", "CD79A", "HLA-DRA"), 13:24] + 20
  obj <- sn_initialize_seurat_object(counts, species = "human")
  obj$cell_type <- rep(c("Tcell", "Bcell"), each = 12)
  Seurat::Idents(obj) <- obj$cell_type
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- sn_find_de(obj, analysis = "markers", group_by = "cell_type",
    layer = "data", min_pct = 0, logfc_threshold = 0,
    store_name = "celltype_markers", return_object = TRUE, verbose = FALSE
  )
  evidence <- sn_prepare_de_evidence(obj, de_name = "celltype_markers", n_genes = 3)
  names(evidence)
}
#> INFO [2026-03-26 18:52:34] Initializing Seurat object for project: Shennong.
#> INFO [2026-03-26 18:52:34] Running QC metrics for human.
#> INFO [2026-03-26 18:52:34] Seurat object initialization complete.
#> [1] "task"             "source_de_name"   "summary"          "top_markers"     
#> [5] "top_marker_table" "caveats"         
```
