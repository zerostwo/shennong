# Prepare differential-expression evidence from a stored DE result

Prepare differential-expression evidence from a stored DE result

## Usage

``` r
sn_prepare_de_evidence(object, de_result_id, n_genes = 15)
```

## Arguments

- object:

  A `Seurat` object.

- de_result_id:

  Identifier of a stored DE result.

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
    result_id = "celltype_markers", return_object = TRUE, verbose = FALSE
  )
  evidence <- sn_prepare_de_evidence(obj, de_result_id = "celltype_markers", n_genes = 3)
  names(evidence)
}
#> INFO [2026-11-02 09:55:45] Initializing Seurat object for project: Shennong.
#> INFO [2026-11-02 09:55:45] Running QC metrics for human.
#> INFO [2026-11-02 09:55:45] Seurat object initialization complete.
#> Warning: restarting interrupted promise evaluation
#> Warning: restarting interrupted promise evaluation
#> Warning: restarting interrupted promise evaluation
#> Warning: restarting interrupted promise evaluation
#> Warning: restarting interrupted promise evaluation
#> Warning: restarting interrupted promise evaluation
#> Warning: restarting interrupted promise evaluation
#> Warning: restarting interrupted promise evaluation
#> Warning: No DE genes identified
#> Warning: The following tests were not performed: 
#> Warning: When testing Tcell versus all:
#>  Failed to get one or more globals that were identified via static code inspection (function (x); {; vals <- data.use[x, ]; return(wilcox.test(x = vals[j], y = vals[-j], ...)$p.value); }). Specifically, the 1 globals ‘...’ give an error. The reason was: argument "ident.1" is missing, with no default
#> Warning: When testing Bcell versus all:
#>  Failed to get one or more globals that were identified via static code inspection (function (x); {; vals <- data.use[x, ]; return(wilcox.test(x = vals[j], y = vals[-j], ...)$p.value); }). Specifically, the 1 globals ‘...’ give an error. The reason was: argument "ident.1" is missing, with no default
#> Warning: Unknown or uninitialised column: `gene`.
#> Error in sn_prepare_de_evidence(obj, de_result_id = "celltype_markers",     n_genes = 3): The stored DE result does not contain a ranking column.
```
