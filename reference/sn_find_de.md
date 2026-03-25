# Run differential expression analysis on a Seurat object

This function provides a single entry point for:

- marker discovery across all groups,

- direct contrasts between two groups, and

- pseudobulk contrasts aggregated by sample.

## Usage

``` r
sn_find_de(
  object,
  analysis = NULL,
  ident_1 = NULL,
  ident_2 = NULL,
  group_by = NULL,
  subset_by = NULL,
  subset_levels = NULL,
  sample_col = NULL,
  assay = "RNA",
  layer = NULL,
  features = NULL,
  method = NULL,
  only_pos = NULL,
  logfc_threshold = 0.1,
  min_pct = 0.25,
  p_val_cutoff = 0.05,
  de_logfc = 0.25,
  min_cells_per_sample = 10,
  store_name = "default",
  return_object = TRUE,
  verbose = TRUE,
  ...
)
```

## Arguments

- object:

  A `Seurat` object.

- analysis:

  One of `"markers"`, `"contrast"`, or `"pseudobulk"`. If `NULL`, the
  function infers the analysis from the other arguments.

- ident_1, ident_2:

  Group labels for direct contrasts. These are required for `"contrast"`
  and `"pseudobulk"` analyses.

- group_by:

  Optional metadata column that defines the groups to compare. When
  omitted, Seurat identities are used.

- subset_by:

  Optional metadata column used to repeat a contrast within each subset,
  for example per cell type.

- subset_levels:

  Optional character vector of subset values to analyze. Defaults to all
  observed values in `subset_by`.

- sample_col:

  Metadata column containing sample IDs. Required for `"pseudobulk"`
  analyses.

- assay:

  Assay used for DE analysis. Defaults to `"RNA"`.

- layer:

  Assay layer used for DE analysis. Defaults to `"data"` for marker and
  contrast analyses and to `"counts"` for pseudobulk analyses.

- features:

  Optional feature subset to test.

- method:

  Statistical method. For `"markers"` and `"contrast"`, this can be any
  Seurat `test.use` value or `"COSGR"` for marker discovery. For
  `"pseudobulk"`, choose one of `"DESeq2"`, `"edgeR"`, or `"limma"`.

- only_pos:

  Whether to return only positive markers. Defaults to `TRUE` for
  `"markers"` and `FALSE` otherwise.

- logfc_threshold, min_pct:

  Standard Seurat marker filtering arguments.

- p_val_cutoff:

  Adjusted p-value threshold used when storing result metadata.

- de_logfc:

  Absolute log fold-change threshold used when storing result metadata.

- min_cells_per_sample:

  Minimum cells required for a sample/group pseudobulk profile to be
  retained.

- store_name:

  Key used under `object@misc$de_results`.

- return_object:

  If `TRUE`, return the updated Seurat object with stored DE results.
  Otherwise return the result table.

- verbose:

  Whether to emit progress information.

- ...:

  Additional arguments passed through to the selected DE method.

## Value

Either a DE result table or an updated `Seurat` object.

## Details

Results can be stored back into `object@misc$de_results[[store_name]]`,
so downstream helpers such as
[`sn_plot_dot()`](https://songqi.org/shennong/reference/sn_plot_dot.md)
can reuse the same marker table.

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

  marker_tbl <- sn_find_de(
    obj,
    analysis = "markers",
    group_by = "cell_type",
    layer = "data",
    min_pct = 0,
    logfc_threshold = 0,
    return_object = FALSE,
    verbose = FALSE
  )
  head(marker_tbl)

  obj <- sn_find_de(
    obj,
    analysis = "markers",
    group_by = "cell_type",
    layer = "data",
    min_pct = 0,
    logfc_threshold = 0,
    store_name = "celltype_markers",
    return_object = TRUE,
    verbose = FALSE
  )
  names(obj@misc$de_results)
}
#> INFO [2026-03-25 17:02:15] Initializing Seurat object for project: Shennong
#> INFO [2026-03-25 17:02:15] Running QC metrics for human ...
#> INFO [2026-03-25 17:02:16] Seurat object initialization complete.
#> For a (much!) faster implementation of the Wilcoxon Rank Sum Test,
#> (default method for FindMarkers) please install the presto package
#> --------------------------------------------
#> install.packages('devtools')
#> devtools::install_github('immunogenomics/presto')
#> --------------------------------------------
#> After installation of presto, Seurat will automatically use the more 
#> efficient implementation (no further action necessary).
#> This message will be shown once per session
#> [1] "celltype_markers"
```
