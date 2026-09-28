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
  sample_by = NULL,
  assay = "RNA",
  layer = NULL,
  features = NULL,
  method = NULL,
  only_pos = NULL,
  logfc_threshold = 0.1,
  min_pct = 0.25,
  min_cells_per_sample = 10,
  result_id = "default",
  return_object = TRUE,
  verbose = TRUE,
  modality = c("auto", "single_cell", "bulk"),
  metadata = NULL,
  design = ~condition,
  contrast = NULL,
  backend_control = list(),
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

- sample_by:

  Metadata column containing sample IDs. Required for `"pseudobulk"`
  analyses.

- assay:

  Assay used for DE analysis. Defaults to `"RNA"`.

- layer:

  Assay layer used for DE analysis. Defaults to `"data"` for marker and
  contrast analyses and to `"counts"` for pseudobulk analyses.

- features:

  Optional feature subset to test. Requested and available candidates
  are stored separately from `input$tested_features`, the features with
  finite backend statistics before marker significance filtering.
  `input$tested_features_by_comparison` retains each comparison's
  background for downstream ORA. COSG records scored features.

- method:

  Statistical method. For `"markers"` and `"contrast"`, this can be any
  Seurat `test.use` value or `"COSGR"` for marker discovery. For
  `"pseudobulk"`, choose one of `"DESeq2"`, `"edgeR"`, or `"limma"`. For
  bulk input, choose `"auto"`, `"edger"`, `"deseq2"`, `"limma"`, or
  `"dream"`.

- only_pos:

  Whether to return only positive markers. Defaults to `TRUE` for
  `"markers"` and `FALSE` otherwise.

- logfc_threshold, min_pct:

  Standard Seurat marker filtering arguments.

- min_cells_per_sample:

  Minimum cells required for a sample/group pseudobulk profile to be
  retained.

- result_id:

  Stable identifier used as the stored-result key and canonical result
  name.

- return_object:

  If `TRUE`, return the updated Seurat object with stored DE results.
  Otherwise return the unified result; retrieve its table with
  [`sn_get_de_result()`](https://zerostwo.github.io/shennong/dev/reference/sn_get_de_result.md)
  or `result$tables$primary`.

- verbose:

  Whether to emit progress information.

- modality:

  Input modality. `"auto"` selects single-cell analysis for Seurat
  objects and bulk analysis for matrices, lists, and
  `SummarizedExperiment` objects.

- metadata:

  Optional sample metadata for bulk analysis.

- design:

  A fixed- or mixed-effects formula for bulk analysis. For pseudobulk
  analysis, an optional fixed-effects formula over cell metadata;
  include `sample_by` when the same biological unit contributes both
  contrast groups. When omitted, independent and completely paired
  designs are detected automatically.

- contrast:

  Character triple giving the contrast as variable, numerator, and
  denominator. For pseudobulk analysis this must agree with
  `c(group_by, ident_1, ident_2)`.

- backend_control:

  Bulk backend controls or a custom bulk `runner`/precomputed `result`.

- ...:

  Additional arguments passed through to the selected DE method.

## Value

For single-cell input, either a unified DE result or an updated `Seurat`
object. For bulk input, a validated Shennong bulk-DE result.

## Details

Results are stored in the canonical Shennong result registry, so
downstream helpers such as
[`sn_plot_dot()`](https://zerostwo.github.io/shennong/dev/reference/sn_plot_dot.md)
can reuse the same marker table.

Testing retains the complete backend result. Select significant genes
afterward with
`sn_get_de_result(p_adjusted_cutoff = ..., logfc_threshold = ...)`.

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

  marker_result <- sn_find_de(
    obj,
    analysis = "markers",
    group_by = "cell_type",
    layer = "data",
    min_pct = 0,
    logfc_threshold = 0,
    return_object = FALSE,
    verbose = FALSE
  )
  head(sn_get_de_result(marker_result))

  obj <- sn_find_de(
    obj,
    analysis = "markers",
    group_by = "cell_type",
    layer = "data",
    min_pct = 0,
    logfc_threshold = 0,
    result_id = "celltype_markers",
    return_object = TRUE,
    verbose = FALSE
  )
  sn_list_results(obj, type = "de")
  sn_get_result(obj, type = "de", result_id = "celltype_markers")
}
#> INFO [2026-09-28 08:49:13] Initializing Seurat object for project: Shennong.
#> INFO [2026-09-28 08:49:14] Running QC metrics for human.
#> INFO [2026-09-28 08:49:14] Seurat object initialization complete.
#> For a (much!) faster implementation of the Wilcoxon Rank Sum Test,
#> (default method for FindMarkers) please install the presto package
#> --------------------------------------------
#> install.packages('devtools')
#> devtools::install_github('immunogenomics/presto')
#> --------------------------------------------
#> After installation of presto, Seurat will automatically use the more 
#> efficient implementation (no further action necessary).
#> This message will be shown once per session
#> $schema_version
#> [1] "2.0.0"
#> 
#> $analysis
#> [1] "markers"
#> 
#> $method
#> [1] "wilcox"
#> 
#> $group_by
#> [1] "cell_type"
#> 
#> $group_col
#> [1] "cluster"
#> 
#> $ident_1
#> NULL
#> 
#> $ident_2
#> NULL
#> 
#> $subset_by
#> NULL
#> 
#> $sample_col
#> NULL
#> 
#> $design
#> NULL
#> 
#> $contrast
#> NULL
#> 
#> $parameters
#> list()
#> 
#> $input
#> $input$assay
#> [1] "RNA"
#> 
#> $input$layer
#> [1] "data"
#> 
#> $input$requested_features
#> NULL
#> 
#> $input$candidate_features
#>  [1] "GENE1"   "GENE2"   "GENE3"   "GENE4"   "GENE5"   "GENE6"   "GENE7"  
#>  [8] "GENE8"   "GENE9"   "GENE10"  "GENE11"  "GENE12"  "GENE13"  "GENE14" 
#> [15] "CD3D"    "CD3E"    "TRAC"    "MS4A1"   "CD79A"   "HLA-DRA"
#> 
#> $input$tested_features
#>  [1] "GENE1"   "GENE2"   "GENE3"   "GENE4"   "GENE5"   "GENE6"   "GENE7"  
#>  [8] "GENE8"   "GENE9"   "GENE10"  "GENE11"  "GENE12"  "GENE13"  "GENE14" 
#> [15] "CD3D"    "CD3E"    "TRAC"    "MS4A1"   "CD79A"   "HLA-DRA"
#> 
#> $input$tested_features_count
#> [1] 20
#> 
#> $input$tested_features_source
#> [1] "backend_statistics"
#> 
#> $input$tested_feature_groups
#> [1] "cluster"
#> 
#> $input$tested_features_by_comparison
#> # A tibble: 20 × 2
#>    gene    cluster
#>    <chr>   <chr>  
#>  1 CD3E    Tcell  
#>  2 CD3D    Tcell  
#>  3 TRAC    Tcell  
#>  4 GENE10  Tcell  
#>  5 GENE14  Tcell  
#>  6 GENE6   Tcell  
#>  7 GENE1   Tcell  
#>  8 GENE2   Tcell  
#>  9 GENE11  Tcell  
#> 10 GENE7   Tcell  
#> 11 GENE3   Tcell  
#> 12 MS4A1   Bcell  
#> 13 HLA-DRA Bcell  
#> 14 CD79A   Bcell  
#> 15 GENE4   Bcell  
#> 16 GENE8   Bcell  
#> 17 GENE12  Bcell  
#> 18 GENE13  Bcell  
#> 19 GENE5   Bcell  
#> 20 GENE9   Bcell  
#> 
#> 
#> $assay
#> [1] "RNA"
#> 
#> $layer
#> [1] "data"
#> 
#> $rank_col
#> [1] "avg_log2FC"
#> 
#> $p_col
#> [1] "p_val_adj"
#> 
#> $min_pct
#> [1] 0
#> 
#> $logfc_threshold
#> [1] 0
#> 
#> $n_genes
#> [1] 6
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
#> [1] "2026-09-28 08:49:14 UTC"
#> 
#> $provenance$acceleration
#> $provenance$acceleration$patches
#> [1] "seurat"
#> 
#> $provenance$acceleration$used_patches
#> character(0)
#> 
#> $provenance$acceleration$suppressed_patches
#> [1] "seurat"
#> 
#> 
#> $provenance$result_id
#> [1] "celltype_markers"
#> 
#> $provenance$analysis_type
#> [1] "de"
#> 
#> 
#> $analysis_type
#> [1] "de"
#> 
#> $result_id
#> [1] "celltype_markers"
#> 
#> $backend
#> [1] "wilcox"
#> 
#> $tables
#> $tables$primary
#> # A tibble: 6 × 7
#>         p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene   
#>         <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>  
#> 1 0.000000740       4.52     1 0.583 0.0000148 Tcell   CD3E   
#> 2 0.000000740       4.32     1 0.667 0.0000148 Tcell   CD3D   
#> 3 0.000000740       4.15     1 0.667 0.0000148 Tcell   TRAC   
#> 4 0.000000740       4.07     1 0.75  0.0000148 Bcell   MS4A1  
#> 5 0.000000740       3.61     1 0.833 0.0000148 Bcell   HLA-DRA
#> 6 0.000000740       4.72     1 0.583 0.0000148 Bcell   CD79A  
#> 
#> 
#> $embeddings
#> list()
#> 
#> $graphs
#> list()
#> 
#> $models
#> list()
#> 
#> $diagnostics
#> list()
#> 
#> $warnings
#> character(0)
#> 

bulk_counts <- matrix(
  stats::rpois(40 * 6, 20), nrow = 40,
  dimnames = list(paste0("gene_", 1:40), paste0("sample_", 1:6))
)
bulk_metadata <- data.frame(
  condition = factor(rep(c("control", "treated"), each = 3)),
  row.names = colnames(bulk_counts)
)
bulk_de <- sn_find_de(
  bulk_counts,
  metadata = bulk_metadata,
  design = ~condition,
  contrast = c("condition", "treated", "control"),
  backend_control = list(result = data.frame(
    gene = rownames(bulk_counts), logFC = 0, PValue = 1, FDR = 1
  ))
)
```
