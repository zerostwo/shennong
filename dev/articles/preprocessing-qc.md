# Preprocessing and QC with PBMC3k

Preprocessing should answer four questions before clustering:

1.  What object did we create from the count matrix?
2.  Which low-quality cells were flagged, and why?
3.  Which genes are retained for downstream analysis?
4.  What QC evidence can be stored and revisited later?

Shennong keeps those decisions explicit while keeping the code short.
The examples use PBMC3k and are evaluated only when
`SHENNONG_RUN_VIGNETTES=true`.

## Load PBMC3k and inspect species-aware QC fields

`sn_load_data("pbmc3k")` downloads or reuses the cached 10x H5 file,
reads it, and initializes a Seurat object with human QC metrics.

``` r

library(Shennong)
library(Seurat)
library(dplyr)

pbmc <- sn_load_data("pbmc3k")
#> INFO [2026-05-06 00:00:42] Initializing Seurat object for project: pbmc3k.
#> INFO [2026-05-06 00:00:42] Running QC metrics for human.
#> INFO [2026-05-06 00:00:42] Seurat object initialization complete.
pbmc$sample <- "pbmc3k"

sn_get_species(pbmc)
#> [1] "human"

head(pbmc[[]][, c("sample", "nCount_RNA", "nFeature_RNA", "percent.mt")])
#>                  sample nCount_RNA nFeature_RNA percent.mt
#> AAACATACAACCAC-1 pbmc3k       2844         1046  2.7777778
#> AAACATTGAGCTAC-1 pbmc3k       5717         1722  3.2534546
#> AAACATTGATCAGC-1 pbmc3k       3766         1500  0.6903877
#> AAACCGTGCTTCCG-1 pbmc3k       2981         1153  1.5095606
#> AAACCGTGTATGCG-1 pbmc3k       1229          698  0.8950366
#> AAACGCACTGGTAC-1 pbmc3k       2536         1032  1.4984227
```

If you start from a raw matrix instead, the equivalent explicit form is:

``` r

counts <- sn_read("pbmc3k_filtered_feature_bc_matrix.h5")

pbmc <- sn_initialize_seurat_object(
  x = counts,
  project = "pbmc3k_demo",
  sample_name = "pbmc3k",
  study = "10x_pbmc",
  species = "human"
)
```

When gene identifiers or mixed symbol styles are a concern, standardize
at initialization or call
[`sn_standardize_gene_symbols()`](https://songqi.org/shennong/dev/reference/sn_standardize_gene_symbols.md)
before downstream analysis. If `HGNChelper` cannot provide an
unambiguous replacement for an otherwise valid symbol, Shennong
preserves the original symbol instead of returning an `NA` feature name.

``` r

pbmc <- sn_standardize_gene_symbols(
  pbmc,
  species = "human",
  is_gene_id = FALSE
)
```

## Filter cells with inspectable QC flags

[`sn_filter_cells()`](https://songqi.org/shennong/dev/reference/sn_filter_cells.md)
uses median absolute deviation thresholds by default. The key point is
that filtering and flagging are separated: set `filter = FALSE` to see
which cells would be removed, then set `filter = TRUE` when the rule is
acceptable.

``` r

pbmc_flagged <- sn_filter_cells(
  x = pbmc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  method = "mad",
  n = c(5, 5, 3),
  plot = FALSE,
  filter = FALSE
)

qc_flag_cols <- grep("^qc_", colnames(pbmc_flagged[[]]), value = TRUE)
head(pbmc_flagged[[]][, qc_flag_cols, drop = FALSE])
#> data frame with 0 columns and 6 rows
```

Once the flags look sensible, apply the same rule:

``` r

pbmc_cells <- sn_filter_cells(
  x = pbmc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  method = "mad",
  n = c(5, 5, 3),
  plot = FALSE,
  filter = TRUE
)

c(before = ncol(pbmc), after = ncol(pbmc_cells))
#> before  after 
#>   2753   2597
```

## Filter genes by expression and annotation

[`sn_filter_genes()`](https://songqi.org/shennong/dev/reference/sn_filter_genes.md)
keeps the threshold and optional annotation rule in one place. The
expression threshold protects sparse downstream steps, while
`gene_class = "coding"` can remove noncoding genes when your workflow is
protein-coding focused.

``` r

pbmc_genes <- sn_filter_genes(
  x = pbmc_cells,
  min_cells = 5,
  plot = FALSE,
  filter = TRUE,
  species = "human",
  gene_class = "coding"
)
#> WARN [2026-05-06 00:00:45] Annotation-based gene filtering could not match 17 features for species 'human'. Those unmatched features will be dropped. Examples: LINC01115.1, PCBP1-AS1.1, LSP1P5.1, DDX11L2.1, LINC01618.1, CAST.1, RAET1E-AS1.1, LINC03021.1, LINC03023.1, BMS1P14.1.

c(before = nrow(pbmc_cells), after = nrow(pbmc_genes))
#> before  after 
#>  54872  12403
```

## Normalize and add cell-cycle scores

[`sn_normalize_data()`](https://songqi.org/shennong/dev/reference/sn_normalize_data.md)
and
[`sn_score_cell_cycle()`](https://songqi.org/shennong/dev/reference/sn_score_cell_cycle.md)
keep preprocessing inside the same Shennong vocabulary used later by
clustering and plotting.

``` r

pbmc_norm <- sn_normalize_data(
  object = pbmc_genes,
  normalization_method = "seurat",
  verbose = FALSE
)

pbmc_norm <- sn_score_cell_cycle(pbmc_norm, species = "human")

table(pbmc_norm$Phase)
#> 
#>  G1 G2M   S 
#> 872 896 829
```

## Optional doublet and ambient-RNA correction

Doublet detection and ambient correction depend on optional backends.
Shennong keeps them as explicit steps rather than hiding them in
clustering.

``` r

pbmc_doublets <- sn_find_doublets(pbmc_norm)
table(pbmc_doublets$doublet_class)
```

For ambient RNA, PBMC3k raw counts can be loaded from the same example
registry. The corrected matrix is written to a new layer by default, so
the original counts remain inspectable.

``` r

pbmc_raw <- sn_load_data("pbmc3k", matrix_type = "raw")

pbmc_decont <- sn_remove_ambient_contamination(
  x = pbmc_norm,
  raw = pbmc_raw,
  method = "decontx",
  layer = "decontaminated_counts",
  verbose = FALSE
)
```

## Store a QC assessment

[`sn_assess_qc()`](https://songqi.org/shennong/dev/reference/sn_assess_qc.md)
summarizes the object after filtering and can compare it to the
pre-filter object. Store the result when you want later reports to
recover the QC decision without rerunning preprocessing.

``` r

qc_report <- sn_assess_qc(
  object = pbmc_norm,
  reference = pbmc,
  sample_by = "sample",
  verbose = FALSE
)

qc_report$overall
#>   n_samples n_cells qc_score qc_label retention_fraction
#> 1         1    2597 84.92677     good          0.9433345
#>   low_quality_removed_fraction doublet_removed_fraction
#> 1                          NaN                      NaN

pbmc_norm <- sn_assess_qc(
  object = pbmc_norm,
  reference = pbmc,
  sample_by = "sample",
  store_name = "pbmc3k_preprocessing",
  return_object = TRUE,
  verbose = FALSE
)

names(pbmc_norm@misc$qc_assessments)
#> [1] "pbmc3k_preprocessing"
```

The object is now ready for clustering. The important part is not just
that cells and genes were filtered, but that the rationale and summary
are stored in a place downstream code can find.
