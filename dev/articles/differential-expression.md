# Differential expression and enrichment

Choose the statistical question before choosing a backend. Finding
cluster markers and testing a replicated treatment effect use different
sampling units.

| Question | Call | Sampling unit |
|----|----|----|
| Which genes distinguish each cluster? | `sn_find_de(analysis = "markers")` | Cells; exploratory cluster description |
| Which genes differ between two cell groups? | `sn_find_de(analysis = "contrast")` | Cells unless a sample-aware method is selected |
| Which genes change between replicated conditions? | `sn_find_de(analysis = "pseudobulk", sample_by = ...)` | Biological samples |
| Which genes change in a bulk count matrix? | `sn_find_de(counts, metadata = ..., design = ..., contrast = ...)` | Matrix columns / samples |

## Run a complete marker example

This example runs without downloading data. The bundled PBMC object
contains precomputed clusters, which are used only to demonstrate marker
extraction.

``` r

library(Shennong)
data("pbmc_small", package = "SeuratObject")
pbmc <- SeuratObject::CreateSeuratObject(
  SeuratObject::GetAssayData(pbmc_small, assay = "RNA", layer = "counts")
)
pbmc$cluster <- as.character(SeuratObject::Idents(pbmc_small))
pbmc <- sn_normalize_data(pbmc, method = "seurat", verbose = FALSE)

pbmc <- sn_find_de(
  pbmc, analysis = "markers", group_by = "cluster", method = "wilcox",
  assay = "RNA", layer = "data", min_pct = 0.1, logfc_threshold = 0,
  result_id = "markers", verbose = FALSE
)
sn_list_results(pbmc, type = "de")
#> # A tibble: 1 × 8
#>   collection       type  result_id analysis method created_at      n_rows source
#>   <chr>            <chr> <chr>     <chr>    <chr>  <chr>            <int> <chr> 
#> 1 shennong.results de    markers   de       wilcox 2026-09-24 00:…    134 NA
```

Each cluster is compared with the remaining cells. Keep the full result
for reproducibility, and filter only the table used for reporting:

``` r

full <- sn_get_result(pbmc, type = "de", result_id = "markers")
markers <- sn_get_de_result(
  pbmc, "markers", direction = "up", p_adjusted_cutoff = 0.05,
  logfc_threshold = 0.25, top_n = 5, top_scope = "group"
)
head(markers)
#> # A tibble: 6 × 7
#>        p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene     
#>        <dbl>      <dbl> <dbl> <dbl>     <dbl> <fct>   <chr>    
#> 1 0.0000151        7.60 0.444 0.068  0.00347  0       GNLY     
#> 2 0.0000201        4.28 0.444 0.068  0.00462  0       LAMP1    
#> 3 0.00000397       3.88 0.528 0.091  0.000913 0       CD7      
#> 4 0.000169         3.19 0.389 0.068  0.0388   0       CD247    
#> 5 0.0000997        2.72 0.5   0.114  0.0229   0       LCK      
#> 6 0.0000421       10.9  0.263 0      0.00968  2       LINC00926
full$input[c("assay", "layer", "tested_features_count", "tested_features_source")]
#> $assay
#> [1] "RNA"
#> 
#> $layer
#> [1] "data"
#> 
#> $tested_features_count
#> [1] 222
#> 
#> $tested_features_source
#> [1] "backend_statistics"
```

Testing filters (`min_pct`, `logfc_threshold` in `sn_find_de`) affect
the set of hypotheses and the tested-gene background. Getter filters
affect the output view only. For enrichment, preserve the tested
background rather than using all genes in the genome or only the
significant genes.

## Compare two groups explicitly

`ident_1` is the numerator/target group and `ident_2` is the reference
group. The example compares two existing clusters; it is not a treatment
comparison.

``` r

levels <- sort(unique(pbmc$cluster))
contrast_result <- sn_find_de(
  pbmc, analysis = "contrast", group_by = "cluster",
  ident_1 = levels[[1]], ident_2 = levels[[2]],
  method = "wilcox", layer = "data", min_pct = 0, logfc_threshold = 0,
  result_id = "cluster_contrast", return_object = FALSE, verbose = FALSE
)
head(sn_get_de_result(contrast_result, direction = "up", top_n = 5))
#> # A tibble: 5 × 6
#>   gene       p_val avg_log2FC pct.1 pct.2 p_val_adj
#>   <chr>      <dbl>      <dbl> <dbl> <dbl>     <dbl>
#> 1 PF4   0.0109           12.2 0.25      0   1      
#> 2 CD7   0.00000579       11.2 0.528     0   0.00133
#> 3 TUBB1 0.0109           10.9 0.25      0   1      
#> 4 GZMH  0.0109           10.8 0.25      0   1      
#> 5 GZMM  0.000158         10.4 0.417     0   0.0363
```

Because `return_object = FALSE`, this returns a unified result without
adding it to `pbmc`. Pass that result to the same table getter, or use
[`sn_store_result()`](https://zerostwo.github.io/shennong/dev/reference/sn_store_result.md)
to save it later.

## Test replicated conditions with pseudobulk

Use your own multi-sample object for this recipe. It needs raw RNA
counts and metadata columns `sample_id`, `condition`, and `cell_type`.
Each `sample_id` must identify a biological sample with a single
condition, and each condition needs biological replication. Do not
invent sample IDs by splitting cells.

``` r

# study is your multi-sample Seurat object, with counts in RNA/counts.
sample_design <- unique(study[[]][, c("sample_id", "condition")])
stopifnot(!anyDuplicated(sample_design$sample_id))
table(sample_design$condition)

study <- sn_find_de(
  study, analysis = "pseudobulk", method = "DESeq2",
  assay = "RNA", layer = "counts",
  group_by = "condition", ident_1 = "treated", ident_2 = "control",
  sample_by = "sample_id", subset_by = "cell_type", subset_levels = "T cell",
  min_cells_per_sample = 20,
  result_id = "T_cell_treated_vs_control"
)
sn_list_results(study, type = "de")
head(sn_get_de_result(study, "T_cell_treated_vs_control", p_adjusted_cutoff = 0.05))
```

The recipe uses DESeq2 through the pseudobulk path. Inspect retained
samples and the design, especially after subsetting rare cell types. For
a paired or covariate-adjusted design, see [bulk
transcriptomics](https://zerostwo.github.io/shennong/dev/articles/bulk-transcriptomics.md)
and explicitly retain donor/sample metadata. A cell-level Wilcoxon test
does not substitute for a donor-level condition test.

## Choose ORA or GSEA

| Analysis | Input | Interpretation |
|----|----|----|
| ORA | Selected genes and a tested-gene universe | Over-representation among selected genes |
| GSEA | A named numeric ranking or a gene/score table | Enrichment along a ranked list |

`mapping = gene ~ group` describes grouped ORA when `group` is
categorical. `mapping = gene ~ score` selects GSEA when the right side
is numeric. `mapping = gene ~ score | group` describes grouped GSEA. If
group labels are numeric codes, explicitly set `analysis = "ora"` to
avoid treating them as scores. Duplicate genes are rejected by default;
choose a documented collapse rule only when it matches the analysis.

### ORA from stored markers

The following recipes use the `pbmc` created above. They require
**clusterProfiler** and **msigdbr**, and MSigDB may need a first-use
download. They are shown but not executed during ordinary documentation
builds. The small marker example may yield no significant pathways.

``` r

pbmc <- sn_run_enrichment(
  pbmc, source_de_result_id = "markers", analysis = "ora",
  species = "human", database = "H",
  de_direction = "up", de_p_adjusted_cutoff = 0.05, de_logfc_threshold = 0.25,
  min_gs_size = 10,
  result_id = "marker_hallmark_ora"
)
sn_list_results(pbmc, type = "enrichment")
terms <- sn_get_enrichment_result(
  pbmc, "marker_hallmark_ora", p_adjusted_cutoff = 0.05, top_n = 10
)
head(terms)
```

`de_*` options select genes from the stored DE analysis before ORA. The
tested gene universe retained by DE is used per comparison when
available. An explicit `universe` replaces that background; supply one
only when justified. Enrichment `pvalue_cutoff` and `qvalue_cutoff` are
backend filters, while the getter applies an adjusted-p-value filter to
returned terms.

### GSEA from a complete ranking

Use all tested genes with a meaningful signed ranking statistic, not
just the top markers selected for a figure. This demonstration uses
log-fold change; choose and document the ranking statistic for your
study.

``` r

ranking <- contrast_result$tables$primary
ranked <- sn_run_enrichment(
  ranking, mapping = gene ~ avg_log2FC, analysis = "gsea",
  species = "human", database = "H",
  result_id = "contrast_hallmark_gsea", return_object = FALSE
)
head(sn_get_enrichment_result(ranked, top_n = 10))
ranked$parameters
```

For a marker table with one ranking per cluster:

``` r

grouped <- sn_run_enrichment(
  full$tables$primary, mapping = gene ~ avg_log2FC | cluster,
  analysis = "gsea", species = "human", database = "H",
  return_object = FALSE
)
head(sn_get_enrichment_result(grouped, top_n = 5, top_scope = "group"))
```

### Several databases and their return values

``` r

combined <- sn_run_enrichment(
  ranking, mapping = gene ~ avg_log2FC, analysis = "gsea",
  species = "human", database = c("H", "C2:CP:REACTOME"),
  return_object = FALSE
)
head(combined$tables$primary)       # includes a database column
names(combined$models$database_results)
names(combined$models$backend_results)
```

With `return_object = FALSE`, the combined result contains a primary
table and individual unified results in `models$database_results`;
native backend objects are in `models$backend_results`. With a Seurat
input and `return_object = TRUE`, the workflow stores one result per
database. Discover those IDs with
`sn_list_results(pbmc, type = "enrichment")` before retrieving one.
Never infer a database from list position.

## Common problems

| Symptom | Check |
|----|----|
| No genes or terms returned | Inspect the complete table and tested-gene coverage; relax display filters only for inspection |
| Multiple stored results error | List the IDs, then pass the intended `result_id` |
| Result already exists | Use a new ID, or explicitly remove the reviewed old result before rerunning |
| Missing grouping/ranking column | Inspect `names(result$tables$primary)` and the formula mapping |
| Sparse or split Seurat layers | Select the intended assay/layer; see [layer-aware workflows](https://zerostwo.github.io/shennong/dev/articles/layered-workflows.md) |
| Unexpected enrichment background | Inspect the DE tested-gene provenance and any explicit `universe` |

Continue with [annotation and
pathways](https://zerostwo.github.io/shennong/dev/articles/annotation-pathways.md)
for feature classes, reference annotation, plotting, and scoring on a
larger real dataset.
