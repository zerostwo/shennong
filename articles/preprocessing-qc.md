# Preprocessing and QC workflow

This article covers the first stage of a typical Shennong analysis:

- constructing a Seurat object
- inferring species and standardizing feature naming
- filtering poor-quality cells
- filtering genes by expression and bundled GENCODE annotations
- normalizing data
- adding cell-cycle scores for downstream modeling

The examples use the built-in `pbmc_small` and `pbmc_small_raw` objects
so the workflow stays check-safe and does not depend on network
downloads.

``` r
library(Shennong)
library(dplyr)
library(knitr)
library(Seurat)

if (!exists("pbmc_small", inherits = FALSE)) {
  try(data("pbmc_small", package = "Shennong", envir = environment()), silent = TRUE)
}
if (!exists("pbmc_small", inherits = FALSE) && file.exists(file.path("data", "pbmc_small.rda"))) {
  load(file.path("data", "pbmc_small.rda"))
}

if (!exists("pbmc_small_raw", inherits = FALSE)) {
  try(data("pbmc_small_raw", package = "Shennong", envir = environment()), silent = TRUE)
}
if (!exists("pbmc_small_raw", inherits = FALSE) && file.exists(file.path("data", "pbmc_small_raw.rda"))) {
  load(file.path("data", "pbmc_small_raw.rda"))
}
```

## 1. Initialize a Seurat object

Start from a matrix-like object and let Shennong create the Seurat
object with species-aware QC metadata.

``` r
knitr::kable(feature_summary)
```

| metric            | value |
|:------------------|:------|
| species           | human |
| median_nFeature   | 1280  |
| median_nCount     | 3353  |
| median_percent_mt | 2.07  |

## 2. Filter poor-quality cells

[`sn_filter_cells()`](https://songqi.org/shennong/reference/sn_filter_cells.md)
works on standard QC columns such as `nFeature_RNA`, `nCount_RNA`, and
`percent.mt`. The function can also draw diagnostics, but the article
disables plotting to keep the rendered workflow compact.

``` r
knitr::kable(qc_summary)
```

| stage         | cells | genes |
|:--------------|------:|------:|
| initialized   |   200 | 54872 |
| cell_filtered |   185 | 54872 |
| gene_filtered |   185 | 10707 |
| normalized    |   185 | 10707 |

## 3. Filter genes with expression and annotation rules

Shennong ships bundled human and mouse GENCODE snapshots. That means
[`sn_filter_genes()`](https://songqi.org/shennong/reference/sn_filter_genes.md)
can combine expression thresholds with annotation-aware filters such as
`gene_class = "coding"` or exact `gene_type` values.

``` r
pbmc_genes <- sn_filter_genes(
  pbmc_cells,
  min_cells = 3,
  gene_class = "coding",
  plot = FALSE
)
```

For finer control you can target exact biotypes:

``` r
pbmc_lnc_only <- sn_filter_genes(
  pbmc_cells,
  min_cells = 3,
  gene_type = c("lncRNA", "antisense"),
  plot = FALSE
)
```

## 4. Normalize and add cell-cycle scores

After QC, normalize the object and add cell-cycle scores if they are
useful for regression, cluster interpretation, or doublet diagnostics.

``` r
knitr::kable(cycle_summary)
```

| Phase | cells |
|:------|------:|
| G1    |    61 |
| G2M   |    62 |
| S     |    62 |

## 5. Optional QC extensions

Shennong also exposes heavier QC helpers that are usually run only when
the analysis requires them:

``` r
# Doublet detection
pbmc_doublets <- sn_find_doublets(pbmc_norm)

# Ambient RNA correction from a raw count matrix
pbmc_decont <- sn_remove_ambient_contamination(
  pbmc_norm,
  raw_counts = pbmc_small_raw
)
```

## Where this stage leads next

Once preprocessing is complete, the usual next step is
[`sn_run_cluster()`](https://songqi.org/shennong/reference/sn_run_cluster.md)
for single-sample clustering or Harmony-based integration. The dedicated
clustering article continues from this point.
