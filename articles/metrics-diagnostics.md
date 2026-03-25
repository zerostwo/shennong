# Metrics and diagnostics workflow

This article focuses on the diagnostic layer that sits after clustering
or integration. In practice, Shennong metrics fall into three scientific
groups:

- batch removal: LISI, batch silhouette, PCR batch score
- biological conservation: graph connectivity, isolated-label
  preservation
- structure diagnostics: challenging groups, cluster entropy, cluster
  purity

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
```

## Run the aggregate assessment first

[`sn_assess_integration()`](https://songqi.org/shennong/reference/sn_assess_integration.md)
is the recommended entry point because it bundles the main metrics into
one object with summary, per-cell, and per-group outputs.

``` r
knitr::kable(summary_tbl, digits = 3)
```

| metric                       | category      | score | scaled_score | n_cells | source                    | note                  |
|:-----------------------------|:--------------|------:|-------------:|--------:|:--------------------------|:----------------------|
| batch_silhouette             | batch_removal | 0.210 |        0.790 |     200 | harmony + seurat_clusters |                       |
| batch_lisi                   | batch_removal | 1.604 |        0.604 |     200 | harmony                   |                       |
| cluster_batch_entropy        | batch_removal | 0.822 |        0.822 |     200 | seurat_clusters vs sample |                       |
| pcr_batch                    | batch_removal | 0.012 |        0.839 |     200 | harmony vs pca            |                       |
| well_resolved_group_fraction | structure     | 1.000 |        1.000 |     200 | seurat_clusters           | 4 rare groups flagged |
| overall_integration_score    | aggregate     | 0.829 |        0.829 |     200 | 0.4 batch + 0.6 biology   |                       |

## Inspect local sample mixing with LISI

Higher sample LISI indicates better local batch mixing. When `lisi` is
installed, Shennong stores the per-cell values so you can examine the
full distribution rather than relying on a single mean.

``` r
if (!is.null(batch_lisi)) {
  knitr::kable(head(batch_lisi, 10), digits = 3)
}
```

| cell_id                   | sample |
|:--------------------------|-------:|
| pbmc1k_ACGTTCCGTGGGTCAA-1 |  1.943 |
| pbmc1k_GGACGTCGTTCAACGT-1 |  1.802 |
| pbmc1k_CCAATTTCATTCGATG-1 |  1.795 |
| pbmc1k_TTTCCTCTCCTACACC-1 |  1.521 |
| pbmc1k_GGGTTATCAGCCATTA-1 |  1.101 |
| pbmc1k_CCGCAAGCATTCAGGT-1 |  1.936 |
| pbmc1k_AGGACGAAGATTAGTG-1 |  1.470 |
| pbmc1k_GTAGAGGCAACTTCTT-1 |  1.351 |
| pbmc1k_TTCACGCGTTAAGTCC-1 |  1.993 |
| pbmc1k_GTCATGAAGACTCATC-1 |  1.607 |

## Surface rare or difficult groups

The challenging-group summary is useful for small populations that do
not form an obvious island in UMAP but still look unstable in neighbor
structure or silhouette space.

``` r
knitr::kable(challenging_tbl, digits = 3)
```

| seurat_clusters | n_cells | fraction_cells | median_neighbor_purity | mean_neighbor_purity | graph_connectivity | mean_silhouette | separation_score | challenge_score | rare_group | challenging_group |
|:----------------|--------:|---------------:|-----------------------:|---------------------:|-------------------:|----------------:|-----------------:|----------------:|:-----------|:------------------|
| 2               |      45 |          0.225 |                  0.879 |                0.843 |                  1 |          -0.013 |            0.791 |           0.209 | TRUE       | FALSE             |
| 4               |      26 |          0.130 |                  0.871 |                0.811 |                  1 |           0.180 |            0.820 |           0.180 | TRUE       | FALSE             |
| 1               |      46 |          0.230 |                  1.000 |                0.995 |                  1 |           0.178 |            0.863 |           0.137 | TRUE       | FALSE             |
| 0               |      55 |          0.275 |                  0.909 |                0.859 |                  1 |           0.393 |            0.869 |           0.131 | FALSE      | FALSE             |
| 3               |      28 |          0.140 |                  1.000 |                0.980 |                  1 |           0.367 |            0.895 |           0.105 | TRUE       | FALSE             |

## Score candidate rare cells directly

[`sn_detect_rare_cells()`](https://songqi.org/shennong/reference/sn_detect_rare_cells.md)
provides a cell-level complement to the group-level diagnostics. The
native `gini` backend uses rare-gene enrichment to score cells, while
optional backends such as FiRE or Python-based scCAD can be used when
installed locally.

``` r
knitr::kable(head(rare_tbl[order(rare_tbl$rare_score, decreasing = TRUE), ], 10), digits = 3)
```

|     | cell_id                   | method | rare_score | rare_cell |
|:----|:--------------------------|:-------|-----------:|:----------|
| 23  | pbmc1k_CGAGAAGTCCAATGCA-1 | gini   |      6.185 | TRUE      |
| 51  | pbmc1k_TCATATCCATGCTGCG-1 | gini   |      5.830 | TRUE      |
| 14  | pbmc1k_TACCGGGTCCTCGATC-1 | gini   |      5.448 | TRUE      |
| 19  | pbmc1k_GTGTGGCGTAAGTTGA-1 | gini   |      5.082 | TRUE      |
| 134 | pbmc3k_TGGACCCTCATGGT-1   | gini   |      4.814 | TRUE      |
| 128 | pbmc3k_GAGTTGTGGTAGCT-1   | gini   |      4.222 | TRUE      |
| 181 | pbmc3k_TGACGATGCAAAGA-1   | gini   |      4.219 | TRUE      |
| 61  | pbmc1k_CATTGAGAGGGACCAT-1 | gini   |      4.186 | TRUE      |
| 156 | pbmc3k_GATATCCTCCCGTT-1   | gini   |      4.148 | TRUE      |
| 191 | pbmc3k_TCAGACGACGTTAG-1   | gini   |      3.914 | TRUE      |

## Inspect isolated-label preservation

[`sn_calculate_isolated_label_score()`](https://songqi.org/shennong/reference/sn_calculate_isolated_label_score.md)
and the corresponding table inside
[`sn_assess_integration()`](https://songqi.org/shennong/reference/sn_assess_integration.md)
highlight low-frequency labels or clusters that remain well-separated
after integration.

``` r
knitr::kable(isolated_tbl, digits = 3)
```

| seurat_clusters | n_cells | fraction_cells | mean_silhouette | isolated_score | isolated_label |
|:----------------|--------:|---------------:|----------------:|---------------:|:---------------|
| 4               |      26 |          0.130 |           0.180 |          0.590 | TRUE           |
| 3               |      28 |          0.140 |           0.367 |          0.684 | TRUE           |
| 2               |      45 |          0.225 |          -0.013 |          0.494 | TRUE           |
| 1               |      46 |          0.230 |           0.178 |          0.589 | TRUE           |
| 0               |      55 |          0.275 |           0.393 |          0.697 | FALSE          |

## Quantify batch mixing inside each cluster

Cluster entropy is a practical complement to LISI because it answers a
simpler question: within a given cluster, are the batches actually
mixed?

``` r
knitr::kable(entropy_tbl, digits = 3)
```

| seurat_clusters | n_cells | n_labels | dominant_label | entropy | normalized_entropy |
|:----------------|--------:|---------:|:---------------|--------:|-------------------:|
| 1               |      46 |        2 | pbmc3k         |   0.678 |              0.978 |
| 4               |      26 |        2 | pbmc3k         |   0.690 |              0.996 |
| 2               |      45 |        2 | pbmc1k         |   0.580 |              0.837 |
| 0               |      55 |        2 | pbmc3k         |   0.212 |              0.305 |
| 3               |      28 |        2 | pbmc3k         |   0.691 |              0.996 |

## Add supervised metrics when labels are available

If you already have a trusted annotation column such as `cell_type`, you
can extend the same workflow with supervised conservation metrics:

``` r
assessment <- sn_assess_integration(
  pbmc_integrated,
  batch = "sample",
  label = "cell_type",
  cluster = "seurat_clusters",
  reduction = "harmony",
  baseline_reduction = "pca"
)
```

That enables `label_lisi`, ARI/NMI, cluster purity against `cell_type`,
and a more interpretable isolated-label summary.
