# Integration metrics and diagnostics

Good integration is not just a pretty UMAP. A useful diagnostic layer
should separate three questions:

- Are technical labels mixed?
- Are biological labels still separated?
- Which clusters or groups are difficult enough to inspect manually?

Shennong exposes both one-line aggregate diagnostics and individual
metrics. This article uses PBMC3k and a teaching-only pseudo-batch
column.

## Prepare a diagnostic object

``` r

library(Shennong)
library(Seurat)
library(dplyr)

pbmc <- sn_load_data("pbmc3k")
#> INFO [2026-05-05 21:31:19] Initializing Seurat object for project: pbmc3k.
#> INFO [2026-05-05 21:31:19] Running QC metrics for human.
#> INFO [2026-05-05 21:31:20] Seurat object initialization complete.
pbmc$library <- rep(c("library_a", "library_b"), length.out = ncol(pbmc))

pbmc <- sn_run_cluster(
  object = pbmc,
  batch = "library",
  normalization_method = "seurat",
  hvg_group_by = "library",
  nfeatures = 1500,
  dims = 1:15,
  resolution = 0.6,
  species = "human",
  verbose = FALSE
)
```

## Start with the aggregate assessment

[`sn_assess_integration()`](https://songqi.org/shennong/dev/reference/sn_assess_integration.md)
runs a compact panel of metrics and returns the results in predictable
sections: `summary`, `per_cell`, `per_group`, and `parameters`.

``` r

assessment <- sn_assess_integration(
  x = pbmc,
  batch = "library",
  label_by = "seurat_clusters",
  cluster_by = "seurat_clusters",
  reduction = "harmony",
  baseline_reduction = "pca",
  dims = 1:15,
  max_cells = 1500
)

assessment$summary
#>                          metric             category        score
#> 1              batch_silhouette        batch_removal 0.0291797570
#> 2              label_silhouette biology_conservation 0.2866008297
#> 3                    batch_lisi        batch_removal 1.8889673188
#> 4                    label_lisi biology_conservation 1.1758400747
#> 5            graph_connectivity biology_conservation 1.0000000000
#> 6                           ari biology_conservation 1.0000000000
#> 7                           nmi biology_conservation 1.0000000000
#> 8          isolated_label_score biology_conservation 0.7662856393
#> 9          cluster_label_purity biology_conservation 1.0000000000
#> 10        cluster_batch_entropy        batch_removal 0.9886574048
#> 11                    pcr_batch        batch_removal 0.0002694834
#> 12 well_resolved_group_fraction            structure 1.0000000000
#> 13           batch_mixing_score            aggregate 0.7294300338
#> 14   biology_conservation_score            aggregate 0.9128640383
#> 15              structure_score            aggregate 1.0000000000
#> 16    overall_integration_score            aggregate 0.8394904365
#>    scaled_score n_cells                             source
#> 1    0.97082024    1500          harmony + seurat_clusters
#> 2    0.64330041    1500                            harmony
#> 3    0.88896732    1500                            harmony
#> 4    0.98046221    1500                            harmony
#> 5    1.00000000    1500                             RNA_nn
#> 6    1.00000000    1500 seurat_clusters vs seurat_clusters
#> 7    1.00000000    1500 seurat_clusters vs seurat_clusters
#> 8    0.76628564      50                    seurat_clusters
#> 9    1.00000000    1500 seurat_clusters vs seurat_clusters
#> 10   0.98865740    1500         seurat_clusters vs library
#> 11   0.06927517    1500                     harmony vs pca
#> 12   1.00000000    1500                    seurat_clusters
#> 13   0.72943003    1500                          aggregate
#> 14   0.91286404    1500                          aggregate
#> 15   1.00000000    1500                          aggregate
#> 16   0.83949044    1500         0.4 batch_by + 0.6 biology
#>                     note
#> 1                       
#> 2                       
#> 3                       
#> 4                       
#> 5                       
#> 6                       
#> 7                       
#> 8      3 isolated labels
#> 9                       
#> 10                      
#> 11                      
#> 12 3 rare groups flagged
#> 13                      
#> 14                      
#> 15                      
#> 16
names(assessment$per_group)
#> [1] "batch_silhouette"     "graph_connectivity"  
#> [3] "clustering_agreement" "isolated_label_score"
#> [5] "cluster_purity"       "cluster_entropy"     
#> [7] "pcr_batch"            "challenging_groups"  
#> [9] "composition"
```

Use the wrapper for routine reporting, then call individual metrics when
one part of the assessment needs more detail.

## Inspect mixing and separation directly

LISI asks how diverse each cell’s local neighborhood is with respect to
a label. Here the same function can measure library mixing or cluster
mixing by changing `label`.

``` r

library_lisi <- sn_calculate_lisi(
  pbmc,
  reduction = "harmony",
  label_by = "library",
  dims = 1:15,
  max_cells = 1500
)

cluster_lisi <- sn_calculate_lisi(
  pbmc,
  reduction = "harmony",
  label_by = "seurat_clusters",
  dims = 1:15,
  max_cells = 1500
)

head(library_lisi)
#>            cell_id  library
#> 1 AAACATACAACCAC-1 1.851349
#> 2 AAACATTGATCAGC-1 1.763703
#> 3 AAACCGTGCTTCCG-1 1.693069
#> 4 AAACGCACTGGTAC-1 1.599863
#> 5 AAACGCTGACCAGT-1 1.818432
#> 6 AAACGCTGGTTCTT-1 1.998104
head(cluster_lisi)
#>            cell_id seurat_clusters
#> 1 AAACATTGAGCTAC-1        1.000000
#> 2 AAACATTGATCAGC-1        1.041738
#> 3 AAACGCTGACCAGT-1        1.085267
#> 4 AAACGCTGTAGCCA-1        2.651311
#> 5 AAACGCTGTTTCTG-1        1.115321
#> 6 AAACTTGAAAAACG-1        1.000000
```

Silhouette asks whether a label is compact in the embedding. For
biological labels, high silhouette can be good; for batch labels, high
silhouette often means residual batch separation.

``` r

cluster_sil <- sn_calculate_silhouette(
  pbmc,
  label_by = "seurat_clusters",
  reduction = "harmony",
  dims = 1:15,
  max_cells = 1500
)

head(cluster_sil)
#>            cell_id seurat_clusters silhouette_width
#> 1 AAACATTGAGCTAC-1               4       0.42128412
#> 2 AAACATTGATCAGC-1               1       0.23931000
#> 3 AAACGCTGACCAGT-1               3       0.10464782
#> 4 AAACGCTGTAGCCA-1               3      -0.04240411
#> 5 AAACGCTGTTTCTG-1               5       0.22616282
#> 6 AAACTTGAAAAACG-1               4       0.37935553
```

Graph connectivity asks whether cells with the same label remain
connected in the neighbor graph.

``` r

connectivity <- sn_calculate_graph_connectivity(
  pbmc,
  label_by = "seurat_clusters",
  reduction = "harmony",
  dims = 1:15,
  max_cells = 1500
)

connectivity
#>    seurat_clusters n_cells largest_component connectivity_score
#> 1                0     253               253                  1
#> 2                1     224               224                  1
#> 3                2     237               237                  1
#> 4                3     221               221                  1
#> 5                4     209               209                  1
#> 6                5     136               136                  1
#> 7                6     131               131                  1
#> 8                7      36                36                  1
#> 9                8      36                36                  1
#> 10               9      17                17                  1
```

## Compare clusters to labels when references exist

When you have a trusted label column, purity, entropy, and clustering
agreement show how the computed clusters relate to that reference.

``` r

pbmc$coarse_label <- paste0("cluster_", pbmc$seurat_clusters)

sn_calculate_cluster_purity(
  pbmc,
  cluster_by = "seurat_clusters",
  label_by = "coarse_label"
)
#>    seurat_clusters n_cells dominant_label dominant_label_n
#> 1                3     368      cluster_3              368
#> 2                4     357      cluster_4              357
#> 3                1     527      cluster_1              527
#> 4                2     487      cluster_2              487
#> 5                6     154      cluster_6              154
#> 6                5     165      cluster_5              165
#> 7                0     606      cluster_0              606
#> 8                9      17      cluster_9               17
#> 9                8      36      cluster_8               36
#> 10               7      36      cluster_7               36
#>    purity_score impurity_score
#> 1             1              0
#> 2             1              0
#> 3             1              0
#> 4             1              0
#> 5             1              0
#> 6             1              0
#> 7             1              0
#> 8             1              0
#> 9             1              0
#> 10            1              0

sn_calculate_cluster_entropy(
  pbmc,
  cluster_by = "seurat_clusters",
  label_by = "coarse_label"
)
#>    seurat_clusters n_cells n_labels dominant_label entropy
#> 1                3     368        1      cluster_3       0
#> 2                4     357        1      cluster_4       0
#> 3                1     527        1      cluster_1       0
#> 4                2     487        1      cluster_2       0
#> 5                6     154        1      cluster_6       0
#> 6                5     165        1      cluster_5       0
#> 7                0     606        1      cluster_0       0
#> 8                9      17        1      cluster_9       0
#> 9                8      36        1      cluster_8       0
#> 10               7      36        1      cluster_7       0
#>    normalized_entropy
#> 1                   0
#> 2                   0
#> 3                   0
#> 4                   0
#> 5                   0
#> 6                   0
#> 7                   0
#> 8                   0
#> 9                   0
#> 10                  0

sn_calculate_clustering_agreement(
  pbmc,
  cluster_by = "seurat_clusters",
  label_by = "coarse_label"
)
#>    cluster_column label_column n_cells ari nmi
#> 1 seurat_clusters coarse_label    2753   1   1

sn_calculate_isolated_label_score(
  pbmc,
  label_by = "coarse_label",
  reduction = "harmony",
  dims = 1:15,
  max_cells = 1500
)
#>    coarse_label n_cells fraction_cells mean_silhouette isolated_score
#> 1     cluster_9      17     0.01133333       0.5761626      0.7880813
#> 2     cluster_7      36     0.02400000       0.4048276      0.7024138
#> 3     cluster_8      36     0.02400000       0.4395448      0.7197724
#> 4     cluster_6     131     0.08733333       0.3224888      0.6612444
#> 5     cluster_5     136     0.09066667       0.4203340      0.7101670
#> 6     cluster_4     209     0.13933333       0.4784576      0.7392288
#> 7     cluster_3     221     0.14733333       0.1120699      0.5560349
#> 8     cluster_1     224     0.14933333       0.1718597      0.5859299
#> 9     cluster_2     237     0.15800000       0.3980331      0.6990165
#> 10    cluster_0     253     0.16866667       0.2395322      0.6197661
#>    isolated_label
#> 1            TRUE
#> 2            TRUE
#> 3            TRUE
#> 4           FALSE
#> 5           FALSE
#> 6           FALSE
#> 7           FALSE
#> 8           FALSE
#> 9           FALSE
#> 10          FALSE
```

The artificial label above is deliberately simple. In a real workflow,
replace it with cell-type annotations, reference-mapping labels, or
curated manual labels.

## Quantify residual batch signal

PCR batch scoring measures how much variation in the embedding is
explained by batch. Supplying a baseline reduction lets you compare
correction against the pre-integration PCA space.

``` r

pcr <- sn_calculate_pcr_batch(
  pbmc,
  batch = "library",
  reduction = "harmony",
  baseline_reduction = "pca",
  dims = 1:15,
  max_cells = 1500
)

pcr
#>   reduction baseline_reduction batch_column n_cells batch_variance
#> 1   harmony                pca      library    1500   0.0002694834
#>   baseline_batch_variance pcr_improvement scaled_score
#> 1            0.0002895414    2.005803e-05   0.06927517
```

When several metadata variables might drive the same residual structure,
rank them directly with variance-explained scoring. Use the
single-variable mode for a quick screen, then the partial mode when the
design is not fully confounded.

``` r

pbmc$study <- paste0("study_", pbmc$library)
pbmc$tissue <- "blood"
pbmc$sample <- rep(paste0("donor_", seq_len(4)), length.out = ncol(pbmc))

variance_drivers <- sn_calculate_variance_explained(
  pbmc,
  variables = c("library", "study", "tissue", "sample"),
  reduction = "pca",
  dims = 1:15,
  max_cells = 1500
)

variance_drivers
#>   variable reduction method n_cells n_dims n_levels variance_explained
#> 1   sample       pca single    1500     15        4       0.0011152765
#> 2  library       pca single    1500     15        2       0.0002026187
#> 3    study       pca single    1500     15        2       0.0002026187
#> 4   tissue       pca single    1500     15        1       0.0000000000
#>   mean_dim_variance_explained max_dim_variance_explained
#> 1                 0.002004462                0.006997445
#> 2                 0.000441611                0.002234708
#> 3                 0.000441611                0.002234708
#> 4                 0.000000000                0.000000000
```

## Flag difficult groups

[`sn_identify_challenging_groups()`](https://songqi.org/shennong/dev/reference/sn_identify_challenging_groups.md)
combines rarity, connectivity, and local mixing diagnostics to surface
groups that deserve manual review. This is useful before final
annotation.

``` r

challenging <- sn_identify_challenging_groups(
  pbmc,
  group_by = "seurat_clusters",
  reduction = "harmony",
  dims = 1:15,
  max_cells = 1500
)

challenging
#>    seurat_clusters n_cells fraction_cells median_neighbor_purity
#> 1                3     221     0.14733333              1.0000000
#> 2                1     224     0.14933333              1.0000000
#> 3                0     253     0.16866667              1.0000000
#> 4                6     131     0.08733333              1.0000000
#> 5                2     237     0.15800000              1.0000000
#> 6                7      36     0.02400000              1.0000000
#> 7                5     136     0.09066667              1.0000000
#> 8                8      36     0.02400000              1.0000000
#> 9                9      17     0.01133333              0.9411765
#> 10               4     209     0.13933333              1.0000000
#>    mean_neighbor_purity graph_connectivity mean_silhouette
#> 1             0.9154503                  1       0.1120699
#> 2             0.8730401                  1       0.1718597
#> 3             0.8912202                  1       0.2395322
#> 4             0.9481837                  1       0.3224888
#> 5             0.9660129                  1       0.3980331
#> 6             0.9466419                  1       0.4048276
#> 7             0.9564398                  1       0.4203340
#> 8             0.9617983                  1       0.4395448
#> 9             0.9091375                  1       0.5761626
#> 10            0.9938152                  1       0.4784576
#>    separation_score challenge_score rare_group challenging_group
#> 1         0.8520116      0.14798835      FALSE             FALSE
#> 2         0.8619766      0.13802338      FALSE             FALSE
#> 3         0.8732554      0.12674464      FALSE             FALSE
#> 4         0.8870815      0.11291853      FALSE             FALSE
#> 5         0.8996722      0.10032782      FALSE             FALSE
#> 6         0.9008046      0.09919540       TRUE             FALSE
#> 7         0.9033890      0.09661099      FALSE             FALSE
#> 8         0.9065908      0.09340920       TRUE             FALSE
#> 9         0.9097526      0.09024741       TRUE             FALSE
#> 10        0.9130763      0.08692373      FALSE             FALSE
```

ROGUE is an optional cluster-purity diagnostic that depends on the
upstream `ROGUE` package. Keep it as a targeted check rather than a
mandatory build step.

``` r

rogue_tbl <- sn_calculate_rogue(
  pbmc,
  cluster_by = "seurat_clusters",
  sample_by = "library",
  max_cells = 1500
)

rogue_tbl
```

The intended workflow is: use the aggregate table for a quick read, then
drill into single metrics only when a value changes a downstream
decision.
