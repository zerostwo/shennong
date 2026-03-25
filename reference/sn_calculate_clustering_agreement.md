# Calculate agreement between clusters and reference labels

The returned table includes both adjusted Rand index (ARI) and
normalized mutual information (NMI), which are commonly used to quantify
how well clustering preserves known cell identities.

## Usage

``` r
sn_calculate_clustering_agreement(x, cluster, label)
```

## Arguments

- x:

  A Seurat object or data frame containing the required columns.

- cluster:

  Metadata/data-frame column containing cluster labels.

- label:

  Metadata/data-frame column containing reference labels.

## Value

A one-row data frame with ARI and NMI.

## Examples

``` r
meta <- data.frame(
  cluster = c("T", "T", "B", "B"),
  label = c("T", "T", "B", "B")
)
sn_calculate_clustering_agreement(meta, cluster = "cluster", label = "label")
#>   cluster_column label_column n_cells ari nmi
#> 1        cluster        label       4   1   1
```
