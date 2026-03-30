# Calculate cluster entropy for a categorical label

Cluster entropy measures how mixed a categorical label is within each
cluster. When used with batch labels, higher normalized entropy
indicates stronger within-cluster batch mixing.

## Usage

``` r
sn_calculate_cluster_entropy(x, cluster, label)
```

## Arguments

- x:

  A Seurat object or data frame containing the required columns.

- cluster:

  Metadata/data-frame column containing cluster labels.

- label:

  Metadata/data-frame column containing the label to evaluate within
  each cluster.

## Value

A data frame with one row per cluster, including raw entropy and a
normalized entropy score in `[0, 1]`.

## Examples

``` r
meta <- data.frame(
  cluster = c("0", "0", "1", "1"),
  batch = c("a", "b", "a", "b")
)
sn_calculate_cluster_entropy(meta, cluster = "cluster", label = "batch")
#>   cluster n_cells n_labels dominant_label   entropy normalized_entropy
#> 1       0       2        2              a 0.6931472                  1
#> 2       1       2        2              a 0.6931472                  1
```
