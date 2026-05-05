# Calculate cluster_by entropy for a categorical label

Cluster entropy measures how mixed a categorical label_by is within each
cluster. When used with batch labels, higher normalized entropy
indicates stronger within-cluster batch mixing.

## Usage

``` r
sn_calculate_cluster_entropy(
  x,
  cluster_by = NULL,
  label_by = NULL,
  cluster = NULL,
  label = NULL
)
```

## Arguments

- x:

  A Seurat object or data frame containing the required columns.

- cluster_by:

  Metadata/data-frame column containing cluster_by labels.

- label_by:

  Metadata/data-frame column containing the label_by to evaluate within
  each cluster.

- cluster:

  Deprecated alias for `cluster_by`.

- label:

  Deprecated alias for `label_by`.

## Value

A data frame with one row per cluster, including raw entropy and a
normalized entropy score in `[0, 1]`.

## Examples

``` r
meta <- data.frame(
  cluster = c("0", "0", "1", "1"),
  batch = c("a", "b", "a", "b")
)
sn_calculate_cluster_entropy(meta, cluster_by = "cluster", label_by = "batch")
#>   cluster n_cells n_labels dominant_label   entropy normalized_entropy
#> 1       0       2        2              a 0.6931472                  1
#> 2       1       2        2              a 0.6931472                  1
```
