# Calculate cluster purity against a reference label

Cluster purity summarizes how homogeneous each cluster is with respect
to a reference label. This is useful for checking whether clustering
preserves known cell identities after integration.

## Usage

``` r
sn_calculate_cluster_purity(x, cluster, label)
```

## Arguments

- x:

  A Seurat object or data frame containing the required columns.

- cluster:

  Metadata/data-frame column containing cluster labels.

- label:

  Metadata/data-frame column containing reference labels.

## Value

A data frame with one row per cluster and purity diagnostics in
`[0, 1]`.

## Examples

``` r
meta <- data.frame(
  cluster = c("T", "T", "B", "B"),
  label = c("T", "T", "B", "B")
)
sn_calculate_cluster_purity(meta, cluster = "cluster", label = "label")
#>   cluster n_cells dominant_label dominant_label_n purity_score impurity_score
#> 1       T       2              T                2            1              0
#> 2       B       2              B                2            1              0
```
