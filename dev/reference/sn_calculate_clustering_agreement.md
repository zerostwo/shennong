# Calculate agreement between clusters and reference labels

The returned table includes both adjusted Rand index (ARI) and
normalized mutual information (NMI), which are commonly used to quantify
how well clustering preserves known cell identities.

## Usage

``` r
sn_calculate_clustering_agreement(
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

  Metadata/data-frame column containing reference labels.

- cluster:

  Deprecated alias for `cluster_by`.

- label:

  Deprecated alias for `label_by`.

## Value

A one-row data frame with ARI and NMI.

## Examples

``` r
meta <- data.frame(
  cluster = c("T", "T", "B", "B"),
  label = c("T", "T", "B", "B")
)
sn_calculate_clustering_agreement(meta, cluster_by = "cluster", label_by = "label")
#>   cluster_column label_column n_cells ari nmi
#> 1        cluster        label       4   1   1
```
