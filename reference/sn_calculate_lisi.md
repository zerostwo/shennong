# Calculate LISI score

This function calculates the Local Intrinsic Dimensionality-based
Outlier score (Lisi score) for a Seurat object.

## Usage

``` r
sn_calculate_lisi(x, reduction = "pca", label = "sample")
```

## Arguments

- x:

  A Seurat object.

- reduction:

  The dimensionality reduction method used to generate the embeddings
  (default is "pca").

- label:

  The column name in object@meta.data that specifies the sample labels
  (default is "sample").

## Value

A data frame with the Lisi score for each cell, along with the cell ID.

## Examples

``` r
# Calculate Lisi score for a Seurat object
```
