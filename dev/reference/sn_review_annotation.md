# Review stored annotation evidence and low-confidence labels

Review stored annotation evidence and low-confidence labels

## Usage

``` r
sn_review_annotation(x, store_name = "annotation", low_confidence_only = TRUE)
```

## Arguments

- x:

  A Seurat object or unified annotation result.

- store_name:

  Annotation result name when `x` is a Seurat object.

- low_confidence_only:

  Return only low-confidence rows in the review tables.

## Value

A list containing cell, cluster, evidence, and diagnostic tables.

## Examples

``` r
if (FALSE) sn_review_annotation(object, "annotation") # \dontrun{}
```
