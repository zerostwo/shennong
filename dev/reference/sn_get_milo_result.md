# Retrieve a stored miloR result from a Seurat object

Retrieve a stored miloR result from a Seurat object

## Usage

``` r
sn_get_milo_result(
  object,
  result_id = NULL,
  annotation = NULL,
  spatial_fdr = NULL
)
```

## Arguments

- object:

  A `Seurat` object or a unified result of this analysis type.

- result_id:

  Name of the stored milo result.

- annotation:

  Optional subset of annotation labels to keep.

- spatial_fdr:

  Optional maximum `SpatialFDR` threshold.

## Value

A filtered tibble.

## Details

This getter always returns the selected table. Use
[`sn_get_result()`](https://zerostwo.github.io/shennong/dev/reference/sn_get_result.md)
for the complete stored result and metadata.
