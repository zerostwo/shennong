# Retrieve a stored deconvolution result from a Seurat object

Retrieve a stored deconvolution result from a Seurat object

## Usage

``` r
sn_get_deconvolution_result(
  object,
  result_id = NULL,
  samples = NULL,
  cell_types = NULL
)
```

## Arguments

- object:

  A `Seurat` object or a unified result of this analysis type.

- result_id:

  Name of the stored result.

- samples:

  Optional subset of bulk samples to keep.

- cell_types:

  Optional subset of cell types to keep.

## Value

A filtered tibble.

## Details

This getter always returns the selected table. Use
[`sn_get_result()`](https://zerostwo.github.io/shennong/dev/reference/sn_get_result.md)
for the complete stored result and metadata.
