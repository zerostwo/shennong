# Retrieve a stored cell-cell communication result

Retrieve a stored cell-cell communication result

## Usage

``` r
sn_get_cell_communication_result(
  object,
  result_id = NULL,
  sources = NULL,
  targets = NULL
)
```

## Arguments

- object:

  A Seurat object.

- result_id:

  Name of the stored result.

- sources, targets:

  Optional source/target labels to keep.

## Value

A filtered tibble.

## Details

This getter always returns the selected table. Use
[`sn_get_result()`](https://zerostwo.github.io/shennong/dev/reference/sn_get_result.md)
for the complete stored result and metadata.
