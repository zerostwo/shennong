# Retrieve stored regulatory activity results

Retrieve stored regulatory activity results

## Usage

``` r
sn_get_regulatory_activity_result(
  object,
  result_id = NULL,
  sources = NULL,
  conditions = NULL
)
```

## Arguments

- object:

  A Seurat object.

- result_id:

  Name of the stored result.

- sources:

  Optional TF or pathway names to keep.

- conditions:

  Optional cell or group names to keep.

## Value

A filtered tibble.

## Details

This getter always returns the selected table. Use
[`sn_get_result()`](https://zerostwo.github.io/shennong/dev/reference/sn_get_result.md)
for the complete stored result and metadata.
