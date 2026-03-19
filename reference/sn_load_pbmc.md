# Load PBMC example datasets

Backward-compatible wrapper around
[`sn_load_data()`](https://songqi.org/shennong/reference/sn_load_data.md).

## Usage

``` r
sn_load_pbmc(...)
```

## Arguments

- ...:

  Passed to
  [`sn_load_data()`](https://songqi.org/shennong/reference/sn_load_data.md).

## Value

See
[`sn_load_data()`](https://songqi.org/shennong/reference/sn_load_data.md).

## Examples

``` r
if (FALSE) { # \dontrun{
pbmc <- sn_load_pbmc(dataset = "pbmc1k")
raw_counts <- sn_load_pbmc(dataset = "pbmc1k", matrix_type = "raw")
} # }
```
