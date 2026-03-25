# Small Built-In PBMC Raw Counts Matrix

A compact raw count matrix matched to `pbmc_small`. It includes the
sampled filtered barcodes plus additional raw-only droplets for
ambient-RNA and preprocessing examples.

## Usage

``` r
pbmc_small_raw
```

## Format

A sparse gene-by-barcode matrix with raw barcodes.

## Source

Derived during development from the Zenodo-backed `pbmc1k` and `pbmc3k`
example datasets via `data-raw/build_shennong_example_data.R`.
