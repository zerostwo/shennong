# Small Built-In PBMC Seurat Object

A compact Seurat object sampled from the package PBMC example assets.
The object combines cells from `pbmc1k` and `pbmc3k` and includes
sample-level metadata, so examples can demonstrate both single-sample
and multi-sample workflows without a network download.

## Usage

``` r
pbmc_small
```

## Format

A `Seurat` object with filtered counts and metadata columns including
`sample`, `source_dataset`, and `source_barcode`.

## Source

Derived during development from the Zenodo-backed `pbmc1k` and `pbmc3k`
example datasets via `data-raw/build_shennong_example_data.R`.
