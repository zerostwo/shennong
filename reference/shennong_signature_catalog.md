# Bundled Shennong Signature Catalog

A package-owned, tree-structured signature catalog built from the
upstream SignatuR hierarchy and stored as the current
Shennong-controlled snapshot under `data/`.

## Usage

``` r
shennong_signature_catalog
```

## Format

A named list with two top-level entries:

- `metadata`:

  Build metadata such as source package, source version, and build date.

- `tree`:

  A nested tree for `human` and `mouse`, where each node records its
  `name`, `kind`, optional `genes`, and optional `children`.

## Source

Imported from SignatuR during development and built into package data
with `data-raw/build_shennong_signatures.R`.
