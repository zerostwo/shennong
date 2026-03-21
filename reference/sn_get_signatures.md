# Retrieve bundled Shennong signature genes by preset category names

Shennong ships a package-owned snapshot of the signature categories
required by its core workflows. The snapshot is generated from SignatuR
during package development and stored inside the package so runtime
behavior is stable and does not depend on the user's installed SignatuR
version.

## Usage

``` r
sn_get_signatures(species = "human", category = NULL)
```

## Arguments

- species:

  One of `"human"` or `"mouse"`.

- category:

  A character vector of predefined signature categories.

## Value

A unique character vector of signature gene symbols.

## Examples

``` r
block_genes <- sn_get_signatures(
  species = "human",
  category = c("mito", "ribo", "tcr", "immunoglobulins")
)
```
