# Retrieve blocklist genes from SignatuR by preset gene set names

This function uses the `SignatuR` package to fetch blocklist (or
related) genes (e.g. ribosomal, mitochondrial, heatshock, etc.). It then
checks those symbols via `HGNChelper` to remove invalid or outdated
entries.

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

A unique character vector of valid gene symbols.

## Examples

``` r
if (FALSE) { # \dontrun{
block_genes <- sn_get_signatures(c("mito", "ribo", "tcr", "immunoglobulins", "pseudogenes", "noncoding"))
} # }
```
