# List bundled Shennong signatures

List bundled Shennong signatures

## Usage

``` r
sn_list_signatures(species = NULL, include_groups = FALSE)
```

## Arguments

- species:

  Optional species filter. Use `NULL` to return all species.

- include_groups:

  If `TRUE`, include non-leaf group nodes from the signature tree.

## Value

A tibble with the available signature paths, node kinds, and gene
counts.

## Examples

``` r
sn_list_signatures(species = "human")
sn_list_signatures(species = "human", include_groups = TRUE)
```
