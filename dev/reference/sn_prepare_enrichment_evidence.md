# Prepare enrichment evidence

Prepare enrichment evidence

## Usage

``` r
sn_prepare_enrichment_evidence(
  object = NULL,
  enrichment_name = NULL,
  result = NULL,
  n_terms = 10
)
```

## Arguments

- object:

  Optional `Seurat` object containing stored enrichment results.

- enrichment_name:

  Name of a stored enrichment result.

- result:

  Optional enrichment result object supplied directly.

- n_terms:

  Number of top terms to keep.

## Value

A structured list ready for prompt construction.

## Examples

``` r
enrich_tbl <- tibble::tibble(
  ID = c("GO:0001", "GO:0002"),
  Description = c("immune response", "lymphocyte activation"),
  NES = c(2.1, 1.7),
  p.adjust = c(0.01, 0.03)
)
evidence <- sn_prepare_enrichment_evidence(result = enrich_tbl, n_terms = 1)
evidence$top_terms
#> # A tibble: 1 × 4
#>   ID      Description       NES p.adjust
#>   <chr>   <chr>           <dbl>    <dbl>
#> 1 GO:0001 immune response   2.1     0.01
```
