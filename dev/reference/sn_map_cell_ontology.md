# Map cell labels to Cell Ontology identifiers

Uses a small, versioned Cell Ontology snapshot shipped with Shennong. A
project-specific mapping can be supplied as a JSON file, data frame, or
list with `id`, `label`, and optional `aliases` fields.

## Usage

``` r
sn_map_cell_ontology(labels, ontology = NULL, strict = FALSE)
```

## Arguments

- labels:

  Character vector of cell-type labels.

- ontology:

  Optional custom ontology mapping.

- strict:

  If `TRUE`, fail when any label is unmapped.

## Value

A tibble with input label, ontology identifier, canonical ontology
label, and the alias that matched.

## Examples

``` r
sn_map_cell_ontology(c("B cells", "T cells", "unknown"))
#> # A tibble: 3 × 4
#>   input_label ontology_id ontology_label matched_alias
#>   <chr>       <chr>       <chr>          <chr>        
#> 1 B cells     CL:0000236  B cell         B cells      
#> 2 T cells     CL:0000084  T cell         T cells      
#> 3 unknown     NA          NA             NA           
```
