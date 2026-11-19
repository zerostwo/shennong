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
```
