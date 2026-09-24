# Upgrade stored Shennong analysis results

Normalize stored analysis results to the current canonical result
envelope.

## Usage

``` r
sn_upgrade_results(object, type = NULL, strict = TRUE)
```

## Arguments

- object:

  A `Seurat` object.

- type:

  Optional analysis type or character vector of types to normalize.

- strict:

  If `TRUE`, stop at the first result that cannot be safely normalized.
  If `FALSE`, leave invalid entries unchanged and warn.

## Value

The modified `Seurat` object.

## Examples

``` r
if (FALSE) object <- sn_upgrade_results(object) # \dontrun{}
```
