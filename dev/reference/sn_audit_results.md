# Audit stored Shennong analysis results

Inspect every registered analysis result and artifact, plus unknown
populated top-level `object@misc` entries, without mutating the object.
The audit distinguishes results that satisfy the current contract from
malformed results that can be normalized safely, and reports unknown
payloads as `unregistered`.

## Usage

``` r
sn_audit_results(object, type = NULL, include_artifacts = TRUE)
```

## Arguments

- object:

  A `Seurat` object.

- type:

  Optional analysis type or character vector of types to inspect.

- include_artifacts:

  If `TRUE`, also report registered runtime, cache, and backend
  artifacts that intentionally remain outside the tabular
  analysis-result contract, as well as unknown populated top-level
  `object@misc` entries.

## Value

A tibble with one row per stored result and validation, migration,
schema-version, and canonical-primary-table status.

## Examples

``` r
if (FALSE) sn_audit_results(object) # \dontrun{}
```
