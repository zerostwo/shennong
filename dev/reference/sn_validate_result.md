# Validate a Shennong analysis result

Validate a Shennong analysis result

## Usage

``` r
sn_validate_result(result, error = TRUE)
```

## Arguments

- result:

  A result list following the Shennong analysis-result contract.

- error:

  If `TRUE`, stop when validation fails. If `FALSE`, return the
  validation report without stopping.

## Value

A validation report with `valid`, `errors`, and `warnings` fields.
Successful reports are returned invisibly when `error = TRUE`.

## Examples

``` r
result <- list(
  schema_version = "1.0", analysis_type = "demo", name = "example",
  method = "mean", backend = "base", input = list(), parameters = list(),
  tables = list(primary = data.frame(value = 1)), embeddings = list(),
  graphs = list(), models = list(), diagnostics = list(), warnings = character(),
  provenance = list(package_versions = list(), random_seed = 1L, timestamp = "2026-01-01 UTC")
)
sn_validate_result(result, error = FALSE)
#> $valid
#> [1] TRUE
#> 
#> $errors
#> character(0)
#> 
#> $warnings
#> character(0)
#> 
#> attr(,"class")
#> [1] "sn_result_validation" "list"                
```
