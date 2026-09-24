# Retrieve a stored Shennong analysis result

Retrieve a stored Shennong analysis result

## Usage

``` r
sn_get_result(object, type, result_id = NULL)
```

## Arguments

- object:

  A `Seurat` object.

- type:

  Analysis type.

- result_id:

  Stored result identifier. Omit only when exactly one result of the
  requested type exists; ambiguous choices are reported as an error.

## Value

A validated Shennong analysis-result list.

## Examples

``` r
if (FALSE) sn_get_result(obj, "trajectory", "cd8_slingshot") # \dontrun{}
```
