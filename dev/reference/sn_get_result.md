# Retrieve a stored Shennong analysis result

Retrieve a stored Shennong analysis result

## Usage

``` r
sn_get_result(object, type, name)
```

## Arguments

- object:

  A `Seurat` object.

- type:

  Analysis type.

- name:

  Stored result name.

## Value

A validated Shennong analysis-result list.

## Examples

``` r
if (FALSE) sn_get_result(obj, "trajectory", "cd8_slingshot") # \dontrun{}
```
