# Delete a stored Shennong analysis result

Delete a stored Shennong analysis result

## Usage

``` r
sn_delete_result(object, type, name)
```

## Arguments

- object:

  A `Seurat` object.

- type:

  Analysis type.

- name:

  Stored result name.

## Value

The modified `Seurat` object.

## Examples

``` r
if (FALSE) obj <- sn_delete_result(obj, "trajectory", "cd8_slingshot") # \dontrun{}
```
