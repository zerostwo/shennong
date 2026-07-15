# Store a Shennong analysis result on a Seurat object

Registered legacy result types are stored in their established
`object@misc` collection. New result types use the generic
`object@misc$analysis_results` collection.

## Usage

``` r
sn_store_result(object, type, name, result)
```

## Arguments

- object:

  A `Seurat` object.

- type:

  Analysis type, for example `"trajectory"` or `"de"`.

- name:

  Stable name used to retrieve the result.

- result:

  A result list. Missing contract fields are filled when they can be
  inferred without changing the analytical content.

## Value

The modified `Seurat` object.

## Examples

``` r
if (FALSE) { # \dontrun{
obj <- sn_store_result(obj, "trajectory", "cd8_slingshot", result)
sn_get_result(obj, "trajectory", "cd8_slingshot")
} # }
```
