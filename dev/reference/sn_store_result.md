# Store a Shennong analysis result on a Seurat object

Every result is stored at
`object@misc$shennong$results[[analysis_type]][[result_id]]`.

## Usage

``` r
sn_store_result(object, type, result_id, result, overwrite = FALSE)
```

## Arguments

- object:

  A `Seurat` object.

- type:

  Analysis type, for example `"trajectory"` or `"de"`.

- result_id:

  Stable identifier used to store and retrieve the result.

- result:

  A result list. Missing contract fields are filled when they can be
  inferred without changing the analytical content.

- overwrite:

  Replace an existing result with the same type and ID. Defaults to
  `FALSE`; choose a new ID to retain both analyses.

## Value

The modified `Seurat` object.

## Examples

``` r
if (FALSE) { # \dontrun{
obj <- sn_store_result(obj, "trajectory", "cd8_slingshot", result)
sn_get_result(obj, "trajectory", "cd8_slingshot")
} # }
```
