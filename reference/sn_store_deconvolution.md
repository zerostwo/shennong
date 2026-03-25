# Store a deconvolution result on a Seurat object

Store a deconvolution result on a Seurat object

## Usage

``` r
sn_store_deconvolution(
  object,
  result,
  store_name = "default",
  method = "bayesprism",
  bulk_samples = NULL,
  reference_label = NULL,
  artifacts = NULL,
  return_object = TRUE
)
```

## Arguments

- object:

  A `Seurat` object.

- result:

  A deconvolution table.

- store_name:

  Name used under `object@misc$deconvolution_results`.

- method:

  Deconvolution backend, for example `"bayesprism"`.

- bulk_samples:

  Optional bulk sample identifiers.

- reference_label:

  Metadata column or label set used as the reference.

- artifacts:

  Optional backend-specific artifacts or file paths.

- return_object:

  If `TRUE`, return the updated object.

## Value

A `Seurat` object or stored-result list.
