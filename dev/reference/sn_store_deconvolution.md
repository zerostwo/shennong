# Store a deconvolution result on a Seurat object

Store a deconvolution result on a Seurat object

## Usage

``` r
sn_store_deconvolution(
  object,
  result,
  result_id = "default",
  method = "bayesprism",
  bulk_samples = NULL,
  reference_label = NULL,
  artifacts = NULL,
  random_seed = NULL,
  return_object = TRUE
)
```

## Arguments

- object:

  A `Seurat` object.

- result:

  A deconvolution table.

- result_id:

  Stable identifier for the stored deconvolution result.

- method:

  Deconvolution backend, for example `"bayesprism"`.

- bulk_samples:

  Optional bulk sample identifiers.

- reference_label:

  Metadata column or label_by set used as the reference.

- artifacts:

  Optional backend-specific artifacts or file paths.

- random_seed:

  Optional random seed recorded in the result provenance.

- return_object:

  If `TRUE`, return the updated object.

## Value

A `Seurat` object or stored-result list.
