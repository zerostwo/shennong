# Store regulatory activity results on a Seurat object

Store regulatory activity results on a Seurat object

## Usage

``` r
sn_store_regulatory_activity(
  object,
  result,
  store_name = "default",
  method = "dorothea",
  group_by = NULL,
  species = NULL,
  network = NULL,
  return_object = TRUE
)
```

## Arguments

- object:

  A Seurat object.

- result:

  Regulatory activity table.

- store_name:

  Name used under `object@misc$regulatory_activity_results`.

- method:

  Inference method.

- group_by:

  Optional grouping column.

- species:

  Optional species label.

- network:

  Optional regulatory network used for inference.

- return_object:

  If `TRUE`, return the updated object.

## Value

A Seurat object or stored-result list.
