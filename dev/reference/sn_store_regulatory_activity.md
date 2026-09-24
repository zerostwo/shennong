# Store regulatory activity results on a Seurat object

Store regulatory activity results on a Seurat object

## Usage

``` r
sn_store_regulatory_activity(
  object,
  result,
  result_id = "default",
  method = "dorothea",
  group_by = NULL,
  species = NULL,
  network = NULL,
  random_seed = NULL,
  return_object = TRUE
)
```

## Arguments

- object:

  A Seurat object.

- result:

  Regulatory activity table.

- result_id:

  Name used under the canonical Shennong result registry.

- method:

  Inference method.

- group_by:

  Optional grouping column.

- species:

  Optional species label.

- network:

  Optional regulatory network used for inference.

- random_seed:

  Optional random seed recorded in the result provenance.

- return_object:

  If `TRUE`, return the updated object.

## Value

A Seurat object or stored-result list.
