# Store a cell-cell communication result on a Seurat object

Store a cell-cell communication result on a Seurat object

## Usage

``` r
sn_store_cell_communication(
  object,
  result,
  store_name = "default",
  method = "cellchat",
  group_by = NULL,
  sender = NULL,
  receiver = NULL,
  species = NULL,
  artifacts = NULL,
  return_object = TRUE
)
```

## Arguments

- object:

  A Seurat object.

- result:

  Communication result table.

- store_name:

  Name used under `object@misc$cell_communication_results`.

- method:

  Communication backend.

- group_by:

  Metadata column used for groups.

- sender, receiver:

  Optional sender/receiver labels.

- species:

  Optional species label.

- artifacts:

  Optional backend-specific artifacts.

- return_object:

  If `TRUE`, return the updated object.

## Value

A Seurat object or stored-result list.
