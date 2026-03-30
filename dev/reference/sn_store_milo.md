# Store a miloR differential-abundance result on a Seurat object

Store a miloR differential-abundance result on a Seurat object

## Usage

``` r
sn_store_milo(
  object,
  result,
  store_name = "default",
  sample_col,
  group_col,
  comparison = NULL,
  reduction = "pca",
  dims = NULL,
  annotation_col = NULL,
  return_object = TRUE
)
```

## Arguments

- object:

  A `Seurat` object.

- result:

  A neighborhood-level differential-abundance table.

- store_name:

  Name used under `object@misc$milo_results`.

- sample_col:

  Sample column used for the design.

- group_col:

  Group column used for the design.

- comparison:

  Human-readable comparison label.

- reduction:

  Reduction used to build the milo neighborhoods.

- dims:

  Optional embedding dimension names or indices used for milo.

- annotation_col:

  Optional neighborhood annotation column.

- return_object:

  If `TRUE`, return the updated object.

## Value

A `Seurat` object or stored-result list.
