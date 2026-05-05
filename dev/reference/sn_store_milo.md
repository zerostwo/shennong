# Store a miloR differential-abundance result on a Seurat object

Store a miloR differential-abundance result on a Seurat object

## Usage

``` r
sn_store_milo(
  object,
  result,
  store_name = "default",
  sample_by = NULL,
  group_by = NULL,
  comparison = NULL,
  reduction = "pca",
  dims = NULL,
  annotation_by = NULL,
  return_object = TRUE,
  sample_col = NULL,
  annotation_col = NULL,
  group_col = NULL
)
```

## Arguments

- object:

  A `Seurat` object.

- result:

  A neighborhood-level differential-abundance table.

- store_name:

  Name used under `object@misc$milo_results`.

- sample_by:

  Sample column used for the design.

- group_by:

  Group column used for the design.

- comparison:

  Human-readable comparison label.

- reduction:

  Reduction used to build the milo neighborhoods.

- dims:

  Optional embedding dimension names or indices used for milo.

- annotation_by:

  Optional neighborhood annotation column.

- return_object:

  If `TRUE`, return the updated object.

- sample_col:

  Deprecated alias for `sample_by`.

- annotation_col:

  Deprecated alias for `annotation_by`.

- group_col:

  Deprecated alias for `group_by`.

## Value

A `Seurat` object or stored-result list.
