# Store a miloR differential-abundance result on a Seurat object

Store a miloR differential-abundance result on a Seurat object

## Usage

``` r
sn_store_milo(
  object,
  result,
  result_id = "default",
  sample_by = NULL,
  group_by = NULL,
  comparison = NULL,
  reduction = "pca",
  dims = NULL,
  annotation_by = NULL,
  random_seed = NULL,
  return_object = TRUE,
  overwrite = FALSE
)
```

## Arguments

- object:

  A `Seurat` object.

- result:

  A neighborhood-level differential-abundance table.

- result_id:

  Stable identifier for the stored Milo result.

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

- random_seed:

  Optional random seed recorded in the result provenance.

- return_object:

  If `TRUE`, return the updated object.

- overwrite:

  Explicitly replace a stored result with the same ID.

## Value

A `Seurat` object or stored-result list.
