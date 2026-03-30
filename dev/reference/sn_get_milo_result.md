# Retrieve a stored miloR result from a Seurat object

Retrieve a stored miloR result from a Seurat object

## Usage

``` r
sn_get_milo_result(
  object,
  milo_name = "default",
  annotation = NULL,
  spatial_fdr = NULL,
  with_metadata = FALSE
)
```

## Arguments

- object:

  A `Seurat` object.

- milo_name:

  Name of the stored milo result.

- annotation:

  Optional subset of annotation labels to keep.

- spatial_fdr:

  Optional maximum `SpatialFDR` threshold.

- with_metadata:

  If `TRUE`, return the full stored-result list.

## Value

A tibble or stored-result list.
