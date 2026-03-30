# Retrieve a stored deconvolution result from a Seurat object

Retrieve a stored deconvolution result from a Seurat object

## Usage

``` r
sn_get_deconvolution_result(
  object,
  deconvolution_name = "default",
  samples = NULL,
  cell_types = NULL,
  with_metadata = FALSE
)
```

## Arguments

- object:

  A `Seurat` object.

- deconvolution_name:

  Name of the stored result.

- samples:

  Optional subset of bulk samples to keep.

- cell_types:

  Optional subset of cell types to keep.

- with_metadata:

  If `TRUE`, return the full stored-result list.

## Value

A tibble or stored-result list.
