# Identify spatial domains

Identify spatial domains

## Usage

``` r
sn_find_spatial_domains(
  object,
  method = c("banksy", "stlearn", "bayesspace", "cellcharter"),
  spatial_cols = NULL,
  assay = NULL,
  layer = "counts",
  store_name = "spatial_domains",
  backend_control = list(),
  return_object = TRUE
)
```

## Arguments

- object:

  A Seurat object with coordinates.

- method:

  BANKSY or an explicit stLearn/BayesSpace/CellCharter adapter.

- spatial_cols:

  Coordinate metadata columns.

- assay, layer:

  Expression assay and layer.

- store_name:

  Stored result name.

- backend_control:

  Backend controls or an explicit `runner`/`result`.

- return_object:

  Return the modified object or result.

## Value

A Seurat object or spatial-domain result.
