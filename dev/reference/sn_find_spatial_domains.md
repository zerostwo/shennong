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
  result_id = "spatial_domains",
  backend_control = list(),
  return_object = TRUE,
  seed = NULL,
  verbose = TRUE,
  sample_by = NULL
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

- result_id:

  Stored result name.

- backend_control:

  Backend controls or an explicit `runner`/`result`.

- return_object:

  Return the modified object or result.

- seed:

  Top-level reproducibility seed. Precedence: `seed` \>
  `backend_control$seed` \> task default; stamped into provenance.

- verbose:

  Top-level progress switch forwarded through `backend_control$verbose`
  when explicitly supplied.

- sample_by:

  Optional metadata column defining independent tissue sections. The
  built-in BANKSY path fails closed for pooled sections.

## Value

A Seurat object or spatial-domain result.
