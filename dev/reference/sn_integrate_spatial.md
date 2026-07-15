# Integrate spatial samples with an explicit backend adapter

Integrate spatial samples with an explicit backend adapter

## Usage

``` r
sn_integrate_spatial(
  object,
  method = c("staligner", "harmony", "custom"),
  spatial_cols = NULL,
  store_name = "spatial_integration",
  backend_control = list(),
  return_object = TRUE
)
```

## Arguments

- object:

  A Seurat object.

- method:

  Integration backend label.

- spatial_cols:

  Coordinate columns.

- store_name:

  Stored result name.

- backend_control:

  Required `runner` or `result` returning an embedding table with a
  `cell` column.

- return_object:

  Return the modified object or result.

## Value

A Seurat object or spatial-integration result.
