# Find spatially variable features

Find spatially variable features

## Usage

``` r
sn_find_spatial_features(
  object,
  method = c("morans_i", "nnsvg", "sparkx"),
  spatial_cols = NULL,
  assay = NULL,
  layer = "data",
  features = NULL,
  store_name = "spatial_features",
  backend_control = list(),
  return_object = TRUE
)
```

## Arguments

- object:

  A Seurat object with coordinate metadata.

- method:

  Moran's I, nnSVG, or an explicit SPARK-X adapter.

- spatial_cols:

  Coordinate metadata columns.

- assay, layer:

  Expression assay and layer.

- features:

  Features to test.

- store_name:

  Stored result name.

- backend_control:

  Method controls or an explicit `runner`/`result`.

- return_object:

  Return the modified object or result.

## Value

A Seurat object or unified spatial-feature result.
