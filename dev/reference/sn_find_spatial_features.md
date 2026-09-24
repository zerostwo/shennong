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
  result_id = "spatial_features",
  backend_control = list(),
  return_object = TRUE,
  seed = NULL,
  verbose = TRUE,
  sample_by = NULL
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

- result_id:

  Stored result name.

- backend_control:

  Method controls or an explicit `runner`/`result`.

- return_object:

  Return the modified object or result.

- seed:

  Top-level reproducibility seed. Precedence: `seed` \>
  `backend_control$seed` \> task default; stamped into provenance.

- verbose:

  Top-level progress switch forwarded through `backend_control$verbose`
  when explicitly supplied.

- sample_by:

  Optional metadata column defining independent samples or tissue
  sections. Spatial graphs and permutations are restricted within these
  boundaries.

## Value

A Seurat object or unified spatial-feature result.
