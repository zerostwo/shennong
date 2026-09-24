# Convert Seurat assay layers to BPCells-backed matrices

`sn_convert_bpcells()` is the compatibility wrapper for
`sn_set_layer_backend(backend = "bpcells")`. New code can use
[`sn_set_layer_backend()`](https://zerostwo.github.io/shennong/dev/reference/sn_set_layer_backend.md)
when layers need to move in either direction.

## Usage

``` r
sn_convert_bpcells(
  object,
  directory,
  assays = NULL,
  layers = "counts",
  overwrite = FALSE,
  verbose = TRUE,
  matrix_type = c("auto", "uint32_t", "double", "float")
)
```

## Arguments

- object:

  A Seurat object.

- directory:

  Output directory containing one BPCells matrix directory per selected
  assay/layer. Required for `backend = "bpcells"` and ignored for
  `backend = "memory"`.

- assays:

  Assays to convert. Defaults to all assays.

- layers:

  Layers to convert within each assay. Defaults to `"counts"`. Use
  `NULL` to convert all layers in each selected assay.

- overwrite:

  Logical; replace existing BPCells matrix directories.

- verbose:

  Whether to print progress messages.

- matrix_type:

  BPCells storage type. `"auto"` stores compatible count-like layers as
  `"uint32_t"` and otherwise preserves their numeric representation.
  Explicit alternatives are `"uint32_t"`, `"double"`, and `"float"`.
  Used only for `backend = "bpcells"`.

## Value

A Seurat object with selected layers backed by BPCells matrices.

## See also

[`sn_set_layer_backend()`](https://zerostwo.github.io/shennong/dev/reference/sn_set_layer_backend.md)
