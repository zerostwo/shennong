# Set the storage backend for Seurat assay layers

`sn_set_layer_backend()` switches selected Seurat assay layers between
on-disk BPCells matrices and in-memory sparse matrices without changing
the assay or layer names. Use `backend = "bpcells"` to write and rebind
layers, or `backend = "memory"` to materialize them as `dgCMatrix`
objects.

## Usage

``` r
sn_set_layer_backend(
  object,
  backend = c("bpcells", "memory"),
  directory = NULL,
  assays = NULL,
  layers = "counts",
  matrix_type = c("auto", "uint32_t", "double", "float"),
  overwrite = FALSE,
  verbose = TRUE
)
```

## Arguments

- object:

  A Seurat object.

- backend:

  Target storage backend: `"bpcells"` or `"memory"`.

- directory:

  Output directory containing one BPCells matrix directory per selected
  assay/layer. Required for `backend = "bpcells"` and ignored for
  `backend = "memory"`.

- assays:

  Assays to convert. Defaults to all assays.

- layers:

  Layers to convert within each assay. Defaults to `"counts"`. Use
  `NULL` to convert all layers in each selected assay.

- matrix_type:

  BPCells storage type. `"auto"` stores compatible count-like layers as
  `"uint32_t"` and otherwise preserves their numeric representation.
  Explicit alternatives are `"uint32_t"`, `"double"`, and `"float"`.
  Used only for `backend = "bpcells"`.

- overwrite:

  Logical; replace existing BPCells matrix directories.

- verbose:

  Whether to print progress messages.

## Value

A Seurat object with the selected layers rebound to the requested
storage backend.

## Details

BPCells matrix directories remain external to a serialized Seurat
object. Keep them alongside the object when moving an analysis.
Materializing a layer with `backend = "memory"` does not delete its
BPCells directory.

## Examples

``` r
if (FALSE) { # \dontrun{
pbmc <- sn_set_layer_backend(
  pbmc,
  backend = "bpcells",
  directory = "data/processed/pbmc_bpcells",
  layers = c("counts", "data")
)
pbmc <- sn_set_layer_backend(
  pbmc,
  backend = "memory",
  layers = "counts"
)
} # }
```
