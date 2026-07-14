# Convert Seurat assay layers to BPCells-backed matrices

`sn_convert_bpcells()` writes selected Seurat assay layers to BPCells
matrix directories and rebinds those layers in the returned Seurat
object. This keeps large count or normalized-expression layers on disk
while preserving the usual Seurat object interface.

## Usage

``` r
sn_convert_bpcells(
  object,
  directory,
  assays = NULL,
  layers = "counts",
  overwrite = FALSE,
  verbose = TRUE
)
```

## Arguments

- object:

  A Seurat object.

- directory:

  Output directory that will contain one BPCells matrix directory per
  selected assay/layer.

- assays:

  Assays to convert. Defaults to all assays.

- layers:

  Layers to convert within each assay. Defaults to `"counts"`. Use
  `NULL` to convert all layers in each selected assay.

- overwrite:

  Logical; overwrite existing BPCells matrix directories.

- verbose:

  Whether to print progress messages.

## Value

A Seurat object with selected layers backed by BPCells matrices.

## Details

BPCells stores matrix directories outside the serialized Seurat object.
If the object is moved to another machine, move the BPCells directory
alongside it and rebind the layers with `BPCells::open_matrix_dir()`
when needed.

## Examples

``` r
if (FALSE) { # \dontrun{
pbmc <- sn_convert_bpcells(
  pbmc,
  directory = "data/processed/pbmc_bpcells",
  layers = c("counts", "data"),
  overwrite = TRUE
)
sn_write(pbmc, "data/processed/pbmc_bpcells_bound.qs2")
} # }
```
