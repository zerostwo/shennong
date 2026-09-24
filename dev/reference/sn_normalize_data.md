# Normalize data in a Seurat object

This function provides a unified normalization entry point for
Seurat-style log-normalization, scran normalization, and SCTransform.
BPCells-backed inputs remain attached to the returned object. When scran
cluster assignments are supplied, size factors are computed in bounded
sparse chunks and the normalized data layer remains BPCells-backed.
Automatic clustering still requires one sparse in-memory materialization
because `scran::quickCluster()` does not accept BPCells iterable
matrices.

## Usage

``` r
sn_normalize_data(
  object,
  method = c("seurat", "scran", "sctransform", "sct"),
  clusters = NULL,
  assay = "RNA",
  layer = "counts",
  ...
)
```

## Arguments

- object:

  A `Seurat` object.

- method:

  One of `"seurat"`, `"scran"`, or `"sctransform"` (alias `"sct"`).

- clusters:

  Optional cluster assignments for scran. Supply either one value per
  cell or the name of a metadata column in `object`. Supplying
  assignments enables bounded-memory processing for BPCells-backed
  inputs.

- assay:

  Assay used for normalization. Defaults to `"RNA"`.

- layer:

  Layer used as the input count matrix. Defaults to `"counts"`.

- ...:

  Additional method-specific arguments passed to
  [`Seurat::NormalizeData()`](https://satijalab.org/seurat/reference/NormalizeData.html),
  `scran::computeSumFactors()`, or
  [`Seurat::SCTransform()`](https://satijalab.org/seurat/reference/SCTransform.html).

## Value

A `Seurat` object with normalized data stored according to the chosen
method.

## Examples

``` r
if (FALSE) { # \dontrun{
seurat_obj <- sn_normalize_data(seurat_obj, method = "scran")
} # }
```
