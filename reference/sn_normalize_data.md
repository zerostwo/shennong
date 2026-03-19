# Normalize data in a Seurat object

This function provides a unified normalization entry point for
Seurat-style log-normalization, scran normalization, and SCTransform.

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

  Optional cluster assignments for
  [`scran::quickCluster`](https://rdrr.io/pkg/scran/man/quickCluster.html).

- assay:

  Assay used for normalization. Defaults to `"RNA"`.

- layer:

  Layer used as the input count matrix. Defaults to `"counts"`.

- ...:

  Additional method-specific arguments passed to
  [`Seurat::NormalizeData()`](https://satijalab.org/seurat/reference/NormalizeData.html),
  [`scran::computeSumFactors()`](https://rdrr.io/pkg/scran/man/computeSumFactors.html),
  or
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
