# Remove ambient RNA contamination from counts

This function provides a unified interface for ambient RNA correction
using either `SoupX` or `decontX`. The input can be a Seurat object, a
count matrix-like object, or a path that
[`sn_read()`](https://songqi.org/shennong/reference/sn_read.md) can
import.

## Usage

``` r
sn_remove_ambient_contamination(
  x,
  raw = NULL,
  method = c("decontx", "soupx"),
  cluster = NULL,
  remove_zero_count_cells = FALSE,
  layer = "decontaminated_counts",
  return_object = TRUE,
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  A Seurat object, count matrix-like object, or path to filtered data.

- raw:

  Optional raw/background counts. Required for `method = "soupx"`. If
  supplied for `decontx`, it is used as the background matrix.

- method:

  One of `"decontx"` or `"soupx"`.

- cluster:

  Optional cluster labels. This can be a vector with one value per cell,
  or a metadata column name when `x` is a Seurat object. If `NULL`,
  clusters are inferred with
  [`sn_run_cluster()`](https://songqi.org/shennong/reference/sn_run_cluster.md).

- remove_zero_count_cells:

  Logical; if `TRUE`, remove cells whose decontX-corrected counts sum to
  zero. If `FALSE`, keep those cells by restoring their original counts
  and emit a warning.

- layer:

  Layer name used when writing corrected counts back to a Seurat object.
  Defaults to `"decontaminated_counts"`. Use `"counts"` to overwrite the
  original counts layer explicitly.

- return_object:

  If `TRUE` and `x` is a Seurat object, return the updated Seurat
  object. Otherwise return the corrected counts matrix.

- verbose:

  Logical; whether to print progress from helper clustering.

- ...:

  Additional method-specific arguments passed to
  [`celda::decontX()`](https://rdrr.io/pkg/celda/man/decontX.html) when
  `method = "decontx"`, or to `SoupX::autoEstCont()` when
  `method = "soupx"`.

## Value

A corrected counts matrix, or an updated Seurat object when
`return_object = TRUE` and `x` is a Seurat object. For decontX-based
Seurat returns, the `decontX_contamination` and `decontX_clusters`
columns are added to `meta.data`.

## Examples

``` r
if (FALSE) { # \dontrun{
data("pbmc_small", package = "Shennong")
data("pbmc_small_raw", package = "Shennong")
pbmc <- sn_remove_ambient_contamination(pbmc, method = "decontx")

corrected <- sn_remove_ambient_contamination(
  x = SeuratObject::LayerData(pbmc, layer = "counts"),
  raw = pbmc_small_raw,
  method = "soupx",
  return_object = FALSE
)
} # }
```
