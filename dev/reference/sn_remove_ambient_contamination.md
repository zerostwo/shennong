# Remove ambient RNA contamination from counts

This function provides a unified interface for ambient RNA or CITE-seq
contamination correction using `SoupX`, `decontX::decontX()`, or
`decontX::decontPro()`. The input can be a Seurat object, a count
matrix-like object, or a path that
[`sn_read()`](https://zerostwo.github.io/shennong/dev/reference/sn_read.md)
can import. When a Seurat object was initialized from a detected 10x
`outs` directory, stored raw matrix metadata is reused automatically if
`raw = NULL`. The `auto` method routes a Seurat object with one
CITE-seq-like ADT, protein, or CITE assay to `decontX::decontPro()`;
otherwise it uses `decontX::decontX()`. Set
`options(shennong.acceleration = FALSE)` or
`SHENNONG_ACCELERATION_DISABLED=true` to prevent automatic acceleration;
a manually active patch remains active until
[`sn_disable_acceleration()`](https://zerostwo.github.io/shennong/dev/reference/sn_disable_acceleration.md)
is called.

## Usage

``` r
sn_remove_ambient_contamination(
  x,
  raw = NULL,
  method = c("auto", "decontx", "decontpro", "soupx"),
  cluster = NULL,
  cluster_backend = NULL,
  remove_zero_count_cells = FALSE,
  layer = "decontaminated_counts",
  return_object = TRUE,
  verbose = FALSE,
  seed = 717L,
  ...,
  assay = NULL
)
```

## Arguments

- x:

  A Seurat object, count matrix-like object, or path to filtered data.

- raw:

  Optional raw/background counts. Required for `method = "soupx"` unless
  recoverable from stored initialization metadata. If supplied for
  `decontx`, it is used as the background matrix.

- method:

  One of `"auto"`, `"decontx"`, `"decontpro"`, or `"soupx"`. The default
  `"auto"` selects `"decontpro"` for a detected CITE-seq-like assay and
  `"decontx"` otherwise.

- cluster:

  Optional cluster labels. This can be a vector with one value per cell,
  or a metadata column name when `x` is a Seurat object. Explicit labels
  take precedence over `cluster_backend`.

- cluster_backend:

  Clustering implementation used when `cluster` is `NULL`. The
  method-specific default is `"native"` for decontX, which lets decontX
  perform its own clustering, and `"shennong"` for SoupX and decontPro,
  which require assignments and obtain them from
  [`sn_run_cluster()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_cluster.md).
  SoupX and decontPro do not support `"native"`.

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

- seed:

  Random seed used only for reproducible stochastic rounding of
  SoupX-adjusted fractional counts. The caller's random-number-generator
  state is restored after rounding.

- ...:

  Additional method-specific arguments passed to `decontX::decontX()`,
  `decontX::decontPro()`, or `SoupX::autoEstCont()` and
  `SoupX::adjustCounts()`.

- assay:

  Optional Seurat assay containing the counts to correct. For
  `method = "auto"`, an assay whose name contains `ADT`, `protein`, or
  `CITE` selects `"decontpro"`; supply it explicitly when more than one
  such assay is present. For matrices and paths, `assay` is not
  applicable.

## Value

A corrected counts matrix, or an updated Seurat object when
`return_object = TRUE` and `x` is a Seurat object. For Seurat returns,
`nCount_<assay>_corrected` and `nFeature_<assay>_corrected` are added to
`meta.data`. For decontX-based Seurat returns, the
`decontX_contamination` and `decontX_clusters` columns are added.
decontPro-based returns also include the inferred cell type and
ambient/background fraction columns. The resolved method, assay,
clustering backend, input-source summary, method-specific `...`
arguments, and backend version are recorded in
`object@commands$sn_remove_ambient_contamination@params` so automatic
and wrapper-controlled choices remain inspectable.

## Examples

``` r
if (FALSE) { # \dontrun{
pbmc <- qs2::qs_read(file.path(
  Sys.getenv("SHENNONG_REAL_DATA_DIR"), "single-cell", "kotliarov_pbmc.qs2"
))
pbmc <- sn_remove_ambient_contamination(pbmc, method = "decontx")
} # }
```
