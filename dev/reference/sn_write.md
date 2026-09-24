# Write tabular and bioinformatics file formats

A Shennong wrapper around
[`rio::export()`](http://gesistsa.github.io/rio/reference/export.md)
with support for common single-cell formats such as BPCells, `.h5ad`,
and 10x `.h5`. Native RDS, RData, and QS2 serialization preserves matrix
classes and dimnames. Tabular formats convert base matrices to data
frames. Dense and sparse matrices are converted to BPCells iterable
matrices for matrix writers.

This adapter is exported so `rio` can dispatch to it.

## Usage

``` r
sn_write(
  x,
  path = NULL,
  to = NULL,
  auto_install = TRUE,
  install_repos = getOption("repos"),
  install_ask = FALSE,
  ...
)

.export.rio_bpcells(file, x, overwrite = FALSE, ...)

.export.rio_h5ad(file, x, mode = "w", ...)

.export.rio_h5(file, x, ...)

.export.rio_qs2(file, x, ...)
```

## Arguments

- x:

  Object to write.

- path:

  Output path. Missing parent directories are created automatically
  before dispatching the selected writer.

- to:

  Optional format override when it cannot be inferred from `path`.

- auto_install:

  Logical; when `TRUE`, install missing writer dependencies before
  writing. This includes rio plus optional Shennong custom writer
  dependencies such as `.qs2`, `.h5ad`, `.h5`, or BPCells.

- install_repos:

  CRAN-like repositories used when `auto_install` needs to install CRAN
  packages.

- install_ask:

  Passed to `BiocManager::install()` when `auto_install` installs
  Bioconductor packages.

- ...:

  Additional arguments forwarded to the underlying exporter.

- file:

  Output path used by the exported `rio` adapter methods.

- overwrite:

  Logical; overwrite an existing BPCells directory.

- mode:

  File mode passed through to `anndataR::write_h5ad()`.

## Value

Invisibly returns the output path.

## See also

[`sn_read()`](https://zerostwo.github.io/shennong/dev/reference/sn_read.md)

## Examples

``` r
tmp <- tempfile(fileext = ".csv")
sn_write(mtcars[1:3, 1:3], tmp)
```
