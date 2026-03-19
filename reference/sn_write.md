# Write tabular and bioinformatics file formats

A Shennong wrapper around
[`rio::export()`](http://gesistsa.github.io/rio/reference/export.md)
with support for common single-cell formats such as BPCells, `.h5ad`,
and 10x `.h5`.

This adapter is exported so `rio` can dispatch to it.

## Usage

``` r
sn_write(x, path = NULL, to = NULL, ...)

.export.rio_bpcells(file, x, overwrite = FALSE, ...)

.export.rio_h5ad(file, x, mode = "w", ...)

.export.rio_h5(file, x, ...)
```

## Arguments

- x:

  Object to write.

- path:

  Output path.

- to:

  Optional format override when it cannot be inferred from `path`.

- ...:

  Additional arguments forwarded to the underlying exporter.

- file:

  Output path used by the exported `rio` adapter methods.

- overwrite:

  Logical; overwrite an existing BPCells directory.

- mode:

  File mode passed through to
  [`anndataR::write_h5ad()`](https://anndataR.scverse.org/reference/write_h5ad.html).

## Value

Invisibly returns the output path.

## See also

[`sn_read()`](https://songqi.org/shennong/reference/sn_read.md)

## Examples

``` r
tmp <- tempfile(fileext = ".csv")
sn_write(mtcars[1:3, 1:3], tmp)
```
