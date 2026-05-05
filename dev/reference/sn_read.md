# Read tabular and bioinformatics file formats

A Shennong wrapper around
[`rio::import()`](http://gesistsa.github.io/rio/reference/import.md)
with support for common single-cell and omics formats such as 10x
directories, `.h5`, `.h5ad`, BPCells directories, and GMT files.

This adapter is exported so `rio` can dispatch to it.

This function imports Gene Matrix Transposed (GMT) files into a
two-column data.frame of term–gene pairs.

Adapted from
[`clusterProfiler::read.gmt()`](https://rdrr.io/pkg/gson/man/read-gmt.html)
(Y. Yu et al., 2012–2024), which is licensed under the Artistic License
2.0.

## Usage

``` r
sn_read(path, format, to = NULL, which, row_names = NULL, ...)

.import.rio_bpcells(file, ...)

.import.rio_10x(file, ...)

.import.rio_10x_spatial(file, ...)

.import.rio_starsolo(file, ...)

.import.rio_h5ad(file, ...)

.import.rio_h5(file, ...)

.import.rio_qs(file, ...)

.import.rio_qs2(file, ...)

.import.rio_gmt(file, ...)
```

## Arguments

- path:

  Path, URL, or supported directory to import.

- format:

  Optional explicit format override.

- to:

  Optional target class forwarded to `rio`.

- which:

  Optional archive member index for compressed archives.

- row_names:

  Optional column index or name to promote to row names when importing a
  data frame.

- ...:

  Additional arguments forwarded to the underlying importer.

- file:

  Input path used by the exported `rio` adapter methods.

## Value

The imported object. Tabular inputs are typically returned as data
frames, while supported matrix-like formats are returned in their native
classes.

A data.frame with two columns: `term` and `gene`.

## Details

The original implementation is available in the clusterProfiler package:
https://github.com/YuLab-SMU/clusterProfiler

## See also

[`sn_write()`](https://songqi.org/shennong/dev/reference/sn_write.md)

The original `read.gmt()` implementation in the `clusterProfiler`
package.

## Author

Yu Guangchuang (original implementation)

## Examples

``` r
tmp <- tempfile(fileext = ".csv")
write.csv(mtcars[1:3, 1:3], tmp, row.names = FALSE)
sn_read(tmp)
#>    mpg cyl disp
#> 1 21.0   6  160
#> 2 21.0   6  160
#> 3 22.8   4  108
```
