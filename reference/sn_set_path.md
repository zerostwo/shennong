# Create a directory path if needed

This helper expands glue expressions in a path string, creates the
directory if it does not already exist, and returns the resulting path.

## Usage

``` r
sn_set_path(path)
```

## Arguments

- path:

  A character scalar giving the directory path to create.

## Value

The resulting path string.

## Examples

``` r
tmp_dir <- tempfile("shennong-")
sn_set_path(tmp_dir)
#> /tmp/RtmpD96fQl/shennong-1fbb6a10a81a
```
