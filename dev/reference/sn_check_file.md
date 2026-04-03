# Check if files exist

This function takes a vector of file paths as input and checks if each
file exists. If a file does not exist, the function returns a message
indicating which file(s) do(es) not exist.

## Usage

``` r
sn_check_file(x, stop = TRUE)
```

## Arguments

- x:

  A vector of file paths.

- stop:

  A logical value indicating whether to stop the function if a file does
  not exist. Default is TRUE.

## Value

If stop is TRUE, the function stops and returns a message indicating
which file(s) do(es) not exist. If stop is FALSE, the function returns a
vector of file paths that do not exist.

## See also

[`file.exists`](https://rdrr.io/r/base/files.html),
[`stop`](https://rdrr.io/r/base/stop.html)

## Examples

``` r
existing <- tempfile("existing-")
file.create(existing)
#> [1] TRUE
missing <- tempfile("missing-")
sn_check_file(c(existing, missing), stop = FALSE)
#> [1] "/tmp/RtmpQTl3iH/missing-1f8b728607c5"
```
