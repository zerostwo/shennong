# Check the current acceleration status

Reports which ShennongOpt patches are registered, installable, and
active in the current session.

## Usage

``` r
sn_check_acceleration()
```

## Value

A named character vector mapping patch name to `"active"` or
`"inactive"`, or an empty vector when ShennongOpt is not installed.

## Examples

``` r
if (FALSE) sn_check_acceleration() # \dontrun{}
```
