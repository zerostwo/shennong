# Check whether Shennong is up to date

This helper compares the installed package version against the latest
available CRAN or GitHub development version. When `channel = "auto"`,
it prefers CRAN if a release exists there and otherwise falls back to
GitHub.

## Usage

``` r
sn_check_version(
  channel = c("auto", "cran", "github"),
  package = "Shennong",
  github_repo = "zerostwo/shennong",
  github_ref = "main",
  repos = getOption("repos"),
  quiet = FALSE
)
```

## Arguments

- channel:

  One of `"auto"`, `"cran"`, or `"github"`.

- package:

  Package name to check. Defaults to `"Shennong"`.

- github_repo:

  GitHub repository in `"owner/repo"` format.

- github_ref:

  GitHub ref to inspect. Defaults to `"main"`.

- repos:

  CRAN-like repositories used for version lookup.

- quiet:

  Logical; if `TRUE`, suppress the summary message.

## Value

A named list with the installed version, remote version, selected
channel, whether the package is up to date, and the recommended install
command.

## Examples

``` r
if (FALSE) { # \dontrun{
sn_check_version()
sn_check_version(channel = "github")
} # }
```
