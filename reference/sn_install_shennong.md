# Install Shennong from CRAN or GitHub

This helper installs the stable CRAN release when available, or the
GitHub development version when requested. When `channel = "auto"`, it
prefers CRAN and falls back to GitHub if no CRAN release is available.

## Usage

``` r
sn_install_shennong(
  channel = c("auto", "cran", "github"),
  package = "Shennong",
  github_repo = "zerostwo/shennong",
  github_ref = "main",
  repos = getOption("repos"),
  ...
)
```

## Arguments

- channel:

  One of `"auto"`, `"cran"`, or `"github"`.

- package:

  Package name. Defaults to `"Shennong"`.

- github_repo:

  GitHub repository in `"owner/repo"` format.

- github_ref:

  GitHub ref to install from. Defaults to `"main"`.

- repos:

  CRAN-like repositories used by
  [`install.packages()`](https://rdrr.io/r/utils/install.packages.html).

- ...:

  Additional arguments passed to
  [`utils::install.packages()`](https://rdrr.io/r/utils/install.packages.html)
  or `remotes::install_github()`.

## Value

Invisibly returns the chosen installation channel.

## Examples

``` r
if (FALSE) { # \dontrun{
sn_install_shennong(channel = "github")
} # }
```
