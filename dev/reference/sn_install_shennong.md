# Install Shennong from CRAN, GitHub, or a local source

This helper installs the stable CRAN release when available, or the
GitHub development version when requested, or installs from a local
source tree or tarball. When `channel = "auto"`, it prefers CRAN and
falls back to GitHub if no CRAN release is available.

## Usage

``` r
sn_install_shennong(
  channel = c("auto", "cran", "github", "local"),
  package = "Shennong",
  source = NULL,
  ref = "main",
  repos = getOption("repos"),
  ...
)
```

## Arguments

- channel:

  One of `"auto"`, `"cran"`, `"github"`, or `"local"`.

- package:

  Package name. Defaults to `"Shennong"`.

- source:

  Installation source. For `channel = "github"`, supply an
  `"owner/repo"` string; when omitted, Shennong uses
  `"zerostwo/shennong"`. For `channel = "local"`, supply a local package
  directory or source tarball path.

- ref:

  GitHub ref used for `channel = "github"`. Defaults to `"main"`.

- repos:

  CRAN-like repositories used by
  [`install.packages()`](https://rdrr.io/r/utils/install.packages.html).

- ...:

  Additional arguments passed to
  [`utils::install.packages()`](https://rdrr.io/r/utils/install.packages.html)
  or `remotes::install_github()` / `remotes::install_local()`. For
  GitHub installs, Shennong defaults to `dependencies = FALSE` and
  `upgrade = "never"` unless you override them explicitly.

## Value

Invisibly returns the chosen installation channel.

## Examples

``` r
if (FALSE) { # \dontrun{
sn_install_shennong(channel = "github")
sn_install_shennong(channel = "github", source = "zerostwo/shennong", ref = "main")
sn_install_shennong(channel = "local", source = "~/personal/packages/shennong")
} # }
```
