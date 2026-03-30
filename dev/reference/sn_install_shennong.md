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
  ref = NULL,
  github_repo = "zerostwo/shennong",
  github_ref = "main",
  local_path = NULL,
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

  Generic installation source. For `channel = "github"`, this should be
  an `"owner/repo"` string. For `channel = "local"`, this should be a
  local package directory or source tarball path. This is the preferred
  source argument for non-CRAN installs.

- ref:

  Generic ref argument. Used for `channel = "github"` and preferred over
  `github_ref`.

- github_repo:

  GitHub repository in `"owner/repo"` format. Deprecated in favor of
  `source`.

- github_ref:

  GitHub ref to install from. Defaults to `"main"`. Deprecated in favor
  of `ref`.

- local_path:

  Local package directory or source tarball used when
  `channel = "local"`. Deprecated in favor of `source`.

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
