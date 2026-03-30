# Install missing Shennong dependencies in one step

This helper installs missing required and/or recommended packages using
the declared source for each dependency. CRAN packages are installed
with
[`install.packages()`](https://rdrr.io/r/utils/install.packages.html),
Bioconductor packages with `BiocManager::install()`, and GitHub packages
with `remotes::install_github()`.

## Usage

``` r
sn_install_dependencies(
  scope = c("all", "required", "recommended"),
  packages = NULL,
  missing_only = TRUE,
  repos = getOption("repos"),
  ask = interactive(),
  upgrade = FALSE,
  ...
)
```

## Arguments

- scope:

  One of `"all"`, `"required"`, or `"recommended"`.

- packages:

  Optional character vector restricting installation to a subset of
  packages returned by
  [`sn_list_dependencies()`](https://songqi.org/shennong/dev/reference/sn_list_dependencies.md).

- missing_only:

  Logical; when `TRUE` (default), install only missing packages.

- repos:

  CRAN-like repositories used for CRAN installs and bootstrap
  installation of helper installers such as BiocManager and remotes.

- ask:

  Passed to `BiocManager::install()`. Defaults to interactive behavior.

- upgrade:

  Logical; when `TRUE`, allow updating already installed GitHub and
  Bioconductor packages during installation.

- ...:

  Additional arguments forwarded to the underlying installer calls.

## Value

Invisibly returns the refreshed dependency table from
[`sn_list_dependencies()`](https://songqi.org/shennong/dev/reference/sn_list_dependencies.md)
after installation.

## Examples

``` r
if (FALSE) { # \dontrun{
sn_install_dependencies(scope = "required")
sn_install_dependencies(packages = c("Seurat", "clusterProfiler"))
} # }
```
