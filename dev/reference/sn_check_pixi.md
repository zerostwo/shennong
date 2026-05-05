# Check, install, and configure pixi for Shennong Python backends

Shennong uses pixi for optional Python backends such as scVI/scANVI. The
R package ships runner scripts, but the Python environments themselves
are created under the user-level Shennong runtime directory
`~/.shennong/`, not inside the current analysis project and not inside
the installed R package.

## Usage

``` r
sn_check_pixi(pixi = NULL, quiet = FALSE)

sn_install_pixi(
  version = "latest",
  pixi_home = "~/.pixi",
  bin_dir = NULL,
  no_path_update = TRUE,
  download_url = NULL,
  force = FALSE,
  quiet = FALSE
)

sn_ensure_pixi(
  pixi = NULL,
  install = TRUE,
  version = "latest",
  pixi_home = "~/.pixi",
  bin_dir = NULL,
  no_path_update = TRUE,
  download_url = NULL,
  quiet = FALSE
)
```

## Arguments

- pixi:

  Optional pixi executable path. When `NULL`, Shennong checks
  `options("shennong.pixi")`, `PATH`, and `~/.pixi/bin/pixi`.

- quiet:

  Logical; suppress status messages where possible.

- version:

  Pixi version passed to the official installer script. Defaults to
  `"latest"`.

- pixi_home:

  Pixi home directory. The installer defaults to `"~/.pixi"`; Shennong
  runtime workflows set `PIXI_HOME` to a user-level path such as
  `"~/.shennong/pixi/home"`.

- bin_dir:

  Optional directory where the standalone pixi binary should be
  installed.

- no_path_update:

  Logical; when `TRUE`, ask the installer not to edit shell startup
  files.

- download_url:

  Optional custom pixi binary download URL. This is useful for
  institutional mirrors.

- force:

  Logical; reinstall even if pixi is already available.

- install:

  Logical used by `sn_ensure_pixi()`; install pixi when it is not
  already available.

## Value

`sn_check_pixi()` and `sn_ensure_pixi()` return a named list with
install status, executable path, and version. `sn_install_pixi()`
invisibly returns the refreshed check result.

## Examples

``` r
info <- sn_check_pixi(quiet = TRUE)
info$installed
#> [1] FALSE
if (FALSE) { # \dontrun{
sn_ensure_pixi()
sn_install_pixi(pixi_home = "~/.pixi", no_path_update = TRUE)
} # }
```
