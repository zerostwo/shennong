# Prepare or call a Shennong pixi environment

`sn_prepare_pixi_environment()` materializes a package-bundled
`inst/pixi/<family>/pixi.toml` template into the user-level
`~/.shennong/pixi/<family>/` workspace. `sn_call_pixi_environment()`
runs a command inside one of these environments.

## Usage

``` r
sn_prepare_pixi_environment(
  environment = NULL,
  pixi_environment = c("auto", "default", "cpu", "gpu"),
  runtime_dir = NULL,
  project_dir = NULL,
  manifest_path = NULL,
  overwrite = FALSE,
  cuda_version = NULL,
  platforms = NULL,
  mirror = c("default", "auto", "china", "tuna", "ustc", "bfsu"),
  install_pixi = FALSE,
  install_environment = FALSE,
  pixi = NULL,
  pixi_version = "latest",
  pixi_download_url = NULL,
  quiet = FALSE
)

sn_call_pixi_environment(
  environment = NULL,
  command,
  args = character(),
  pixi_environment = c("auto", "default", "cpu", "gpu"),
  runtime_dir = NULL,
  project_dir = NULL,
  manifest_path = NULL,
  overwrite = FALSE,
  cuda_version = NULL,
  platforms = NULL,
  mirror = c("default", "auto", "china", "tuna", "ustc", "bfsu"),
  install_pixi = TRUE,
  pixi = NULL,
  pixi_version = "latest",
  pixi_download_url = NULL,
  quiet = FALSE
)

sn_call_scvi(command, args = character(), ...)

sn_call_scanvi(command, args = character(), ...)

sn_call_scarches(command, args = character(), ...)

sn_call_scpoli(command, args = character(), ...)

sn_call_infercnvpy(command, args = character(), ...)

sn_call_cellphonedb(command, args = character(), ...)

sn_call_cell2location(command, args = character(), ...)

sn_call_tangram(command, args = character(), ...)

sn_call_squidpy(command, args = character(), ...)

sn_call_spatialdata(command, args = character(), ...)

sn_call_stlearn(command, args = character(), ...)
```

## Arguments

- environment:

  Python environment name.

- pixi_environment:

  Pixi environment inside the manifest, for example `"cpu"`, `"gpu"`, or
  `"default"`. `"auto"` uses CUDA when available for GPU-aware configs
  and otherwise CPU/default.

- runtime_dir:

  Optional Shennong runtime directory.

- project_dir:

  Optional explicit pixi workspace directory.

- manifest_path:

  Optional explicit materialized manifest path.

- overwrite:

  Whether to overwrite an existing materialized manifest.

- cuda_version:

  CUDA runtime version used when rendering templates.

- platforms:

  Pixi platform vector. Defaults to the current platform.

- mirror:

  Mirror setting passed to
  [`sn_configure_pixi_mirror()`](https://songqi.org/shennong/dev/reference/sn_configure_pixi_mirror.md).

- install_pixi:

  Ensure the standalone pixi binary is available.

- install_environment:

  Run `pixi install` for the selected environment after materializing
  the manifest.

- pixi:

  Optional pixi executable path.

- pixi_version:

  Pixi version used if installation is needed.

- pixi_download_url:

  Optional custom pixi binary download URL.

- quiet:

  Logical; suppress status messages where possible.

- command:

  Command to run inside the pixi environment.

- args:

  Character vector of command arguments.

- ...:

  Additional arguments passed from environment-specific helpers to
  `sn_call_pixi_environment()`.

## Value

`sn_prepare_pixi_environment()` returns a named list of paths and
selected environment metadata. `sn_call_pixi_environment()` invisibly
returns command output.

## Examples

``` r
sn_prepare_pixi_environment("scvi", runtime_dir = tempfile("shennong-home-"))
if (FALSE) { # \dontrun{
sn_call_pixi_environment("scvi", command = "python", args = "--version")
} # }
```
