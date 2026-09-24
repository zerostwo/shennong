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
  pixi_version = "0.69.0",
  pixi_download_url = NULL,
  pixi_sha256 = NULL,
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
  pixi_version = "0.69.0",
  pixi_download_url = NULL,
  pixi_sha256 = NULL,
  quiet = FALSE
)

sn_call_scvi(command, args = character(), ...)

sn_call_scanvi(command, args = character(), ...)

sn_call_mmochi(command, args = character(), ...)

sn_call_scarches(command, args = character(), ...)

sn_call_scpoli(command, args = character(), ...)

sn_call_infercnvpy(command, args = character(), ...)

sn_call_trajectory(command, args = character(), ...)

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

  Platforms that the bundled lock must cover. This does not narrow the
  platform inventory declared by the bundled manifest, because doing so
  would make its multi-platform lock stale. Defaults to the current
  platform.

- mirror:

  Mirror setting passed to
  [`sn_configure_pixi_mirror()`](https://zerostwo.github.io/shennong/dev/reference/sn_configure_pixi_mirror.md).

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

- pixi_sha256:

  Optional expected SHA-256 digest for `pixi_download_url`. Custom URLs
  require this value; official pinned release downloads obtain their
  checksum sidecar automatically.

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

`sn_prepare_pixi_environment()` returns a named list of paths (including
`manifest_path` and `lock_path`) and selected environment metadata.
`sn_call_pixi_environment()` invisibly returns command output.

## Details

The environment-specific aliases `sn_call_scvi()`, `sn_call_scanvi()`,
`sn_call_mmochi()`, `sn_call_scarches()`, `sn_call_scpoli()`,
`sn_call_infercnvpy()`, `sn_call_trajectory()`, `sn_call_cellphonedb()`,
`sn_call_cell2location()`, `sn_call_tangram()`, `sn_call_squidpy()`,
`sn_call_spatialdata()`, and `sn_call_stlearn()` are deprecated
compatibility wrappers that only forward to
`sn_call_pixi_environment()`. Call
`sn_call_pixi_environment("<environment>", command = ..., args = ...)`
directly; the aliases emit a deprecation warning and will be removed in
a future major release.

## Examples

``` r
sn_prepare_pixi_environment("scvi", runtime_dir = tempfile("shennong-home-"))
if (FALSE) { # \dontrun{
sn_call_pixi_environment("scvi", command = "python", args = "--version")
} # }
```
