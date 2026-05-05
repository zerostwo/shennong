# Configure pixi mirrors for Shennong runtime environments

Writes a pixi `config.toml` under the selected `PIXI_HOME`. By default
Shennong uses a user-level pixi home such as `~/.shennong/pixi/home`, so
this does not mutate pixi's global `~/.pixi/config.toml` unless
`pixi_home` points there explicitly.

## Usage

``` r
sn_configure_pixi_mirror(
  mirror = c("default", "auto", "china", "tuna", "ustc", "bfsu"),
  pixi_home = NULL,
  runtime_dir = NULL,
  config_path = NULL,
  append_original = TRUE
)
```

## Arguments

- mirror:

  One of `"default"`, `"auto"`, `"china"`, `"tuna"`, `"ustc"`, or
  `"bfsu"`.

- pixi_home:

  Pixi home directory. When `NULL`, it is derived from `runtime_dir`.

- runtime_dir:

  Optional Shennong runtime directory used to derive `pixi_home`.

- config_path:

  Optional explicit config path.

- append_original:

  Whether to keep the original conda-forge URL as a fallback after
  mirror URLs.

## Value

Invisibly returns the config path, or `NA_character_` when no mirror
config was written.

## Examples

``` r
tmp <- tempfile("pixi-home-")
sn_configure_pixi_mirror("tuna", pixi_home = tmp)
```
