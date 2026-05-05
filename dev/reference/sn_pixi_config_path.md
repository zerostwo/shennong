# Locate a bundled pixi config

Locate a bundled pixi config

## Usage

``` r
sn_pixi_config_path(environment = NULL)
```

## Arguments

- environment:

  Python environment name.

## Value

A path to the package-bundled `pixi.toml` template.

## Examples

``` r
sn_pixi_config_path("scvi")
#> [1] "/home/runner/work/_temp/Library/Shennong/pixi/scvi/pixi.toml"
```
