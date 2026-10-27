# Inspect Shennong pixi runtime paths

Returns the user-level paths used for optional pixi-managed Python
environments. Shennong follows the same convention as downloaded example
data: runtime files are generated under `~/.shennong/` by default, not
under the current analysis project and not inside the installed R
package.

## Usage

``` r
sn_get_pixi_paths(environment = NULL, runtime_dir = NULL)
```

## Arguments

- environment:

  Python environment name. Use
  [`sn_list_pixi_environments()`](https://zerostwo.github.io/shennong/dev/reference/sn_list_pixi_environments.md)
  to see bundled configs.

- runtime_dir:

  Optional explicit Shennong runtime directory. Defaults to
  `getOption("shennong.runtime_dir")`, `SHENNONG_RUNTIME_DIR`,
  `SHENNONG_HOME`, then `"~/.shennong"`.

## Value

A named list of runtime paths.

## Examples

``` r
sn_get_pixi_paths("scvi", runtime_dir = tempfile("shennong-home-"))
```
