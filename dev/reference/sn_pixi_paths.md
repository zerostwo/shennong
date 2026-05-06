# Inspect Shennong pixi runtime paths

Returns the user-level paths used for optional pixi-managed Python
environments. Shennong follows the same convention as downloaded example
data: runtime files are generated under `~/.shennong/` by default, not
under the current analysis project and not inside the installed R
package.

## Usage

``` r
sn_pixi_paths(environment = NULL, path = ".", runtime_dir = NULL)
```

## Arguments

- environment:

  Python environment name. Use
  [`sn_list_pixi_environments()`](https://songqi.org/shennong/dev/reference/sn_list_pixi_environments.md)
  to see bundled configs.

- path:

  Deprecated; retained for compatibility and ignored unless
  `runtime_dir` is supplied explicitly.

- runtime_dir:

  Optional explicit Shennong runtime directory. Defaults to
  `getOption("shennong.runtime_dir")`, `SHENNONG_RUNTIME_DIR`,
  `SHENNONG_HOME`, then `"~/.shennong"`.

## Value

A named list of runtime paths.

## Examples

``` r
sn_pixi_paths("scvi", runtime_dir = tempfile("shennong-home-"))
#> $environment
#> [1] "scvi"
#> 
#> $family
#> [1] "scvi"
#> 
#> $runtime_dir
#> [1] "/tmp/RtmpvOjCeK/shennong-home-1c69285f08a"
#> 
#> $pixi_root
#> [1] "/tmp/RtmpvOjCeK/shennong-home-1c69285f08a/pixi"
#> 
#> $pixi_home
#> [1] "/tmp/RtmpvOjCeK/shennong-home-1c69285f08a/pixi/home"
#> 
#> $project_dir
#> [1] "/tmp/RtmpvOjCeK/shennong-home-1c69285f08a/pixi/scvi"
#> 
#> $source_config_path
#> [1] "/home/runner/work/_temp/Library/Shennong/pixi/scvi/pixi.toml"
#> 
#> $manifest_path
#> [1] "/tmp/RtmpvOjCeK/shennong-home-1c69285f08a/pixi/scvi/pixi.toml"
#> 
#> $workspace_env_dir
#> [1] "/tmp/RtmpvOjCeK/shennong-home-1c69285f08a/pixi/scvi/.pixi/envs"
#> 
#> $runs_dir
#> [1] "/tmp/RtmpvOjCeK/shennong-home-1c69285f08a/runs"
#> 
```
