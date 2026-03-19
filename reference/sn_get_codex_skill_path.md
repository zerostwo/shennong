# Return the installed Shennong Codex skill path

This helper returns the skill directory bundled with the installed
Shennong package.

## Usage

``` r
sn_get_codex_skill_path()
```

## Value

A character scalar path to the bundled skill directory.

## Examples

``` r
if (requireNamespace("Shennong", quietly = TRUE)) {
  path <- sn_get_codex_skill_path()
  dir.exists(path)
}
#> [1] TRUE
```
