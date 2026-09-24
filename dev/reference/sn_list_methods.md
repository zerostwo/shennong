# List registered Shennong analysis methods

Reads the shipped method registry and reports whether each backend is
implemented and currently available. Registered roadmap methods remain
discoverable even when their optional dependency or pixi runtime has not
been installed.

## Usage

``` r
sn_list_methods(task = NULL, available = NULL)
```

## Arguments

- task:

  Optional task name such as `"trajectory"`, `"annotation"`, or
  `"bulk"`.

- available:

  Optional logical filter. Use `TRUE` for methods that can run in the
  current session and `FALSE` for unavailable methods.

## Value

A tibble with one row per task/method pair.

## Examples

``` r
sn_list_methods("trajectory")
sn_list_methods(available = TRUE)
```
