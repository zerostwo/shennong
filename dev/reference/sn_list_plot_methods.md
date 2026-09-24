# List result-aware plot methods

Returns the canonical analysis-type and view vocabulary accepted by
[`sn_plot_result()`](https://zerostwo.github.io/shennong/dev/reference/sn_plot_result.md).
The table is intended for discovery and programmatic UI generation; one
row represents one valid analysis/view pair.

## Usage

``` r
sn_list_plot_methods(analysis_type = NULL)
```

## Arguments

- analysis_type:

  Optional analysis type or documented alias used to filter the returned
  table.

## Value

A tibble with `analysis_type`, `view`, `default`, `required_parameters`,
and `accepted_input` columns.

## Examples

``` r
sn_list_plot_methods("de")
```
