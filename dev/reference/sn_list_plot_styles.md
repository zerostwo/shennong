# List publication-inspired themes and structural plotting styles

Themes control appearance across ggplot renderers. Structural styles
select different geometry in the listed function, with explicit input
requirements and paper-panel references. They do not fit models or infer
biological modules.

## Usage

``` r
sn_list_plot_styles(function_name = NULL)
```

## Arguments

- function_name:

  Optional function name, e.g. `"sn_plot_matrix"`.

## Value

A tibble retaining `style`, `description`, `reference`,
`changes_analysis`, plus `kind`, `function_name`, and
`data_requirements`. `changes_analysis = FALSE` means no biological
inference; documented summaries, density estimation, intersections and
display clustering remain.
