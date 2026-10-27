# Plot bulk survival associations

Plot bulk survival associations

## Usage

``` r
sn_plot_survival(
  x,
  adjusted_p_value = NULL,
  view = c("forest", "km", "risk_table", "ph", "ph_test", "cumulative_hazard"),
  feature = NULL,
  object = NULL,
  style = "classic"
)
```

## Arguments

- x:

  A bulk-survival result.

- adjusted_p_value:

  Optional adjusted p-value cutoff.

- view:

  Survival view: hazard-ratio forest, Kaplan-Meier curve, risk table,
  scaled Schoenfeld residual diagnostics, proportional-hazards test
  p-values, or cumulative hazard.

- feature:

  Optional feature subset for the selected view.

- object:

  Alias for `x`; supply only one of `x` and `object`.

- style:

  Presentation preset; see
  [`sn_list_plot_styles()`](https://zerostwo.github.io/shennong/dev/reference/sn_list_plot_styles.md).

## Value

A survival `ggplot` object.
