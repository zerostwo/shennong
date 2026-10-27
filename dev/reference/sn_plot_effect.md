# Plot supplied effect estimates, confidence intervals or signed scores

Renders a result table without fitting models or computing confidence
intervals. Coefficients use a zero reference; ratios use one and a log
axis.

## Usage

``` r
sn_plot_effect(
  data,
  style = c("forest", "diverging_bar", "lollipop"),
  scale = c("identity", "ratio"),
  group_by = NULL,
  palette = NULL,
  value_label = "Effect",
  title = NULL
)
```

## Arguments

- data:

  Data frame containing `label`, `estimate` and, for `forest`, `lower`,
  `upper`. Intervals must bracket the estimates. Row order is retained.

- style:

  Structural style: `forest`, `diverging_bar`, or `lollipop`.

- scale:

  `identity` for coefficients/scores or `ratio` for positive HR/OR.

- group_by:

  Optional column for categorical color.

- palette:

  Named palette or explicit colors.

- value_label:

  Axis title specifying the statistic and units.

- title:

  Optional plot title.

## Value

A ggplot with the original effect table attached as source data.

## Examples

``` r
d <- data.frame(label = c("NK", "B"), estimate = c(-0.2, 0.4),
                lower = c(-0.3, 0.1), upper = c(-0.1, 0.7))
sn_plot_effect(d, value_label = "Linear coefficient (95% CI)")
```
