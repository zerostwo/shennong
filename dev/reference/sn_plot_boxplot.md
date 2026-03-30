# Create a boxplot from a data frame

Create a boxplot from a data frame

## Usage

``` r
sn_plot_boxplot(
  data,
  x,
  y,
  sort = FALSE,
  aspect_ratio = NULL,
  panel_widths = NULL,
  panel_heights = NULL,
  x_label = NULL,
  y_label = NULL,
  title = NULL,
  angle_x = 0
)
```

## Arguments

- data:

  A data frame.

- x, y:

  Columns mapped to the x- and y-axes.

- sort:

  Currently reserved for future sorting support.

- aspect_ratio:

  Optional panel aspect ratio. When used together with `panel_widths` or
  `panel_heights`, Shennong derives the missing panel dimension
  automatically.

- panel_widths, panel_heights:

  Optional panel size arguments forwarded to
  [`catplot::theme_cat()`](https://rdrr.io/pkg/catplot/man/theme_cat.html)
  when available.

- x_label, y_label, title:

  Optional plot labels.

- angle_x:

  Rotation angle for x-axis text.

## Value

A ggplot object.

## Examples

``` r
sn_plot_boxplot(mtcars, x = cyl, y = mpg)
#> Warning: Orientation is not uniquely specified when both the x and y aesthetics
#> are continuous. Picking default orientation 'x'.
#> Warning: Continuous x aesthetic
#> ℹ did you forget `aes(group = ...)`?

```
