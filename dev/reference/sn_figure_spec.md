# Inspect or calculate a Shennong figure specification

Inspect or calculate a Shennong figure specification

## Usage

``` r
sn_figure_spec(
  plot = NULL,
  plot_type = NULL,
  data_summary = NULL,
  profile = NULL,
  ...
)
```

## Arguments

- plot:

  A ggplot, patchwork, ComplexHeatmap, or `NULL` when calculating from
  synthetic metadata only.

- plot_type:

  Optional plot type used by the sizing rules.

- data_summary:

  Optional named list describing points, groups, panels, features,
  labels, network size, or spatial aspect ratio.

- profile:

  Generic output profile.

- ...:

  Explicit recommendation overrides such as `width_mm` or `rasterize`.

## Value

A figure specification list.
