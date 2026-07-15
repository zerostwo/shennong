# Recommend output dimensions for a figure

Recommend output dimensions for a figure

## Usage

``` r
sn_recommend_figure_size(
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

The `recommended` component of
[`sn_figure_spec()`](https://songqi.org/shennong/dev/reference/sn_figure_spec.md).
