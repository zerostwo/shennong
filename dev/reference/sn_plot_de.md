# Plot a differential-expression result

Plot a differential-expression result

## Usage

``` r
sn_plot_de(
  result,
  type = c("volcano", "ma", "effect", "heatmap"),
  adjusted_p_value = 0.05,
  log2_fold_change = 1,
  n = 40L
)
```

## Arguments

- result:

  A Shennong result or differential-expression table.

- type:

  Volcano, MA, ranked effect, or comparison heatmap.

- adjusted_p_value, log2_fold_change:

  Significance thresholds.

- n:

  Maximum effects displayed for effect/heatmap views.

## Value

A result-aware `ggplot` with an attached figure specification.
