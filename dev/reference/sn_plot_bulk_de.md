# Plot bulk differential expression

Plot bulk differential expression

## Usage

``` r
sn_plot_bulk_de(x, adjusted_p_value = 0.05, log2_fold_change = 1, labels = 10L)
```

## Arguments

- x:

  A bulk-DE result.

- adjusted_p_value:

  Adjusted p-value threshold.

- log2_fold_change:

  Absolute fold-change threshold.

- labels:

  Number of top genes to label when ggrepel is installed.

## Value

A volcano `ggplot` object.
