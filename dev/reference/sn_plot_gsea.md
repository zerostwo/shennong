# Plot a GSEA running-score curve or summary

Plot a GSEA running-score curve or summary

## Usage

``` r
sn_plot_gsea(result, pathway = NULL)
```

## Arguments

- result:

  A result containing `tables$running_score`, or a GSEA table.

- pathway:

  Optional pathway ID/term used to subset running-score data.

## Value

A `ggplot` object with a figure specification.
