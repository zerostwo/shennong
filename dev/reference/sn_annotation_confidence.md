# Calibrate confidence from annotation evidence

Scores are normalized within each entity/backend before weighted
consensus. The output retains the best and second-best labels, their
margin, method support, and a conservative low-confidence flag.

## Usage

``` r
sn_annotation_confidence(
  evidence,
  weights = NULL,
  low_confidence_threshold = 0.55,
  margin_threshold = 0.1
)
```

## Arguments

- evidence:

  Long data frame with `entity`, `label`, `score`, and `method` columns.

- weights:

  Optional named numeric vector of backend weights.

- low_confidence_threshold:

  Minimum calibrated prediction score.

- margin_threshold:

  Minimum best-versus-second-best margin.

## Value

A tibble with one calibrated prediction per entity.

## Examples

``` r
evidence <- data.frame(
  entity = rep(c("0", "1"), each = 2),
  label = rep(c("B cells", "T cells"), 2),
  score = c(2, 0.2, 0.1, 3), method = "markers"
)
sn_annotation_confidence(evidence)
#> # A tibble: 2 × 8
#>   entity prediction second_best_label prediction_score margin
#>   <chr>  <chr>      <chr>                        <dbl>  <dbl>
#> 1 0      B cells    T cells                          1      1
#> 2 1      T cells    B cells                          1      1
#> # ℹ 3 more variables: method_count <int>, methods <chr>,
#> #   low_confidence <lgl>
```
