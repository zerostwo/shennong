# Build consensus annotation labels from evidence

Build consensus annotation labels from evidence

## Usage

``` r
sn_annotation_consensus(
  evidence,
  weights = NULL,
  ontology = TRUE,
  ontology_mapping = NULL,
  low_confidence_threshold = 0.55,
  margin_threshold = 0.1
)
```

## Arguments

- evidence:

  Long evidence table accepted by
  [`sn_annotation_confidence()`](https://songqi.org/shennong/dev/reference/sn_annotation_confidence.md).
  Optional columns include `parent_label`, `supporting_markers`,
  `conflicting_markers`, and `reference_coverage`.

- weights:

  Optional named backend weights.

- ontology:

  Logical; map final labels to the bundled Cell Ontology snapshot.

- ontology_mapping:

  Optional custom ontology mapping.

- low_confidence_threshold, margin_threshold:

  Confidence thresholds.

## Value

A traceable annotation tibble with hierarchical labels, evidence,
confidence, and ontology identifiers.

## Examples

``` r
evidence <- data.frame(
  entity = rep(c("0", "1"), each = 2),
  label = rep(c("B cells", "T cells"), 2),
  score = c(3, 0.1, 0.2, 2.5), method = "markers",
  parent_label = rep(c("B cells", "T cells"), 2)
)
sn_annotation_consensus(evidence)
#> # A tibble: 2 × 16
#>   entity prediction second_best_label prediction_score margin
#>   <chr>  <chr>      <chr>                        <dbl>  <dbl>
#> 1 0      B cells    T cells                          1      1
#> 2 1      T cells    B cells                          1      1
#> # ℹ 11 more variables: method_count <int>, methods <chr>,
#> #   low_confidence <lgl>, level_1 <chr>, level_2 <chr>, level_3 <chr>,
#> #   supporting_markers <chr>, conflicting_markers <chr>,
#> #   reference_coverage <dbl>, ontology_id <chr>, ontology_label <chr>
```
