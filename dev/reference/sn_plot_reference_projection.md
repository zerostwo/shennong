# Plot reference annotation projection

Plot reference annotation projection

## Usage

``` r
sn_plot_reference_projection(
  object,
  result_id = "annotation",
  reduction = NULL,
  color_by = c("prediction", "prediction_score")
)
```

## Arguments

- object:

  A Seurat object with a stored annotation result.

- result_id:

  Annotation result name.

- reduction:

  Reduction used for coordinates.

- color_by:

  Prediction or confidence.

## Value

A reference-projection embedding.
