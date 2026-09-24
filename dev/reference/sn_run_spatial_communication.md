# Add spatial distance evidence to a communication result

Add spatial distance evidence to a communication result

## Usage

``` r
sn_run_spatial_communication(
  object,
  source_result_id = "communication",
  communication = NULL,
  group_by,
  spatial_cols = NULL,
  max_distance = NULL,
  result_id = "spatial_communication",
  return_object = TRUE,
  sample_by = NULL
)
```

## Arguments

- object:

  A Seurat object.

- source_result_id:

  Stored communication result name.

- communication:

  Optional communication result supplied directly.

- group_by:

  Metadata column matching communication source/target labels.

- spatial_cols:

  Coordinate metadata columns.

- max_distance:

  Optional finite non-negative maximum mean nearest-group distance.

- result_id:

  Stored result name.

- return_object:

  Return the modified object or result.

- sample_by:

  Optional metadata column defining independent samples or tissue
  sections. Distances are calculated within sections and then
  aggregated, never between sections. When communication rows contain a
  non-missing `sample` column, distances are matched by source, target,
  and sample instead of using the cross-sample aggregate.

## Value

A Seurat object or spatial-communication result.
