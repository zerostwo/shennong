# Add spatial distance evidence to a communication result

Add spatial distance evidence to a communication result

## Usage

``` r
sn_run_spatial_communication(
  object,
  communication_name = "communication",
  communication = NULL,
  group_by,
  spatial_cols = NULL,
  max_distance = NULL,
  store_name = "spatial_communication",
  return_object = TRUE
)
```

## Arguments

- object:

  A Seurat object.

- communication_name:

  Stored communication result name.

- communication:

  Optional communication result supplied directly.

- group_by:

  Metadata column matching communication source/target labels.

- spatial_cols:

  Coordinate metadata columns.

- max_distance:

  Optional maximum mean nearest-group distance.

- store_name:

  Stored result name.

- return_object:

  Return the modified object or result.

## Value

A Seurat object or spatial-communication result.
