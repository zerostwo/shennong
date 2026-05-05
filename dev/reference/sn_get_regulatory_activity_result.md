# Retrieve stored regulatory activity results

Retrieve stored regulatory activity results

## Usage

``` r
sn_get_regulatory_activity_result(
  object,
  activity_name = "default",
  sources = NULL,
  conditions = NULL,
  with_metadata = FALSE
)
```

## Arguments

- object:

  A Seurat object.

- activity_name:

  Name of the stored result.

- sources:

  Optional TF or pathway names to keep.

- conditions:

  Optional cell or group names to keep.

- with_metadata:

  If `TRUE`, return the full stored-result list.

## Value

A tibble or stored-result list.
