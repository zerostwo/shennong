# Retrieve a stored cell-cell communication result

Retrieve a stored cell-cell communication result

## Usage

``` r
sn_get_cell_communication_result(
  object,
  communication_name = "default",
  sources = NULL,
  targets = NULL,
  with_metadata = FALSE
)
```

## Arguments

- object:

  A Seurat object.

- communication_name:

  Name of the stored result.

- sources, targets:

  Optional source/target labels to keep.

- with_metadata:

  If `TRUE`, return the full stored-result list.

## Value

A tibble or stored-result list.
