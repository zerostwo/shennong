# Retrieve a stored DE result from a Seurat object

Retrieve a stored DE result from a Seurat object

## Usage

``` r
sn_get_de_result(
  object,
  de_name = "default",
  top_n = NULL,
  direction = c("all", "up", "down"),
  groups = NULL,
  with_metadata = FALSE
)
```

## Arguments

- object:

  A `Seurat` object.

- de_name:

  Name of the stored DE result.

- top_n:

  Optional number of rows to keep. When supplied together with a ranking
  column, results are reduced to the top rows overall or per group.

- direction:

  One of `"all"`, `"up"`, or `"down"`.

- groups:

  Optional subset of group labels to keep.

- with_metadata:

  If `TRUE`, return the full stored result list instead of just the
  result table.

## Value

A tibble or stored-result list.

## Examples

``` r
if (FALSE) { # \dontrun{
markers <- sn_get_de_result(seurat_obj, de_name = "cluster_markers", top_n = 5)
} # }
```
