# Retrieve a DE result table

Retrieve a DE result table

## Usage

``` r
sn_get_de_result(
  object,
  result_id = NULL,
  top_n = NULL,
  direction = c("all", "up", "down"),
  groups = NULL,
  top_scope = c("group", "all"),
  p_adjusted_cutoff = NULL,
  logfc_threshold = NULL
)
```

## Arguments

- object:

  A `Seurat` object or unified DE result.

- result_id:

  Identifier of the stored DE result. May be omitted when exactly one DE
  result is stored.

- top_n:

  Optional number of rows to keep. When supplied together with a ranking
  column, results are reduced to the top rows overall or per group.

- direction:

  One of `"all"`, `"up"`, or `"down"`.

- groups:

  Optional subset of group labels to keep.

- top_scope:

  Whether `top_n` applies within each group (the default) or across all
  selected rows. An ungrouped table always uses all rows.

- p_adjusted_cutoff:

  Optional maximum adjusted p-value, from zero to one.

- logfc_threshold:

  Optional minimum absolute log fold change.

## Value

A filtered tibble. Use `sn_get_result(object, "de")` to retrieve the
complete result, including metadata and unfiltered tables.

## Examples

``` r
if (FALSE) { # \dontrun{
markers <- sn_get_de_result(seurat_obj, result_id = "cluster_markers", top_n = 5)
} # }
```
