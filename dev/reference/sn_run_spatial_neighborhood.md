# Analyze spatial neighborhoods

Analyze spatial neighborhoods

## Usage

``` r
sn_run_spatial_neighborhood(
  object,
  method = c("knn", "squidpy"),
  group_by,
  spatial_cols = NULL,
  result_id = "spatial_neighborhood",
  backend_control = list(),
  return_object = TRUE,
  sample_by = NULL
)
```

## Arguments

- object:

  A Seurat object with coordinates and labels.

- method:

  Local k-nearest-neighbor analysis or an explicit Squidpy adapter.

- group_by:

  Metadata labels used for enrichment and co-occurrence.

- spatial_cols:

  Coordinate metadata columns.

- result_id:

  Stored result name.

- backend_control:

  Graph/permutation controls or `runner`/`result`.

- return_object:

  Return the modified object or result.

- sample_by:

  Optional metadata column defining independent samples or tissue
  sections. Neighbors and label permutations never cross boundaries.

## Value

A Seurat object or spatial-neighborhood result.
