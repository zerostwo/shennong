# Calculate graph connectivity for a grouping label

Graph connectivity quantifies whether cells from the same group remain
connected in a neighbor graph. It is widely used to evaluate biological
conservation after integration.

## Usage

``` r
sn_calculate_graph_connectivity(
  x,
  label_by = NULL,
  graph = NULL,
  reduction = "pca",
  dims = NULL,
  cells = NULL,
  k = 20,
  neighbor_method = c("auto", "graph", "annoy", "exact"),
  max_cells = NULL,
  stratify_by = NULL,
  seed = 717,
  n_trees = 50,
  label = NULL
)
```

## Arguments

- x:

  A Seurat object.

- label_by:

  Metadata column used to define groups.

- graph:

  Optional graph name stored in `x@graphs`. If `NULL`, the function
  tries to reuse an existing nearest-neighbor graph and falls back to a
  kNN graph built from the selected embedding.

- reduction:

  Reduction name used when a graph must be built from embeddings.
  Defaults to `"pca"`.

- dims:

  Optional integer vector of embedding dimensions to retain.

- cells:

  Optional character vector of cell names to include.

- k:

  Number of neighbors used when building a graph from embeddings.

- neighbor_method:

  Strategy used when a graph must be built. One of `"auto"`, `"graph"`,
  `"annoy"`, or `"exact"`.

- max_cells:

  Optional integer cap used to subsample cells before running the
  metric.

- stratify_by:

  Optional metadata column used to preserve representation during
  subsampling. Defaults to `label`.

- seed:

  Random seed used when `max_cells` triggers subsampling.

- n_trees:

  Number of Annoy trees when `neighbor_method = "annoy"`.

- label:

  Deprecated alias for `label_by`.

## Value

A data frame with one row per group and a `connectivity_score` column in
`[0, 1]`.

## Examples

``` r
if (FALSE) { # \dontrun{
data("pbmc_small", package = "Shennong")
pbmc <- sn_run_cluster(
  pbmc_small,
  batch_by = "sample",
  species = "human",
  verbose = FALSE
)
connectivity_tbl <- sn_calculate_graph_connectivity(
  pbmc,
  label_by = "seurat_clusters",
  reduction = "harmony"
)
connectivity_tbl
} # }
```
