# Identify rare or difficult-to-separate groups

This helper summarizes group size, local neighbor purity, graph
connectivity, and silhouette width for a grouping column. It is intended
to surface rare populations or poorly separated groups that may be
obscured in standard UMAP visual inspection.

## Usage

``` r
sn_identify_challenging_groups(
  x,
  group_by = NULL,
  graph = NULL,
  reduction = .sn_default_metric_reduction(x),
  dims = NULL,
  cells = NULL,
  k = 20,
  neighbor_method = c("auto", "graph", "annoy", "exact"),
  max_cells = 3000,
  stratify_by = NULL,
  rare_fraction = 0.02,
  rare_n = 50,
  challenge_threshold = 0.5,
  seed = 717,
  n_trees = 50,
  group = NULL
)
```

## Arguments

- x:

  A Seurat object.

- group_by:

  Metadata column defining the groups to inspect.

- graph:

  Optional graph name stored in `x@graphs`. If `NULL`, the function
  tries to reuse an existing nearest-neighbor graph and falls back to a
  kNN graph built from the selected embedding.

- reduction:

  Reduction name used when a graph must be built from embeddings.
  Defaults to `"harmony"` when present, otherwise `"pca"`.

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
  subsampling. Defaults to `group_by`.

- rare_fraction:

  Fraction-of-cells threshold for flagging rare groups.

- rare_n:

  Absolute cell-count threshold for flagging rare groups.

- challenge_threshold:

  Threshold on the derived `challenge_score` used to flag difficult
  groups.

- seed:

  Random seed used when `max_cells` triggers subsampling.

- n_trees:

  Number of Annoy trees when `neighbor_method = "annoy"`.

- group:

  Deprecated alias for `group_by`.

## Value

A data frame with one row per group and rarity/separation flags.

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
difficult_tbl <- sn_identify_challenging_groups(
  pbmc,
  group_by = "seurat_clusters",
  reduction = "harmony"
)
difficult_tbl
} # }
```
