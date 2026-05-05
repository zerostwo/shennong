# Sweep clustering resolutions and summarize cluster-quality diagnostics

This helper reruns graph clustering across a grid of Seurat resolution
values and evaluates each partition with one or more quality
diagnostics. It is intended to support empirical resolution selection
rather than relying on a fixed default cluster_by count.

## Usage

``` r
sn_sweep_cluster_resolution(
  x,
  resolutions = seq(0.2, 1.2, by = 0.2),
  reduction = .sn_default_metric_reduction(x),
  dims = NULL,
  cluster_name = "seurat_clusters",
  metrics = c("n_clusters", "mean_silhouette", "graph_connectivity"),
  label_by = NULL,
  max_cells = 3000,
  rogue_max_cells = max_cells,
  assay = "RNA",
  layer = "counts",
  neighbor_method = c("auto", "graph", "annoy", "exact"),
  k = 20,
  seed = 717,
  n_trees = 50,
  label = NULL
)
```

## Arguments

- x:

  A Seurat object containing the reduction used for clustering.

- resolutions:

  Numeric vector of Seurat resolution values to evaluate.

- reduction:

  Reduction name used for neighbor construction and optional silhouette
  scoring. Defaults to `"harmony"` when present, otherwise `"pca"`.

- dims:

  Optional integer vector of embedding dimensions used by
  `FindNeighbors()`, `FindClusters()`, and silhouette scoring.

- cluster_name:

  Metadata column used to store the temporary cluster assignments.
  Defaults to `"seurat_clusters"`.

- metrics:

  Metrics to calculate for each resolution. Supported values are
  `"n_clusters"`, `"mean_silhouette"`, `"graph_connectivity"`,
  `"cluster_purity"`, `"clustering_agreement"`, and `"rogue"`.

- label_by:

  Optional metadata column containing reference labels. Required for
  `"cluster_purity"` and `"clustering_agreement"`.

- max_cells:

  Optional integer cap used by expensive embedding-based metrics such as
  silhouette.

- rogue_max_cells:

  Optional integer cap used by `"rogue"`.

- assay:

  Assay used when `"rogue"` is requested.

- layer:

  Layer used when `"rogue"` is requested.

- neighbor_method:

  Strategy used when graph connectivity must rebuild a graph from
  embeddings.

- k:

  Number of neighbors used for graph connectivity when rebuilding a
  graph from embeddings.

- seed:

  Random seed forwarded to Seurat clustering and subsampling.

- n_trees:

  Number of Annoy trees used by graph connectivity when
  `neighbor_method = "annoy"`.

- label:

  Deprecated alias for `label_by`.

## Value

A list with a per-resolution summary table and the recommended
resolution based on the mean of the available scaled quality metrics.

## Examples

``` r
if (FALSE) { # \dontrun{
sweep <- sn_sweep_cluster_resolution(
  pbmc,
  resolutions = seq(0.2, 1.0, by = 0.2),
  reduction = "pca",
  dims = 1:20
)
sweep$summary
sweep$recommended_resolution
} # }
```
