# Assess integration quality across multiple metrics

This wrapper combines fast, broadly useful metrics for batch_by mixing,
biological conservation, and cluster-level diagnostics. It reuses
existing neighbor graphs when available and otherwise builds a kNN graph
with Annoy or an exact fallback.

## Usage

``` r
sn_assess_integration(
  x,
  batch_by = NULL,
  label_by = NULL,
  cluster_by = NULL,
  reduction = .sn_default_metric_reduction(x),
  baseline = NULL,
  baseline_reduction = .sn_default_baseline_reduction(x, reduction),
  graph = NULL,
  dims = NULL,
  k = 20,
  metrics = NULL,
  neighbor_method = c("auto", "graph", "annoy", "exact"),
  max_cells = 5000,
  stratify_by = NULL,
  rare_fraction = 0.02,
  rare_n = 50,
  challenge_threshold = 0.5,
  seed = 717,
  n_trees = 50,
  batch = NULL,
  label = NULL,
  cluster = NULL
)
```

## Arguments

- x:

  A Seurat object.

- batch_by:

  Metadata column containing batch labels.

- label_by:

  Optional metadata column containing biological labels such as cell
  type annotations.

- cluster_by:

  Optional metadata column containing the current clustering.

- reduction:

  Reduction name used for the primary integrated embedding. Defaults to
  `"harmony"` when present, otherwise `"pca"`.

- baseline:

  Optional Seurat object used as the baseline reference for PCR batch_by
  scoring. If `NULL`, the baseline is taken from `x`.

- baseline_reduction:

  Optional reduction name used as the baseline.

- graph:

  Optional graph name stored in `x@graphs`.

- dims:

  Optional integer vector of embedding dimensions to retain.

- k:

  Number of neighbors used when a graph must be built from embeddings.

- metrics:

  Optional character vector selecting which metrics to compute.
  Supported values are `"batch_lisi"`, `"label_lisi"`,
  `"batch_silhouette"`, `"label_silhouette"`, `"graph_connectivity"`,
  `"clustering_agreement"`, `"isolated_label_score"`,
  `"cluster_purity"`, `"cluster_entropy"`, `"pcr_batch"`, and
  `"challenging_groups"`.

- neighbor_method:

  Strategy used when a graph must be built. One of `"auto"`, `"graph"`,
  `"annoy"`, or `"exact"`.

- max_cells:

  Optional integer cap used to subsample cells before running the
  metrics.

- stratify_by:

  Optional metadata column used to preserve representation during
  subsampling. Defaults to `batch`.

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

- batch:

  Deprecated alias for `batch_by`.

- label:

  Deprecated alias for `label_by`.

- cluster:

  Deprecated alias for `cluster_by`.

## Value

A list with four top-level elements:

- `summary`: one row per summary metric, plus aggregate scores.

- `per_cell`: per-cell outputs such as LISI and silhouette.

- `per_group`: per-group outputs such as connectivity, composition, and
  difficult-group diagnostics.

- `parameters`: the effective runtime configuration.

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
metrics <- sn_assess_integration(
  pbmc,
  batch_by = "sample",
  cluster_by = "seurat_clusters",
  reduction = "harmony",
  baseline_reduction = "pca"
)
metrics$summary
} # }
```
