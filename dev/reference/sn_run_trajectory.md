# Infer trajectories and test dynamic genes

Runs a cluster-aware Slingshot trajectory, stores per-cell pseudotime
and lineage probabilities, and optionally fits tradeSeq
negative-binomial GAMs. The returned result follows the unified Shennong
analysis-result contract.

## Usage

``` r
sn_run_trajectory(
  object,
  method = c("slingshot", "monocle3", "palantir"),
  reduction = NULL,
  cluster_by = NULL,
  start = NULL,
  end = NULL,
  lineages = NULL,
  store_name = "trajectory",
  dims = NULL,
  assay = NULL,
  counts_layer = "counts",
  test_dynamic = TRUE,
  dynamic_features = NULL,
  max_dynamic_features = 2000L,
  nknots = 6L,
  trend_features = NULL,
  trend_points = 100L,
  backend_control = list(),
  return_object = TRUE
)
```

## Arguments

- object:

  A Seurat object with a dimensional reduction and cluster labels.

- method:

  Trajectory backend. Slingshot and Monocle 3 run directly; Palantir
  accepts a standardized external runner or result.

- reduction:

  Reduction used for inference. Defaults to PCA, then UMAP or the first
  available reduction.

- cluster_by:

  Metadata column containing cluster labels. Defaults to
  `seurat_clusters`, then active identities.

- start, end:

  Optional start and terminal cluster labels.

- lineages:

  Optional named list of expected cluster paths. Their endpoints
  constrain Slingshot and inferred paths are checked against them.

- store_name:

  Name used under the `trajectory` result type.

- dims:

  Reduction dimensions used for inference.

- assay:

  Assay used for dynamic-gene counts.

- counts_layer:

  Raw-count layer used by tradeSeq.

- test_dynamic:

  Fit tradeSeq dynamic and branch tests.

- dynamic_features:

  Optional features tested by tradeSeq.

- max_dynamic_features:

  Maximum number of dynamic features fitted by default.

- nknots:

  Number of tradeSeq spline knots.

- trend_features:

  Features for which fitted trends are retained.

- trend_points:

  Number of fitted points per lineage and feature.

- backend_control:

  Backend controls. Use named `slingshot`, `monocle3`, and `tradeSeq`
  argument lists for direct backends. Palantir accepts `runner` or
  `result`; the same explicit adapter boundary can override Monocle 3
  for externally managed execution.

- return_object:

  Return the updated object instead of the result.

## Value

A Seurat object or a unified trajectory result.

## Examples

``` r
if (FALSE) { # \dontrun{
object <- sn_run_trajectory(
  object, reduction = "pca", cluster_by = "seurat_clusters",
  start = "0", store_name = "development"
)
result <- sn_get_result(object, "trajectory", "development")
result$tables$cells
} # }
```
