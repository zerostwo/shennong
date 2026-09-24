# Compare integration methods with scib-metrics

Benchmarks the integration reductions registered by a multi-method
[`sn_run_cluster()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_cluster.md)
call. UMAP and t-SNE coordinates are never used as metric inputs. The
normalized expression matrix is exported sparsely and every method is
evaluated on the same cells. Parameter-grid objects are grouped by
preprocessing identity, each group uses its matching unintegrated
baseline, and each unique native embedding is scored once before its
score is mapped back to all cluster-resolution runs that share it.
Stored clustering runtime and peak-memory provenance is joined into the
scIB summary, metrics, and ranking tables so quality and compute cost
can be compared without a separate object lookup.

## Usage

``` r
sn_compare_integrations(
  object,
  batch_by = NULL,
  label_by,
  methods = NULL,
  assay = "RNA",
  normalized_layer = "data",
  features = NULL,
  max_cells = 100000L,
  accelerator = c("auto", "cpu", "gpu"),
  n_workers = 1L,
  result_id = "integration_benchmark",
  return_object = TRUE,
  backend_control = list(),
  seed = 717L,
  verbose = TRUE
)
```

## Arguments

- object:

  A Seurat object containing `misc$integration_comparison`.

- batch_by:

  Metadata column containing batch identities.

- label_by:

  Metadata column containing biological labels.

- methods:

  Optional method names or run IDs to evaluate. A method name selects
  all of its registered parameter-grid runs.

- assay:

  RNA assay containing normalized expression.

- normalized_layer:

  Normalized, unintegrated expression layer used by scib-metrics for its
  baseline preparation.

- features:

  Optional common features. Defaults to variable features.

- max_cells:

  Optional shared stratified cell cap.

- accelerator:

  One of `"auto"`, `"cpu"`, or `"gpu"`.

- n_workers:

  Number of parallel neighbor-search workers.

- result_id:

  Stable identifier for the stored benchmark result.

- return_object:

  Return the updated Seurat object when `TRUE`, otherwise return the
  analysis-result object.

- backend_control:

  Advanced pixi, runner, metric, and filesystem controls.

- seed:

  Sampling seed.

- verbose:

  Show backend progress.

## Value

A Seurat object with an `integration_benchmark` result, or the result.
