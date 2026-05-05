# Calculate PCR batch_by effect scores

This function estimates how much variance in an embedding is still
explained by batch. If a baseline reduction or baseline object is
supplied, it also reports the improvement relative to the unintegrated
state.

## Usage

``` r
sn_calculate_pcr_batch(
  x,
  batch_by = NULL,
  reduction = .sn_default_metric_reduction(x),
  dims = NULL,
  baseline = NULL,
  baseline_reduction = NULL,
  cells = NULL,
  max_cells = NULL,
  stratify_by = NULL,
  seed = 717,
  batch = NULL
)
```

## Arguments

- x:

  A Seurat object.

- batch_by:

  Metadata column containing batch labels.

- reduction:

  Reduction name used for the primary score. Defaults to `"harmony"`
  when present, otherwise `"pca"`.

- dims:

  Optional integer vector of embedding dimensions to retain.

- baseline:

  Optional Seurat object used as the baseline reference. If `NULL`, the
  baseline is taken from `x`.

- baseline_reduction:

  Optional reduction name used as the baseline.

- cells:

  Optional character vector of cell names to include.

- max_cells:

  Optional integer cap used to subsample cells before running the
  metric.

- stratify_by:

  Optional metadata column used to preserve representation during
  subsampling. Defaults to `batch`.

- seed:

  Random seed used when `max_cells` triggers subsampling.

- batch:

  Deprecated alias for `batch_by`.

## Value

A one-row data frame containing the weighted batch_by variance explained
by the selected reduction and, when available, the baseline comparison.

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
sn_calculate_pcr_batch(
  pbmc,
  batch_by = "sample",
  reduction = "harmony",
  baseline_reduction = "pca"
)
} # }
```
