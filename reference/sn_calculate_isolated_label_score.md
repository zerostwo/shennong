# Calculate isolated-label preservation scores

This helper focuses on rare or low-frequency labels and summarizes how
well they remain separated in the selected embedding. The score is based
on the mean silhouette width of isolated labels and is scaled to
`[0, 1]` where larger values indicate better preservation.

## Usage

``` r
sn_calculate_isolated_label_score(
  x,
  label,
  reduction = .sn_default_metric_reduction(x),
  dims = NULL,
  cells = NULL,
  max_cells = 3000,
  stratify_by = label,
  isolated_fraction = 0.05,
  isolated_n = 100,
  seed = 717
)
```

## Arguments

- x:

  A Seurat object.

- label:

  Metadata column containing biological labels.

- reduction:

  Reduction name used to extract embeddings. Defaults to `"harmony"`
  when present, otherwise `"pca"`.

- dims:

  Optional integer vector of embedding dimensions to retain.

- cells:

  Optional character vector of cell names to score.

- max_cells:

  Optional integer cap used to subsample cells before running the
  metric. Defaults to `3000` because silhouette needs a full distance
  matrix.

- stratify_by:

  Optional metadata column used to preserve representation during
  subsampling. Defaults to `label`.

- isolated_fraction:

  Fraction-of-cells threshold used to flag isolated labels.

- isolated_n:

  Absolute cell-count threshold used to flag isolated labels.

- seed:

  Random seed used when `max_cells` triggers subsampling.

## Value

A data frame with one row per label and columns describing label
abundance, silhouette separation, and whether the label is considered
isolated. The attributes `overall_score` and `isolated_labels` summarize
the isolated-label subset.

## Examples

``` r
if (FALSE) { # \dontrun{
data("pbmc_small", package = "Shennong")
pbmc <- sn_run_cluster(
  pbmc_small,
  batch = "sample",
  species = "human",
  verbose = FALSE
)
isolated_tbl <- sn_calculate_isolated_label_score(
  pbmc,
  label = "seurat_clusters",
  reduction = "harmony"
)
isolated_tbl
} # }
```
