# Calculate silhouette widths from a Seurat embedding

Silhouette widths summarize how well cells are separated by a
categorical metadata label in the selected embedding.

## Usage

``` r
sn_calculate_silhouette(
  x,
  label,
  reduction = "pca",
  dims = NULL,
  cells = NULL,
  max_cells = 3000,
  stratify_by = label,
  seed = 717
)
```

## Arguments

- x:

  A Seurat object.

- label:

  Metadata column used as the grouping label.

- reduction:

  Reduction name used to extract embeddings. Defaults to `"pca"`.

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

- seed:

  Random seed used when `max_cells` triggers subsampling.

## Value

A data frame with per-cell silhouette widths.

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
sil_tbl <- sn_calculate_silhouette(
  pbmc,
  label = "seurat_clusters",
  reduction = "harmony"
)
head(sil_tbl)
} # }
```
