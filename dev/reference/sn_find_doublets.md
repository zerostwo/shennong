# Find doublets with scDblFinder or Scrublet

This function identifies potential doublets in a Seurat object. The
default `method = "scdblfinder"` converts the object to a
SingleCellExperiment and runs the `scDblFinder` package; for the
validated scDblFinder default workflow on an in-memory `dgCMatrix` with
at most 33,000 cells, Shennong lazily registers and uses its guarded
ShennongOpt fast path when available, falling back to the captured
upstream implementation otherwise (compatibility comes from the targeted
function fingerprints rather than the package version label).
`method = "scrublet"` runs Scrublet through scanpy's native wrapper
inside the managed `scrublet` pixi environment and always consumes raw
counts.

## Usage

``` r
sn_find_doublets(
  object,
  method = c("scdblfinder", "scrublet"),
  clusters = NULL,
  cluster_backend = c("native", "shennong"),
  group_by = NULL,
  dbr_sd = NULL,
  n_workers = 1,
  assay = "RNA",
  layer = "counts",
  min_features = 200,
  backend_control = list()
)
```

## Arguments

- object:

  A `Seurat` object.

- method:

  Doublet-detection backend. One of `"scdblfinder"` (default) or
  `"scrublet"`.

- clusters:

  Optional cluster assignments. A metadata column name or one value per
  cell. Explicit assignments take precedence over `cluster_backend`.
  scDblFinder-specific; ignored by scrublet.

- cluster_backend:

  Clustering implementation used when `clusters` is `NULL`. `"native"`
  (the default) lets scDblFinder perform its own automatic clustering.
  `"shennong"` first obtains assignments from
  [`sn_run_cluster()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_cluster.md)
  using the selected assay and layer. scDblFinder-specific; ignored by
  scrublet.

- group_by:

  An optional metadata column used as the donor or sample grouping.
  scDblFinder-specific; ignored by scrublet.

- dbr_sd:

  A numeric value for adjusting the doublet rate; see `scDblFinder`
  documentation. scDblFinder-specific; ignored by scrublet.

- n_workers:

  Number of sample groups to process concurrently. For a BPCells
  backend, `n_workers = 1` materializes one sample-sized sparse matrix
  at a time and gives the lowest peak memory use. Higher values can hold
  up to `n_workers` sample matrices in memory concurrently.
  scDblFinder-specific; ignored by scrublet.

- assay:

  Assay used for doublet detection. Defaults to `"RNA"`.

- layer:

  Layer used as the input count matrix. Defaults to `"counts"`. Scrublet
  requires raw counts: a non-`"counts"` layer stops with an error
  instead of silently scoring normalized values.

- min_features:

  Minimum number of detected features required for a cell to enter
  doublet detection. Defaults to `200`. Cells below this threshold are
  skipped and retain `NA` in the output columns.

- backend_control:

  Optional named list of backend-specific settings. For scrublet:
  `seed`, `n_prin_comps`, `expected_doublet_rate`,
  `synthetic_doublet_umi_subsampling`, `min_counts`, `min_cells`,
  `min_gene_variability_pctl`, `threshold_min`, `threshold_max`,
  `output_dir`, `runtime_dir`, `keep_run_dir`, and `quiet`. Default
  Scrublet runs use a unique package-owned temporary directory and
  remove it after import. An explicit `output_dir` is retained only when
  `keep_run_dir = TRUE` (the default for an explicit path); with
  `keep_run_dir = FALSE`, it is a parent for an isolated child and is
  never recursively deleted. Ignored by scDblFinder.

## Value

The input Seurat object with two new columns in `meta.data`:
`scDblFinder.class` and `scDblFinder.score` for `method = "scdblfinder"`
(with a `_corrected` suffix on non-default corrected layers), or
`scrublet.class` and `scrublet.score` for `method = "scrublet"`. Cells
whose selected layer sums to zero or whose detected-feature count is
below `min_features` are skipped and retain `NA` in the corresponding
output columns. Because `scDblFinder()` requires an in-memory sparse
matrix, BPCells-backed inputs must supply `group_by`; Shennong then
materializes each sample independently instead of the complete matrix.
Scrublet materializes the retained cells once and scores them through
scanpy's native `sc.pp.scrublet()` wrapper.

## Examples

``` r
if (FALSE) { # \dontrun{
seurat_obj <- sn_find_doublets(seurat_obj, n_workers = 4)

# Scrublet through the managed pixi environment:
sn_prepare_pixi_environment("scrublet", install_environment = TRUE)
seurat_obj <- sn_find_doublets(
  seurat_obj,
  method = "scrublet",
  backend_control = list(seed = 717)
)
} # }
```
