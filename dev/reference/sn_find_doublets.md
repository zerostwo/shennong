# Find doublets using scDblFinder

This function identifies potential doublets in a Seurat object by
converting it to a SingleCellExperiment and using the `scDblFinder`
package.

## Usage

``` r
sn_find_doublets(
  object,
  clusters = NULL,
  group_by = NULL,
  dbr_sd = NULL,
  ncores = 1,
  assay = "RNA",
  layer = "counts",
  min_features = 200
)
```

## Arguments

- object:

  A `Seurat` object.

- clusters:

  Optional cluster assignments. If not provided, scDblFinder will
  attempt automatic clustering.

- group_by:

  An optional metadata column used as the donor or sample grouping.

- dbr_sd:

  A numeric value for adjusting the doublet rate; see `scDblFinder`
  documentation.

- ncores:

  Number of cores to use (for parallel processing).

- assay:

  Assay used for doublet detection. Defaults to `"RNA"`.

- layer:

  Layer used as the input count matrix. Defaults to `"counts"`.

- min_features:

  Minimum number of detected features required for a cell to be passed
  to `scDblFinder()`. Defaults to `200`. Cells below this threshold are
  skipped and retain `NA` in the output columns.

## Value

The input Seurat object with two new columns in `meta.data`:
`scDblFinder.class` and `scDblFinder.score` when `layer = "counts"`, or
`scDblFinder.class_corrected` and `scDblFinder.score_corrected` for
non-default corrected layers. Cells whose selected layer sums to zero or
whose detected-feature count is below `min_features` are skipped and
retain `NA` in the corresponding output columns.

## Examples

``` r
if (FALSE) { # \dontrun{
seurat_obj <- sn_find_doublets(
  seurat_obj,
  clusters = NULL,
  group_by = NULL,
  dbr_sd = NULL,
  ncores = 4
)
} # }
```
