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
  layer = "counts"
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

## Value

The input Seurat object with two new columns in `meta.data`:
`scDblFinder.class` and `scDblFinder.score`.

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
