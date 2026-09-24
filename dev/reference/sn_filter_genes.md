# Filter genes based on the number of expressing cells

This function filters genes in a Seurat object based on the number of
cells in which they are expressed. Optionally, it can visualize the
effect of different filtering thresholds and retain only genes matching
bundled GENCODE-based annotation classes.

## Usage

``` r
sn_filter_genes(
  x,
  min_cells = 3,
  plot = TRUE,
  filter = TRUE,
  assay = "RNA",
  layer = "counts",
  species = NULL,
  gene_class = NULL,
  gene_type = NULL,
  object = NULL
)
```

## Arguments

- x:

  A Seurat object.

- min_cells:

  An integer specifying the minimum number of cells in which a gene must
  be expressed to be retained. Default is 3.

- plot:

  Logical; if TRUE, a bar plot is generated showing the number of
  remaining genes at different filtering thresholds. Default is TRUE.

- filter:

  Logical; if TRUE, returns the filtered Seurat object. If FALSE,
  returns the original object. Default is TRUE.

- assay:

  Assay used when extracting expression values. Defaults to `"RNA"`.

- layer:

  Layer used for gene filtering. Defaults to `"counts"`.

- species:

  Optional species label used for annotation-aware filtering. Required
  only when it cannot be inferred from the object.

- gene_class:

  Optional coarse annotation class to retain. One of `"coding"` or
  `"noncoding"`.

- gene_type:

  Optional character vector of exact GENCODE `gene_type` values to
  retain.

- object:

  Alias for `x`; supply only one of `x` and `object`.

## Value

A filtered Seurat object if `filter = TRUE`, otherwise the original
object.

## Details

The function computes the number of cells expressing each gene in the
Seurat object and filters out genes expressed in fewer than `min_cells`
cells. If `plot = TRUE`, it visualizes the effect of filtering using
`ggplot2`. When `gene_class` or `gene_type` is supplied, the function
also matches the feature set against the bundled
`shennong_gene_annotations` data and keeps only the requested annotation
subset.

## Examples

``` r
if (FALSE) { # \dontrun{
  pbmc <- qs2::qs_read(file.path(
    Sys.getenv("SHENNONG_REAL_DATA_DIR"), "single-cell", "kotliarov_pbmc.qs2"
  ))
  pbmc_filtered <- sn_filter_genes(
    pbmc,
    min_cells = 5,
    plot = FALSE,
    filter = TRUE
  )
  pbmc_coding <- sn_filter_genes(
    pbmc,
    min_cells = 1,
    plot = FALSE,
    filter = TRUE,
    species = "human",
    gene_class = "coding"
  )
} # }
```
