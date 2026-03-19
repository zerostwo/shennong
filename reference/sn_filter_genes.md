# Filter genes based on the number of expressing cells

This function filters genes in a Seurat object based on the number of
cells in which they are expressed. Optionally, it can visualize the
effect of different filtering thresholds.

## Usage

``` r
sn_filter_genes(
  x,
  min_cells = 3,
  plot = TRUE,
  filter = TRUE,
  assay = "RNA",
  layer = "counts"
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

## Value

A filtered Seurat object if `filter = TRUE`, otherwise the original
object.

## Details

The function computes the number of cells expressing each gene in the
Seurat object and filters out genes expressed in fewer than `min_cells`
cells. If `plot = TRUE`, it visualizes the effect of filtering using
`ggplot2`.

## Examples

``` r
library(Seurat)
#> Loading required package: SeuratObject
#> Loading required package: sp
#> 
#> Attaching package: ‘SeuratObject’
#> The following objects are masked from ‘package:base’:
#> 
#>     intersect, t

# Load example Seurat object
pbmc_small_filtered <- sn_filter_genes(pbmc_small, min_cells = 5, plot = TRUE, filter = TRUE)

```
