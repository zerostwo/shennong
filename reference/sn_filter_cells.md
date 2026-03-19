# Filter cells in a Seurat object based on QC metrics

This function filters cells in a Seurat object using Median Absolute
Deviation (MAD) to identify outliers across specified metadata features.
Supports grouped analysis and provides visual diagnostics.

## Usage

``` r
sn_filter_cells(
  x,
  features,
  group_by = NULL,
  method = "mad",
  n = 5,
  plot = TRUE,
  filter = TRUE
)
```

## Arguments

- x:

  A Seurat object

- features:

  Character vector of metadata column names to use for filtering

- group_by:

  (Optional) Metadata column name to group by for group-wise
  calculations

- method:

  Outlier detection method (currently only "mad" supported)

- n:

  Numeric threshold(s) for MAD multiplier. Can be single value or vector
  matching features.

- plot:

  Logical indicating whether to generate QC diagnostic plots

- filter:

  Logical indicating whether to filter out flagged cells

## Value

Seurat object with QC flags in metadata. If filter=TRUE, returns
subsetted object.

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage
pbmc <- sn_filter_cells(pbmc, features = c("nCount_RNA", "nFeature_RNA"))

# Grouped analysis
pbmc <- sn_filter_cells(pbmc, features = "percent.mt", group_by = "sample")

# Custom threshold with plotting
pbmc <- sn_filter_cells(pbmc, features = "nCount_RNA", n = 3, plot = TRUE)
} # }
```
