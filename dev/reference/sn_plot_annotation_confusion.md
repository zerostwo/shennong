# Plot annotation confusion against known labels

Plot annotation confusion against known labels

## Usage

``` r
sn_plot_annotation_confusion(
  x,
  truth,
  store_name = "annotation",
  normalize = TRUE
)
```

## Arguments

- x:

  A Seurat object or annotation result.

- truth:

  Known labels. For a Seurat object this can be a metadata column name;
  otherwise it must be a vector aligned to the cell table.

- store_name:

  Stored annotation name.

- normalize:

  Normalize rows to proportions.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) sn_plot_annotation_confusion(object, truth = "known_cell_type") # \dontrun{}
```
