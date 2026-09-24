# Plot a program activity heatmap

Plot a program activity heatmap

## Usage

``` r
sn_plot_program_heatmap(
  x,
  result_id,
  programs = NULL,
  group_by = NULL,
  scale_rows = TRUE,
  object = NULL
)
```

## Arguments

- x:

  A Seurat object or program-scoring result.

- result_id:

  Stored result result_id.

- programs:

  Optional programs to keep.

- group_by:

  Optional Seurat metadata column used to aggregate cell-level scores
  before plotting.

- scale_rows:

  Standardize each program across displayed groups.

- object:

  Alias for `x`; supply only one of `x` and `object`.

## Value

A `ggplot` heatmap.

## Examples

``` r
if (FALSE) sn_plot_program_heatmap(object, "immune_programs", group_by = "cell_type") # \dontrun{}
```
