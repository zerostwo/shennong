# Plot program activity distributions

Plot program activity distributions

## Usage

``` r
sn_plot_program_activity(
  x,
  result_id,
  programs = NULL,
  group_by = NULL,
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

  Optional Seurat metadata column used on the x-axis.

- object:

  Alias for `x`; supply only one of `x` and `object`.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) sn_plot_program_activity(object, "immune_programs", group_by = "cell_type") # \dontrun{}
```
