# Plot program activity distributions

Plot program activity distributions

## Usage

``` r
sn_plot_program_activity(x, name, programs = NULL, group_by = NULL)
```

## Arguments

- x:

  A Seurat object or program-scoring result.

- name:

  Stored result name.

- programs:

  Optional programs to keep.

- group_by:

  Optional Seurat metadata column used on the x-axis.

## Value

A `ggplot` object.

## Examples

``` r
if (FALSE) sn_plot_program_activity(object, "immune_programs", group_by = "cell_type") # \dontrun{}
```
