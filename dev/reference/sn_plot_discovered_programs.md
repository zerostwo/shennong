# Plot discovered gene programs

Plot discovered gene programs

## Usage

``` r
sn_plot_discovered_programs(
  x,
  name = NULL,
  type = c("weights", "activity", "stability"),
  programs = NULL,
  n = 20L
)
```

## Arguments

- x:

  A Seurat object or program-discovery result.

- name:

  Stored result name when `x` is a Seurat object.

- type:

  Gene weights, per-cell activity, or run stability.

- programs:

  Optional programs to retain.

- n:

  Maximum genes per program.

## Value

A `ggplot` object.
