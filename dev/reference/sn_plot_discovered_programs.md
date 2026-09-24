# Plot discovered gene programs

Plot discovered gene programs

## Usage

``` r
sn_plot_discovered_programs(
  x,
  result_id = NULL,
  type = c("weights", "activity", "stability"),
  programs = NULL,
  n = 20L,
  object = NULL
)
```

## Arguments

- x:

  A Seurat object or program-discovery result.

- result_id:

  Stored result result_id when `x` is a Seurat object.

- type:

  Gene weights, per-cell activity, or run stability.

- programs:

  Optional programs to retain.

- n:

  Maximum genes per program.

- object:

  Alias for `x`; supply only one of `x` and `object`.

## Value

A `ggplot` object.
