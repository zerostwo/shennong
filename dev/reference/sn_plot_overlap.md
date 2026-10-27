# Plot exact set intersections as an UpSet panel

Counts exclusive intersections of supplied identifiers, with
intersection bars above a connected membership matrix. Each identifier
contributes to exactly one intersection. Empty intersections are
omitted, not imputed.

## Usage

``` r
sn_plot_overlap(
  sets,
  style = "upset",
  sort_by = c("size", "degree"),
  color = "#BD5050",
  title = NULL
)
```

## Arguments

- sets:

  Named list of character vectors; entries must be unique within each
  set. Empty sets are permitted. No background universe is inferred.

- style:

  `upset`, for exclusive-intersection bars and membership dots.

- sort_by:

  `size` (descending intersection count) or `degree` (descending number
  of participating sets, then count).

- color:

  Color of bars and active membership dots.

- title:

  Optional title.

## Value

A patchwork retaining membership and exact intersection tables.

## Examples

``` r
sn_plot_overlap(list(B = c("S1", "S2"), T = c("S2", "S3"),
                     Myeloid = c("S2", "S4")))
```
