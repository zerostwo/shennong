# Resolve a palette into explicit colors

Resolve a palette into explicit colors

## Usage

``` r
sn_get_palette(
  palette = "Paired",
  n = NULL,
  palette_type = c("auto", "discrete", "continuous"),
  direction = 1
)
```

## Arguments

- palette:

  Palette name or explicit character vector of colors.

- n:

  Number of colors to return. When omitted, the palette's native length
  is returned for discrete use and `256` colors are returned for
  continuous use.

- palette_type:

  One of `"auto"`, `"discrete"`, or `"continuous"`. Defaults to
  `"auto"`.

- direction:

  Direction for ordered palettes. Use `1` for the default order and `-1`
  to reverse it.

## Value

A character vector of hex colors.

## Examples

``` r
sn_get_palette("Paired", n = 14)
sn_get_palette("RdBu", palette_type = "continuous", direction = -1)
```
