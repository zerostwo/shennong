# List available color palettes

Returns the built-in Shennong palettes together with the available
`RColorBrewer` palettes. By default it prints a swatch-style plot and
invisibly returns the underlying metadata data frame.

## Usage

``` r
sn_list_palettes(
  display = c("plot", "preview", "table", "none"),
  source = NULL,
  palette_type = NULL
)
```

## Arguments

- display:

  One of `"plot"`, `"preview"`, `"table"`, or `"none"`. Defaults to
  `"plot"`.

- source:

  Optional source filter such as `"shennong"`, `"ggokabeito"`, or
  `"RColorBrewer"`.

- palette_type:

  Optional palette-type filter.

## Value

Invisibly returns a data frame with palette names, source, palette type,
maximum native size, preview colors, and whether each palette supports
discrete and continuous use.

## Examples

``` r
sn_list_palettes()
sn_list_palettes(source = "ggokabeito", display = "table")
```
