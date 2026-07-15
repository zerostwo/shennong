# Export a publication figure

Alias of
[`sn_save_figure()`](https://songqi.org/shennong/dev/reference/sn_save_figure.md)
for workflows that use export terminology.

## Usage

``` r
sn_export_figure(
  plot,
  filename,
  profile = NULL,
  width = "auto",
  height = "auto",
  units = c("mm", "cm", "in"),
  dpi = "auto",
  background = "white",
  embed_fonts = TRUE,
  validate = TRUE
)
```

## Arguments

- plot:

  A ggplot, patchwork, or ComplexHeatmap object.

- filename:

  Output PDF, SVG, TIFF, or PNG path.

- profile:

  Optional output profile.

- width, height:

  Numeric dimensions or `"auto"`.

- units:

  Units for explicit dimensions.

- dpi:

  Numeric DPI or `"auto"`.

- background:

  Explicit output background.

- embed_fonts:

  Use Cairo PDF output for embedded/subsettable fonts.

- validate:

  Run figure QA before saving.
