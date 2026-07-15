# Save a publication figure with deterministic dimensions

Save a publication figure with deterministic dimensions

## Usage

``` r
sn_save_figure(
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

## Value

The normalized output path, invisibly.
