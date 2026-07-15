# Export a figure, source data, specification, and manifest bundle

Export a figure, source data, specification, and manifest bundle

## Usage

``` r
sn_export_figure_bundle(
  plot,
  path,
  formats = c("pdf", "png"),
  include_data = TRUE,
  include_spec = TRUE,
  include_session = TRUE,
  profile = NULL,
  ...
)
```

## Arguments

- plot:

  A publication figure.

- path:

  Bundle directory, whose base name is used for output files.

- formats:

  Vector containing PDF, SVG, TIFF, or PNG.

- include_data:

  Include available source data as CSV.

- include_spec:

  Include the calculated specification.

- include_session:

  Include [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html)
  output.

- profile:

  Optional output profile.

- ...:

  Arguments passed to
  [`sn_save_figure()`](https://songqi.org/shennong/dev/reference/sn_save_figure.md).

## Value

The manifest list, invisibly.
