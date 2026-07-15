# Publication figure specifications and export bundles

Shennong figures remain native `ggplot`, patchwork, or ComplexHeatmap
objects. The publication engine attaches a lightweight specification
instead of introducing a new plotting class, so ordinary `+ theme(...)`
composition keeps working.

## Profiles and automatic sizing

Generic profiles describe screen, column, page, and slide constraints.
They do not claim to replace current journal artwork instructions;
verify those before submission.

``` r

sn_list_figure_profiles()
#> # A tibble: 6 × 7
#>   profile     width_mm height_mm max_height_mm base_font_pt min_font_pt
#>   <chr>          <dbl>     <dbl>         <dbl>        <dbl>       <dbl>
#> 1 screen           160       100           180           10           8
#> 2 single_col…       85        NA           120            7           6
#> 3 one_half_c…      120        NA           160            7           6
#> 4 double_col…      180        NA           220            7           6
#> 5 supplement…      180        NA           240            7           6
#> 6 slide_16_9       254       143           143           16          12
#> # ℹ 1 more variable: raster_dpi <dbl>

large_spec <- sn_figure_spec(
  plot_type = "embedding",
  data_summary = list(
    n_points = 5000000,
    n_groups = 18,
    n_panels = 1,
    max_label_chars = 16
  ),
  profile = "single_column"
)
large_spec$recommended
#> $width_mm
#> [1] 85
#> 
#> $height_mm
#> [1] 87
#> 
#> $point_size
#> [1] 0.061
#> 
#> $alpha
#> [1] 0.55
#> 
#> $font_size_pt
#> [1] 7
#> 
#> $line_width
#> [1] 0.35
#> 
#> $legend_position
#> [1] "bottom"
#> 
#> $rasterize
#> [1] TRUE
#> 
#> $raster_dpi
#> [1] 600
#> 
#> $panel_layout
#> $panel_layout$columns
#> [1] 1
#> 
#> $panel_layout$rows
#> [1] 1
#> 
#> 
#> $pagination
#> $pagination$required
#> [1] FALSE
#> 
#> $pagination$pages
#> [1] 1
#> 
#> $pagination$per_page
#> [1] NA
#> 
#> $pagination$axis
#> [1] NA
```

The five-million-point example calculates metadata only; it does not
allocate or render five million rows. Point size follows a clamped power
rule, and high-density layers are recommended for rasterization while
text and geometry remain vector in PDF/SVG output.

## Inspect and apply a specification

Core Shennong plots attach `shennong_figure_spec` automatically. Any
native plot can also be profiled after construction.

``` r

plot <- ggplot2::ggplot(
  mtcars,
  ggplot2::aes(wt, mpg, color = factor(cyl))
) +
  ggplot2::geom_point()

plot <- sn_apply_figure_profile(plot, "single_column")
sn_figure_spec(plot)$recommended
#> $width_mm
#> [1] 85
#> 
#> $height_mm
#> [1] 85
#> 
#> $point_size
#> [1] 1.6
#> 
#> $alpha
#> [1] 0.85
#> 
#> $font_size_pt
#> [1] 7
#> 
#> $line_width
#> [1] 0.35
#> 
#> $legend_position
#> [1] "right"
#> 
#> $rasterize
#> [1] FALSE
#> 
#> $raster_dpi
#> [1] 600
#> 
#> $panel_layout
#> $panel_layout$columns
#> [1] 1
#> 
#> $panel_layout$rows
#> [1] 1
#> 
#> 
#> $pagination
#> $pagination$required
#> [1] FALSE
#> 
#> $pagination$pages
#> [1] 1
#> 
#> $pagination$per_page
#> [1] NA
#> 
#> $pagination$axis
#> [1] NA
sn_validate_figure(plot)
#> Figure validation: PASS
#> - width: 85 mm
#> - height: 85 mm
#> - font_size: 7 pt
#> - raster_dpi: 600 dpi
#> - panel_width: 85 mm
#> - panel_height: 85 mm
```

Validation reports calculated dimensions, text and raster settings,
small panels, category/legend overload, long labels, pagination, extreme
aspect ratios, and dense networks. It predicts likely problems; final
visual review is still required.

## Deterministic export

[`sn_save_figure()`](https://songqi.org/shennong/dev/reference/sn_save_figure.md)
supports PDF, SVG, TIFF, and PNG without depending on the current
RStudio graphics device. Explicit dimensions override automatic values.

``` r

sn_save_figure(
  plot,
  "figures/Figure_1A.pdf",
  profile = "single_column"
)
```

For a reproducible handoff, export the rendering, source table,
calculated specification, session information, checksums, and validation
result together.

``` r

manifest <- sn_export_figure_bundle(
  plot,
  path = "figures/Figure_1A",
  formats = c("pdf", "png"),
  include_data = TRUE,
  include_spec = TRUE,
  include_session = TRUE,
  profile = "single_column"
)
manifest$files
```

Result-aware
[`sn_plot_de()`](https://songqi.org/shennong/dev/reference/sn_plot_de.md),
[`sn_plot_enrichment()`](https://songqi.org/shennong/dev/reference/sn_plot_enrichment.md),
and
[`sn_plot_gsea()`](https://songqi.org/shennong/dev/reference/sn_plot_gsea.md)
attach their standardized plotting tables as bundle source data. The QC,
clustering, integration, annotation projection, bulk, dimensional,
feature, dot, heatmap, violin, box, bar, composition, and Milo helpers
also expose figure metadata for automatic sizing and QA.
