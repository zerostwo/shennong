# Plotting stored results and interpretation summaries

Use only the stages needed by the requested analysis. Each section is an API
reference, not a requirement to run all listed methods.

## 8

Prefer `sn_plot_result()` for stored/downstream analysis results; inspect
   valid analysis/view pairs with `sn_list_plot_methods()`. Use
   `sn_plot_distribution()` for numeric violin/box/histogram/density/ridge
   views and `sn_plot_association()` for cell- or sample-level correlations.
   Use `sn_plot_dim()`, `sn_plot_feature()`, `sn_plot_heatmap()`,
   `sn_plot_violin()`, `sn_plot_dot()`, `sn_plot_boxplot()`,
   `sn_plot_barplot()`, `sn_plot_composition()`, and `sn_plot_milo()` for
   package-style plots; resolve reusable colors with `sn_list_palettes()` and
   `sn_get_palette()`.
   `sn_plot_composition()` accepts a Seurat object directly. Use ordinary
   `type = "bar"` for descriptive cell counts/proportions; use
   `type = "sample_bar"` or `"sample_boxplot"` with an explicit `sample_by`
   when biological samples are the replicates; use `type = "alluvial"` with
   two or more `flow_by` metadata columns for source/annotation flows. Use
   `unit_by` for donor metadata repeated across cells so donor distributions do
   not count cells as independent donors.

## 9

Build prompts or stored-result summaries with the interpretation helpers
   when a narrative or report-ready output is needed.

## Three-dimensional nebula and glass maps

`sn_plot_dim()` and `sn_plot_feature()` accept `style = "nebula"` (bright rims)
or `style = "glass"` (translucent envelopes). `style = "classic"` preserves the
existing plots. Both new styles require **three real embedding dimensions**;
plotting never adds a random Z coordinate or recomputes UMAP. Use a separate
reduction to preserve your original two-dimensional map:

```r
# Optional packages for this workflow:
# install.packages(c("misc3d", "htmlwidgets", "chromote", "png"))
obj <- Seurat::RunUMAP(
  obj, reduction = "pca", dims = 1:20, n.components = 3L,
  reduction.name = "umap3d", reduction.key = "UMAP3D_", seed.use = 717
)

viewer <- sn_plot_dim(
  obj, reduction = "umap3d", dims = 1:3, group_by = "cell_type",
  style = "glass", interactive = TRUE, label = TRUE
)
htmlwidgets::saveWidget(viewer, "umap-viewer.html", selfcontained = TRUE)
utils::browseURL(normalizePath("umap-viewer.html"))

# Drag to rotate; shift-drag to pan; scroll to zoom.
# Click Download camera JSON, then read the downloaded file:
camera <- sn_get_plot_camera("~/Downloads/shennong-camera.json")
# Alternatively, paste the list from Copy R camera.

p <- sn_plot_dim(
  obj, reduction = "umap3d", dims = 1:3, group_by = "cell_type",
  style = "glass", camera = camera, label = TRUE, raster_dpi = 600
)
ggplot2::ggsave("umap-glass.pdf", p, width = 6, height = 5, dpi = 600)

# Reuse exactly the same view for gene expression. Surface geometry is defined
# by group_by; point colors come from RNA/data, not smoothed expression.
p_gene <- sn_plot_feature(
  obj, features = "NKG7", assay = "RNA", layer = "data",
  reduction = "umap3d", dims = 1:3, group_by = "cell_type",
  style = "nebula", camera = camera, raster_dpi = 600,
  palette = c("#27334D", "#35CFC4", "#FFF4A3")
)
ggplot2::ggsave("NKG7-nebula.pdf", p_gene, width = 6, height = 5, dpi = 600)
```

The browser is a self-contained WebGL viewer with locally bundled shaders
and no CDN dependency. Smooth normals, per-pixel rim lighting, soft particles
and bloom distinguish glass from nebula. The original coordinates and number
of cells are preserved; no synthetic particles are inserted. Static export
requires Chrome/Chromium (or `CHROMOTE_CHROME` pointing to its executable),
`chromote` and `png`. It uses local software WebGL rendering, so no dedicated
GPU is required. A standalone browser cannot change a running R variable: the JSON
or copied R list is the explicit return path. `sn_get_plot_camera(viewer)` returns
its **initial** camera, not subsequent browser interactions. Downloading/copying
pauses automatic rotation. Preserve the reduction, cell subset and grouping
when transferring a view; camera pan/zoom use scene-normalized coordinates.

Static output is a regular ggplot. Its **entire point/surface layer is rasterized
at 600 dpi**, including when saved to PDF. Labels and legends remain vector.
`ggsave(dpi = 600)` alone does not configure a ggplot raster layer; use
`raster_dpi = 600` as above (also the new styles' default). The actual raster size
follows the physical panel size at draw time. Browser and PDF use the same geometry and WebGL shader code; export renders
at the requested pixel size instead of upsampling a browser screenshot.
GPU/CPU antialiasing and vector text layout may differ. Saving a widget with `ggsave()` is unsupported; regenerate the static
plot with `interactive = FALSE` (the default).

`style_control` accepts `surface_alpha`, `point_alpha`, `glow` (all 0--1),
`surface_mass` (0.5--0.99, default 0.95), `bandwidth` (0.2--5, default 0.75),
`grid_size` (integer 16--64, default 48), `background`, and `auto_rotate`.
The surface is a Gaussian-smoothed, binned three-dimensional density envelope
of each group's coordinates. Its smoothing scale is estimated from local
neighbor distances, so distant islands do not inflate a global covariance
ellipsoid. At most 128 evenly spaced cells probe distances to all group cells
for this bandwidth estimate; all selected cells enter the KDE grid. It is a visual aid, **not a measured tissue
boundary or feature-expression isosurface**. Disconnected components are
retained when resolved by the grid; smoothing can still merge nearby islands.
Groups with fewer than five cells or rank-deficient coordinates show points
without a surface. All selected cells enter the density and point layer; no
hidden subsampling is applied. Browser responsiveness depends on cell count
and grid size; use an explicit `cells` subset when necessary.

Static plots support multiple features and `split_by`; browser viewing currently
supports one panel. Three-dimensional styles use a square borderless panel,
white labels with dark backplates and colored anchors, and numeric feature legends.
The default 3D categorical palette uses muted luminous colors; an explicit
`palette` or `cols` still takes precedence. Default point size is 0.20 mm. They reject shape/highlight/repel,
expression blending, `mode = "density"`, and extra Seurat arguments in `...`.
These remain available through `style = "classic"`. `group_by` in feature plots
is only for the new 3D styles; `keep_scale` controls shared expression limits.
