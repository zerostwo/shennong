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
   Group-averaged `sn_plot_heatmap()` can return ggplot2, pheatmap, or
   ComplexHeatmap output. In its compact group-by-feature layout,
   `cluster_rows` clusters groups and `cluster_cols` clusters features.
   `sn_plot_composition()` accepts a Seurat object directly. Use ordinary
   `type = "bar"` for descriptive cell counts/proportions; use
   `type = "sample_bar"` or `"sample_boxplot"` with an explicit `sample_by`
   when biological samples are the replicates; use `type = "alluvial"` with
   two or more `flow_by` metadata columns for source/annotation flows. Use
   `unit_by` for donor metadata repeated across cells so donor distributions do
   not count cells as independent donors.

For dense expression maps, use `sn_plot_feature(object, features = "score",
reduction = "umap", pt_size = 0.05)`. High values are drawn last; NULL
`stroke_size` means zero stroke, preserving small points with or without
rasterization. Set a nonzero stroke explicitly only when wanted.

## 9

Build prompts or stored-result summaries with the interpretation helpers
   when a narrative or report-ready output is needed.

## Two- and three-dimensional nebula and glass maps

`sn_plot_dim()` and `sn_plot_feature()` accept `style = "nebula"` (bright rims)
or `style = "glass"` (translucent envelopes). `style = "classic"` preserves the
existing plots. Omitting `dims` uses **UMAP 1 and UMAP 2**, with planar density
contours and a front-facing camera. Points and contours stay in the same plane;
there is no inferred depth. For example:

```r
p <- sn_plot_dim(obj, reduction = "umap", group_by = "cell_type",
                 style = "nebula", raster_dpi = 600)
ggplot2::ggsave("umap-2d.pdf", p, width = 6, height = 6, dpi = 600)
```

An explicit `dims = 1:3` uses three real embedding dimensions. Plotting never
adds a random Z coordinate or recomputes UMAP. To compute a separate 3D map:

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
The surface is a Gaussian-smoothed, binned density contour (2D) or density envelope (3D)
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
The default styled categorical palette uses muted luminous colors; an explicit
`palette` or `cols` still takes precedence. Default point size is 0.20 mm. They reject shape/highlight/repel,
expression blending, `mode = "density"`, and extra Seurat arguments in `...`.
These remain available through `style = "classic"`. `group_by` in feature plots
is only for the new styles; `keep_scale` controls shared expression limits.


## Spatial coordinate maps

Use `sn_plot_spatial(obj, group_by = "region", split_by = "section",
pt_size = 0.8, label = TRUE)` for categorical maps. Named `cols` keep labels
matched to colors; palette, legend, axis, border, alpha, and title controls
follow dimension plots. Keep `point_size` compatibility but never pass both
size arguments. Metadata coordinates are preferred; supply `images` to select
native Seurat image/FOV centroids. Native images always occupy separate panels.
These maps do not draw histology or segmentation outlines. Preserve coordinate
units and use `flip_y = FALSE` for Cartesian coordinates.

`sn_plot_spatial_feature(obj, c("CD3D", "MS4A1"), assay = "Spatial",
layer = "data", split_by = "section")` uses the exact layer, with no counts
fallback. Numeric metadata is accepted and has name precedence over genes.
Panels share a common color scale; `limits` clamps colors without changing
source values. `raster = TRUE` requires ggrastr and leaves labels vector.
For domains, discover `sn_list_results(obj, type = "spatial_domains")`, retrieve
`sn_get_result(obj, "spatial_domains", result_id)`, and use
`sn_plot_spatial_domain(obj, result_id = result_id, label = TRUE)`; the retrieved
result is also accepted. Stored section boundaries are retained.


## Atlas figure styles

- Discover presets with `sn_list_plot_styles()`. Use `style = "atlas"`,
  `"atlas_outline"`, or `"atlas_module"` on supported ggplot-based embedding,
  expression, composition, association, DE/enrichment, survival, communication,
  and spatial views. Keep `classic` for previous defaults.
- Read the `atlas-figure-styles` article and installed
  `extdata/zhang-atlas-{papers,figures}.csv` for the six-paper figure inventory.
  T-cell Science entries are caption-only; supplementary figures are not an
  exhaustive inventory. Presets are visual adaptations, not pixel replicas.
- Use `sn_plot_matrix(x, type = "bubble", size = size_matrix)` when color and
  area describe different quantities. Keep aligned dimnames, missing values,
  units, statistic names and explicit common limits. `sn_plot_dot()` retains
  percent-expression area; it must not stand in for a significance bubble.
- `sn_plot_network(edges, nodes, layout = "supplied")` renders provided scores
  and fixed coordinates. It does not infer CoVarNet CMs, correlation, or LR
  communication. Preserve isolates and original weights.
- Atlas heatmaps require `backend = "ggplot2"`; split-dot color encodings are
  not supported by atlas styles. Source matrices and node/edge tables can be
  exported with `sn_export_figure_bundle()` for readback.

### Structural atlas layouts

`sn_list_plot_styles(function_name = "sn_plot_matrix")` lists structural styles
and the exact paper panels/input requirements. Themes (`atlas`,
`atlas_outline`, `atlas_module`) are separate from the 23 geometry choices.
Use matrix `split_triangle`, `circular`, `annotated`, `dendrogram`,
`grouped_bubble`, `heatmap_bar`, `heatmap_density`; distribution `violin_box`,
`raincloud`, `paired_box`, `paired_slope`, `split_violin`, `stacked_violin`;
`sn_plot_effect()` for supplied coefficient/ratio forests, signed bars and
lollipops; `sn_plot_overlap()` for exclusive UpSet intersections; and
`sn_plot_network()` for bipartite, supplied clone-tree, weighted chord and
shared-coordinate rewiring structures. Existing `sn_plot_heatmap(mode="average")`
and `sn_plot_violin()` expose the matching matrix/distribution layouts.

Require explicit complete subject pairs (`pair_by`, optional specimen
`sample_by`), an explicit coefficient versus ratio scale, original interval
bounds, and supplied network/tree topology. No automatic P-values, BCR trees,
CoVarNet discovery or imputation. Reusable source tables survive PDF/PNG bundle
export. Consult `extdata/zhang-atlas-layouts.csv` and the rendered atlas article;
synthetic layout examples must not be described as original cancer results.
