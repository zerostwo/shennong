#!/usr/bin/env Rscript
# Source-data reproduction and real-PBMC style review; no inference or relabeling.
args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(getwd(), mustWork = TRUE)
stopifnot(file.exists(file.path(root, "DESCRIPTION")))
pkgload::load_all(root, quiet = TRUE)
source_dir <- file.path(root, "data-local", "zhang-atlas-figures")
out <- file.path(root, "data-local", "atlas-style-review")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
filename <- "41586_2025_9053_MOESM6_ESM.xlsx"
url <- paste0("https://static-content.springer.com/esm/",
  "art%3A10.1038%2Fs41586-025-09053-4/MediaObjects/", filename)
file <- file.path(source_dir, filename)
if (!file.exists(file) && "--download" %in% args) {
  utils::download.file(url, file, mode = "wb", quiet = TRUE)
}
if (!file.exists(file)) stop("Missing CM source workbook; rerun with --download.")
stopifnot(requireNamespace("readxl", quietly = TRUE))
d <- readxl::read_excel(file, sheet = "Fig. 2c", skip = 1)
roe <- vapply(d[-1], as.numeric, numeric(nrow(d)))
rownames(roe) <- d[[1]]
stopifnot(identical(dim(roe), c(12L, 34L)), all(is.finite(roe)))
p <- sn_plot_matrix(roe, style = "atlas_module", palette = c("#FFFFFF", "#FCAE91", "#CB181D"),
  limits = c(0, 3), value_label = "Ro/e", title = "CM tissue prevalence: original Fig. 2c numbers",
  panel_widths = 600, panel_heights = 180)
stopifnot(identical(attr(p, "shennong_figure_data")$matrix, roe))
sn_export_figure_bundle(p, file.path(out, "nature-cm-fig2c"),
  width = 280, height = 145, dpi = 180)
readback <- utils::read.csv(file.path(out, "nature-cm-fig2c", "nature-cm-fig2c_data.csv"),
  check.names = FALSE, row.names = 1)
stopifnot(isTRUE(all.equal(as.matrix(readback), roe, tolerance = 1e-14)))
jsonlite::write_json(list(doi = "10.1038/s41586-025-09053-4",
  source_url = url, source_sha256 = digest::digest(file = file, algo = "sha256"),
  sheet = "Fig. 2c", skipped_rows = 1L, dimensions = dim(roe),
  transformation = "none; numeric text parsed as double", display_limits = c(0, 3), order = "source workbook; no reclustering",
  readback_all_values_equal = TRUE,
  fidelity = "Original numerical values; approximate sequential red palette; display capped at 3; no original dendrogram"),
  file.path(out, "source-provenance.json"), pretty = TRUE, auto_unbox = TRUE)

pbmc <- SeuratObject::pbmc_small
pbmc$atlas_identity <- SeuratObject::Idents(pbmc)
markers <- list(B = c("MS4A1", "CD79A"), T = c("CD3D", "CD2"))
styles <- c("classic", "atlas", "atlas_outline", "atlas_module")
dot_plots <- lapply(styles, function(style) {
  sn_plot_dot(pbmc, features = markers, style = style,
    panel_widths = 95, panel_heights = 150) + ggplot2::ggtitle(style)
})
dot_grid <- patchwork::wrap_plots(dot_plots, ncol = 2) +
  patchwork::plot_annotation(title = "Same PBMC marker values, four styles",
    theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 8)))
sn_save_figure(dot_grid, file.path(out, "marker-style-comparison.pdf"), width = 230, height = 230)
sn_save_figure(dot_grid, file.path(out, "marker-style-comparison.png"), width = 230, height = 230, dpi = 160)
panels <- list(
  embedding = sn_plot_dim(pbmc, reduction = "tsne", style = "atlas",
    raster = FALSE, pt_size = 1.4, panel_widths = 150, panel_heights = 150,
    title = "PBMC: original t-SNE and identities"),
  feature = sn_plot_feature(pbmc, features = "MS4A1", reduction = "tsne",
    mode = "expression", style = "atlas", raster = FALSE, pt_size = 1.4,
    panel_widths = 150, panel_heights = 150),
  matrix = sn_plot_heatmap(pbmc, unlist(markers), layer = "data", mode = "average",
    style = "atlas_outline", panel_widths = 150, panel_heights = 100),
  distribution = sn_plot_distribution(pbmc, "nFeature_RNA", view = "box", style = "atlas_outline",
    group_by = "atlas_identity", show_points = TRUE, title = "Detected genes per cell",
    y_label = "Detected genes", x_label = "Original identity",
    panel_widths = 150, panel_heights = 150)
)
gallery <- patchwork::wrap_plots(panels, ncol = 2) +
  patchwork::plot_annotation(title = "Atlas visual grammar on real teaching data",
    subtitle = "SeuratObject::pbmc_small; not the cancer atlas datasets",
    theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 8),
      plot.subtitle = ggplot2::element_text(size = 8)))
sn_save_figure(gallery, file.path(out, "atlas-gallery.pdf"), width = 240, height = 205)
sn_save_figure(gallery, file.path(out, "atlas-gallery.png"), width = 240, height = 205, dpi = 160)
edges <- data.frame(source = c("A", "A", "B", "C"), target = c("B", "C", "D", "D"),
  weight = c(0.3, 0.8, 0.4, 0.6))
nodes <- data.frame(node = c("A", "B", "C", "D", "Isolate"),
  lineage = c("Group 1", "Group 2", "Group 1", "Group 3", "Group 2"))
network <- sn_plot_network(edges, nodes, group_by = "lineage", style = "atlas_module",
  weight_label = "Supplied score", title = "Artificial topology demo; no CM inference")
sn_export_figure_bundle(network, file.path(out, "network-interface-demo"), width = 160, height = 120)
writeLines(capture.output(sessionInfo()), file.path(out, "session.txt"))
cat("Verified 408 original CM values; review artifacts: ", out, "\n", sep = "")
