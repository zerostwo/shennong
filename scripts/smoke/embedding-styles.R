#!/usr/bin/env Rscript
# Local, immutable-input PBMC demonstration. No downloads; ~2,000 cells.
devtools::load_all(quiet = TRUE)
input <- "data-local/pkgdown-real/single-cell/kotliarov_pbmc.qs2"
out <- "dev/umap/style-demo"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
obj <- qs2::qs_read(input)
obj <- Seurat::NormalizeData(obj, verbose = FALSE)
obj <- Seurat::FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
obj <- Seurat::ScaleData(obj, verbose = FALSE)
obj <- Seurat::RunPCA(obj, npcs = 30, verbose = FALSE, seed.use = 717)
obj <- Seurat::FindNeighbors(obj, dims = 1:20, verbose = FALSE)
obj <- Seurat::FindClusters(obj, resolution = .5, random.seed = 717, verbose = FALSE)
obj <- Seurat::RunUMAP(obj, dims = 1:20, n.components = 3L, reduction.name = "umap3d",
                       reduction.key = "UMAP3D_", seed.use = 717, verbose = FALSE)
qs2::qs_save(obj, file.path(out, "pbmc-3d.qs2"))
for (style in c("nebula", "glass")) {
  p <- sn_plot_dim(obj, reduction = "umap3d", dims = 1:3, group_by = "seurat_clusters",
                   style = style, label = TRUE, title = paste("PBMC clusters /", style), raster_dpi = 600)
  ggplot2::ggsave(file.path(out, paste0(style, ".pdf")), p, width = 6, height = 5, dpi = 600)
  ggplot2::ggsave(file.path(out, paste0(style, ".png")), p, width = 6, height = 5, dpi = 150)
  w <- sn_plot_dim(obj, reduction = "umap3d", dims = 1:3, group_by = "seurat_clusters",
                   style = style, label = TRUE, interactive = TRUE,
                   title = paste("PBMC clusters /", style))
  htmlwidgets::saveWidget(w, file.path(out, paste0(style, ".html")), selfcontained = TRUE)
}
p <- sn_plot_feature(obj, "NKG7", reduction = "umap3d", dims = 1:3, group_by = "seurat_clusters",
                     style = "glass", raster_dpi = 600, palette = c("#27334D", "#35CFC4", "#FFF4A3"))
ggplot2::ggsave(file.path(out, "NKG7.pdf"), p, width = 6, height = 5, dpi = 600)
ggplot2::ggsave(file.path(out, "NKG7.png"), p, width = 6, height = 5, dpi = 150)
w <- sn_plot_feature(obj, "NKG7", reduction = "umap3d", dims = 1:3, group_by = "seurat_clusters",
                     style = "glass", interactive = TRUE, palette = c("#27334D", "#35CFC4", "#FFF4A3"))
htmlwidgets::saveWidget(w, file.path(out, "NKG7.html"), selfcontained = TRUE)
jsonlite::write_json(list(input = input, input_sha256 = digest::digest(file = input, algo = "sha256"),
                         cells = ncol(obj), features = nrow(obj), seed = 717, raster_dpi = 600,
                         session = capture.output(sessionInfo())),
                    file.path(out, "manifest.json"), auto_unbox = TRUE, pretty = TRUE)
