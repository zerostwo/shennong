#!/usr/bin/env Rscript
# Local, offline visual check using the real Visium fixture; never edits it.
args <- commandArgs(trailingOnly = TRUE)
fixture <- if (length(args)) args[[1]] else "data-local/pkgdown-real/spatial/visium_lymph_node.qs2"
destination <- if (length(args) > 1L) args[[2]] else "data-local/spatial-plot-review"
if (!file.exists(fixture)) stop("Supply the existing Visium fixture path.")
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
pkgload::load_all(quiet = TRUE)
artifact <- qs2::qs_read(fixture)
stopifnot(is.list(artifact$sections), length(artifact$sections) > 0L)
plots <- list()
for (section in names(artifact$sections)) {
  object <- artifact$sections[[section]]
  object <- Seurat::NormalizeData(object, assay = "Spatial", verbose = FALSE)
  # Descriptive QC strata, not biological labels or inferred tissue domains.
  depth <- object$nCount_Spatial
  object$library_size_quartile <- cut(depth,
    breaks = unique(stats::quantile(depth, seq(0, 1, 0.25))), include.lowest = TRUE)
  categorical <- sn_plot_spatial(object, "library_size_quartile", palette = "Paired",
    pt_size = 1.4, title = paste(section, "- library size (QC)"))
  feature <- sn_plot_spatial_feature(object, c("CD3D", "MS4A1"), assay = "Spatial",
    layer = "data", pt_size = 1.4, title = paste(section, "- normalized expression"))
  for (name in c("categorical", "feature")) {
    plot <- get(name)
    stopifnot(nrow(plot$data) == ncol(object) * if (name == "feature") 2L else 1L)
    native <- SeuratObject::GetTissueCoordinates(object[[SeuratObject::Images(object)[[1]]]], scale = NULL)
    index <- match(plot$data$cell, native$cell)
    stopifnot(identical(as.numeric(plot$data$spatial_x), native$x[index]),
              identical(as.numeric(plot$data$spatial_y), native$y[index]))
    stem <- paste(section, name, sep = "-")
    ggplot2::ggsave(file.path(destination, paste0(stem, ".png")), plot,
      width = if (name == "feature") 7 else 4.5, height = 4, dpi = 180)
    ggplot2::ggsave(file.path(destination, paste0(stem, ".pdf")), plot,
      width = if (name == "feature") 7 else 4.5, height = 4)
    utils::write.csv(plot$data, file.path(destination, paste0(stem, ".csv")), row.names = FALSE)
    plots[[stem]] <- plot
  }
}
saveRDS(plots, file.path(destination, "plots.rds"))
cat("Validated native coordinate identity for", length(plots), "real Visium maps.\n")
cat(normalizePath(destination), "\n")
