test_that("feature drawing order is preserved by default and sorting is opt-in", {
  skip_if_not_installed("Seurat")
  object <- SeuratObject::pbmc_small
  object$score <- rev(seq_len(ncol(object)))
  for (raster in c(FALSE, TRUE)) {
    args <- list(object = object, features = "score", reduction = "tsne",
                 raster = raster, palette = NULL)
    for (setting in list(NULL, FALSE, TRUE)) {
      call_args <- args
      if (!is.null(setting)) call_args$order <- setting
      plot <- suppressWarnings(do.call(sn_plot_feature, call_args))
      leaf <- Shennong:::.sn_patchwork_leaf_plots(plot)[[1]]
      expected <- if (isTRUE(setting)) sort(object$score) else object$score
      expect_equal(leaf$data$score, unname(expected))
      points <- ggplot2::ggplot_build(leaf)$data[[1]]
      expect_equal(nrow(points), ncol(object))
      expect_equal(points$x, leaf$data$tSNE_1)
    }
  }
})
