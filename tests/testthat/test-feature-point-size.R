test_that("ordered feature points retain small sizes without a default stroke", {
  skip_if_not_installed("Seurat")
  object <- SeuratObject::pbmc_small
  object$score <- rev(seq_len(ncol(object)))
  plots <- lapply(c(FALSE, TRUE), function(raster) {
    suppressWarnings(sn_plot_feature(
      object, features = "score", reduction = "tsne",
      pt_size = 0.01, raster = raster, palette = NULL
    ))
  })
  for (plot in plots) {
    leaf <- Shennong:::.sn_patchwork_leaf_plots(plot)[[1]]
    points <- ggplot2::ggplot_build(leaf)$data[[1]]
    expect_equal(unique(points$size), 0.01)
    expect_equal(unique(points$stroke), 0)
    expect_equal(leaf$data$score, unname(sort(object$score)))
    expect_s3_class(ggplot2::ggplotGrob(leaf), "gtable")
  }
  explicit <- suppressWarnings(sn_plot_feature(
    object, features = "score", reduction = "tsne",
    pt_size = 0.01, stroke_size = 0.2, raster = FALSE, palette = NULL
  ))
  leaf <- Shennong:::.sn_patchwork_leaf_plots(explicit)[[1]]
  expect_equal(unique(ggplot2::ggplot_build(leaf)$data[[1]]$stroke), 0.2)
})
