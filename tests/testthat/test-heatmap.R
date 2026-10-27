make_heatmap_test_object <- function() {
  counts <- matrix(
    c(
      10, 9, 1, 2, 6, 5,
      1, 2, 9, 10, 4, 5,
      3, 4, 8, 7, 9, 8
    ),
    nrow = 3,
    byrow = TRUE,
    dimnames = list(c("GENE1", "GENE2", "GENE3"), paste0("cell", seq_len(6)))
  )
  object <- Seurat::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  object$group <- factor(rep(c("control", "treated"), each = 3), levels = c("control", "treated"))
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  Seurat::ScaleData(object, features = rownames(object), verbose = FALSE)
}

test_that("average heatmaps retain every group and honor the selected layer", {
  skip_if_not_installed("Seurat")
  object <- make_heatmap_test_object()
  plot <- sn_plot_heatmap(
    object,
    features = rownames(object),
    mode = "average",
    group_by = "group"
  )

  expect_s3_class(plot, "ggplot")
  expect_identical(levels(plot$data$group), c("control", "treated"))
  expect_equal(nrow(plot$data), 6L)
  expect_true(all(table(plot$data$group) == 3L))

  scaled <- SeuratObject::LayerData(object, assay = "RNA", layer = "scale.data")
  expected <- rowMeans(scaled["GENE1", object$group == "control", drop = FALSE])
  observed <- plot$data$expression[
    plot$data$feature == "GENE1" & plot$data$group == "control"
  ]
  expect_equal(observed, unname(expected))
  expect_gt(length(unique(plot$data$value)), 1L)
  expect_no_error(ggplot2::ggplotGrob(plot))
})

test_that("average heatmaps support native pheatmap and ComplexHeatmap backends", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("pheatmap")
  skip_if_not_installed("ComplexHeatmap")
  object <- make_heatmap_test_object()

  pheatmap_plot <- sn_plot_heatmap(
    object,
    features = rownames(object),
    mode = "average",
    group_by = "group",
    backend = "pheatmap",
    cluster_rows = TRUE,
    cluster_cols = TRUE
  )
  complex_plot <- sn_plot_heatmap(
    object,
    features = rownames(object),
    mode = "average",
    group_by = "group",
    backend = "ComplexHeatmap",
    cluster_rows = TRUE,
    cluster_cols = TRUE
  )

  expect_s3_class(pheatmap_plot, "pheatmap")
  expect_s3_class(pheatmap_plot$tree_row, "hclust")
  expect_s3_class(pheatmap_plot$tree_col, "hclust")
  expect_s4_class(complex_plot, "Heatmap")
  expect_no_error(grid::grid.draw(pheatmap_plot$gtable))
  expect_no_error(ComplexHeatmap::draw(complex_plot))
})

test_that("ggplot average heatmaps expose row and column clustering", {
  skip_if_not_installed("Seurat")
  object <- make_heatmap_test_object()
  plot <- sn_plot_heatmap(
    object,
    features = rownames(object),
    mode = "average",
    group_by = "group",
    cluster_rows = TRUE,
    cluster_cols = TRUE
  )

  expect_setequal(levels(plot$data$group), c("control", "treated"))
  expect_setequal(levels(plot$data$feature), rownames(object))
  expect_no_error(ggplot2::ggplot_build(plot))
})
