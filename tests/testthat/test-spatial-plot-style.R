.spatial_style_object <- function() {
  counts <- matrix(seq_len(24), 3, dimnames = list(c("G1", "G2", "G3"), paste0("c", 1:8)))
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  object$x <- rep(1:4, 2)
  object$y <- rep(c(1, 2), each = 4)
  object$section <- factor(rep(c("B", "A"), each = 4), levels = c("B", "A"))
  object$region <- factor(rep(c("right", "left"), 4), levels = c("right", "left"))
  object$score <- c(NA, 2:8)
  Seurat::NormalizeData(object, verbose = FALSE)
}

test_that("spatial facets preserve cells, factor order, colors, and geometry", {
  object <- .spatial_style_object()
  original <- object
  colors <- c(left = "#0000FF", right = "#FF0000")
  plot <- sn_plot_spatial(object, "region", split_by = "section", cols = colors,
                          pt_size = 0.4, label = TRUE, show_border = FALSE)
  built <- ggplot2::ggplot_build(plot)
  expect_equal(nrow(plot$data), 8L)
  expect_identical(levels(plot$data$value), c("right", "left"))
  expect_equal(as.character(built$layout$layout$spatial_sample), c("B", "A"))
  expect_equal(built$data[[1]]$colour, unname(colors[as.character(plot$data$value)]))
  expect_true(all(built$data[[1]]$size == 0.4))
  expect_true(all(built$data[[1]]$stroke == 0))
  expect_equal(plot$coordinates$ratio, 1)
  expect_equal(built$data[[1]]$y, -plot$data$spatial_y)
  expect_equal(nrow(plot$layers[[2]]$data), 4L)
  expect_equal(plot$theme$text$size, 8)
  expect_silent(ggplot2::ggplotGrob(plot))
  expect_identical(object, original)
  expect_equal(sn_get_figure_spec(plot)$data_summary$n_panels, 2)
  selected <- sn_plot_spatial(object, "region", cells = c("c4", "c2"), cols = colors)
  expect_equal(selected$data$cell, c("c4", "c2"))
  expect_error(sn_plot_spatial(object, cells = "unknown"), "cells")
  expect_error(sn_plot_spatial(object, "region", cols = c(left = "red")), "cols")
  expect_error(sn_plot_spatial(object, point_size = 2, pt_size = 1), "only one")
  expect_error(sn_plot_spatial(object, split_by = "absent"), "split_by")
  plot <- sn_plot_spatial(object, "region", cols = colors, label = TRUE, label_color = "group")
  labels <- ggplot2::ggplot_build(plot)$data[[2]]
  expect_equal(labels$colour, unname(colors[as.character(labels$label)]))
  expect_equal(sn_get_figure_spec(plot)$data_summary$aspect_ratio, 3)
})

test_that("spatial expression uses exact layers and preserves raw values and feature order", {
  object <- .spatial_style_object()
  plot <- sn_plot_spatial_feature(object, c("G2", "G1"), layer = "counts",
                                  split_by = "section", limits = c(0, 10))
  expect_identical(levels(plot$data$feature), c("G2", "G1"))
  expect_equal(nrow(plot$data), 16L)
  raw <- SeuratObject::LayerData(object, layer = "counts")
  expect_equal(plot$data$expression, unname(raw[cbind(match(plot$data$feature, rownames(raw)),
                                                     match(plot$data$cell, colnames(raw)))]))
  built <- ggplot2::ggplot_build(plot)
  expect_equal(nrow(built$layout$layout), 4L)
  expect_equal(plot$scales$get_scales("colour")$limits, c(0, 10))
  expect_silent(ggplot2::ggplotGrob(plot))
  expect_error(sn_plot_spatial_feature(object, "G1", layer = "missing"), "layer|Layer")
  expect_warning(partial <- sn_plot_spatial_feature(object, c("missing", "G1")), "missing")
  expect_equal(unique(as.character(partial$data$feature)), "G1")
  score <- sn_plot_spatial_feature(object, "score", na_value = "pink", order = TRUE)
  expect_true(is.na(score$data$expression[[1]]))
  expect_equal(score$data$expression[-1], 2:8)
  expect_equal(ggplot2::ggplot_build(score)$data[[1]]$colour[[1]], "pink")
  expect_error(sn_plot_spatial_feature(object, "region"), "numeric")
})

test_that("native image coordinates work without metadata and sections never overlay", {
  object <- .spatial_style_object()
  object$x <- NULL
  object$y <- NULL
  for (section in c("one", "two")) {
    cells <- if (section == "one") colnames(object)[1:4] else colnames(object)[5:8]
    coords <- data.frame(x = 1:4, y = 5:8, row.names = cells)
    object[[section]] <- SeuratObject::CreateFOV(coords, type = "centroids", assay = "RNA", key = paste0(section, "_"))
  }
  plot <- sn_plot_spatial(object, "region")
  expect_equal(nrow(plot$data), 8L)
  expect_equal(plot$data$spatial_x, rep(1:4, 2))
  expect_equal(plot$data$spatial_y, rep(5:8, 2))
  expect_equal(nrow(ggplot2::ggplot_build(plot)$layout$layout), 2L)
  expect_silent(ggplot2::ggplotGrob(plot))
  expect_equal(nrow(sn_plot_spatial(object, images = "two")$data), 4L)
  expect_error(sn_plot_spatial(object, images = "missing"), "images")
  expect_error(sn_plot_spatial(object, images = "one", cells = "c8"), "coordinates")
})

test_that("domain and deconvolution maps share styling and section-aware facets", {
  object <- .spatial_style_object()
  object <- sn_find_spatial_domains(object, method = "cellcharter", sample_by = "section",
    result_id = "domains", backend_control = list(result = list(domains = data.frame(
      cell = colnames(object), domain = object$region))))
  plot <- sn_plot_spatial_domain(object = object, result_id = "domains", label = TRUE,
                                  palette = "OkabeIto", pt_size = 0.5)
  expect_equal(nrow(ggplot2::ggplot_build(plot)$layout$layout), 2L)
  expect_silent(ggplot2::ggplotGrob(plot))
  table <- data.frame(cell = rep(colnames(object), 2), cell_type = rep(c("T", "B"), each = 8),
                      spatial_x = rep(object$x, 2), spatial_y = rep(object$y, 2),
                      section = rep(object$section, 2), proportion = seq_len(16) / 16)
  plot <- sn_plot_spatial_deconvolution(table, split_by = "section", pt_size = 0.3)
  expect_equal(nrow(ggplot2::ggplot_build(plot)$layout$layout), 4L)
  expect_silent(ggplot2::ggplotGrob(plot))
  expect_equal(sort(plot$data$expression), table$proportion)
})

test_that("VisiumV1 uses full-resolution column/row coordinates in barcode order", {
  object <- .spatial_style_object()
  coordinates <- data.frame(imagerow = (1:8) * 100, imagecol = (8:1) * 200,
                             row.names = rev(colnames(object)))
  expect_warning(object[["slice"]] <- methods::new("VisiumV1", image = array(1, c(10, 10, 3)),
    scale.factors = Seurat::scalefactors(lowres = 0.1, hires = 0.5),
    coordinates = coordinates, spot.radius = 0.01, assay = "RNA", key = "slice_"), "unordered cells")
  plot <- sn_plot_spatial(object, "region", images = "slice", cells = c("c1", "c7"))
  expect_equal(plot$data$spatial_x, coordinates[c("c1", "c7"), "imagecol"])
  expect_equal(plot$data$spatial_y, coordinates[c("c1", "c7"), "imagerow"])
  expect_equal(as.character(plot$data$value), as.character(object$region[c("c1", "c7")]))
  expect_error(sn_plot_spatial(object, images = "slice", spatial_cols = c("x", "y")), "only one")
  normal <- sn_plot_spatial(object, flip_y = FALSE)
  expect_equal(ggplot2::ggplot_build(normal)$data[[1]]$y, unname(object$y))
})

test_that("missing labels, invalid coordinates and abundance values remain explicit", {
  object <- .spatial_style_object()
  object$region <- rep(NA_character_, ncol(object))
  plot <- sn_plot_spatial(object, "region", label = TRUE, na_value = "pink")
  expect_true(all(ggplot2::ggplot_build(plot)$data[[1]]$colour == "pink"))
  expect_silent(ggplot2::ggplotGrob(plot))
  object$x[[1L]] <- Inf
  expect_error(sn_plot_spatial(object), "finite")
  table <- data.frame(cell = "c1", cell_type = "T", spatial_x = 1, spatial_y = 2, proportion = 12)
  expect_equal(sn_plot_spatial_deconvolution(table)$data$expression, 12)
  expect_error(sn_plot_spatial_deconvolution(table, limits = c(1, 0)), "limits")
})

test_that("rasterized maps keep vector labels and explicit point sizing", {
  skip_if_not_installed("ggrastr")
  plot <- sn_plot_spatial(.spatial_style_object(), "region", raster = TRUE,
                          raster_dpi = 150, pt_size = 0.15, label = TRUE)
  expect_true(all(ggplot2::ggplot_build(plot)$data[[1]]$size == 0.15))
  file <- tempfile(fileext = ".pdf")
  on.exit(unlink(file))
  expect_silent(ggplot2::ggsave(file, plot, width = 4, height = 3))
  expect_gt(file.info(file)$size, 1000)
})
