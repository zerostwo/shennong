make_plot_scale_fixture <- function() {
  set.seed(42)
  counts <- matrix(rpois(8 * 36, 4), nrow = 8,
    dimnames = list(paste0("g", 1:8), paste0("c", 1:36)))
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  object$group <- factor(rep(c("A", "B"), each = 18))
  object$condition <- rep(c("x", "y", "z"), 12)
  SeuratObject::Idents(object) <- object$group
  Seurat::NormalizeData(object, verbose = FALSE)
}

test_that("unscaled dot colors preserve average expression including a single identity", {
  skip_if_not_installed("Seurat")
  object <- make_plot_scale_fixture()
  for (idents in list(NULL, "A")) {
    for (style in c("classic", "atlas", "atlas_outline", "atlas_module")) {
      plot <- suppressWarnings(sn_plot_dot(object, features = c("g1", "g2"),
        group_by = "group", idents = idents, scale = FALSE, style = style))
      reference <- suppressWarnings(Seurat::DotPlot(object, features = c("g1", "g2"),
        group.by = "group", idents = idents, scale = FALSE))
      expect_equal(plot$data, reference$data)
      built <- ggplot2::ggplot_build(plot)
      colors <- if (style == "atlas_outline") built$data[[1]]$fill else built$data[[1]]$colour
      aesthetic <- if (style == "atlas_outline") "fill" else "colour"
      expect_false(any(colors == "grey50"))
      expect_gt(length(unique(colors)), 1L)
      expect_equal(plot$labels[[aesthetic]], "Average expression")
      expect_true("Average expression" %in% vapply(built$plot$guides$params, `[[`, character(1), "title"))
      expect_equal(plot$scales$get_scales(aesthetic)$limits, range(plot$data$avg.exp.scaled))
      expect_s3_class(ggplot2::ggplotGrob(plot), "gtable")
    }
  }
  single_object <- object[, object$group == "A"]
  single_object$group <- as.character(single_object$group)
  single <- suppressWarnings(sn_plot_dot(single_object,
    features = c("g1", "g2"), group_by = "group"))
  expect_equal(single$labels$colour, "Average expression")
  expect_false(any(ggplot2::ggplot_build(single)$data[[1]]$colour == "grey50"))
})

test_that("dot colors retain explicit limits without censoring valid expression", {
  skip_if_not_installed("Seurat")
  object <- make_plot_scale_fixture()
  plot <- sn_plot_dot(object, features = c("g1", "g2"), group_by = "group",
    scale = FALSE, col_min = 0, col_max = 1)
  expect_gt(min(plot$data$avg.exp.scaled), 1)
  expect_equal(plot$scales$get_scales("colour")$limits, c(0, 1))
  expect_false(any(ggplot2::ggplot_build(plot)$data[[1]]$colour == "grey50"))
  capped <- sn_plot_dot(object, features = c("g1", "g2"), group_by = "group",
    scale = FALSE, col_max = 1)
  expect_equal(capped$scales$get_scales("colour")$limits, c(1, 1))
  expect_false(any(ggplot2::ggplot_build(capped)$data[[1]]$colour == "grey50"))
  scaled <- suppressWarnings(sn_plot_dot(object, features = c("g1", "g2"),
    group_by = "group"))
  expect_equal(scaled$labels$colour, "Z score")
  expect_true("Z score" %in% vapply(ggplot2::ggplot_build(scaled)$plot$guides$params, `[[`, character(1), "title"))
  expect_equal(scaled$scales$get_scales("colour")$limits, c(-2.5, 2.5))
})

test_that("split dot plots preserve upstream identity colors for two or more conditions", {
  skip_if_not_installed("Seurat")
  object <- make_plot_scale_fixture()
  palette <- c(z = "#008000", x = "#FF0000", y = "#0000FF")
  for (conditions in list(c("x", "y"), c("x", "y", "z"))) {
    subset <- object[, object$condition %in% conditions]
    for (scaled in c(TRUE, FALSE)) {
      plot <- suppressWarnings(sn_plot_dot(subset, features = c("g1", "g2"),
        group_by = "group", split_by = "condition", scale = scaled, palette = palette))
      reference <- suppressWarnings(Seurat::DotPlot(subset, features = c("g1", "g2"),
        group.by = "group", split.by = "condition", scale = scaled,
        cols = unname(palette[unique(subset$condition)])))
      expect_equal(plot$data, reference$data)
      expect_equal(ggplot2::ggplot_build(plot)$data[[1]]$colour, plot$data$colors)
      expect_s3_class(ggplot2::ggplotGrob(plot), "gtable")
    }
  }
  automatic <- sn_plot_dot(object, features = c("g1", "g2"),
    group_by = "group", split_by = "condition")
  expect_equal(ggplot2::ggplot_build(automatic)$data[[1]]$colour, automatic$data$colors)
  expect_s3_class(ggplot2::ggplotGrob(automatic), "gtable")
  expect_error(sn_plot_dot(object, features = "g1", split_by = "condition",
    palette = c(x = "red", y = "blue")), "palette.*z")
})

test_that("violin palettes cover the actual fill levels in every feature panel", {
  skip_if_not_installed("Seurat")
  object <- make_plot_scale_fixture()
  for (split in list(NULL, "condition")) {
    palette <- if (is.null(split)) c(A = "#FF0000", B = "#0000FF") else
      c(z = "#008000", x = "#FF0000", y = "#0000FF")
    plot <- suppressMessages(sn_plot_violin(object, features = c("g1", "g2"),
      group_by = "group", split_by = split, palette = palette))
    for (i in seq_len(2L)) {
      expect_setequal(unique(ggplot2::ggplot_build(plot[[i]])$data[[1]]$fill), unname(palette))
    }
    expect_s3_class(patchwork::patchworkGrob(plot), "gtable")
  }
  automatic <- suppressMessages(sn_plot_violin(object, features = c("g1", "g2"),
    group_by = "group", split_by = "condition"))
  for (i in seq_len(2L)) {
    expect_length(unique(ggplot2::ggplot_build(automatic[[i]])$data[[1]]$fill), 3L)
  }
  expect_s3_class(patchwork::patchworkGrob(automatic), "gtable")
})

test_that("violin feature titles survive multi-panel assembly and an overall title", {
  skip_if_not_installed("Seurat")
  object <- make_plot_scale_fixture()
  for (title in list(NULL, "Marker expression")) {
    plot <- suppressMessages(sn_plot_violin(object, features = c("g1", "g2"),
      group_by = "group", title = title))
    expect_equal(vapply(seq_len(2L), function(i) plot[[i]]$labels$title, character(1)), c("g1", "g2"))
    if (!is.null(title)) expect_equal(plot$patches$annotation$title, title)
    expect_s3_class(patchwork::patchworkGrob(plot), "gtable")
  }
})

test_that("split dot labels with regex characters retain colors and original labels", {
  skip_if_not_installed("Seurat")
  object <- make_plot_scale_fixture()
  palette <- c("x+y" = "#FF0000", "z[1]" = "#0000FF", "snSplit2" = "#008000")
  object$condition <- rep(names(palette), 12)
  metadata <- object[[]]
  plot <- sn_plot_dot(object, features = c("g1", "g2"),
    group_by = "group", split_by = "condition", palette = palette)
  reference_object <- object
  reference_object$condition <- rep(c("x", "y", "z"), 12)
  reference <- Seurat::DotPlot(reference_object, features = c("g1", "g2"),
    group.by = "group", split.by = "condition", cols = unname(palette))
  expect_equal(plot$data$avg.exp, reference$data$avg.exp)
  expect_equal(plot$data$avg.exp.scaled, reference$data$avg.exp.scaled)
  expect_identical(plot$data$colors, reference$data$colors)
  expect_identical(ggplot2::ggplot_build(plot)$data[[1]]$colour, reference$data$colors)
  expected_ids <- as.vector(t(outer(c("A", "B"), names(palette), paste, sep = "_")))
  expect_setequal(as.character(plot$data$id), expected_ids)
  expect_identical(object[[]], metadata)
  expect_s3_class(ggplot2::ggplotGrob(plot), "gtable")
})

test_that("split dot plots reject group and split labels that would merge identities", {
  skip_if_not_installed("Seurat")
  object <- make_plot_scale_fixture()
  object$group <- rep(c("A_x", "A"), each = 18)
  object$condition <- rep(c("y", "x_y"), each = 18)
  expect_error(sn_plot_dot(object, features = c("g1", "g2"),
    group_by = "group", split_by = "condition"), "ambiguous.*identities")
})
