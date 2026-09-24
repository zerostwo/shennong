.embedding_style_fixture <- function() {
  skip_if_not_installed("Seurat")
  suppressWarnings(skip_if_not_installed("misc3d"))
  skip_if_not_installed("htmlwidgets")
  skip_if_not_installed("chromote")
  skip_if_not_installed("png")
  set.seed(78)
  counts <- matrix(rpois(6 * 90, 3), 6, dimnames = list(paste0("gene", 1:6), paste0("cell", 1:90)))
  obj <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  xyz <- matrix(rnorm(270, sd = .35), ncol = 3)
  xyz[, 1] <- xyz[, 1] + rep(c(-1.5, 0, 1.5), each = 30)
  rownames(xyz) <- colnames(obj)
  colnames(xyz) <- paste0("UMAP3D_", 1:3)
  obj[["umap3d"]] <- SeuratObject::CreateDimReducObject(xyz, key = "UMAP3D_", assay = "RNA")
  obj$cell_type <- rep(c("A", "B", "C"), each = 30)
  obj$batch <- rep(c("one", "two"), 45)
  obj
}

test_that("public embedding styles preserve cells and rasterize the complete scene at 600 dpi", {
  obj <- .embedding_style_fixture()
  for (style in c("nebula", "glass")) {
    p <- sn_plot_dim(obj, reduction = "umap3d", dims = 1:3,
                     group_by = "cell_type", style = style, label = TRUE)
    expect_s3_class(p, "ggplot")
    scene <- attr(p, "shennong_embedding_scene")
    expect_identical(scene$cells, colnames(obj))
    expect_equal(scene$xyz, unname(SeuratObject::Embeddings(obj[["umap3d"]])))
    expect_gt(sum(vapply(scene$surfaces, function(s) nrow(s$vertices), integer(1))), 0)
    grobs <- p$layers[[1]]$draw_geom(ggplot2::ggplot_build(p)$data[[1]], ggplot2::ggplot_build(p)$layout)
    expect_equal(grobs[[1]]$dpi, 600)
    expect_s3_class(grobs[[1]], "sn_embedding_webgl")
    if (!nzchar(chromote::find_chrome() %||% "")) next
    file <- tempfile(fileext = ".pdf")
    ggplot2::ggsave(file, p, width = 4, height = 4, dpi = 600)
    expect_gt(file.info(file)$size, 1000)
    unlink(file)
  }
})

test_that("browser camera JSON round-trips through public static plotting", {
  obj <- .embedding_style_fixture()
  camera <- list(azimuth = 120, elevation = -25, roll = 10, zoom = 1.3, pan = c(.1, -.2))
  path <- tempfile(fileext = ".json")
  jsonlite::write_json(camera, path, auto_unbox = TRUE)
  on.exit(unlink(path))
  recovered <- sn_get_plot_camera(path)
  p <- sn_plot_dim(obj, reduction = "umap3d", dims = 1:3, style = "glass", camera = recovered)
  expect_equal(sn_get_plot_camera(p), camera)
  skip_if_not_installed("htmlwidgets")
  w <- sn_plot_dim(obj, reduction = "umap3d", dims = 1:3, style = "glass", camera = recovered, interactive = TRUE)
  expect_s3_class(w, "htmlwidget")
  expect_equal(w$x$camera, recovered)
  expect_equal(w$x$scene$xyz, attr(p, "shennong_embedding_scene")$xyz)
  expect_error(sn_get_plot_camera(list(zoom = 0)), "zoom")
  expect_error(sn_get_plot_camera(list(azmuth = 0)), "Unknown")
})

test_that("feature styles honor assay layer, selected cell order, cutoffs and shared camera", {
  obj <- .embedding_style_fixture()
  selected <- rev(colnames(obj)[1:60])
  p <- sn_plot_feature(obj, "gene2", reduction = "umap3d", dims = 1:3,
                       cells = selected, group_by = "cell_type", style = "nebula",
                       layer = "counts", min_cutoff = 2, max_cutoff = 4)
  scene <- attr(p, "shennong_embedding_scene")
  expect_identical(scene$cells, selected)
  expected <- as.numeric(SeuratObject::LayerData(obj, layer = "counts")["gene2", selected])
  expect_equal(scene$values, pmin(4, pmax(2, expected)))
  expect_equal(scene$limits, c(2, 4))
  expect_s3_class(sn_plot_feature(obj, c("gene1", "gene2"), reduction = "umap3d",
                                 dims = 1:3, style = "glass"), "patchwork")
  expect_s3_class(sn_plot_dim(obj, reduction = "umap3d", dims = 1:3,
                             style = "glass", split_by = "batch"), "patchwork")
})

test_that("styles reject missing 3D inputs and incompatible options explicitly", {
  obj <- .embedding_style_fixture()
  expect_error(sn_plot_dim(obj, reduction = "umap3d", dims = 1, style = "glass"), "two or three")
  expect_error(sn_plot_dim(obj, reduction = "umap3d", dims = 1:3, style = "glass", raster = FALSE), "raster")
  expect_error(sn_plot_feature(obj, "gene1", reduction = "umap3d", dims = 1:3,
                               style = "glass", mode = "density"), "mode")
  expect_error(sn_plot_dim(obj, reduction = "umap3d", dims = 1:3,
                           style = "glass", style_control = list(glwo = 1)), "Unknown")
  expect_error(sn_plot_dim(obj, reduction = "umap3d", dims = 1:3,
                           style = "glass", cells = c("missing")), "cells")
  expect_s3_class(sn_plot_dim(obj, reduction = "umap3d", style = "classic", raster = FALSE), "ggplot")
})

test_that("single-group widget vectors stay arrays and missing layers cannot silently fall back", {
  obj <- .embedding_style_fixture()
  skip_if_not_installed("htmlwidgets")
  w <- sn_plot_dim(obj, reduction = "umap3d", dims = 1:3, style = "glass", label = TRUE, interactive = TRUE)
  json <- htmlwidgets:::toJSON(w$x)
  restored <- jsonlite::fromJSON(json, simplifyVector = FALSE)
  expect_type(restored$scene$labels, "list")
  expect_length(restored$scene$labels, 1)
  expect_type(restored$scene$group_labels, "list")
  expect_error(sn_plot_feature(obj, "gene1", reduction = "umap3d", dims = 1:3,
                               style = "glass", layer = "absent"), "layer")
  expect_identical(SeuratObject::Reductions(obj), "umap3d")
})

test_that("geometry handles small groups and preserves disconnected islands", {
  skip_if_not_installed("misc3d")
  set.seed(701)
  ctl <- .sn_embedding_control(list(grid_size = 32L), "glass")
  expect_equal(nrow(.sn_embedding_surface(matrix(1:12, 4), ctl, 5)), 0)
  xyz <- rbind(matrix(rnorm(150, sd = .1), ncol = 3), matrix(rnorm(150, sd = .1), ncol = 3) + 5)
  vertices <- .sn_embedding_surface(xyz, ctl, 5)
  expect_gt(nrow(vertices), 0)
  # No fabricated bridge across the empty middle of the two clouds.
  expect_false(any(vertices[, 1] > 1.5 & vertices[, 1] < 3.5))
})

test_that("WebGL export renders real pixels at the requested size and does not upsample", {
  obj <- .embedding_style_fixture()
  skip_if(!nzchar(chromote::find_chrome() %||% ""), "Chrome is needed for WebGL capture.")
  p <- sn_plot_dim(obj, reduction = "umap3d", dims = 1:3, style = "nebula", group_by = "cell_type")
  scene <- attr(p, "shennong_embedding_scene")
  raster <- .sn_embedding_capture(scene, sn_get_plot_camera(p), 600L, 600L, 600)
  expect_equal(dim(raster), c(600, 600))
  expect_s3_class(raster, "nativeRaster")
  expect_gt(length(unique(as.integer(raster))), 100)
  expect_equal(length(scene$cells), ncol(obj))
  # The browser dependency must ship locally, rather than fetching a CDN.
  w <- .sn_embedding_widget(scene, sn_get_plot_camera(p))
  directory <- tempfile("widget-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  htmlwidgets::saveWidget(w, file.path(directory, "viewer.html"), selfcontained = FALSE)
  expect_true(length(list.files(directory, "shennong-webgl[.]js", recursive = TRUE)) == 1L)
})


test_that("omitted dimensions use the first two coordinates without synthetic depth", {
  obj <- .embedding_style_fixture()
  for (style in c("glass", "nebula")) {
    p <- sn_plot_dim(obj, reduction = "umap3d", style = style, group_by = "cell_type")
    scene <- attr(p, "shennong_embedding_scene")
    expect_equal(scene$dimensions, 2)
    expect_equal(scene$xyz[, 1:2], unname(SeuratObject::Embeddings(obj[["umap3d"]])[, 1:2]))
    expect_true(all(scene$xyz[, 3] == 0))
    expect_equal(sn_get_plot_camera(p)$elevation, -90)
    expect_true(all(vapply(scene$surfaces, function(s) all(s$vertices[, 3] == 0), logical(1))))
    expect_true(all(vapply(scene$surfaces, function(s) any(s$density > 1) && any(s$density < 1), logical(1))))
    w <- sn_plot_feature(obj, "gene2", reduction = "umap3d", style = style, interactive = TRUE)
    expect_equal(w$x$scene$dimensions, 2)
    expect_true(all(w$x$scene$xyz[, 3] == 0))
    if (nzchar(chromote::find_chrome() %||% "")) {
      raster <- .sn_embedding_capture(scene, sn_get_plot_camera(p), 600L, 600L, 600)
      expect_equal(dim(raster), c(600L, 600L))
      expect_gt(length(unique(as.vector(raster))), 100)
    }
  }
})
