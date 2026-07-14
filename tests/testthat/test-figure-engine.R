test_that("figure profiles and large-data specs are deterministic", {
  profiles <- sn_list_figure_profiles()
  expect_true(all(c("screen", "single_column", "double_column", "slide_16_9") %in% profiles$profile))
  sizes <- c(500, 5000, 50000, 500000, 5000000)
  specs <- lapply(sizes, function(n) sn_figure_spec(
    plot_type = "embedding", data_summary = list(n_points = n, n_groups = 18)
  ))
  points <- vapply(specs, function(x) x$recommended$point_size, numeric(1))
  raster <- vapply(specs, function(x) x$recommended$rasterize, logical(1))
  expect_true(all(diff(points) < 0))
  expect_identical(raster, c(FALSE, FALSE, FALSE, TRUE, TRUE))
  expect_gte(min(points), 0.03)
})

test_that("auto sizing responds monotonically to panels, genes, and labels", {
  one <- sn_recommend_figure_size(plot_type = "embedding", data_summary = list(n_points = 5000, n_panels = 1))
  four <- sn_recommend_figure_size(plot_type = "embedding", data_summary = list(n_points = 5000, n_panels = 4))
  dot10 <- sn_figure_spec(plot_type = "dot", data_summary = list(n_features = 10, n_groups = 8, max_label_chars = 8))
  dot40 <- sn_figure_spec(plot_type = "dot", data_summary = list(n_features = 40, n_groups = 8, max_label_chars = 32))
  expect_gte(four$height_mm, one$height_mm)
  expect_gt(dot40$recommended$height_mm, dot10$recommended$height_mm)
  expect_gte(dot40$data_summary$max_label_chars, dot10$data_summary$max_label_chars)
  expect_identical(dot40$recommended$legend_position, "bottom")
})

test_that("profiles cap dimensions and explicit overrides win", {
  spec <- sn_figure_spec(
    plot_type = "dot", profile = "single_column",
    data_summary = list(n_features = 100, n_groups = 30),
    width_mm = 72, rasterize = FALSE
  )
  expect_equal(spec$recommended$width_mm, 72)
  expect_lte(spec$recommended$height_mm, 120)
  expect_false(spec$recommended$rasterize)
  expect_true(spec$recommended$pagination$required)
  expect_error(sn_figure_spec(profile = "imaginary"), "Unknown figure profile")
})

test_that("figure profile attachment and validation preserve native ggplot", {
  plot <- ggplot2::ggplot(mtcars, ggplot2::aes(.data$wt, .data$mpg, color = factor(.data$cyl))) + ggplot2::geom_point()
  plot <- Shennong:::.sn_attach_figure_spec(
    plot, "embedding",
    list(n_points = nrow(mtcars), n_groups = 3, labels = c("four", "six", "eight")),
    source_data = mtcars
  )
  applied <- sn_apply_figure_profile(plot, "single_column")
  report <- sn_validate_figure(applied)
  expect_s3_class(applied, "ggplot")
  expect_s3_class(sn_figure_spec(applied), "sn_figure_spec")
  expect_s3_class(report, "sn_figure_validation")
  expect_true(report$valid)
  expect_equal(report$spec$profile, "single_column")
})

test_that("figure validation predicts category, label, and network overload", {
  plot <- ggplot2::ggplot(data.frame(x = 1, y = 1), ggplot2::aes(.data$x, .data$y)) + ggplot2::geom_point()
  plot <- Shennong:::.sn_attach_figure_spec(
    plot, "network",
    list(n_points = 1, n_groups = 35, n_legend_items = 35, max_label_chars = 40, n_nodes = 50, n_edges = 2000)
  )
  report <- sn_validate_figure(plot)
  expect_identical(report$status, "WARNING")
  expect_gte(length(report$warnings), 3L)
  expect_true(any(grepl("1,000 edges", report$warnings, fixed = TRUE)))
})

test_that("publication export supports vector and raster formats", {
  plot <- ggplot2::ggplot(mtcars, ggplot2::aes(.data$wt, .data$mpg)) + ggplot2::geom_point()
  plot <- Shennong:::.sn_attach_figure_spec(plot, "embedding", list(n_points = nrow(mtcars)), source_data = mtcars)
  directory <- withr::local_tempdir()
  files <- file.path(directory, paste0("figure.", c("pdf", "svg", "png", "tiff")))
  invisible(lapply(files, function(file) sn_save_figure(plot, file, profile = "single_column")))
  expect_true(all(file.exists(files)))
  expect_true(all(file.info(files)$size > 100))
  expect_error(sn_save_figure(plot, file.path(directory, "bad.jpg")), "Unsupported figure format")
})

test_that("figure bundles retain plots, source data, specs, sessions, and checksums", {
  plot <- ggplot2::ggplot(mtcars, ggplot2::aes(.data$wt, .data$mpg)) + ggplot2::geom_point()
  plot <- Shennong:::.sn_attach_figure_spec(plot, "embedding", list(n_points = nrow(mtcars)), source_data = mtcars)
  directory <- file.path(withr::local_tempdir(), "Figure_2A")
  manifest <- sn_export_figure_bundle(plot, directory, formats = c("pdf", "png"))
  expected <- c("Figure_2A.pdf", "Figure_2A.png", "Figure_2A_data.csv", "Figure_2A_spec.yml", "Figure_2A_session.txt", "Figure_2A_manifest.json")
  expect_true(all(file.exists(file.path(directory, expected))))
  expect_equal(length(manifest$files), 5L)
  expect_true(all(vapply(manifest$files, function(x) nzchar(x$md5), logical(1))))
})

test_that("DE and enrichment plots consume standardized evidence", {
  set.seed(77)
  de <- data.frame(
    gene = paste0("gene_", 1:40), log2FoldChange = stats::rnorm(40),
    padj = stats::runif(40), baseMean = seq_len(40), comparison = rep(c("A", "B"), each = 20)
  )
  enrichment <- data.frame(
    ID = paste0("P", 1:5), Description = paste("Pathway", 1:5),
    p.adjust = seq(0.01, 0.05, by = 0.01), NES = c(2, 1.5, -1.2, 1, -2),
    core_enrichment = c("G1/G2", "G2/G3", "G4/G5", "G1/G3", "G5/G6")
  )
  plots <- c(
    lapply(c("volcano", "ma", "effect", "heatmap"), function(type) sn_plot_de(de, type = type)),
    lapply(c("dot", "bar", "network", "emap"), function(type) sn_plot_enrichment(enrichment, type = type)),
    list(sn_plot_gsea(enrichment))
  )
  expect_true(all(vapply(plots, inherits, logical(1), "ggplot")))
  expect_true(all(vapply(plots, function(x) !is.null(attr(x, "shennong_figure_spec")), logical(1))))
})

test_that("core diagnostic plot helpers expose figure specs", {
  qc <- list(by_sample = data.frame(sample = c("S1", "S2"), qc_score = c(.8, .6), n_cells = c(100, 80), retention_fraction = c(.9, .7)))
  metadata <- data.frame(nFeature_RNA = 1:20, nCount_RNA = 21:40, percent.mt = seq(1, 5, length.out = 20), sample = rep(c("S1", "S2"), 10), row.names = paste0("C", 1:20))
  before <- matrix(1:30, 6, dimnames = list(paste0("G", 1:6), paste0("C", 1:5)))
  after <- pmax(before - 2, 0); dimnames(after) <- dimnames(before)
  sweep <- list(summary = data.frame(resolution = c(.2, .4, .6), n_clusters = c(3, 5, 7), composite_score = c(.5, .8, .7)))
  integration <- list(summary = data.frame(metric = c("lisi", "silhouette"), category = c("batch_removal", "biology"), scaled_score = c(.8, .7)))
  clusters <- data.frame(res.0.2 = rep(1:2, each = 10), res.0.4 = rep(1:4, each = 5))
  plots <- list(
    sn_plot_qc(qc),
    sn_plot_qc_thresholds(metadata, thresholds = list(percent.mt = c(2, 4)), sample_by = "sample"),
    sn_plot_ambient_correction(before, after),
    sn_plot_resolution_sweep(sweep),
    sn_plot_integration(integration),
    sn_plot_cluster_tree(clusters)
  )
  expect_true(all(vapply(plots, inherits, logical(1), "ggplot")))
  expect_true(all(vapply(plots, function(x) !is.null(sn_figure_spec(x)), logical(1))))
})

test_that("Seurat diagnostic plots and reference projection run on bundled data", {
  skip_if_not_installed("Seurat")
  object <- pbmc_small
  object$doublet_class <- rep(c("singlet", "doublet"), length.out = ncol(object))
  object <- suppressWarnings(Seurat::NormalizeData(object, verbose = FALSE))
  object <- suppressWarnings(Seurat::FindVariableFeatures(object, nfeatures = 20, verbose = FALSE))
  object <- suppressWarnings(Seurat::ScaleData(object, features = Seurat::VariableFeatures(object), verbose = FALSE))
  object <- suppressWarnings(Seurat::RunPCA(object, features = Seurat::VariableFeatures(object), npcs = 5, verbose = FALSE))
  doublets <- sn_plot_doublets(object, class_col = "doublet_class", reduction = "pca")
  hvg <- sn_plot_hvg(object, label_n = 3)
  elbow <- sn_plot_elbow(object, reduction = "pca", ndims = 5)
  cells <- data.frame(cell = colnames(object), prediction = rep(c("T", "B"), length.out = ncol(object)), prediction_score = seq(.5, .9, length.out = ncol(object)))
  result <- list(
    schema_version = "1.0", analysis_type = "annotation", name = "annotation", method = "test", backend = "test",
    input = list(), parameters = list(), tables = list(primary = cells, cells = cells), embeddings = list(), graphs = list(), models = list(),
    diagnostics = list(), warnings = character(), provenance = Shennong:::.sn_analysis_provenance()
  )
  object <- sn_store_result(object, "annotation", "annotation", result)
  projection <- sn_plot_reference_projection(object, reduction = "pca")
  plots <- list(doublets, hvg, elbow, projection)
  expect_true(all(vapply(plots, inherits, logical(1), "ggplot")))
  expect_true(all(vapply(plots, function(x) !is.null(attr(x, "shennong_figure_spec")), logical(1))))
})

test_that("migrated core plots carry automatic figure specs", {
  data <- transform(mtcars, cyl = factor(cyl))
  plots <- list(
    sn_plot_boxplot(data, x = cyl, y = mpg),
    sn_plot_barplot(data, x = cyl, y = mpg),
    sn_plot_composition(data.frame(sample = rep(c("S1", "S2"), each = 2), cell_type = rep(c("T", "B"), 2), proportion = c(60, 40, 30, 70)), x = sample, fill = cell_type)
  )
  expect_true(all(vapply(plots, function(x) inherits(attr(x, "shennong_figure_spec"), "sn_figure_spec"), logical(1))))
})
