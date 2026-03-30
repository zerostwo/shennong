library(testthat)

make_dotplot_test_object <- function() {
  set.seed(717)
  genes <- c(
    paste0("GENE", seq_len(40)),
    "CD3D", "CD3E", "TRAC", "LCK", "LTB",
    "MS4A1", "CD79A", "HLA-DRA", "HLA-DPA1"
  )
  counts <- matrix(rpois(length(genes) * 40, lambda = 2), nrow = length(genes), ncol = 40)
  rownames(counts) <- genes
  colnames(counts) <- paste0("cell", seq_len(40))
  cell_type <- rep(rep(c("Tcell", "Bcell"), each = 5), 4)

  counts[rownames(counts) %in% c("CD3D", "CD3E", "TRAC", "LCK", "LTB"), cell_type == "Tcell"] <-
    counts[rownames(counts) %in% c("CD3D", "CD3E", "TRAC", "LCK", "LTB"), cell_type == "Tcell"] + 8
  counts[rownames(counts) %in% c("MS4A1", "CD79A", "HLA-DRA", "HLA-DPA1"), cell_type == "Bcell"] <-
    counts[rownames(counts) %in% c("MS4A1", "CD79A", "HLA-DRA", "HLA-DPA1"), cell_type == "Bcell"] + 8

  object <- sn_initialize_seurat_object(
    x = Matrix::Matrix(counts, sparse = TRUE),
    project = "dotplot-test"
  )
  object$cell_type <- cell_type
  Seurat::Idents(object) <- object$cell_type
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  sn_find_de(
    object = object,
    analysis = "markers",
    group_by = "cell_type",
    layer = "data",
    min_pct = 0,
    logfc_threshold = 0,
    store_name = "celltype_markers",
    return_object = TRUE,
    verbose = FALSE
  )
}

test_that("sn_write and sn_read round-trip a csv file", {
  input <- data.frame(a = 1:3, b = c("x", "y", "z"))
  path <- tempfile(fileext = ".csv")

  sn_write(input, path)
  output <- sn_read(path)

  expect_s3_class(output, "data.frame")
  expect_equal(output$a, input$a)
  expect_equal(output$b, input$b)
})

test_that("plot helpers return ggplot objects", {
  data <- data.frame(
    group = c("A", "B", "A", "B"),
    value = c(1, 2, 3, 4),
    fill_group = c("x", "x", "y", "y")
  )

  expect_s3_class(sn_plot_boxplot(data, x = group, y = value), "ggplot")
  expect_s3_class(sn_plot_barplot(data, x = group, y = value, fill = fill_group), "ggplot")
})

test_that("sn_plot_barplot sorts the discrete axis and works without fill", {
  data <- data.frame(
    cell_type = factor(c("Tcell", "Bcell", "Mono"), levels = c("Mono", "Tcell", "Bcell")),
    log2_fc = c(1.5, -0.2, 0.7),
    stringsAsFactors = TRUE
  )

  plot <- sn_plot_barplot(
    data,
    x = log2_fc,
    y = cell_type,
    sort_by = "log2_fc"
  )

  expect_s3_class(plot, "ggplot")
  expect_identical(levels(plot$data$cell_type), c("Bcell", "Mono", "Tcell"))
})

test_that("sn_plot_barplot can summarize repeated observations with error bars and points", {
  data <- data.frame(
    group = rep(c("A", "B"), each = 4),
    sample = rep(paste0("S", 1:4), times = 2),
    value = c(1, 2, 3, 4, 2, 3, 4, 5),
    subtype = rep(c("x", "y"), times = 4),
    stringsAsFactors = FALSE
  )

  plot <- sn_plot_barplot(
    data,
    x = group,
    y = value,
    fill = subtype,
    stat = "auto",
    errorbar = "se",
    show_points = TRUE
  )

  expect_s3_class(plot, "ggplot")
})

test_that("sn_plot_barplot maps summarized bars to the summarized value column", {
  data <- data.frame(
    cell_type_level1 = rep(c("Tcell", "Bcell"), each = 3),
    sample = rep(paste0("S", 1:3), times = 2),
    proportion = c(40, 45, 50, 20, 25, 30),
    stringsAsFactors = FALSE
  )

  plot <- sn_plot_barplot(
    data,
    x = cell_type_level1,
    y = proportion
  )

  expect_s3_class(plot, "ggplot")
  expect_true(".sn_center" %in% colnames(plot$data))
})

test_that("sn_plot_composition supports default y selection, ordering, and faceting", {
  data <- data.frame(
    cell_type = factor(c("Mono", "Mono", "cDC1", "cDC1", "B", "B"), levels = c("Mono", "cDC1", "B")),
    Mutation = factor(c("WT", "PPP2R1A", "WT", "PPP2R1A", "WT", "PPP2R1A"), levels = c("WT", "PPP2R1A")),
    study = c("A", "A", "A", "A", "B", "B"),
    proportion = c(80, 20, 60, 40, 10, 90)
  )

  plot <- sn_plot_composition(
    data,
    x = cell_type,
    fill = Mutation,
    facet_col = study,
    order_value = "WT"
  )

  expect_s3_class(plot, "ggplot")
  expect_identical(levels(plot$data$cell_type), c("B", "cDC1", "Mono"))
})

test_that("sn_plot_composition resolves named palettes and keeps axis labels", {
  data <- expand.grid(
    Mutation = c("WT", "PPP2R1A"),
    cell_type_level1 = paste0("type", seq_len(11)),
    stringsAsFactors = FALSE
  )
  data$proportion <- seq_len(nrow(data))

  plot <- sn_plot_composition(
    data,
    x = Mutation,
    y = proportion,
    fill = cell_type_level1,
    y_label = "Proportion (%)",
    palette = "Paired",
    panel_widths = 160,
    panel_heights = 80
  )

  expect_s3_class(plot, "ggplot")
  expect_identical(plot$labels$y, "Proportion (%)")
  expect_true(length(plot$scales$scales) >= 2)
})

test_that("Seurat plotting helpers return ggplot objects", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(matrix(rpois(200, lambda = 3), nrow = 20, ncol = 10), sparse = TRUE)
  rownames(counts) <- paste0("gene", seq_len(20))
  rownames(counts)[1:3] <- c("CD3D", "CD3E", "LYZ")
  colnames(counts) <- paste0("cell", seq_len(10))

  object <- sn_initialize_seurat_object(x = counts, project = "plots")
  object$group <- rep(c("A", "B"), each = 5)
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  object <- suppressWarnings(
    Seurat::FindVariableFeatures(object, nfeatures = 10, verbose = FALSE)
  )
  object <- Seurat::ScaleData(object, features = rownames(object), verbose = FALSE)
  object <- suppressWarnings(
    Seurat::RunPCA(
      object,
      features = Seurat::VariableFeatures(object),
      npcs = 5,
      verbose = FALSE
    )
  )
  object <- suppressWarnings(
    Seurat::RunUMAP(object, dims = 1:5, n.neighbors = 5, verbose = FALSE)
  )

  feature_plot <- sn_plot_feature(object, features = "CD3D", reduction = "umap", palette = "YlOrRd", direction = -1)
  multi_feature_plot <- sn_plot_feature(object, features = c("CD3D", "LYZ"), reduction = "umap", palette = "RdBu")
  dim_plot <- sn_plot_dim(object, reduction = "umap", group_by = "group")
  dim_plot_repel <- sn_plot_dim(object, reduction = "umap", group_by = "group", label = TRUE, repel = TRUE)
  dim_plot_panel_widths <- sn_plot_dim(object, reduction = "umap", group_by = "group", panel_widths = 80)
  dim_plot_no_aspect <- sn_plot_dim(object, reduction = "umap", group_by = "group", panel_widths = 80, aspect_ratio = NULL)
  feature_color_scale <- Filter(function(scale) "colour" %in% scale$aesthetics, feature_plot$scales$scales)[[1]]

  expect_s3_class(feature_plot, "ggplot")
  expect_s3_class(multi_feature_plot, "ggplot")
  expect_s3_class(dim_plot, "ggplot")
  expect_s3_class(dim_plot_repel, "ggplot")
  expect_s3_class(dim_plot_panel_widths, "ggplot")
  expect_s3_class(dim_plot_no_aspect, "ggplot")
  expect_s3_class(feature_plot$theme$axis.text, "element_blank")
  expect_s3_class(feature_plot$theme$axis.title, "element_blank")
  expect_s3_class(dim_plot$theme$axis.text, "element_blank")
  expect_s3_class(dim_plot$theme$axis.title, "element_blank")
  expect_equal(feature_color_scale$labels[[1]], "Min")
  expect_equal(tail(feature_color_scale$labels, 1), "Max")
  expect_gt(length(feature_color_scale$breaks), 2)
  expect_true(is.language(feature_plot$labels$title))
  multi_feature_scales <- lapply(
    multi_feature_plot$patches$plots,
    function(plot) Filter(function(scale) "colour" %in% scale$aesthetics, plot$scales$scales)
  )
  expect_true(all(vapply(multi_feature_scales, length, integer(1)) >= 1L))
  expect_equal(length(multi_feature_plot$guides$guides), 1)
  expected_geom <- if (rlang::is_installed("shadowtext")) "GeomShadowText" else "GeomText"
  expect_identical(class(dim_plot_repel$layers[[2]]$geom)[[1]], expected_geom)
  expect_equal(as.numeric(dim_plot$theme$legend.text$margin[[4]]), -4)
  expect_true(all(vapply(
    c(list(multi_feature_plot), multi_feature_plot$patches$plots),
    function(plot) is.language(plot$labels$title),
    logical(1)
  )))
  expect_gt(Shennong:::.sn_auto_point_size(object, NULL), 0)
  expect_s3_class(
    suppressWarnings(
      sn_plot_violin(object, features = c("CD3D", "LYZ"), group_by = "group")
    ),
    "ggplot"
  )
  if (requireNamespace("Nebulosa", quietly = TRUE)) {
    density_plot <- sn_plot_feature(
      object,
      features = c("CD3D", "LYZ"),
      reduction = "umap",
      mode = "density"
    )
    density_scale <- Filter(function(scale) "fill" %in% scale$aesthetics, density_plot$patches$plots[[1]]$scales$scales)[[1]]
    expect_s3_class(density_plot, "ggplot")
    expect_equal(density_plot$patches$layout$guides, "collect")
    expect_equal(density_scale$labels[[1]], "Min")
    expect_equal(tail(density_scale$labels, 1), "Max")
    expect_identical(density_plot$theme$plot.background$fill, "#030711")
  }
})

test_that("sn_plot_dot can reuse stored top markers from object@misc", {
  skip_if_not_installed("Seurat")

  object <- make_dotplot_test_object()
  feature_info <- Shennong:::.sn_resolve_dotplot_features(
    object = object,
    features = "top_markers",
    de_name = "celltype_markers",
    n = 2,
    marker_groups = "Tcell"
  )
  plot <- suppressWarnings(
    sn_plot_dot(
      x = object,
      features = "top_markers",
      de_name = "celltype_markers",
      n = 2,
      marker_groups = "Tcell",
      palette = "RdBu",
      legend_position = "top",
      zscore_legend_labels = "text"
    )
  )
  color_scale <- Filter(function(scale) "colour" %in% scale$aesthetics, plot$scales$scales)[[1]]

  expect_equal(feature_info$group_by, "cell_type")
  expect_length(feature_info$features, 2)
  expect_s3_class(plot, "ggplot")
  expect_equal(plot$theme$legend.position, "top")
  expect_equal(plot$labels$colour, "Z score")
  expect_equal(plot$labels$size, "Percent (%)")
  expect_equal(color_scale$breaks[[1]], -2.5)
  expect_equal(tail(color_scale$breaks, 1), 2.5)
  expect_equal(color_scale$labels[[1]], "Min")
  expect_equal(tail(color_scale$labels, 1), "Max")
  expect_gt(length(color_scale$breaks), 2)
  expect_equal(color_scale$limits, c(-2.5, 2.5))
  expect_equal(as.numeric(color_scale$guide$params$theme$legend.key.width), 8)
  expect_equal(as.numeric(color_scale$guide$params$theme$legend.key.height), 32)
  expect_error(
    Shennong:::.sn_resolve_dotplot_features(
      object = object,
      features = "top_markers",
      de_name = "missing"
    ),
    "No stored DE result named 'missing'"
  )
})

test_that("palette helpers list and resolve palettes", {
  expect_no_error(sn_list_palettes())
  expect_no_error(sn_list_palettes(display = "plot", source = "ggokabeito"))
  palette_tbl <- sn_list_palettes(display = "none")

  expect_s3_class(palette_tbl, "data.frame")
  expect_true(all(c("name", "source", "max_n", "preview") %in% colnames(palette_tbl)))
  expect_true("Paired" %in% palette_tbl$name)
  expect_true("ZhangJian2024" %in% palette_tbl$name)
  expect_true("OkabeIto" %in% palette_tbl$name)
  expect_true("viridis" %in% palette_tbl$name)
  expect_equal(
    palette_tbl$source[match("OkabeIto", palette_tbl$name)],
    "ggokabeito"
  )
  expect_equal(
    palette_tbl$source[match("viridis", palette_tbl$name)],
    "viridis"
  )

  paired <- sn_get_palette("Paired", n = 14)
  expect_length(paired, 14)
  expect_true(all(nzchar(paired)))

  paired_native <- sn_get_palette("Paired", n = 12)
  expect_identical(toupper(paired_native[[11]]), "#ECD577")
  expect_false("#FFFF99" %in% toupper(paired_native))

  rd_bu <- sn_get_palette("RdBu", palette_type = "continuous", n = 32, direction = -1)
  expect_length(rd_bu, 32)
  expect_true(all(nzchar(rd_bu)))

  okabe <- sn_get_palette("OkabeIto")
  expect_identical(okabe[[1]], "#E69F00")
  expect_identical(okabe[[9]], "#000000")

  viridis <- sn_get_palette("viridis", n = 8)
  expect_length(viridis, 8)
  expect_true(all(nzchar(viridis)))

  custom <- sn_get_palette(c("#111111", "#222222", "#333333"))
  expect_identical(custom, c("#111111", "#222222", "#333333"))
})

test_that("sn_add_data_from_anndata adds metadata and embeddings", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(matrix(rpois(40, lambda = 3), nrow = 10, ncol = 4), sparse = TRUE)
  rownames(counts) <- paste0("gene", seq_len(10))
  colnames(counts) <- paste0("cell", seq_len(4))
  object <- sn_initialize_seurat_object(x = counts, project = "ann")

  metadata <- data.frame(cell_type = c("T", "T", "B", "B"), row.names = colnames(object))
  meta_path <- tempfile(fileext = ".csv")
  write.csv(metadata, meta_path)

  umap <- data.frame(UMAP_1 = c(1, 2, 3, 4), UMAP_2 = c(4, 3, 2, 1), row.names = colnames(object))
  umap_path <- tempfile(fileext = ".csv")
  write.csv(umap, umap_path)

  updated <- sn_add_data_from_anndata(
    object = object,
    metadata_path = meta_path,
    umap_path = umap_path
  )

  expect_true("cell_type" %in% colnames(updated[[]]))
  expect_true("umap" %in% names(updated@reductions))
})
