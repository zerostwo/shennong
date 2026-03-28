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

  expect_s3_class(sn_plot_feature(object, features = "CD3D", reduction = "umap"), "ggplot")
  expect_s3_class(
    suppressWarnings(
      sn_plot_violin(object, features = c("CD3D", "LYZ"), group_by = "group")
    ),
    "ggplot"
  )
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
      marker_groups = "Tcell"
    )
  )

  expect_equal(feature_info$group_by, "cell_type")
  expect_length(feature_info$features, 2)
  expect_s3_class(plot, "ggplot")
  expect_error(
    Shennong:::.sn_resolve_dotplot_features(
      object = object,
      features = "top_markers",
      de_name = "missing"
    ),
    "No stored DE result named 'missing'"
  )
})

test_that("show_all_palettes prints without error", {
  expect_no_error(show_all_palettes())
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
