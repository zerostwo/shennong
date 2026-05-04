library(testthat)

test_that("sn_run_milo returns neighborhood DA results and can annotate neighborhoods", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("miloR")
  skip_if_not_installed("SingleCellExperiment")

  set.seed(1)
  counts <- Matrix::Matrix(matrix(rpois(60 * 24, lambda = 3), nrow = 60, ncol = 24), sparse = TRUE)
  rownames(counts) <- paste0("gene", seq_len(60))
  colnames(counts) <- paste0("cell", seq_len(24))

  object <- sn_initialize_seurat_object(
    x = counts,
    project = "milo-test",
    species = "human"
  )
  object$sample <- rep(paste0("S", 1:6), each = 4)
  object$Mutation <- rep(c("WT", "WT", "WT", "PPP2R1A", "PPP2R1A", "PPP2R1A"), each = 4)
  object$cell_type <- rep(c("Tcell", "Bcell"), length.out = ncol(object))
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  object <- suppressWarnings(Seurat::FindVariableFeatures(object, nfeatures = 30, verbose = FALSE))
  object <- Seurat::ScaleData(object, features = Seurat::VariableFeatures(object), verbose = FALSE)
  object <- suppressWarnings(
    Seurat::RunPCA(
      object,
      features = Seurat::VariableFeatures(object),
      npcs = 10,
      verbose = FALSE
    )
  )

  result <- sn_run_milo(
    object,
    sample_by = "sample",
    group_by = "Mutation",
    contrast = c("PPP2R1A", "WT"),
    reduction = "pca",
    dims = 1:5,
    k = 5,
    prop = 0.5,
    annotation_by = "cell_type"
  )

  expect_s3_class(result, "data.frame")
  expect_true(all(c("logFC", "PValue", "FDR", "Nhood", "comparison", "cell_type") %in% colnames(result)))
  expect_true(all(result$comparison == "PPP2R1A vs WT"))
})

test_that("sn_run_milo can return intermediate milo objects", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("miloR")
  skip_if_not_installed("SingleCellExperiment")

  set.seed(2)
  counts <- Matrix::Matrix(matrix(rpois(40 * 16, lambda = 2), nrow = 40, ncol = 16), sparse = TRUE)
  rownames(counts) <- paste0("gene", seq_len(40))
  colnames(counts) <- paste0("cell", seq_len(16))

  object <- sn_initialize_seurat_object(
    x = counts,
    project = "milo-intermediate",
    species = "human"
  )
  object$sample <- rep(paste0("S", 1:4), each = 4)
  object$group <- rep(c("A", "A", "B", "B"), each = 4)
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  object <- suppressWarnings(Seurat::FindVariableFeatures(object, nfeatures = 20, verbose = FALSE))
  object <- Seurat::ScaleData(object, features = Seurat::VariableFeatures(object), verbose = FALSE)
  object <- suppressWarnings(
    Seurat::RunPCA(
      object,
      features = Seurat::VariableFeatures(object),
      npcs = 8,
      verbose = FALSE
    )
  )

  result <- suppressWarnings(
    sn_run_milo(
      object,
      sample_by = "sample",
      group_by = "group",
      reduction = "pca",
      dims = 1:4,
      k = 4,
      prop = 0.5,
      return_intermediate = TRUE,
      verbose = FALSE
    )
  )

  expect_true(all(c("table", "design_df", "milo") %in% names(result)))
  expect_s3_class(result$table, "data.frame")
  expect_true(inherits(result$milo, "Milo"))
})

test_that("milo results can be stored, retrieved, and plotted", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(matrix(rpois(30 * 6, lambda = 2), nrow = 30, ncol = 6), sparse = TRUE)
  rownames(counts) <- paste0("gene", seq_len(30))
  colnames(counts) <- paste0("cell", seq_len(6))

  object <- sn_initialize_seurat_object(
    x = counts,
    project = "milo-store",
    species = "human"
  )

  milo_tbl <- data.frame(
    logFC = c(1.2, -0.8, 0.4),
    PValue = c(0.01, 0.2, 0.03),
    FDR = c(0.02, 0.3, 0.05),
    SpatialFDR = c(0.03, 0.25, 0.06),
    Nhood = 1:3,
    cell_type = c("Tcell", "Bcell", "Myeloid"),
    stringsAsFactors = FALSE
  )

  object <- sn_store_milo(
    object,
    result = milo_tbl,
    store_name = "demo_milo",
    sample_by = "sample",
    group_by = "group",
    annotation_by = "cell_type",
    return_object = TRUE
  )

  retrieved <- sn_get_milo_result(
    object,
    milo_name = "demo_milo",
    annotation = "Tcell",
    spatial_fdr = 0.05
  )
  plot <- sn_plot_milo(object, milo_name = "demo_milo", annotation_by = "cell_type")

  expect_equal(nrow(retrieved), 1)
  expect_equal(retrieved$cell_type[[1]], "Tcell")
  expect_s3_class(plot, "ggplot")
  expect_true("milo_results" %in% names(object@misc))
})

test_that("sn_run_milo requires a constant group label within each sample", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("miloR")
  skip_if_not_installed("SingleCellExperiment")

  counts <- Matrix::Matrix(matrix(rpois(20 * 8, lambda = 2), nrow = 20, ncol = 8), sparse = TRUE)
  rownames(counts) <- paste0("gene", seq_len(20))
  colnames(counts) <- paste0("cell", seq_len(8))

  object <- sn_initialize_seurat_object(
    x = counts,
    project = "milo-bad-design",
    species = "human"
  )
  object$sample <- rep(c("S1", "S2"), each = 4)
  object$group <- c("A", "A", "B", "B", "B", "B", "B", "B")
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  object <- suppressWarnings(Seurat::FindVariableFeatures(object, nfeatures = 10, verbose = FALSE))
  object <- Seurat::ScaleData(object, features = Seurat::VariableFeatures(object), verbose = FALSE)
  object <- suppressWarnings(
    Seurat::RunPCA(
      object,
      features = Seurat::VariableFeatures(object),
      npcs = 5,
      verbose = FALSE
    )
  )

  expect_error(
    sn_run_milo(
      object,
      sample_by = "sample",
      group_by = "group",
      reduction = "pca",
      dims = 1:3,
      k = 3,
      prop = 0.5,
      verbose = FALSE
    ),
    "not constant within samples"
  )
})
