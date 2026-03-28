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
    sample_col = "sample",
    group_col = "Mutation",
    contrast = c("PPP2R1A", "WT"),
    reduction = "pca",
    dims = 1:5,
    k = 5,
    prop = 0.5,
    annotation_col = "cell_type"
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
      sample_col = "sample",
      group_col = "group",
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
      sample_col = "sample",
      group_col = "group",
      reduction = "pca",
      dims = 1:3,
      k = 3,
      prop = 0.5,
      verbose = FALSE
    ),
    "not constant within samples"
  )
})
