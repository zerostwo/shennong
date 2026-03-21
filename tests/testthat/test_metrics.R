library(testthat)

make_metrics_test_object <- function() {
  set.seed(101)
  counts <- Matrix::Matrix(matrix(rpois(80 * 120, lambda = 3), nrow = 80, ncol = 120), sparse = TRUE)
  rownames(counts) <- paste0("gene", seq_len(80))
  colnames(counts) <- paste0("cell", seq_len(120))

  object <- sn_initialize_seurat_object(x = counts, project = "metrics-test")
  object$sample <- rep(c("s1", "s2", "s3", "s4"), each = 30)
  object$cluster_id <- rep(c("A", "B"), each = 60)
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  object <- suppressWarnings(Seurat::FindVariableFeatures(object, nfeatures = 30, verbose = FALSE))
  object <- Seurat::ScaleData(object, features = Seurat::VariableFeatures(object), verbose = FALSE)
  suppressWarnings(Seurat::RunPCA(
    object,
    features = Seurat::VariableFeatures(object),
    npcs = 5,
    verbose = FALSE
  ))
}

test_that("sn_calculate_lisi returns one score per cell from the requested reduction", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("lisi")

  object <- make_metrics_test_object()
  lisi_scores <- sn_calculate_lisi(object, reduction = "pca", label = "sample")

  expect_s3_class(lisi_scores, "data.frame")
  expect_equal(nrow(lisi_scores), ncol(object))
  expect_true(all(c("cell_id", "sample") %in% colnames(lisi_scores)))
})

test_that("sn_calculate_rogue returns a score and validates metadata columns", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("ROGUE")
  skip_if_not_installed("tibble")

  library(tibble)

  object <- make_metrics_test_object()
  rogue_score <- sn_calculate_rogue(object, cluster = "cluster_id")

  expect_true(is.numeric(rogue_score))
  expect_length(rogue_score, 1)
  expect_error(
    sn_calculate_rogue(object, cluster = "missing_cluster", sample = "sample"),
    "Specified cluster column not found"
  )
  expect_error(
    sn_calculate_rogue(object, cluster = "cluster_id", sample = "missing_sample"),
    "Specified sample column not found"
  )
})
