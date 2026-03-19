library(testthat)

test_that("sn_initialize_seurat_object accepts base matrices and adds metadata", {
  skip_if_not_installed("Seurat")

  counts <- matrix(
    c(
      5, 0, 3, 1,
      2, 4, 1, 0,
      0, 1, 0, 2,
      3, 3, 3, 3
    ),
    nrow = 4,
    byrow = TRUE
  )
  rownames(counts) <- c("MT-CO1", "RPLP0", "HBB", "CD3D")
  colnames(counts) <- paste0("cell", 1:4)
  metadata <- data.frame(batch = c("a", "a", "b", "b"), row.names = colnames(counts))

  object <- sn_initialize_seurat_object(
    x = counts,
    metadata = metadata,
    project = "prep",
    sample_name = "sample_a",
    study = "study_a",
    species = "human"
  )

  expect_s4_class(object, "Seurat")
  expect_equal(unname(object$sample), rep("sample_a", 4))
  expect_equal(unname(object$study), rep("study_a", 4))
  expect_true(all(c("percent.mt", "percent.ribo", "percent.hb", "batch") %in% colnames(object[[]])))
  expect_true("sn_initialize_seurat_object" %in% names(object@commands))
})

test_that("sn_initialize_seurat_object infers species and computes QC metrics", {
  skip_if_not_installed("Seurat")

  counts <- matrix(
    c(
      5, 0, 3, 1,
      2, 4, 1, 0,
      0, 1, 0, 2,
      3, 3, 3, 3
    ),
    nrow = 4,
    byrow = TRUE
  )
  rownames(counts) <- c("MT-CO1", "RPLP0", "HBB", "CD3D")
  colnames(counts) <- paste0("cell", 1:4)

  object <- sn_initialize_seurat_object(
    x = counts,
    project = "prep-inferred"
  )

  expect_equal(sn_get_species(object), "human")
  expect_true(all(c("percent.mt", "percent.ribo", "percent.hb") %in% colnames(object[[]])))
})

test_that("sn_filter_genes drops genes below the minimum cell threshold", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(
    matrix(
      c(
        1, 1, 1, 1,
        1, 0, 0, 0,
        0, 0, 0, 1,
        2, 2, 2, 2
      ),
      nrow = 4,
      byrow = TRUE
    ),
    sparse = TRUE
  )
  rownames(counts) <- paste0("gene", 1:4)
  colnames(counts) <- paste0("cell", 1:4)
  object <- sn_initialize_seurat_object(x = counts, project = "genes")

  filtered <- sn_filter_genes(object, min_cells = 2, plot = FALSE, filter = TRUE)

  expect_equal(rownames(filtered), c("gene1", "gene4"))
})

test_that("sn_filter_genes can use a non-default layer", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(
    matrix(
      c(
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 0
      ),
      nrow = 4,
      byrow = TRUE
    ),
    sparse = TRUE
  )
  rownames(counts) <- paste0("gene", 1:4)
  colnames(counts) <- paste0("cell", 1:4)

  alt_counts <- Matrix::Matrix(
    matrix(
      c(
        1, 1, 1, 1,
        1, 0, 0, 0,
        0, 0, 0, 1,
        2, 2, 2, 2
      ),
      nrow = 4,
      byrow = TRUE
    ),
    sparse = TRUE
  )
  rownames(alt_counts) <- rownames(counts)
  colnames(alt_counts) <- colnames(counts)

  object <- sn_initialize_seurat_object(x = counts, project = "genes-layer")
  SeuratObject::LayerData(object, layer = "decontaminated_counts") <- alt_counts

  filtered <- sn_filter_genes(
    object,
    min_cells = 2,
    plot = FALSE,
    filter = TRUE,
    assay = "RNA",
    layer = "decontaminated_counts"
  )

  expect_equal(rownames(filtered), c("gene1", "gene4"))
})

test_that("sn_filter_cells adds qc flags and can subset outliers", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(matrix(rpois(80, lambda = 3), nrow = 20, ncol = 4), sparse = TRUE)
  rownames(counts) <- paste0("gene", 1:20)
  colnames(counts) <- paste0("cell", 1:4)

  object <- sn_initialize_seurat_object(x = counts, project = "cells")
  object$nCount_RNA[4] <- 1000
  object$nFeature_RNA[4] <- 1000

  flagged <- sn_filter_cells(
    x = object,
    features = c("nCount_RNA", "nFeature_RNA"),
    n = 1,
    plot = FALSE,
    filter = FALSE
  )

  expect_true(all(c("nCount_RNA_qc", "nFeature_RNA_qc") %in% colnames(flagged[[]])))
  expect_equal(unname(flagged$nCount_RNA_qc[4]), "Failed")

  filtered <- sn_filter_cells(
    x = object,
    features = c("nCount_RNA", "nFeature_RNA"),
    n = 1,
    plot = FALSE,
    filter = TRUE
  )

  expect_lt(ncol(filtered), ncol(object))
})

test_that("Seurat-returning helpers record commands in the object history", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(matrix(rpois(80, lambda = 3), nrow = 20, ncol = 4), sparse = TRUE)
  rownames(counts) <- paste0("gene", 1:20)
  colnames(counts) <- paste0("cell", 1:4)

  object <- sn_initialize_seurat_object(x = counts, project = "commands")
  object <- sn_filter_genes(object, min_cells = 1, plot = FALSE, filter = FALSE)

  command_names <- names(object@commands)
  expect_true(any(grepl("^sn_initialize_seurat_object", command_names)))
  expect_true(any(grepl("^sn_filter_genes", command_names)))
})

test_that("sn_normalize_data can normalize from a non-default layer", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("scran")
  skip_if_not_installed("SingleCellExperiment")

  counts <- Matrix::Matrix(0, nrow = 50, ncol = 60, sparse = TRUE)
  rownames(counts) <- paste0("gene", 1:50)
  colnames(counts) <- paste0("cell", 1:60)

  alt_counts <- Matrix::Matrix(matrix(rpois(50 * 60, lambda = 3), nrow = 50), sparse = TRUE)
  rownames(alt_counts) <- rownames(counts)
  colnames(alt_counts) <- colnames(counts)

  object <- sn_initialize_seurat_object(x = counts, project = "normalize-layer")
  SeuratObject::LayerData(object, layer = "decontaminated_counts") <- alt_counts

  normalized <- sn_normalize_data(
    object = object,
    method = "scran",
    assay = "RNA",
    layer = "decontaminated_counts",
    clusters = rep(c("a", "b"), each = 30)
  )

  expect_gt(sum(SeuratObject::LayerData(normalized, assay = "RNA", layer = "data")), 0)
  expect_equal(
    as.matrix(SeuratObject::LayerData(normalized, assay = "RNA", layer = "counts")),
    as.matrix(counts)
  )
})
