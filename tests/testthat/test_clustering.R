library(testthat)

make_test_object <- function(seed, prefix, n_genes = 200, n_cells = 40) {
  set.seed(seed)
  counts <- matrix(rpois(n_genes * n_cells, lambda = 2), nrow = n_genes, ncol = n_cells)
  counts[1:20, 1:(n_cells / 2)] <- counts[1:20, 1:(n_cells / 2)] + 5
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  rownames(counts) <- paste0("gene", seq_len(n_genes))
  colnames(counts) <- paste0(prefix, "_cell", seq_len(n_cells))

  sn_initialize_seurat_object(
    x = counts,
    project = prefix
  )
}

test_that("sn_run_cluster clusters a single dataset with the standard workflow", {
  skip_if_not_installed("Seurat")

  object <- make_test_object(seed = 1, prefix = "single")
  clustered <- sn_run_cluster(
    object = object,
    pipeline = "standard",
    nfeatures = 50,
    block_genes = NULL,
    npcs = 10,
    dims = 1:10,
    verbose = FALSE
  )

  expect_s4_class(clustered, "Seurat")
  expect_true("seurat_clusters" %in% colnames(clustered[[]]))
  expect_true("umap" %in% names(clustered@reductions))
})

test_that("sn_run_cluster supports the SCTransform workflow for a single dataset", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("glmGamPoi")

  object <- make_test_object(seed = 2, prefix = "sct")
  clusters <- sn_run_cluster(
    object = object,
    pipeline = "sctransform",
    nfeatures = 50,
    npcs = 10,
    dims = 1:10,
    return_cluster = TRUE,
    verbose = FALSE
  )

  expect_length(clusters, ncol(object))
})

test_that("sn_run_cluster integrates batches with harmony", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("harmony")

  object1 <- make_test_object(seed = 3, prefix = "batch1")
  object1$sample <- "pbmc1k"
  object2 <- make_test_object(seed = 4, prefix = "batch2")
  object2$sample <- "pbmc3k"

  merged <- merge(x = object1, y = object2, add.cell.ids = c("pbmc1k", "pbmc3k"))

  clustered <- sn_run_cluster(
    object = merged,
    batch = "sample",
    pipeline = "standard",
    nfeatures = 50,
    block_genes = NULL,
    npcs = 10,
    dims = 1:10,
    verbose = FALSE
  )

  expect_s4_class(clustered, "Seurat")
  expect_true("harmony" %in% names(clustered@reductions))
  expect_true("seurat_clusters" %in% colnames(clustered[[]]))
})

test_that("sn_remove_ambient_contamination requires raw counts for SoupX", {
  counts <- matrix(rpois(100 * 20, lambda = 3), nrow = 100, ncol = 20)
  rownames(counts) <- paste0("gene", seq_len(100))
  colnames(counts) <- paste0("cell", seq_len(20))

  expect_error(
    sn_remove_ambient_contamination(
      x = counts,
      method = "soupx",
      return_object = FALSE
    ),
    "`raw` is required"
  )
})

test_that("sn_remove_ambient_contamination supports decontX on matrices", {
  skip_if_not_installed("celda")
  skip_if_not_installed("SingleCellExperiment")

  counts <- matrix(rpois(200 * 20, lambda = 3), nrow = 200, ncol = 20)
  rownames(counts) <- paste0("gene", seq_len(200))
  colnames(counts) <- paste0("cell", seq_len(20))
  cluster <- rep(c("a", "b"), each = 10)

  corrected <- sn_remove_ambient_contamination(
    x = counts,
    method = "decontx",
    cluster = cluster,
    return_object = FALSE,
    verbose = FALSE,
    estimateDelta = FALSE,
    maxIter = 5
  )

  expect_equal(dim(corrected), dim(counts))
  expect_equal(rownames(corrected), rownames(counts))
  expect_equal(colnames(corrected), colnames(counts))
  expect_true(all(corrected == round(corrected)))
})

test_that("sn_remove_ambient_contamination defaults to decontX and writes a new layer for Seurat", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("celda")
  skip_if_not_installed("SingleCellExperiment")

  object <- make_test_object(seed = 5, prefix = "decontx", n_genes = 120, n_cells = 20)
  cluster <- rep(c("a", "b"), each = 10)

  updated <- sn_remove_ambient_contamination(
    x = object,
    cluster = cluster,
    verbose = FALSE,
    estimateDelta = FALSE,
    maxIter = 5
  )

  expect_s4_class(updated, "Seurat")
  expect_true(all(c("decontX_contamination", "decontX_clusters") %in% colnames(updated[[]])))
  corrected <- SeuratObject::LayerData(updated, layer = "decontaminated_counts")
  expect_equal(dim(corrected), dim(SeuratObject::LayerData(object, layer = "counts")))
  expect_true(all(corrected == round(corrected)))
})
