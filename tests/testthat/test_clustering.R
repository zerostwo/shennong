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

make_zero_layer <- function(object) {
  Matrix::Matrix(
    0,
    nrow = nrow(object),
    ncol = ncol(object),
    sparse = TRUE,
    dimnames = dimnames(SeuratObject::LayerData(object, layer = "counts"))
  )
}

test_that("sn_run_cluster clusters a single dataset with the standard workflow", {
  skip_if_not_installed("Seurat")

  object <- make_test_object(seed = 1, prefix = "single")
  clustered <- sn_run_cluster(
    object = object,
    normalization_method = "seurat",
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
    normalization_method = "sctransform",
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
    normalization_method = "seurat",
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

test_that("sn_run_cluster skips cell-cycle scoring cleanly when markers do not overlap", {
  skip_if_not_installed("Seurat")

  object1 <- make_test_object(seed = 13, prefix = "batchx")
  object1$sample <- "pbmc1k"
  object2 <- make_test_object(seed = 14, prefix = "batchy")
  object2$sample <- "pbmc3k"
  merged <- merge(x = object1, y = object2, add.cell.ids = c("pbmc1k", "pbmc3k"))

  expect_no_error(
    clustered <- suppressWarnings(sn_run_cluster(
      object = merged,
      batch = "sample",
      species = "human",
      normalization_method = "seurat",
      nfeatures = 50,
      npcs = 10,
      dims = 1:10,
      verbose = FALSE
    ))
  )

  expect_s4_class(clustered, "Seurat")
  expect_true("seurat_clusters" %in% colnames(clustered[[]]))
  expect_false(all(c("S.Score", "G2M.Score", "Phase", "CC.Difference") %in% colnames(clustered[[]])))
})

test_that("sn_run_cluster can use a non-default layer without overwriting counts", {
  skip_if_not_installed("Seurat")

  object <- make_test_object(seed = 11, prefix = "layered")
  original_counts <- SeuratObject::LayerData(object, layer = "counts")
  zero_counts <- make_zero_layer(object)

  SeuratObject::LayerData(object, layer = "counts") <- zero_counts
  SeuratObject::LayerData(object, layer = "decontaminated_counts") <- original_counts

  clustered <- sn_run_cluster(
    object = object,
    normalization_method = "seurat",
    assay = "RNA",
    layer = "decontaminated_counts",
    nfeatures = 50,
    block_genes = NULL,
    npcs = 10,
    dims = 1:10,
    verbose = FALSE
  )

  expect_true("seurat_clusters" %in% colnames(clustered[[]]))
  expect_equal(
    as.matrix(SeuratObject::LayerData(clustered, layer = "counts")),
    as.matrix(zero_counts)
  )
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
  expect_true(all(c(
    "decontX_contamination",
    "decontX_clusters",
    "nCount_RNA_corrected",
    "nFeature_RNA_corrected"
  ) %in% colnames(updated[[]])))
  corrected <- SeuratObject::LayerData(updated, layer = "decontaminated_counts")
  expect_equal(dim(corrected), dim(SeuratObject::LayerData(object, layer = "counts")))
  expect_true(all(corrected == round(corrected)))
})

test_that("decontX zero-count handling restores original cells by default", {
  original_counts <- Matrix::Matrix(
    matrix(
      c(
        5, 0,
        1, 2,
        0, 3
      ),
      nrow = 3,
      byrow = TRUE
    ),
    sparse = TRUE
  )
  rownames(original_counts) <- paste0("gene", 1:3)
  colnames(original_counts) <- paste0("cell", 1:2)
  corrected_counts <- Matrix::Matrix(
    matrix(
      c(
        0.4, 0,
        0.3, 0,
        0.2, 0
      ),
      nrow = 3,
      byrow = TRUE
    ),
    sparse = TRUE
  )
  rownames(corrected_counts) <- rownames(original_counts)
  colnames(corrected_counts) <- colnames(original_counts)

  handled <- Shennong:::.sn_handle_zero_count_cells(
    original_counts = original_counts,
    corrected_counts = round(corrected_counts),
    remove_zero_count_cells = FALSE
  )

  expect_equal(as.matrix(handled$counts), as.matrix(original_counts))
  expect_equal(handled$zero_cells, colnames(original_counts))
  expect_length(handled$removed_cells, 0)
})

test_that("decontX zero-count handling can remove affected cells", {
  original_counts <- Matrix::Matrix(
    matrix(
      c(
        5, 0,
        1, 2,
        0, 3
      ),
      nrow = 3,
      byrow = TRUE
    ),
    sparse = TRUE
  )
  rownames(original_counts) <- paste0("gene", 1:3)
  colnames(original_counts) <- paste0("cell", 1:2)
  corrected_counts <- Matrix::Matrix(0, nrow = 3, ncol = 2, sparse = TRUE)
  colnames(corrected_counts) <- colnames(original_counts)
  rownames(corrected_counts) <- rownames(original_counts)

  handled <- Shennong:::.sn_handle_zero_count_cells(
    original_counts = original_counts,
    corrected_counts = corrected_counts,
    remove_zero_count_cells = TRUE
  )

  expect_equal(ncol(handled$counts), 0)
  expect_equal(handled$removed_cells, colnames(original_counts))
})

test_that("sn_find_doublets can analyze a non-default layer", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("scDblFinder")
  skip_if_not_installed("SingleCellExperiment")

  object <- make_test_object(seed = 12, prefix = "doublets", n_genes = 250, n_cells = 40)
  original_counts <- SeuratObject::LayerData(object, layer = "counts")
  zero_counts <- make_zero_layer(object)

  object$precluster <- rep(c("a", "b"), each = 20)
  SeuratObject::LayerData(object, layer = "counts") <- zero_counts
  SeuratObject::LayerData(object, layer = "decontaminated_counts") <- original_counts

  updated <- sn_find_doublets(
    object = object,
    clusters = "precluster",
    assay = "RNA",
    layer = "decontaminated_counts",
    ncores = 1
  )

  expect_true(all(c("scDblFinder.class", "scDblFinder.score") %in% colnames(updated[[]])))
  expect_equal(
    as.matrix(SeuratObject::LayerData(updated, layer = "counts")),
    as.matrix(zero_counts)
  )
})
