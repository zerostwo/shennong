make_bpcells_preprocessing_fixture <- function(n_features = 12L, n_cells = 8L) {
  skip_if_not_installed("BPCells")
  skip_if_not_installed("SeuratObject")

  counts <- Matrix::Matrix(
    matrix(
      as.integer((seq_len(n_features * n_cells) %% 5L) > 0L),
      nrow = n_features,
      ncol = n_cells
    ),
    sparse = TRUE
  )
  rownames(counts) <- paste0("gene", seq_len(n_features))
  colnames(counts) <- paste0("cell", seq_len(n_cells))

  path <- tempfile("shennong-bpcells-")
  BPCells::write_matrix_dir(
    BPCells::convert_matrix_type(counts, "uint32_t"),
    dir = path
  )
  list(
    counts = counts,
    path = path,
    matrix = BPCells::open_matrix_dir(path)
  )
}

test_that("sn_initialize_seurat_object preserves BPCells counts on disk", {
  fixture <- make_bpcells_preprocessing_fixture()
  on.exit(unlink(fixture$path, recursive = TRUE), add = TRUE)

  object_from_dir <- sn_initialize_seurat_object(
    x = fixture$matrix,
    project = "bpcells-init"
  )
  rename_dims <- SeuratObject::LayerData(
    object_from_dir,
    assay = "RNA",
    layer = "counts"
  )
  object_from_rename <- sn_initialize_seurat_object(
    x = rename_dims,
    project = "bpcells-rename-dims-init"
  )

  for (object in list(object_from_dir, object_from_rename)) {
    layer <- SeuratObject::LayerData(object, assay = "RNA", layer = "counts")
    expect_true(Shennong:::.sn_is_iterable_matrix(layer))
    expect_false(inherits(layer, "Matrix"))
    expect_equal(
      as.matrix(methods::as(layer, "dgCMatrix")),
      as.matrix(fixture$counts)
    )
  }
})

test_that("BPCells sparse materialization does not use a dense intermediate", {
  fixture <- make_bpcells_preprocessing_fixture()
  on.exit(unlink(fixture$path, recursive = TRUE), add = TRUE)
  object <- SeuratObject::CreateSeuratObject(counts = fixture$matrix)
  layer <- SeuratObject::LayerData(object, assay = "RNA", layer = "counts")

  materialized <- Shennong:::.sn_as_sparse_matrix(layer)

  expect_s4_class(materialized, "dgCMatrix")
  expect_equal(as.matrix(materialized), as.matrix(fixture$counts))
})

test_that("sn_find_doublets materializes BPCells counts one sample at a time", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("scDblFinder")
  fixture <- make_bpcells_preprocessing_fixture()
  on.exit(unlink(fixture$path, recursive = TRUE), add = TRUE)

  object <- sn_initialize_seurat_object(x = fixture$matrix, project = "bpcells-doublets")
  object$sample <- rep(c("sample_a", "sample_b"), each = ncol(object) / 2L)
  expected_scores <- stats::setNames(seq_len(ncol(object)) / 100, colnames(object))
  seen_cells <- list()

  local_mocked_bindings(
    .sn_run_scDblFinder = function(sce, ..., returnType, BPPARAM) {
      sample_counts <- SummarizedExperiment::assay(sce, "counts")
      expect_s4_class(sample_counts, "dgCMatrix")
      expect_identical(returnType, "scores")
      expect_s4_class(BPPARAM, "SerialParam")
      seen_cells[[length(seen_cells) + 1L]] <<- colnames(sce)
      data.frame(
        score = unname(expected_scores[colnames(sce)]),
        class = rep("singlet", ncol(sce)),
        row.names = colnames(sce)
      )
    },
    .env = asNamespace("Shennong")
  )

  updated <- sn_find_doublets(
    object,
    group_by = "sample",
    min_features = 1,
    ncores = 1
  )

  expect_length(seen_cells, 2L)
  expect_true(all(lengths(seen_cells) == ncol(object) / 2L))
  expect_setequal(unlist(seen_cells, use.names = FALSE), colnames(object))
  expect_equal(
    unname(updated$scDblFinder.score),
    unname(expected_scores[colnames(object)])
  )
  expect_true(all(as.character(updated$scDblFinder.class) == "singlet"))
  expect_true(Shennong:::.sn_is_iterable_matrix(
    SeuratObject::LayerData(updated, assay = "RNA", layer = "counts")
  ))
})

test_that("sn_find_doublets requires grouping for BPCells counts", {
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("scDblFinder")
  fixture <- make_bpcells_preprocessing_fixture()
  on.exit(unlink(fixture$path, recursive = TRUE), add = TRUE)
  object <- sn_initialize_seurat_object(x = fixture$matrix, project = "bpcells-ungrouped")

  expect_error(
    sn_find_doublets(object, min_features = 1),
    "requires `group_by`"
  )
})
