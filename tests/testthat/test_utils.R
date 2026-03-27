library(testthat)

make_utils_test_object <- function(seed = 1, prefix = "utils", n_genes = 12, n_cells = 6) {
  set.seed(seed)
  counts <- matrix(rpois(n_genes * n_cells, lambda = 2), nrow = n_genes, ncol = n_cells)
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  rownames(counts) <- paste0("gene", seq_len(n_genes))
  colnames(counts) <- paste0(prefix, "_cell", seq_len(n_cells))

  sn_initialize_seurat_object(
    x = counts,
    project = prefix,
    species = "human"
  )
}

make_split_layer_object <- function() {
  object1 <- make_utils_test_object(seed = 11, prefix = "splita", n_genes = 8, n_cells = 4)
  object2 <- make_utils_test_object(seed = 12, prefix = "splitb", n_genes = 8, n_cells = 4)
  object2 <- object2[c(1:6, 8), ]

  merge(x = object1, y = object2, add.cell.ids = c("splita", "splitb"))
}

test_that("sn_check_file reports missing files without stopping when requested", {
  existing <- tempfile("existing-")
  file.create(existing)
  missing <- tempfile("missing-")

  expect_equal(
    sn_check_file(c(existing, missing), stop = FALSE),
    missing
  )
})

test_that("sn_check_file errors for missing files by default", {
  expect_error(
    sn_check_file(tempfile("missing-")),
    "does not exist"
  )
})

test_that("sn_set_path creates a directory and returns its path", {
  dir_path <- tempfile("shennong-dir-")

  expect_equal(sn_set_path(dir_path), dir_path)
  expect_true(dir.exists(dir_path))
})

test_that("sn_get_species returns the explicit species when provided", {
  expect_equal(sn_get_species(object = NULL, species = "human"), "human")
})

test_that("sn_get_species infers human species from feature names", {
  features <- c("MT-CO1", "RPLP0", "CD3D", "TRAC")
  expect_equal(sn_get_species(object = features), "human")
})

test_that("sn_get_species infers mouse species from feature names", {
  features <- c("mt-Co1", "Rplp0", "Cd3d", "Trac")
  expect_equal(sn_get_species(object = features), "mouse")
})

test_that("sparse-matrix helpers coerce common inputs without changing values", {
  dense <- matrix(
    c(
      1, 4,
      2, 3
    ),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("gene1", "gene2"), c("cell1", "cell2"))
  )
  dense_df <- as.data.frame(dense)
  dense_sparse <- Matrix::Matrix(dense, sparse = TRUE)
  general_sparse <- methods::as(dense_sparse, "generalMatrix")

  from_matrix <- Shennong:::.sn_as_sparse_matrix(dense)
  from_df <- Shennong:::.sn_as_sparse_matrix(dense_df)
  from_sparse <- Shennong:::.sn_as_sparse_matrix(dense_sparse)
  from_general <- Shennong:::.sn_as_sparse_matrix(general_sparse)

  expect_s4_class(from_matrix, "dgCMatrix")
  expect_s4_class(from_df, "dgCMatrix")
  expect_identical(from_sparse, dense_sparse)
  expect_true(inherits(from_general, "CsparseMatrix"))
  expect_equal(as.matrix(from_matrix), dense)
  expect_equal(as.matrix(from_df), dense)
  expect_equal(as.matrix(from_general), dense)
})

test_that("sparse aggregation helpers preserve grouped sums", {
  counts <- Matrix::Matrix(
    matrix(
      c(
        1, 0, 2, 0,
        0, 3, 0, 4,
        5, 0, 6, 0
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(
        c("gene1", "gene2", "gene3"),
        c("cell1", "cell2", "cell3", "cell4")
      )
    ),
    sparse = TRUE
  )

  by_sample <- Shennong:::.sn_aggregate_columns_by_group(
    counts,
    groups = c("s1", "s2", "s1", "s2")
  )
  by_gene <- Shennong:::.sn_aggregate_rows_by_group(
    counts,
    groups = c("set1", "set2", "set1")
  )

  expect_equal(
    as.matrix(by_sample),
    matrix(
      c(
        3, 0,
        0, 7,
        11, 0
      ),
      nrow = 3,
      byrow = TRUE,
      dimnames = list(rownames(counts), c("s1", "s2"))
    )
  )
  expect_equal(
    as.matrix(by_gene),
    matrix(
      c(
        6, 0, 8, 0,
        0, 3, 0, 4
      ),
      nrow = 2,
      byrow = TRUE,
      dimnames = list(c("set1", "set2"), colnames(counts))
    )
  )
})

test_that("exact knn helper matches brute-force neighbors on small embeddings", {
  embeddings <- rbind(
    c(0, 0),
    c(1, 0),
    c(0, 2),
    c(3, 0)
  )
  rownames(embeddings) <- paste0("cell", seq_len(nrow(embeddings)))

  knn <- Shennong:::.sn_exact_knn(
    embeddings = embeddings,
    k = 2,
    include_distance = TRUE,
    block_size = 2
  )

  expect_equal(knn$idx[1, ], c(2, 3))
  expect_equal(round(knn$dist[1, ], 6), round(c(1, 2), 6))
  expect_equal(knn$idx[4, ], c(2, 1))
})

test_that("Seurat command logging stores named command objects on the Seurat object", {
  skip_if_not_installed("Seurat")

  object <- make_utils_test_object()
  logged <- Shennong:::.sn_log_seurat_command(
    object = object,
    assay = "RNA",
    name = "test_command"
  )

  expect_true("test_command" %in% names(logged@commands))
  expect_s4_class(logged@commands$test_command, "SeuratCommand")
  expect_equal(methods::slot(logged@commands$test_command, "name"), "test_command")
  expect_equal(methods::slot(logged@commands$test_command, "assay.used"), "RNA")
  expect_identical(dim(logged), dim(object))
})

test_that("Seurat layer helpers validate assays and combine split layers safely", {
  skip_if_not_installed("Seurat")

  object <- make_split_layer_object()
  matched_layers <- Shennong:::.sn_match_seurat_layers(
    object = object,
    assay = "RNA",
    layer = "counts"
  )
  combined <- Shennong:::.sn_get_seurat_layer_data(
    object = object,
    assay = "RNA",
    layer = "counts"
  )

  expect_true(isTRUE(Shennong:::.sn_validate_seurat_assay_layer(object, assay = "RNA", layer = "counts")))
  expect_setequal(matched_layers, SeuratObject::Layers(object[["RNA"]]))
  expect_s4_class(combined, "dgCMatrix")
  expect_identical(dimnames(combined), list(rownames(object), colnames(object)))
  expect_false("counts" %in% SeuratObject::Layers(object[["RNA"]]))

  expect_error(
    Shennong:::.sn_validate_seurat_assay_layer(matrix(1, nrow = 2, ncol = 2)),
    "Input must be a Seurat object"
  )
  expect_error(
    Shennong:::.sn_validate_seurat_assay_layer(object, assay = "ADT", layer = "counts"),
    "Assay 'ADT' was not found"
  )
  expect_error(
    Shennong:::.sn_validate_seurat_assay_layer(object, assay = "RNA", layer = "data"),
    "Layer 'data' was not found"
  )
})

test_that("temporary analysis counts are restored after using a non-default layer", {
  skip_if_not_installed("Seurat")

  object <- make_utils_test_object()
  original_counts <- SeuratObject::LayerData(object, assay = "RNA", layer = "counts")
  SeuratObject::LayerData(object, assay = "RNA", layer = "counts.alt") <- original_counts * 10

  prepared <- Shennong:::.sn_prepare_seurat_analysis_input(
    object = object,
    assay = "RNA",
    layer = "counts.alt"
  )
  restored <- Shennong:::.sn_restore_seurat_analysis_input(
    object = prepared$object,
    context = prepared$context
  )

  expect_true(isTRUE(prepared$context$needs_temp_counts))
  expect_equal(
    as.matrix(SeuratObject::LayerData(prepared$object, assay = "RNA", layer = "counts")),
    as.matrix(original_counts * 10)
  )
  expect_equal(
    as.matrix(SeuratObject::LayerData(restored, assay = "RNA", layer = "counts")),
    as.matrix(original_counts)
  )
  expect_true("counts.alt" %in% SeuratObject::Layers(restored[["RNA"]]))
})

test_that("temporary combined counts for split layers are removed after restoration", {
  skip_if_not_installed("Seurat")

  object <- make_split_layer_object()
  original_layers <- SeuratObject::Layers(object[["RNA"]])
  combined <- Shennong:::.sn_get_seurat_layer_data(object, assay = "RNA", layer = "counts")

  prepared <- Shennong:::.sn_prepare_seurat_analysis_input(
    object = object,
    assay = "RNA",
    layer = "counts"
  )
  restored <- Shennong:::.sn_restore_seurat_analysis_input(
    object = prepared$object,
    context = prepared$context
  )

  expect_true("counts" %in% SeuratObject::Layers(prepared$object[["RNA"]]))
  expect_equal(
    as.matrix(SeuratObject::LayerData(prepared$object, assay = "RNA", layer = "counts")),
    as.matrix(combined)
  )
  expect_setequal(SeuratObject::Layers(restored[["RNA"]]), original_layers)
  expect_false("counts" %in% SeuratObject::Layers(restored[["RNA"]]))
})

test_that("layer alias helpers restore pre-existing target layers", {
  skip_if_not_installed("Seurat")

  object <- make_utils_test_object()
  original_counts <- SeuratObject::LayerData(object, assay = "RNA", layer = "counts")
  SeuratObject::LayerData(object, assay = "RNA", layer = "data") <- original_counts
  SeuratObject::LayerData(object, assay = "RNA", layer = "data.alt") <- original_counts * 2

  aliased <- Shennong:::.sn_prepare_seurat_layer_alias(
    object = object,
    assay = "RNA",
    source_layer = "data.alt"
  )
  restored <- Shennong:::.sn_restore_seurat_layer_alias(
    object = aliased$object,
    context = aliased$context
  )

  expect_equal(aliased$target_layer, "data")
  expect_equal(Shennong:::.sn_guess_seurat_target_layer("counts.alt"), "counts")
  expect_equal(Shennong:::.sn_guess_seurat_target_layer("scale.data.batch1"), "scale.data")
  expect_equal(Shennong:::.sn_guess_seurat_target_layer("data.batch1"), "data")
  expect_equal(
    as.matrix(SeuratObject::LayerData(aliased$object, assay = "RNA", layer = "data")),
    as.matrix(original_counts * 2)
  )
  expect_equal(
    as.matrix(SeuratObject::LayerData(restored, assay = "RNA", layer = "data")),
    as.matrix(original_counts)
  )
  expect_true("data.alt" %in% SeuratObject::Layers(restored[["RNA"]]))
})
