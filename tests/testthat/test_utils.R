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

  getFromNamespace("merge.Seurat", ns = "SeuratObject")(
    x = object1,
    y = object2,
    add.cell.ids = c("splita", "splitb")
  )
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

test_that("future globals max size is resolved from memory and existing options", {
  old_override <- getOption("shennong.future.globals.maxSize", NULL)
  on.exit(options(shennong.future.globals.maxSize = old_override), add = TRUE)
  options(shennong.future.globals.maxSize = NULL)

  gib <- 1024^3

  expect_equal(
    Shennong:::.sn_resolve_future_globals_max_size(
      current = NULL,
      memory_limit = 64 * gib,
      min_size = 8 * gib,
      cap = 128 * gib
    ),
    8 * gib
  )

  expect_equal(
    Shennong:::.sn_resolve_future_globals_max_size(
      current = NULL,
      memory_limit = 4 * gib,
      min_size = 8 * gib,
      memory_fraction = 0.75,
      cap = 128 * gib
    ),
    ceiling(3 * gib)
  )

  expect_equal(
    Shennong:::.sn_resolve_future_globals_max_size(
      current = 16 * gib,
      memory_limit = 64 * gib,
      min_size = 8 * gib,
      cap = 128 * gib
    ),
    16 * gib
  )

  options(shennong.future.globals.maxSize = 12 * gib)
  expect_equal(
    Shennong:::.sn_resolve_future_globals_max_size(
      current = NULL,
      memory_limit = 64 * gib,
      min_size = 8 * gib,
      cap = 128 * gib
    ),
    12 * gib
  )
})

test_that("future globals helper restores the caller option", {
  old_future <- getOption("future.globals.maxSize", NULL)
  old_override <- getOption("shennong.future.globals.maxSize", NULL)
  on.exit({
    options(future.globals.maxSize = old_future)
    options(shennong.future.globals.maxSize = old_override)
  }, add = TRUE)

  options(future.globals.maxSize = 100)
  options(shennong.future.globals.maxSize = 1024)
  observed <- Shennong:::.sn_with_auto_future_globals(
    getOption("future.globals.maxSize"),
    context = "test future operation",
    verbose = FALSE
  )

  expect_equal(observed, 1024)
  expect_equal(getOption("future.globals.maxSize"), 100)

  options(future.globals.maxSize = NULL)
  options(shennong.future.globals.maxSize = 1024^3)
  observed <- Shennong:::.sn_with_auto_future_globals(
    getOption("future.globals.maxSize"),
    context = "test future operation",
    verbose = FALSE
  )

  expect_equal(observed, 1024^3)
  expect_null(getOption("future.globals.maxSize", NULL))
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

  explicit <- Shennong:::.sn_log_seurat_command(
    object = object,
    assay = "RNA",
    name = "explicit_command",
    call_string = quote(sn_example(object, method = "resolved")),
    params = list(method = "resolved", hidden_default = 717L)
  )
  command <- explicit@commands$explicit_command
  expect_identical(methods::slot(command, "call.string"), "sn_example(object, method = \"resolved\")")
  expect_identical(
    methods::slot(command, "params"),
    list(method = "resolved", hidden_default = 717L)
  )
})

test_that("logged commands with nested list parameters remain printable", {
  skip_if_not_installed("Seurat")

  object <- make_utils_test_object()
  logged <- Shennong:::.sn_log_seurat_command(
    object = object,
    assay = "RNA",
    name = "nested_command",
    params = list(
      schema_version = "1.0.0",
      requested = list(method = "auto", assay = NULL),
      backend_args = list(maxIter = 7L),
      input = list(x = list(type = "Seurat", assay = "RNA")),
      hidden = function(x) x
    )
  )
  command <- logged@commands$nested_command

  expect_s4_class(command, "sn_seurat_command")
  expect_true(methods::is(command, "SeuratCommand"))
  params <- methods::slot(command, "params")
  expect_identical(params$schema_version, "1.0.0")
  expect_identical(params$requested$method, "auto")
  expect_null(params$requested$assay)
  expect_identical(params$backend_args, list(maxIter = 7L))
  expect_identical(params$input$x, list(type = "Seurat", assay = "RNA"))

  output <- expect_no_error(capture.output(print(command)))
  rendered <- paste(output, collapse = "\n")
  expect_match(rendered, "Command: ", fixed = TRUE)
  expect_match(rendered, "requested : list(method = \"auto\", assay = NULL)", fixed = TRUE)
  expect_match(rendered, "backend_args : list(maxIter = 7L)", fixed = TRUE)
  expect_match(rendered, "input : list(x = list(type = \"Seurat\", assay = \"RNA\"))", fixed = TRUE)
  expect_no_match(rendered, "hidden :", fixed = TRUE)
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

test_that("standard counts avoid unnecessary temporary layer snapshots", {
  skip_if_not_installed("Seurat")

  object <- make_utils_test_object()
  before <- SeuratObject::LayerData(object, assay = "RNA", layer = "counts")
  prepared <- Shennong:::.sn_prepare_seurat_analysis_input(
    object = object,
    assay = "RNA",
    layer = "counts"
  )
  restored <- Shennong:::.sn_restore_seurat_analysis_input(
    object = prepared$object,
    context = prepared$context
  )

  expect_false(prepared$context$needs_temp_counts)
  expect_null(prepared$context$original_counts)
  expect_length(prepared$context$original_analysis_layers, 0L)
  expect_identical(
    SeuratObject::LayerData(restored, assay = "RNA", layer = "counts"),
    before
  )
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

test_that("scoped seeds restore the caller RNG state", {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) original_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", original_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(90210)
  caller_seed <- get(".Random.seed", envir = .GlobalEnv)
  first <- Shennong:::.sn_with_seed(17, stats::runif(4))
  expect_identical(get(".Random.seed", envir = .GlobalEnv), caller_seed)
  second <- Shennong:::.sn_with_seed(17, stats::runif(4))
  expect_identical(first, second)
  expect_identical(get(".Random.seed", envir = .GlobalEnv), caller_seed)

  rm(".Random.seed", envir = .GlobalEnv)
  Shennong:::.sn_with_seed(17, stats::runif(1))
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
})

test_that("owned run cleanup never removes a user parent directory", {
  parent <- withr::local_tempdir()
  sentinel <- file.path(parent, "keep-me.txt")
  writeLines("user data", sentinel)

  run_dir <- Shennong:::.sn_create_owned_run_dir(parent, "backend-")
  writeLines("temporary input", file.path(run_dir, "matrix.csv"))
  expect_true(dir.exists(run_dir))
  expect_invisible(Shennong:::.sn_cleanup_owned_run_dir(run_dir))
  expect_false(dir.exists(run_dir))
  expect_true(file.exists(sentinel))
  expect_error(
    Shennong:::.sn_cleanup_owned_run_dir(parent),
    "without a Shennong ownership marker"
  )
})
