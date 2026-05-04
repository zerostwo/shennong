library(testthat)

test_that("sn_load_data returns cached file path when return_object is FALSE", {
  save_dir <- tempfile("shennong-cache-")
  dir.create(save_dir)
  local_h5 <- file.path(save_dir, "pbmc1k_filtered_feature_bc_matrix.h5")
  file.create(local_h5)

  expect_equal(
    sn_load_data(
      dataset = "pbmc1k",
      save_dir = save_dir,
      return_object = FALSE
    ),
    local_h5
  )
})

test_that("example-data registry exposes the expected bundled dataset metadata", {
  catalog <- Shennong:::.sn_example_data_catalog()

  expect_s3_class(catalog, "data.frame")
  expect_true(all(c("dataset", "zenodo_record", "species") %in% colnames(catalog)))
  expect_setequal(catalog$dataset, c("pbmc1k", "pbmc3k", "pbmc4k", "pbmc8k"))
  expect_true(all(catalog$species == "human"))
})

test_that("sn_load_data returns cached raw matrix paths and validates dataset names", {
  save_dir <- tempfile("shennong-cache-")
  dir.create(save_dir)
  local_h5 <- file.path(save_dir, "pbmc4k_raw_feature_bc_matrix.h5")
  file.create(local_h5)

  expect_equal(
    sn_load_data(
      dataset = "pbmc4k",
      matrix_type = "raw",
      save_dir = save_dir,
      return_object = FALSE
    ),
    local_h5
  )

  expect_error(
    sn_load_data(dataset = "not-a-dataset", save_dir = save_dir, return_object = FALSE),
    "must contain only known"
  )
})

test_that("sn_load_data downloads missing files into the cache and returns the path", {
  save_dir <- tempfile("shennong-cache-")
  dir.create(save_dir)

  downloaded <- NULL
  path <- with_mocked_bindings(
    sn_load_data(
      dataset = "pbmc3k",
      matrix_type = "filtered",
      save_dir = save_dir,
      return_object = FALSE
    ),
    sn_download_zenodo = function(record_id, files, save_dir, token = NULL, overwrite = FALSE, quiet = FALSE) {
      destfile <- file.path(save_dir, files)
      downloaded <<- list(
        record_id = record_id,
        files = files,
        save_dir = save_dir,
        token = token,
        overwrite = overwrite,
        quiet = quiet,
        destfile = destfile
      )
      file.create(destfile)
      stats::setNames(destfile, files)
    },
    .package = "Shennong"
  )

  expect_true(file.exists(path))
  expect_equal(downloaded$record_id, "14884845")
  expect_equal(downloaded$files, "pbmc3k_filtered_feature_bc_matrix.h5")
  expect_equal(downloaded$destfile, path)
  expect_null(downloaded$token)
  expect_false(downloaded$overwrite)
  expect_false(downloaded$quiet)
})

test_that("sn_load_data returns raw matrices through sn_read without Seurat initialization", {
  save_dir <- tempfile("shennong-cache-")
  dir.create(save_dir)

  counts <- Matrix::Matrix(
    matrix(1:12, nrow = 3, dimnames = list(paste0("gene", 1:3), paste0("cell", 1:4))),
    sparse = TRUE
  )
  local_h5 <- file.path(save_dir, "pbmc1k_raw_feature_bc_matrix.h5")
  file.create(local_h5)

  returned <- with_mocked_bindings(
    sn_load_data(
      dataset = "pbmc1k",
      matrix_type = "raw",
      save_dir = save_dir,
      return_object = TRUE
    ),
    sn_read = function(path, ...) {
      expect_equal(path, local_h5)
      counts
    },
    sn_initialize_seurat_object = function(...) {
      stop("filtered initialization should not be used for raw matrices")
    },
    .package = "Shennong"
  )

  expect_s4_class(returned, "dgCMatrix")
  expect_equal(as.matrix(returned), as.matrix(counts))
})

test_that("sn_load_data initializes filtered matrices with the resolved species", {
  save_dir <- tempfile("shennong-cache-")
  dir.create(save_dir)

  counts <- Matrix::Matrix(
    matrix(1:12, nrow = 3, dimnames = list(paste0("gene", 1:3), paste0("cell", 1:4))),
    sparse = TRUE
  )
  local_h5 <- file.path(save_dir, "pbmc8k_filtered_feature_bc_matrix.h5")
  file.create(local_h5)

  captured <- NULL
  sentinel <- list(kind = "seurat-object")
  returned <- with_mocked_bindings(
    sn_load_data(
      dataset = "pbmc8k",
      save_dir = save_dir,
      return_object = TRUE
    ),
    sn_read = function(path, ...) {
      expect_equal(path, local_h5)
      counts
    },
    sn_initialize_seurat_object = function(x, species, ...) {
      captured <<- list(x = x, species = species)
      sentinel
    },
    .package = "Shennong"
  )

  expect_identical(returned, sentinel)
  expect_equal(as.matrix(captured$x), as.matrix(counts))
  expect_equal(captured$species, "human")
})

test_that("sn_load_data returns named paths for multiple example datasets", {
  save_dir <- tempfile("shennong-cache-")
  dir.create(save_dir)

  paths <- with_mocked_bindings(
    sn_load_data(
      dataset = c("pbmc1k", "pbmc3k"),
      matrix_type = "filtered",
      save_dir = save_dir,
      return_object = FALSE,
      token = "restricted-token",
      overwrite = TRUE,
      quiet = TRUE
    ),
    sn_download_zenodo = function(record_id, files, save_dir, token = NULL, overwrite = FALSE, quiet = FALSE) {
      expect_equal(record_id, "14884845")
      expect_equal(token, "restricted-token")
      expect_true(overwrite)
      expect_true(quiet)
      destfile <- file.path(save_dir, files)
      file.create(destfile)
      stats::setNames(destfile, files)
    },
    .package = "Shennong"
  )

  expect_named(paths, c("pbmc1k", "pbmc3k"))
  expect_equal(
    unname(paths),
    file.path(
      save_dir,
      c("pbmc1k_filtered_feature_bc_matrix.h5", "pbmc3k_filtered_feature_bc_matrix.h5")
    )
  )
})

test_that("sn_load_data returns a named list for multiple raw datasets", {
  save_dir <- tempfile("shennong-cache-")
  dir.create(save_dir)
  file.create(file.path(save_dir, "pbmc1k_raw_feature_bc_matrix.h5"))
  file.create(file.path(save_dir, "pbmc3k_raw_feature_bc_matrix.h5"))

  returned <- with_mocked_bindings(
    sn_load_data(
      dataset = c("pbmc1k", "pbmc3k"),
      matrix_type = "raw",
      save_dir = save_dir,
      return_object = TRUE
    ),
    sn_read = function(path, ...) {
      Matrix::Matrix(
        matrix(seq_len(6), nrow = 2, dimnames = list(c("gene1", "gene2"), paste0(basename(path), "_cell", 1:3))),
        sparse = TRUE
      )
    },
    sn_initialize_seurat_object = function(...) {
      stop("filtered initialization should not be used for raw matrices")
    },
    .package = "Shennong"
  )

  expect_type(returned, "list")
  expect_named(returned, c("pbmc1k", "pbmc3k"))
  expect_s4_class(returned$pbmc1k, "dgCMatrix")
  expect_s4_class(returned$pbmc3k, "dgCMatrix")
})

test_that("sn_load_data merges multiple filtered example datasets", {
  save_dir <- tempfile("shennong-cache-")
  dir.create(save_dir)
  file.create(file.path(save_dir, "pbmc1k_filtered_feature_bc_matrix.h5"))
  file.create(file.path(save_dir, "pbmc3k_filtered_feature_bc_matrix.h5"))

  initialized <- list()
  sentinel <- list(kind = "merged")
  returned <- with_mocked_bindings(
    sn_load_data(
      dataset = c("pbmc1k", "pbmc3k"),
      save_dir = save_dir,
      return_object = TRUE,
      species = "human"
    ),
    sn_read = function(path, ...) basename(path),
    sn_initialize_seurat_object = function(x, species, sample_name, project, ...) {
      initialized[[sample_name]] <<- list(x = x, species = species, sample_name = sample_name, project = project)
      list(sample_name = sample_name)
    },
    .sn_merge_example_objects = function(objects) {
      expect_named(objects, c("pbmc1k", "pbmc3k"))
      sentinel
    },
    .package = "Shennong"
  )

  expect_identical(returned, sentinel)
  expect_named(initialized, c("pbmc1k", "pbmc3k"))
  expect_equal(initialized$pbmc1k$sample_name, "pbmc1k")
  expect_equal(initialized$pbmc3k$project, "pbmc3k")
})

test_that("sn_load_data validates vectorized dataset and species inputs", {
  expect_error(
    sn_load_data(dataset = c("pbmc1k", "pbmc1k"), return_object = FALSE),
    "duplicated"
  )
  expect_error(
    sn_load_data(dataset = c("pbmc1k", "not-a-dataset"), return_object = FALSE),
    "Invalid value"
  )
  expect_error(
    sn_load_data(dataset = c("pbmc1k", "pbmc3k"), species = c("human", "mouse", "human"), return_object = FALSE),
    "same length"
  )
})
