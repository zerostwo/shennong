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
    "must be one of"
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
    .sn_download_example_file = function(url, destfile) {
      downloaded <<- list(url = url, destfile = destfile)
      file.create(destfile)
      invisible(destfile)
    },
    .package = "Shennong"
  )

  expect_true(file.exists(path))
  expect_match(downloaded$url, "pbmc3k_filtered_feature_bc_matrix\\.h5$")
  expect_equal(downloaded$destfile, path)
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
