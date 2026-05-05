library(testthat)

make_public_index <- function() {
  list(
    schema_version = "1.0.0",
    project = list(project_id = "shennong-public-001"),
    zenodo = list(record_id = "20044788", doi = "10.5281/zenodo.20044788"),
    studies = list(
      list(
        study_id = "StudyA",
        display_name = "Study A",
        zip_file = "StudyA.zip",
        samples = list(
          list(
            sample_id = "Sample1",
            display_name = "Sample 1",
            organism = "Homo sapiens",
            taxon_id = 9606,
            assay = "single_cell_transcriptomics",
            technology = "10x_genomics",
            source_accessions = list(geo = list("GSE1"), bioproject = list("PRJ1")),
            processing = list(
              pipeline = "cellranger",
              pipeline_version = "9.0.1",
              reference = "GENCODE:48:GRCh38.p14"
            ),
            files = list(
              filtered_h5 = "StudyA/Sample1/filtered_feature_bc_matrix.h5",
              raw_h5 = "StudyA/Sample1/raw_feature_bc_matrix.h5",
              metrics_summary = "StudyA/Sample1/metrics_summary.csv"
            )
          ),
          list(
            sample_id = "Sample2",
            display_name = "Sample 2",
            organism = "Homo sapiens",
            taxon_id = 9606,
            assay = "single_cell_transcriptomics",
            technology = "cite_seq",
            source_accessions = list(geo = list("GSE1")),
            processing = list(
              pipeline = "cellranger",
              pipeline_version = "9.0.1",
              reference = "GENCODE:48:GRCh38.p14"
            ),
            files = list(
              filtered_h5 = "StudyA/Sample2/filtered_feature_bc_matrix.h5",
              raw_h5 = "StudyA/Sample2/raw_feature_bc_matrix.h5",
              metrics_summary = "StudyA/Sample2/metrics_summary.csv"
            )
          )
        )
      )
    )
  )
}

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

test_that("sn_list_datasets exposes public collection samples from the Zenodo index", {
  catalog <- with_mocked_bindings(
    sn_list_datasets(source = "public", save_dir = tempfile("shennong-cache-")),
    .sn_read_public_data_index = function(...) make_public_index(),
    .package = "Shennong"
  )

  expect_s3_class(catalog, "data.frame")
  expect_equal(catalog$dataset, c("Sample1", "Sample2"))
  expect_equal(catalog$study_id, c("StudyA", "StudyA"))
  expect_equal(catalog$zenodo_record, c("20044788", "20044788"))
  expect_equal(catalog$filtered_h5[[1]], "StudyA/Sample1/filtered_feature_bc_matrix.h5")
  expect_equal(catalog$geo_accession[[1]], "GSE1")

  combined <- with_mocked_bindings(
    sn_list_datasets(source = "all", save_dir = tempfile("shennong-cache-")),
    .sn_read_public_data_index = function(...) make_public_index(),
    .package = "Shennong"
  )
  expect_true(all(c("Sample1", "pbmc3k") %in% combined$dataset))
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
    with_mocked_bindings(
      sn_load_data(dataset = "not-a-dataset", save_dir = save_dir, return_object = FALSE),
      .sn_read_public_data_index = function(...) make_public_index(),
      .package = "Shennong"
    ),
    "Invalid value|known sample IDs"
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

test_that("sn_load_data resolves public collection samples by sample or study/sample", {
  save_dir <- tempfile("shennong-cache-")
  dir.create(save_dir)
  local_h5 <- file.path(save_dir, "public-filtered.h5")
  file.create(local_h5)

  counts <- Matrix::Matrix(
    matrix(1:12, nrow = 3, dimnames = list(paste0("gene", 1:3), paste0("cell", 1:4))),
    sparse = TRUE
  )
  captured <- list()
  sentinel <- list(kind = "public-seurat")

  returned <- with_mocked_bindings(
    sn_load_data(
      dataset = "Sample1",
      save_dir = save_dir,
      return_object = TRUE
    ),
    .sn_read_public_data_index = function(...) make_public_index(),
    .sn_prepare_public_sample_file = function(row, matrix_type, ...) {
      captured$prepared <<- list(row = row, matrix_type = matrix_type)
      local_h5
    },
    sn_read = function(path, ...) {
      expect_equal(path, local_h5)
      counts
    },
    sn_initialize_seurat_object = function(x, species, sample_name, project, ...) {
      captured$initialized <<- list(x = x, species = species, sample_name = sample_name, project = project)
      sentinel
    },
    .package = "Shennong"
  )

  expect_identical(returned, sentinel)
  expect_equal(captured$prepared$row$sample_id, "Sample1")
  expect_equal(captured$prepared$matrix_type, "filtered")
  expect_equal(captured$initialized$species, "human")
  expect_equal(captured$initialized$sample_name, "Sample1")
  expect_equal(captured$initialized$project, "StudyA")

  path <- with_mocked_bindings(
    sn_load_data(
      dataset = "StudyA",
      sample_id = "Sample2",
      matrix_type = "raw",
      save_dir = save_dir,
      return_object = FALSE
    ),
    .sn_read_public_data_index = function(...) make_public_index(),
    .sn_prepare_public_sample_file = function(row, matrix_type, ...) {
      expect_equal(row$sample_id, "Sample2")
      expect_equal(matrix_type, "raw")
      local_h5
    },
    .package = "Shennong"
  )
  expect_equal(path, local_h5)
})

test_that("sn_load_data returns public metrics tables", {
  save_dir <- tempfile("shennong-cache-")
  dir.create(save_dir)
  metrics_path <- file.path(save_dir, "metrics_summary.csv")
  write.csv(data.frame(Metric = "Estimated Number of Cells", Value = 10), metrics_path, row.names = FALSE)

  metrics <- with_mocked_bindings(
    sn_load_data(
      dataset = "StudyA",
      sample_id = "Sample1",
      matrix_type = "metrics",
      save_dir = save_dir
    ),
    .sn_read_public_data_index = function(...) make_public_index(),
    .sn_prepare_public_sample_file = function(row, matrix_type, ...) {
      expect_equal(matrix_type, "metrics")
      metrics_path
    },
    .package = "Shennong"
  )

  expect_s3_class(metrics, "data.frame")
  expect_equal(metrics$Metric[[1]], "Estimated Number of Cells")
})

test_that("sn_load_data asks for sample_id when a public study has multiple samples", {
  expect_error(
    with_mocked_bindings(
      sn_load_data(dataset = "StudyA", return_object = FALSE),
      .sn_read_public_data_index = function(...) make_public_index(),
      .package = "Shennong"
    ),
    "`sample_id` must be supplied"
  )
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
    "cannot be mixed"
  )
  expect_error(
    sn_load_data(dataset = c("pbmc1k", "pbmc3k"), species = c("human", "mouse", "human"), return_object = FALSE),
    "same length"
  )
})
