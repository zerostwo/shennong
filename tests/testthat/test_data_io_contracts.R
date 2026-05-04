library(testthat)

test_that("io helper functions detect urls and infer formats from files and directories", {
  tenx_dir <- tempfile("tenx-")
  dir.create(tenx_dir)
  file.create(file.path(tenx_dir, c("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz")))

  starsolo_dir <- tempfile("starsolo-")
  dir.create(starsolo_dir)
  file.create(file.path(starsolo_dir, c("matrix.mtx", "barcodes.tsv", "features.tsv")))

  bpcells_dir <- tempfile("bpcells-")
  dir.create(bpcells_dir)
  file.create(file.path(bpcells_dir, "index.val"))

  spatial_dir <- tempfile("spatial-")
  dir.create(spatial_dir)
  file.create(file.path(spatial_dir, "filtered_feature_bc_matrix.h5"))

  unknown_dir <- tempfile("unknown-")
  dir.create(unknown_dir)
  file.create(file.path(unknown_dir, "notes.txt"))

  gz_path <- tempfile(fileext = ".csv.gz")
  file.create(gz_path)

  expect_true(Shennong:::.sn_is_url("https://example.org/data.tsv?download=1"))
  expect_false(Shennong:::.sn_is_url("/tmp/data.tsv"))
  expect_true(Shennong:::.sn_is_custom_format("gmt"))
  expect_true(Shennong:::.sn_is_custom_format("gtf"))
  expect_true(Shennong:::.sn_is_custom_format("10x_spatial"))
  expect_false(Shennong:::.sn_is_custom_format("csv"))
  expect_equal(Shennong:::.sn_standardize_format("CSV"), "csv")
  expect_equal(Shennong:::.sn_standardize_format(","), ",")
  expect_equal(Shennong:::.sn_infer_format("clipboard"), "clipboard")
  expect_equal(Shennong:::.sn_infer_format("https://example.org/data.tsv?download=1"), "tsv")
  expect_equal(Shennong:::.sn_infer_format(gz_path), "csv")
  expect_equal(Shennong:::.guess_dir_format(tenx_dir), "10x")
  expect_equal(Shennong:::.guess_dir_format(starsolo_dir), "starsolo")
  expect_equal(Shennong:::.guess_dir_format(bpcells_dir), "bpcells")
  expect_equal(Shennong:::.guess_dir_format(spatial_dir), "10x_spatial")
  expect_null(Shennong:::.guess_dir_format(unknown_dir))
})

test_that("sn_list_10x_paths finds detected outs and selected 10x assets", {
  root <- tempfile("inputs-")
  dir.create(root)

  outs_dir <- file.path(root, "sample1", "outs")
  dir.create(file.path(outs_dir, "filtered_feature_bc_matrix"), recursive = TRUE)
  dir.create(file.path(outs_dir, "raw_feature_bc_matrix"), recursive = TRUE)
  file.create(file.path(outs_dir, "metrics_summary.csv"))

  h5_outs_dir <- file.path(root, "sample2", "outs")
  dir.create(h5_outs_dir, recursive = TRUE)
  file.create(file.path(h5_outs_dir, "filtered_feature_bc_matrix.h5"))
  file.create(file.path(h5_outs_dir, "raw_feature_bc_matrix.h5"))

  outs_paths <- sn_list_10x_paths(root)

  expect_equal(
    sort(normalizePath(unname(outs_paths), winslash = "/", mustWork = TRUE)),
    sort(normalizePath(c(outs_dir, h5_outs_dir), winslash = "/", mustWork = TRUE))
  )
  expect_equal(sort(names(outs_paths)), c("sample1", "sample2"))

  expect_equal(
    sort(normalizePath(
      unname(sn_list_10x_paths(root, path_type = "filtered")),
      winslash = "/",
      mustWork = TRUE
    )),
    sort(normalizePath(
      c(
        file.path(outs_dir, "filtered_feature_bc_matrix"),
        file.path(h5_outs_dir, "filtered_feature_bc_matrix.h5")
      ),
      winslash = "/",
      mustWork = TRUE
    ))
  )

  expect_equal(
    sort(normalizePath(
      unname(sn_list_10x_paths(root, path_type = "raw")),
      winslash = "/",
      mustWork = TRUE
    )),
    sort(normalizePath(
      c(
        file.path(outs_dir, "raw_feature_bc_matrix"),
        file.path(h5_outs_dir, "raw_feature_bc_matrix.h5")
      ),
      winslash = "/",
      mustWork = TRUE
    ))
  )

  expect_equal(
    sort(normalizePath(
      unname(sn_list_10x_paths(root, recursive = FALSE, include_root = FALSE)),
      winslash = "/",
      mustWork = TRUE
    )),
    sort(normalizePath(c(outs_dir, h5_outs_dir), winslash = "/", mustWork = TRUE))
  )

  expect_equal(
    normalizePath(
      unname(sn_list_10x_paths(root, path_type = "metrics")),
      winslash = "/",
      mustWork = TRUE
    ),
    normalizePath(file.path(outs_dir, "metrics_summary.csv"), winslash = "/", mustWork = TRUE)
  )
})

test_that("sn_list_10x_paths deduplicates root and nested outs candidates", {
  root <- tempfile("single-sample-root-")
  dir.create(file.path(root, "outs", "filtered_feature_bc_matrix"), recursive = TRUE)

  paths <- sn_list_10x_paths(root, recursive = TRUE, include_root = TRUE)

  expect_length(paths, 1)
  expect_equal(
    normalizePath(unname(paths), winslash = "/", mustWork = TRUE),
    normalizePath(file.path(root, "outs"), winslash = "/", mustWork = TRUE)
  )
})

test_that("10x path helpers select h5 assets correctly", {
  h5_root <- tempfile("tenx-h5-select-")
  dir.create(file.path(h5_root, "outs"), recursive = TRUE)
  filtered_h5 <- file.path(h5_root, "outs", "filtered_feature_bc_matrix.h5")
  raw_h5 <- file.path(h5_root, "outs", "raw_feature_bc_matrix.h5")
  file.create(filtered_h5)
  file.create(raw_h5)

  info <- list(
    outs_path = file.path(h5_root, "outs"),
    filtered_path = filtered_h5,
    raw_path = raw_h5,
    metrics_path = NULL
  )

  expect_equal(Shennong:::.sn_select_10x_path(info, "outs"), file.path(h5_root, "outs"))
  expect_equal(Shennong:::.sn_select_10x_path(info, "metrics"), NULL)
  expect_equal(Shennong:::.sn_select_10x_path(info, "filtered_h5"), filtered_h5)
  expect_equal(Shennong:::.sn_select_10x_path(info, "raw_h5"), raw_h5)
})

test_that("10x candidate discovery prefers outs directories during recursive scans", {
  root <- tempfile("inputs-candidates-")
  dir.create(root)

  outs_dir <- file.path(root, "sampleA", "outs")
  dir.create(file.path(outs_dir, "filtered_feature_bc_matrix"), recursive = TRUE)

  nested_noise <- file.path(root, "misc", "deep", "tree")
  dir.create(nested_noise, recursive = TRUE)

  candidates <- Shennong:::.sn_find_10x_candidates(root, recursive = TRUE, include_root = TRUE)

  expect_true(normalizePath(root, winslash = "/", mustWork = TRUE) %in% normalizePath(candidates, winslash = "/", mustWork = TRUE))
  expect_true(normalizePath(outs_dir, winslash = "/", mustWork = TRUE) %in% normalizePath(candidates, winslash = "/", mustWork = TRUE))
  expect_false(normalizePath(nested_noise, winslash = "/", mustWork = TRUE) %in% normalizePath(candidates, winslash = "/", mustWork = TRUE))
})

test_that("GMT import and custom dispatchers follow the expected IO contracts", {
  gmt_path <- tempfile(fileext = ".gmt")
  writeLines(
    c(
      "TERM_A\tdescription\tGENE1\tGENE2",
      "TERM_B\tdescription\tGENE3"
    ),
    gmt_path
  )

  imported <- Shennong:::.import.rio_gmt(gmt_path)

  expect_equal(colnames(imported), c("term", "gene"))
  expect_equal(as.character(imported$term), c("TERM_A", "TERM_A", "TERM_B"))
  expect_equal(as.character(imported$gene), c("GENE1", "GENE2", "GENE3"))
  expect_error(
    Shennong:::.sn_dispatch_custom_reader(gmt_path, format = "unsupported"),
    "Unsupported import format"
  )
  expect_error(
    Shennong:::.sn_dispatch_custom_writer(
      path = tempfile(fileext = ".foo"),
      x = data.frame(a = 1),
      format = "unsupported"
    ),
    "Unsupported export format"
  )
})

test_that("sn_read validates paths and delegates custom readers for inferred formats", {
  skip_if_not_installed("rio")

  expect_error(
    sn_read(file.path(tempdir(), "definitely-missing.csv")),
    "No such path"
  )

  tenx_dir <- tempfile("tenx-read-")
  dir.create(tenx_dir)
  file.create(file.path(tenx_dir, c("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz")))
  spatial_dir <- tempfile("spatial-read-")
  dir.create(spatial_dir)
  file.create(file.path(spatial_dir, "filtered_feature_bc_matrix.h5"))

  local_mocked_bindings(
    .import.rio_10x = function(file, ...) {
      list(kind = "tenx", path = file)
    },
    .import.rio_10x_spatial = function(file, ...) {
      list(kind = "spatial", path = file)
    },
    .env = asNamespace("Shennong")
  )

  imported <- sn_read(tenx_dir)
  expect_equal(imported$kind, "tenx")
  expect_equal(normalizePath(imported$path, winslash = "/", mustWork = TRUE), normalizePath(tenx_dir, winslash = "/", mustWork = TRUE))
  imported_spatial <- sn_read(spatial_dir)
  expect_equal(imported_spatial$kind, "spatial")
  expect_equal(normalizePath(imported_spatial$path, winslash = "/", mustWork = TRUE), normalizePath(spatial_dir, winslash = "/", mustWork = TRUE))
})

test_that("sn_read and sn_write preserve row names and matrix-like tabular inputs", {
  skip_if_not_installed("rio")

  input <- data.frame(id = c("row1", "row2"), value = c(1, 2), stringsAsFactors = FALSE)
  csv_path <- tempfile(fileext = ".csv")
  utils::write.csv(input, csv_path, row.names = FALSE)

  read_back <- sn_read(csv_path, row_names = 1)
  expect_equal(rownames(read_back), c("row1", "row2"))
  expect_equal(read_back$value, c(1, 2))
  read_back_by_name <- sn_read(csv_path, row_names = "id")
  expect_equal(rownames(read_back_by_name), c("row1", "row2"))
  expect_equal(read_back_by_name$value, c(1, 2))
  expect_error(
    sn_read(csv_path, row_names = "missing_id"),
    "`row_names` must identify a column"
  )

  matrix_path <- tempfile(fileext = ".csv")
  matrix_input <- matrix(c(1, 2, 3, 4), nrow = 2, dimnames = list(c("g1", "g2"), c("c1", "c2")))
  sn_write(matrix_input, matrix_path)
  matrix_output <- sn_read(matrix_path)

  expect_s3_class(matrix_output, "data.frame")
  expect_equal(unname(as.matrix(matrix_output)), unname(matrix_input))
  expect_equal(colnames(matrix_output), colnames(matrix_input))

  list_path <- tempfile(fileext = ".csv")
  sn_write(list(data.frame(a = 1:2, b = c("x", "y"))), list_path)
  list_output <- sn_read(list_path)

  expect_equal(list_output$a, 1:2)
  expect_equal(list_output$b, c("x", "y"))
  gmt_path <- tempfile(fileext = ".gmt")
  writeLines("TERM_A\tdescription\tGENE1\tGENE2", gmt_path)
  gmt_read <- sn_read(gmt_path, format = "gmt")
  expect_s3_class(gmt_read, "data.frame")
  expect_error(
    sn_write(list(a = 1, b = 2), tempfile(fileext = ".csv")),
    "'x' is not a data.frame or matrix"
  )
  expect_error(
    sn_write(data.frame(a = 1)),
    "Must specify 'path' and/or 'to'"
  )
})

test_that("sn_add_data_from_anndata respects custom reduction keys and logs the import", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(matrix(rpois(40, lambda = 3), nrow = 10, ncol = 4), sparse = TRUE)
  rownames(counts) <- paste0("gene", seq_len(10))
  colnames(counts) <- paste0("cell", seq_len(4))
  object <- sn_initialize_seurat_object(x = counts, project = "ann-key")

  umap <- data.frame(
    UMAP_1 = c(1, 2, 3, 4),
    UMAP_2 = c(4, 3, 2, 1),
    row.names = colnames(object)
  )
  umap_path <- tempfile(fileext = ".csv")
  utils::write.csv(umap, umap_path)

  updated <- sn_add_data_from_anndata(
    object = object,
    umap_path = umap_path,
    key = "latent"
  )

  expect_true("latent" %in% names(updated@reductions))
  expect_equal(updated[["latent"]]@key, "latent_")
  expect_true("sn_add_data_from_anndata" %in% names(updated@commands))
})

test_that("sn_add_data_from_anndata can add metadata without embeddings", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(matrix(rpois(40, lambda = 3), nrow = 10, ncol = 4), sparse = TRUE)
  rownames(counts) <- paste0("gene", seq_len(10))
  colnames(counts) <- paste0("cell", seq_len(4))
  object <- sn_initialize_seurat_object(x = counts, project = "ann-meta")

  metadata <- data.frame(
    sample = c("A", "A", "B", "B"),
    score = c(1, 2, 3, 4),
    row.names = colnames(object)
  )
  metadata_path <- tempfile(fileext = ".csv")
  utils::write.csv(metadata, metadata_path)

  updated <- sn_add_data_from_anndata(
    object = object,
    metadata_path = metadata_path
  )

  expect_equal(as.character(updated$sample), metadata$sample)
  expect_equal(as.numeric(updated$score), metadata$score)
})

test_that("sn_write supports explicit format selection and custom writer dispatch", {
  skip_if_not_installed("rio")

  tmpdir <- tempfile("sn-write-")
  dir.create(tmpdir)
  oldwd <- setwd(tmpdir)
  on.exit(setwd(oldwd), add = TRUE)

  tbl <- data.frame(a = 1:2, b = c("x", "y"))
  auto_path <- sn_write(tbl, to = "csv")

  expect_equal(auto_path, "tbl.csv")
  expect_true(file.exists(file.path(tmpdir, "tbl.csv")))

  skip_if_not_installed("BPCells")

  bpcells_path <- file.path(tmpdir, "matrix-dir")
  counts <- methods::as(
    Matrix::Matrix(matrix(c(1L, 0L, 2L, 3L), nrow = 2), sparse = TRUE),
    "dgCMatrix"
  )
  sn_write(counts, path = bpcells_path, to = "bpcells", overwrite = TRUE)

  expect_true(dir.exists(bpcells_path))
  expect_true(any(grepl("val$", list.files(bpcells_path))))
})

test_that("io helpers cover custom source passthrough and non-csv AnnData imports", {
  skip_if_not_installed("Seurat")
  skip_if_not_installed("rio")

  local_path <- tempfile(fileext = ".h5ad")
  file.create(local_path)
  expect_equal(Shennong:::.sn_download_custom_source(local_path, format = "h5ad"), local_path)

  counts <- Matrix::Matrix(matrix(rpois(40, lambda = 3), nrow = 10, ncol = 4), sparse = TRUE)
  rownames(counts) <- paste0("gene", seq_len(10))
  colnames(counts) <- paste0("cell", seq_len(4))
  object <- sn_initialize_seurat_object(x = counts, project = "ann-tsv")

  metadata_path <- tempfile(fileext = ".tsv")
  writeLines(
    c(
      "cell\tsample\tscore",
      "cell1\tA\t1",
      "cell2\tA\t2",
      "cell3\tB\t3",
      "cell4\tB\t4"
    ),
    metadata_path
  )

  updated <- sn_add_data_from_anndata(
    object = object,
    metadata_path = metadata_path
  )

  expect_equal(as.character(updated$sample), c("A", "A", "B", "B"))
  expect_equal(as.numeric(updated$score), c(1, 2, 3, 4))
})

test_that("sn_read and sn_write support qs and qs2 serialized objects", {
  skip_if_not_installed("rio")

  expect_true(Shennong:::.sn_is_custom_format("qs"))
  expect_true(Shennong:::.sn_is_custom_format("qs2"))

  if (requireNamespace("qs", quietly = TRUE)) {
    qs_path <- tempfile(fileext = ".qs")
    payload <- list(a = 1:3, b = data.frame(label = c("x", "y"), value = c(2, 4)))
    sn_write(payload, qs_path)
    restored <- sn_read(qs_path)
    expect_equal(restored, payload)
  }

  if (requireNamespace("qs2", quietly = TRUE)) {
    qs2_path <- tempfile(fileext = ".qs2")
    payload2 <- list(alpha = letters[1:4], beta = matrix(1:6, nrow = 2))
    sn_write(payload2, qs2_path)
    restored2 <- sn_read(qs2_path)
    expect_equal(restored2, payload2)
  } else {
    expect_error(
      Shennong:::.import.rio_qs2(tempfile(fileext = ".qs2")),
      "qs2"
    )
  }
})

test_that("h5ad writer accepts existing SingleCellExperiment objects", {
  skip_if_not_installed("anndataR")
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("rhdf5")

  sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = matrix(1:6, nrow = 2))
  )
  h5ad_path <- tempfile(fileext = ".h5ad")

  sn_write(sce, h5ad_path)

  expect_true(file.exists(h5ad_path))
})
