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

test_that("sn_read and sn_write preserve row names and matrix-like tabular inputs", {
  skip_if_not_installed("rio")

  input <- data.frame(id = c("row1", "row2"), value = c(1, 2), stringsAsFactors = FALSE)
  csv_path <- tempfile(fileext = ".csv")
  utils::write.csv(input, csv_path, row.names = FALSE)

  read_back <- sn_read(csv_path, row_names = 1)
  expect_equal(rownames(read_back), c("row1", "row2"))
  expect_equal(read_back$value, c(1, 2))

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
