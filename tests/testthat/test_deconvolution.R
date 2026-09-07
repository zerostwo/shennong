library(testthat)

make_deconvolution_object <- function() {
  set.seed(123)
  counts <- matrix(rpois(30 * 18, lambda = 3), nrow = 30, ncol = 18)
  rownames(counts) <- paste0("gene", seq_len(30))
  colnames(counts) <- paste0("cell", seq_len(18))
  object <- SeuratObject::CreateSeuratObject(
    counts = Matrix::Matrix(counts, sparse = TRUE),
    project = "deconv"
  )
  object <- Seurat::AddMetaData(
    object = object,
    metadata = data.frame(
      cell_type = rep(c("Tcell", "Bcell", "Mono"), each = 6),
      cell_state = rep(c("T1", "T2", "B1", "B2", "M1", "M2"), each = 3),
      row.names = colnames(object)
    )
  )
  object
}

make_bulk_matrix <- function(object) {
  counts <- Shennong:::.sn_get_seurat_layer_data(object, assay = "RNA", layer = "counts")
  cbind(
    sample_a = Matrix::rowSums(counts[, 1:9, drop = FALSE]),
    sample_b = Matrix::rowSums(counts[, 10:18, drop = FALSE])
  )
}

test_that("sn_run_bulk_deconvolution prepares local CIBERSORTx commands from a Seurat reference", {
  skip_if_not_installed("Seurat")

  object <- make_deconvolution_object()
  bulk <- make_bulk_matrix(object)
  outdir <- tempfile("cibersortx-")

  bundle <- sn_run_bulk_deconvolution(
    object,
    bulk = bulk,
    method = "cibersortx",
    cell_type_by = "cell_type",
    layer = "counts",
    outdir = outdir,
    prefix = "demo",
    cibersortx_email = "demo@example.org",
    cibersortx_token = "fake-token",
    cibersortx_dry_run = TRUE,
    return_object = FALSE
  )

  expect_identical(bundle$method, "cibersortx")
  expect_true(file.exists(bundle$files$single_cell_reference))
  expect_true(file.exists(bundle$files$mixture))
  expect_match(bundle$artifacts$commands$create_signature, "cibersortx/fractions|CIBERSORTxFractions")
  expect_match(bundle$artifacts$commands$deconvolve, "--label")
  expect_false(grepl("fake-token", paste(bundle$artifacts$commands, collapse = " "), fixed = TRUE))
  expect_false(grepl("demo@example.org", paste(bundle$artifacts$commands, collapse = " "), fixed = TRUE))
  expect_true(bundle$artifacts$commands_redacted)
  expect_false(bundle$artifacts$host_argv_contains_credentials)
  expect_identical(bundle$artifacts$credential_transport, "ephemeral_read_only_bind_file")
  expect_true(bundle$artifacts$run_dir_retained)
  expect_true(Shennong:::.sn_is_owned_run_dir(bundle$artifacts$run_dir))
  expect_true(startsWith(
    dirname(bundle$files$single_cell_reference),
    bundle$artifacts$run_dir
  ))
  expect_identical(bundle$scale_provenance$reference, "counts")
  expect_identical(bundle$scale_provenance$bulk, "counts")
  expect_true(bundle$scale_provenance$matched)
  expect_identical(bundle$scale_provenance$reference_layer, "counts")
  expect_length(
    list.files(tempdir(), pattern = "^shennong-cibersortx-credentials-"),
    0L
  )
})

test_that("CIBERSORTx host commands never contain credential values", {
  io_dir <- tempfile("cibersortx-command-")
  dir.create(io_dir)
  credentials <- Shennong:::.sn_write_cibersortx_credentials(
    "private@example.org",
    "top-secret-token"
  )
  on.exit(unlink(credentials, force = TRUE), add = TRUE)
  command <- Shennong:::.sn_create_cibersortx_command(
    input_dir = io_dir,
    output_dir = io_dir,
    credentials_path = credentials,
    container = "docker",
    method = "create_sig"
  )
  host_args <- paste(command$args, collapse = " ")
  expect_false(grepl("private@example.org", host_args, fixed = TRUE))
  expect_false(grepl("top-secret-token", host_args, fixed = TRUE))
  expect_match(host_args, "/run/secrets/shennong_cibersortx")
})

test_that("deconvolution refuses unbounded dense expression materialization", {
  object <- make_deconvolution_object()
  expect_error(
    sn_run_bulk_deconvolution(
      object,
      bulk = make_bulk_matrix(object),
      method = "cibersortx",
      cell_type_by = "cell_type",
      max_dense_gb = 1e-12,
      cibersortx_email = "demo@example.org",
      cibersortx_token = "fake-token",
      cibersortx_dry_run = TRUE,
      return_object = FALSE
    ),
    "exceeding `max_dense_gb"
  )
})

test_that("CIBERSORTx matrix export streams sparse input without character coercion", {
  output_dir <- withr::local_tempdir()
  matrix <- Matrix::sparseMatrix(
    i = c(1L, 2L, 3L), j = c(1L, 2L, 1L), x = c(1, 2, 3),
    dims = c(3L, 2L),
    dimnames = list(c("G1", "G2", "G3"), c("C1", "C2"))
  )
  path <- Shennong:::.sn_transform_and_save_cibersortx_single_cell(
    matrix,
    cell_type_labels = c("T", "B"),
    path = output_dir,
    max_dense_gb = 0.001
  )
  lines <- readLines(path)
  expect_identical(lines[[1L]], "GeneSymbol\tT\tB")
  expect_length(lines, 4L)
  expect_match(lines[[2L]], "^G1\\t1\\t0$")

  expect_error(
    Shennong:::.sn_transform_and_save_cibersortx_single_cell(
      matrix,
      cell_type_labels = c("T", "B"),
      path = output_dir,
      max_dense_gb = 1e-12
    ),
    "exceed `max_dense_gb"
  )
})

test_that("sn_run_bulk_deconvolution can import a CIBERSORTx fractions table and store it", {
  skip_if_not_installed("Seurat")

  object <- make_deconvolution_object()
  bulk <- make_bulk_matrix(object)
  result_path <- tempfile(fileext = ".tsv")
  fractions <- data.frame(
    Mixture = c("sample_a", "sample_b"),
    Tcell = c(0.4, 0.2),
    Bcell = c(0.3, 0.5),
    Mono = c(0.3, 0.3),
    `P-value` = c(0.01, 0.01),
    Correlation = c(0.9, 0.85),
    RMSE = c(0.05, 0.08),
    check.names = FALSE
  )
  utils::write.table(fractions, file = result_path, sep = "\t", quote = FALSE, row.names = FALSE)

  object <- sn_run_bulk_deconvolution(
    object,
    bulk = bulk,
    method = "cibersortx",
    cell_type_by = "cell_type",
    layer = "counts",
    outdir = tempdir(),
    prefix = "import",
    cibersortx_result = result_path,
    result_id = "cibersortx_demo",
    return_object = TRUE
  )

  expect_true("cibersortx_demo" %in% names(object@misc$shennong$results$deconvolution))
  imported <- sn_get_deconvolution_result(object, "cibersortx_demo")
  expect_true(all(c("sample", "cell_type", "fraction") %in% colnames(imported)))
  expect_equal(nrow(imported), 6)
  listed <- sn_list_results(object)
  expect_true("deconvolution" %in% listed$type)
})

test_that("sn_run_bulk_deconvolution validates BayesPrism availability cleanly", {
  skip_if_not_installed("Seurat")

  object <- make_deconvolution_object()
  bulk <- make_bulk_matrix(object)

  expect_error(
    testthat::with_mocked_bindings(
      sn_run_bulk_deconvolution(
        object,
        bulk = bulk,
        method = "bayesprism",
        cell_type_by = "cell_type",
        cell_state_by = "cell_state",
        layer = "counts",
        return_object = FALSE
      ),
      check_installed_github = function(pkg, repo, reason = NULL) {
        stop(
          paste0(
            "Package '", pkg, "' is required ",
            reason,
            "\nInstall it with:\n  remotes::install_github('", repo, "')"
          ),
          call. = FALSE
        )
      },
      .package = "Shennong"
    ),
    "BayesPrism"
  )
})

test_that("sn_run_bulk_deconvolution validates local CIBERSORTx credentials", {
  skip_if_not_installed("Seurat")

  object <- make_deconvolution_object()
  bulk <- make_bulk_matrix(object)

  expect_error(
    sn_run_bulk_deconvolution(
      object,
      bulk = bulk,
      method = "cibersortx",
      cell_type_by = "cell_type",
      layer = "counts",
      outdir = tempdir(),
      prefix = "local",
      cibersortx_dry_run = TRUE,
      return_object = FALSE
    ),
    "CIBERSORTx credentials"
  )
})

test_that("BayesPrism rejects non-count matrices, invalid genes, and empty labels at its boundary", {
  reference <- matrix(
    c(1, 0, 2, 3, 1, 4),
    nrow = 2,
    dimnames = list(c("cell1", "cell2"), c("G1", "G2", "G3"))
  )
  bulk <- matrix(
    c(4, 2, 1, 3, 5, 2),
    nrow = 2,
    dimnames = list(c("sample1", "sample2"), c("G1", "G2", "G3"))
  )
  run_invalid <- function(reference_input = reference,
                          bulk_input = bulk,
                          type_labels = c("T", "B"),
                          state_labels = c("T1", "B1")) {
    Shennong:::.sn_run_bayesprism(
      reference_cells_by_gene = reference_input,
      cell_type_labels = type_labels,
      cell_state_labels = state_labels,
      bulk_samples_by_gene = bulk_input
    )
  }

  fractional <- reference
  fractional[1L, 1L] <- fractional[1L, 1L] + 0.25
  expect_error(run_invalid(reference_input = fractional), "integer-like raw counts")
  negative <- reference
  negative[1L, 1L] <- -1
  expect_error(run_invalid(reference_input = negative), "finite, non-negative")
  missing_bulk <- bulk
  missing_bulk[1L, 1L] <- NA_real_
  expect_error(run_invalid(bulk_input = missing_bulk), "finite, non-negative")
  duplicate_reference <- reference
  colnames(duplicate_reference)[[2L]] <- colnames(duplicate_reference)[[1L]]
  expect_error(run_invalid(reference_input = duplicate_reference), "gene identifiers must be unique")
  duplicate_bulk <- bulk
  colnames(duplicate_bulk)[[2L]] <- colnames(duplicate_bulk)[[1L]]
  expect_error(run_invalid(bulk_input = duplicate_bulk), "gene identifiers must be unique")
  expect_error(run_invalid(type_labels = c("T", " ")), "non-missing, non-empty")
  expect_error(run_invalid(state_labels = c("T1", NA_character_)), "non-missing, non-empty")

  object <- make_deconvolution_object()
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  expect_error(
    sn_run_bulk_deconvolution(
      object,
      bulk = make_bulk_matrix(object),
      method = "bayesprism",
      cell_type_by = "cell_type",
      layer = "data",
      return_object = FALSE
    ),
    "raw/count-like reference layer"
  )
})

test_that("CIBERSORTx accepts matched linear input and rejects unsafe or mixed scales", {
  reference <- matrix(
    c(1, 0, 2, 3, 1, 4),
    nrow = 3,
    dimnames = list(c("G1", "G2", "G3"), c("cell1", "cell2"))
  )
  bulk <- matrix(
    c(4, 2, 1, 3, 5, 2),
    nrow = 3,
    dimnames = list(c("G1", "G2", "G3"), c("sample1", "sample2"))
  )
  reference_linear <- sweep(reference + 0.25, 2L, colSums(reference + 0.25), "/")
  bulk_linear <- sweep(bulk + 0.25, 2L, colSums(bulk + 0.25), "/")
  parent <- withr::local_tempdir(pattern = "cibersortx-linear-")
  linear <- Shennong:::.sn_run_cibersortx_local(
    reference_genes_by_cells = reference_linear,
    cell_type_labels = c("T", "B"),
    bulk_genes_by_samples = bulk_linear,
    outdir = parent,
    prefix = "linear",
    email = "demo@example.org",
    token = "fake-token",
    dry_run = TRUE
  )
  expect_identical(linear$scale_provenance$reference, "linear_nonlog")
  expect_identical(linear$scale_provenance$bulk, "linear_nonlog")
  expect_true(linear$scale_provenance$matched)

  expect_error(
    Shennong:::.sn_run_cibersortx_local(
      reference, c("T", "B"), bulk_linear,
      prefix = "mixed", email = "demo@example.org", token = "fake-token", dry_run = TRUE
    ),
    "same expression scale"
  )
  negative <- reference
  negative[1L, 1L] <- -1
  expect_error(
    Shennong:::.sn_run_cibersortx_local(
      negative, c("T", "B"), bulk,
      prefix = "negative", email = "demo@example.org", token = "fake-token", dry_run = TRUE
    ),
    "finite, non-negative"
  )
  missing_bulk <- bulk
  missing_bulk[1L, 1L] <- NA_real_
  expect_error(
    Shennong:::.sn_run_cibersortx_local(
      reference, c("T", "B"), missing_bulk,
      prefix = "missing", email = "demo@example.org", token = "fake-token", dry_run = TRUE
    ),
    "finite, non-negative"
  )
  duplicate <- reference
  rownames(duplicate)[[2L]] <- rownames(duplicate)[[1L]]
  expect_error(
    Shennong:::.sn_run_cibersortx_local(
      duplicate, c("T", "B"), bulk,
      prefix = "duplicate", email = "demo@example.org", token = "fake-token", dry_run = TRUE
    ),
    "gene identifiers must be unique"
  )

  object <- Seurat::NormalizeData(make_deconvolution_object(), verbose = FALSE)
  expect_error(
    sn_run_bulk_deconvolution(
      object,
      bulk = make_bulk_matrix(object),
      method = "cibersortx",
      cell_type_by = "cell_type",
      layer = "data",
      cibersortx_email = "demo@example.org",
      cibersortx_token = "fake-token",
      cibersortx_dry_run = TRUE,
      return_object = FALSE
    ),
    "raw counts or non-log linear expression"
  )
})

test_that("CIBERSORTx uses owned unique children and cleans default successful runs", {
  object <- make_deconvolution_object()
  bulk <- make_bulk_matrix(object)
  parent <- withr::local_tempdir(pattern = "cibersortx-parent-")
  sentinel <- file.path(parent, "user-sentinel.txt")
  writeLines("keep", sentinel)

  retained <- sn_run_bulk_deconvolution(
    object,
    bulk = bulk,
    method = "cibersortx",
    cell_type_by = "cell_type",
    outdir = parent,
    prefix = "retained",
    cibersortx_email = "demo@example.org",
    cibersortx_token = "fake-token",
    cibersortx_dry_run = TRUE,
    return_object = FALSE
  )
  retained_again <- sn_run_bulk_deconvolution(
    object,
    bulk = bulk,
    method = "cibersortx",
    cell_type_by = "cell_type",
    outdir = parent,
    prefix = "retained2",
    cibersortx_email = "demo@example.org",
    cibersortx_token = "fake-token",
    cibersortx_dry_run = TRUE,
    return_object = FALSE
  )
  run_dir <- retained$artifacts$run_dir
  second_run_dir <- retained_again$artifacts$run_dir
  expect_false(identical(run_dir, second_run_dir))
  expect_true(startsWith(run_dir, normalizePath(parent)))
  expect_true(Shennong:::.sn_is_owned_run_dir(run_dir))
  expect_identical(readLines(sentinel), "keep")

  captured_run <- NULL
  create_owned <- Shennong:::.sn_create_owned_run_dir
  cleaned <- testthat::with_mocked_bindings(
    sn_run_bulk_deconvolution(
      object,
      bulk = bulk,
      method = "cibersortx",
      cell_type_by = "cell_type",
      prefix = "temporary",
      cibersortx_email = "demo@example.org",
      cibersortx_token = "fake-token",
      cibersortx_dry_run = TRUE,
      return_object = FALSE
    ),
    .sn_create_owned_run_dir = function(parent = NULL, prefix = "shennong-run-") {
      captured_run <<- create_owned(parent = parent, prefix = prefix)
      captured_run
    },
    .package = "Shennong"
  )
  expect_false(cleaned$artifacts$run_dir_retained)
  expect_length(cleaned$files, 0L)
  expect_false(dir.exists(captured_run))
  expect_false(file.exists(paste0(captured_run, ".shennong-owned-run")))
  expect_identical(readLines(sentinel), "keep")
})

test_that("CIBERSORTx sanitizes failed runs and import-only mode never exports", {
  object <- make_deconvolution_object()
  bulk <- make_bulk_matrix(object)
  parent <- withr::local_tempdir(pattern = "cibersortx-failure-")
  sentinel <- file.path(parent, "user-sentinel.txt")
  writeLines("keep", sentinel)

  expect_error(
    testthat::with_mocked_bindings(
      sn_run_bulk_deconvolution(
        object,
        bulk = bulk,
        method = "cibersortx",
        cell_type_by = "cell_type",
        outdir = parent,
        prefix = "failure",
        cibersortx_email = "demo@example.org",
        cibersortx_token = "fake-token",
        return_object = FALSE
      ),
      .sn_check_container_binary = function(...) "/bin/true",
      .sn_run_command = function(...) list(status = 1L, command = "redacted"),
      .package = "Shennong"
    ),
    "signature-matrix creation failed.*Sanitized diagnostics"
  )
  failed_runs <- list.dirs(parent, recursive = FALSE, full.names = TRUE)
  expect_length(failed_runs, 1L)
  expect_true(file.exists(file.path(failed_runs, "failure.json")))
  expect_false(dir.exists(file.path(failed_runs, "input")))
  expect_identical(readLines(sentinel), "keep")

  result_path <- file.path(parent, "completed.tsv")
  utils::write.table(
    data.frame(Mixture = "sample1", T = 0.6, B = 0.4, check.names = FALSE),
    result_path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  imported <- testthat::with_mocked_bindings(
    Shennong:::.sn_run_cibersortx_local(
      reference_genes_by_cells = matrix(NA_real_, 1L, 1L),
      cell_type_labels = NA_character_,
      bulk_genes_by_samples = matrix(NA_real_, 1L, 1L),
      outdir = parent,
      prefix = "import",
      result_path = result_path
    ),
    .sn_transform_and_save_cibersortx_single_cell = function(...) stop("export called"),
    .sn_transform_and_save_cibersortx_bulk = function(...) stop("export called"),
    .package = "Shennong"
  )
  expect_true(imported$artifacts$import_only)
  expect_identical(imported$scale_provenance$contract, "import_only")
  expect_identical(names(imported$files), "result")
})

test_that("sn_store_deconvolution stores and subsets stored results", {
  skip_if_not_installed("Seurat")

  object <- make_deconvolution_object()
  tbl <- tibble::tibble(
    sample = c("sample_a", "sample_a", "sample_b", "sample_b"),
    cell_type = c("Tcell", "Bcell", "Tcell", "Bcell"),
    fraction = c(0.6, 0.4, 0.3, 0.7)
  )

  object <- sn_store_deconvolution(
    object,
    result = tbl,
    result_id = "manual"
  )

  filtered <- sn_get_deconvolution_result(
    object,
    result_id = "manual",
    samples = "sample_b",
    cell_types = "Bcell"
  )
  expect_equal(nrow(filtered), 1)
  expect_equal(filtered$fraction[[1]], 0.7)
})
