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

test_that("sn_deconvolve_bulk prepares local CIBERSORTx commands from a Seurat reference", {
  skip_if_not_installed("Seurat")

  object <- make_deconvolution_object()
  bulk <- make_bulk_matrix(object)
  outdir <- tempfile("cibersortx-")

  bundle <- sn_deconvolve_bulk(
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
})

test_that("sn_deconvolve_bulk can import a CIBERSORTx fractions table and store it", {
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

  object <- sn_deconvolve_bulk(
    object,
    bulk = bulk,
    method = "cibersortx",
    cell_type_by = "cell_type",
    layer = "counts",
    outdir = tempdir(),
    prefix = "import",
    cibersortx_result = result_path,
    store_name = "cibersortx_demo",
    return_object = TRUE
  )

  expect_true("deconvolution_results" %in% names(object@misc))
  imported <- sn_get_deconvolution_result(object, "cibersortx_demo")
  expect_true(all(c("sample", "cell_type", "fraction") %in% colnames(imported)))
  expect_equal(nrow(imported), 6)
  listed <- sn_list_results(object)
  expect_true("deconvolution" %in% listed$type)
})

test_that("sn_deconvolve_bulk validates BayesPrism availability cleanly", {
  skip_if_not_installed("Seurat")

  object <- make_deconvolution_object()
  bulk <- make_bulk_matrix(object)

  expect_error(
    testthat::with_mocked_bindings(
      sn_deconvolve_bulk(
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

test_that("sn_deconvolve_bulk validates local CIBERSORTx credentials", {
  skip_if_not_installed("Seurat")

  object <- make_deconvolution_object()
  bulk <- make_bulk_matrix(object)

  expect_error(
    sn_deconvolve_bulk(
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
    store_name = "manual"
  )

  filtered <- sn_get_deconvolution_result(
    object,
    deconvolution_name = "manual",
    samples = "sample_b",
    cell_types = "Bcell"
  )
  expect_equal(nrow(filtered), 1)
  expect_equal(filtered$fraction[[1]], 0.7)
})
