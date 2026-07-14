library(testthat)

make_analysis_result_test_object <- function() {
  skip_if_not_installed("SeuratObject")
  counts <- Matrix::Matrix(
    matrix(c(1, 0, 2, 1, 0, 3, 2, 1), nrow = 2),
    sparse = TRUE
  )
  rownames(counts) <- c("gene1", "gene2")
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  SeuratObject::CreateSeuratObject(counts = counts, project = "analysis-result-test")
}

test_that("method registry is unique, queryable, and reports roadmap backends", {
  methods <- sn_list_methods()

  expect_s3_class(methods, "data.frame")
  expect_true(all(c(
    "task", "name", "runtime", "implemented", "default", "available",
    "input_requirements", "outputs", "install_action", "citation"
  ) %in% colnames(methods)))
  expect_equal(anyDuplicated(paste(methods$task, methods$name, sep = "::")), 0L)
  expect_true(any(methods$task == "trajectory" & methods$name == "slingshot" & methods$default))
  expect_true(any(methods$task == "annotation" & methods$name == "singleR"))
  expect_true(any(methods$task == "bulk_de" & methods$name == "auto"))

  status <- sn_method_status("slingshot", task = "trajectory")
  expect_equal(status$runtime, "r")
  expect_equal(status$package, "slingshot")
  expect_true(status$implemented)
  expect_match(status$install_action, "slingshot", ignore.case = TRUE)

  expect_error(sn_method_status("not-a-method"), "Unknown Shennong method")
  expect_error(sn_list_methods(available = "yes"), "available")
})

test_that("analysis result contract validates required structure", {
  result <- list(
    schema_version = "1.0",
    analysis_type = "trajectory",
    name = "demo",
    method = "slingshot",
    backend = "slingshot",
    input = list(reduction = "pca", cells = 4L),
    parameters = list(start = "0"),
    tables = list(primary = tibble::tibble(cell = paste0("cell", 1:4), pseudotime = 1:4)),
    embeddings = list(),
    graphs = list(),
    models = list(),
    diagnostics = list(),
    warnings = character(),
    provenance = list(package_versions = list(Shennong = "0.1.4"), random_seed = 1L, timestamp = "2026-07-15 UTC")
  )

  expect_true(sn_validate_result(result, error = FALSE)$valid)
  invalid <- result
  invalid$provenance <- NULL
  report <- sn_validate_result(invalid, error = FALSE)
  expect_false(report$valid)
  expect_match(paste(report$errors, collapse = " "), "provenance")
  expect_error(sn_validate_result(invalid), "Invalid Shennong analysis result")
})

test_that("generic results round-trip, list, filter, and delete", {
  object <- make_analysis_result_test_object()
  compact <- list(
    method = "slingshot",
    backend = "slingshot",
    input = list(cells = ncol(object)),
    parameters = list(start = "0"),
    tables = list(primary = tibble::tibble(cell = colnames(object), pseudotime = seq_len(ncol(object)))),
    diagnostics = list(converged = TRUE),
    warnings = character(),
    provenance = list(random_seed = 1L)
  )

  object <- sn_store_result(object, "trajectory", "cd8", compact)
  stored <- sn_get_result(object, "trajectory", "cd8")
  expect_equal(stored$analysis_type, "trajectory")
  expect_equal(stored$name, "cd8")
  expect_equal(stored$tables$primary$pseudotime, seq_len(ncol(object)))
  expect_true(sn_validate_result(stored, error = FALSE)$valid)

  listing <- sn_list_results(object, type = "trajectory")
  expect_equal(listing$name, "cd8")
  expect_equal(listing$type, "trajectory")
  expect_equal(listing$n_rows, ncol(object))
  expect_equal(nrow(sn_list_results(object, type = "de")), 0L)

  object <- sn_delete_result(object, "trajectory", "cd8")
  expect_error(sn_get_result(object, "trajectory", "cd8"), "No stored result")
})

test_that("legacy stored results are upgraded to the unified contract", {
  object <- make_analysis_result_test_object()
  legacy <- list(
    schema_version = "1.0.0",
    package_version = "0.1.4",
    created_at = "2026-07-15 UTC",
    table = tibble::tibble(gene = "gene1", score = 1),
    analysis = "markers",
    method = "seurat"
  )
  object <- Shennong:::.sn_store_misc_result(object, "de_results", "markers", legacy)

  generic <- sn_get_result(object, "de", "markers")
  expect_true(sn_validate_result(generic, error = FALSE)$valid)
  expect_equal(generic$tables$primary, legacy$table)
  expect_equal(sn_get_de_result(object, "markers"), legacy$table)
})
