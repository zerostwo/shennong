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

  available_methods <- sn_list_methods(available = TRUE)
  unavailable_methods <- sn_list_methods(available = FALSE)
  expect_true(all(available_methods$available))
  expect_true(all(!unavailable_methods$available))
  expect_equal(
    nrow(available_methods) + nrow(unavailable_methods),
    nrow(methods)
  )

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
    schema_version = "1.0.0",
    analysis_type = "trajectory",
    name = "demo",
    method = "slingshot",
    backend = "slingshot",
    input = list(reduction = "pca", cells = 4L),
    parameters = list(start = "0"),
    tables = list(primary = tibble::tibble(
      cell = paste0("cell", 1:4),
      primary_pseudotime = 1:4
    )),
    embeddings = list(),
    graphs = list(),
    models = list(),
    diagnostics = list(),
    warnings = character(),
    provenance = list(package_versions = list(Shennong = "0.1.4"), random_seed = 1L, timestamp = "2026-07-15 UTC")
  )

  expect_true(sn_validate_result(result, error = FALSE)$valid)
  legacy_version <- result
  legacy_version$schema_version <- "1.0"
  legacy_report <- sn_validate_result(legacy_version, error = FALSE)
  expect_true(legacy_report$valid)
  expect_match(paste(legacy_report$warnings, collapse = " "), "upgraded")
  malformed_primary <- result
  malformed_primary$tables$primary <- tibble::tibble(cell = paste0("cell", 1:4))
  expect_false(sn_validate_result(malformed_primary, error = FALSE)$valid)
  invalid <- result
  invalid$provenance <- NULL
  report <- sn_validate_result(invalid, error = FALSE)
  expect_false(report$valid)
  expect_match(paste(report$errors, collapse = " "), "provenance")
  expect_error(sn_validate_result(invalid), "Invalid Shennong analysis result")

  malformed_versions <- result
  malformed_versions$provenance$package_versions <- "bad"
  expect_false(sn_validate_result(malformed_versions, error = FALSE)$valid)
  malformed_seed <- result
  malformed_seed$provenance$random_seed <- "bad"
  expect_false(sn_validate_result(malformed_seed, error = FALSE)$valid)
  malformed_timestamp <- result
  malformed_timestamp$provenance$timestamp <- c("first", "second")
  expect_false(sn_validate_result(malformed_timestamp, error = FALSE)$valid)
})

test_that("generic results round-trip, list, filter, and delete", {
  object <- make_analysis_result_test_object()
  compact <- list(
    method = "slingshot",
    backend = "slingshot",
    input = list(cells = ncol(object)),
    parameters = list(start = "0"),
    tables = list(primary = tibble::tibble(
      cell = colnames(object),
      primary_pseudotime = seq_len(ncol(object))
    )),
    diagnostics = list(converged = TRUE),
    warnings = character(),
    provenance = list(random_seed = 1L)
  )

  object <- sn_store_result(object, "trajectory", "cd8", compact)
  stored <- sn_get_result(object, "trajectory", "cd8")
  expect_equal(stored$analysis_type, "trajectory")
  expect_equal(stored$name, "cd8")
  expect_equal(stored$tables$primary$primary_pseudotime, seq_len(ncol(object)))
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

test_that("stored-result audit identifies and upgrades legacy envelopes", {
  object <- make_analysis_result_test_object()
  legacy <- list(
    schema_version = "1.0",
    package_version = "0.1.4",
    created_at = "2026-07-15 UTC",
    table = tibble::tibble(gene = "gene1", score = 1),
    analysis = "markers",
    method = "seurat"
  )
  object@misc$de_results <- list(legacy_markers = legacy)
  object@misc$input_source <- list(path = "/data/example", format = "10x")

  audit <- sn_audit_results(object)
  legacy_audit <- audit[audit$name == "legacy_markers", , drop = FALSE]
  artifact_audit <- audit[audit$contract_scope == "artifact", , drop = FALSE]
  expect_equal(legacy_audit$status, "legacy")
  expect_false(legacy_audit$unified)
  expect_true(legacy_audit$repairable)
  expect_equal(legacy_audit$primary_rows, 1L)
  expect_equal(artifact_audit$type, "input_source")
  expect_equal(artifact_audit$status, "artifact")
  expect_true(is.na(artifact_audit$unified))

  upgraded <- sn_upgrade_results(object)
  upgraded_audit <- sn_audit_results(upgraded)
  upgraded_result_audit <- upgraded_audit[
    upgraded_audit$name == "legacy_markers",
    ,
    drop = FALSE
  ]
  expect_equal(upgraded_result_audit$status, "valid")
  expect_true(upgraded_result_audit$unified)
  stored <- sn_get_result(upgraded, "de", "legacy_markers")
  expect_equal(stored$schema_version, "1.0.0")
  expect_identical(stored$table, stored$tables$primary)
})

test_that("audit never labels an unpreparable result as unified", {
  object <- make_analysis_result_test_object()
  malformed <- list(
    schema_version = "1.0.0", analysis_type = "trajectory", name = "malformed",
    method = "slingshot", backend = "slingshot", input = list(), parameters = list(),
    tables = list(primary = tibble::tibble(cell = "cell1", primary_pseudotime = 1)),
    embeddings = list(), graphs = list(), models = list(), diagnostics = list(),
    warnings = character(),
    provenance = list(package_versions = "bad", random_seed = "bad", timestamp = c("a", "b"))
  )
  object@misc$analysis_results <- list(trajectory = list(malformed = malformed))

  audit <- sn_audit_results(object, type = "trajectory")
  expect_identical(audit$status, "invalid")
  expect_false(audit$unified)
  expect_false(audit$repairable)
  expect_match(paste(audit$errors, audit$upgrade_errors), "provenance")
})

test_that("audit covers canonical label transfer, artifacts, and unknown misc", {
  object <- make_analysis_result_test_object()
  object$transfer_label <- c("T cell", "T cell", "B cell", "B cell")
  object$transfer_score <- c(0.9, 0.8, 0.7, 0.6)
  object@misc$label_transfer <- list(
    transfer = list(
      method = "seurat",
      prediction_columns = c("transfer_label", "transfer_score")
    )
  )
  object@misc$qc <- list(retained_cells = ncol(object))
  object@misc$user_payload <- list(note = "belongs to the user")
  object <- Shennong:::.sn_store_label_transfer_result(
    object = object,
    prediction_prefix = "transfer",
    method = "seurat",
    label_by = "cell_type",
    prediction_columns = c("transfer_label", "transfer_score")
  )

  audit <- sn_audit_results(object)
  canonical <- audit[
    audit$contract_scope == "analysis_result" & audit$type == "annotation",
    ,
    drop = FALSE
  ]
  expect_equal(canonical$name, "transfer")
  expect_true(canonical$unified)
  expect_equal(canonical$primary_rows, ncol(object))
  expect_true(any(audit$collection == "label_transfer" & audit$status == "artifact"))
  expect_true(any(audit$collection == "qc" & audit$status == "artifact"))
  expect_true(any(
    audit$collection == "user_payload" &
      audit$contract_scope == "unregistered" &
      audit$status == "unregistered" &
      is.na(audit$unified)
  ))

  artifacts_before <- object@misc[c("label_transfer", "qc", "user_payload")]
  upgraded <- sn_upgrade_results(object)
  expect_identical(
    upgraded@misc[c("label_transfer", "qc", "user_payload")],
    artifacts_before
  )
})

test_that("legacy storage helpers return the unified envelope additively", {
  object <- make_analysis_result_test_object()
  table <- tibble::tibble(feature = "gene1", score = 1)

  returned <- list(
    sn_store_deconvolution(object, table, return_object = FALSE),
    sn_store_regulatory_activity(object, table, return_object = FALSE),
    sn_store_milo(
      object,
      table,
      sample_by = "sample",
      group_by = "condition",
      return_object = FALSE
    ),
    sn_store_enrichment(object, table, return_object = FALSE)
  )
  expect_true(all(vapply(returned, function(result) {
    sn_validate_result(result, error = FALSE)$valid &&
      identical(result$table, result$tables$primary)
  }, logical(1))))
})

test_that("future schemas and partial field names are never treated as legacy", {
  object <- make_analysis_result_test_object()
  future <- list(
    schema_version = "2.0.0",
    analysis_type = "trajectory",
    name = "future",
    method = "future_backend",
    backend = "future_backend",
    input = list(cells = ncol(object)),
    parameters = list(),
    tables = list(primary = tibble::tibble(
      cell = colnames(object),
      primary_pseudotime = seq_len(ncol(object))
    )),
    embeddings = list(),
    graphs = list(),
    models = list(),
    diagnostics = list(),
    warnings = character(),
    provenance = list(
      package_versions = list(Shennong = "future"),
      random_seed = NA_integer_,
      timestamp = "future"
    )
  )
  object@misc$analysis_results <- list(trajectory = list(future = future))

  audit <- sn_audit_results(object, type = "trajectory")
  expect_identical(audit$status, "invalid")
  expect_false(audit$repairable)
  expect_error(
    sn_upgrade_results(object, type = "trajectory"),
    "unsupported `schema_version` '2.0.0'"
  )
  expect_warning(
    unchanged <- sn_upgrade_results(object, type = "trajectory", strict = FALSE),
    "unsupported `schema_version` '2.0.0'"
  )
  expect_identical(
    unchanged@misc$analysis_results$trajectory$future$schema_version,
    "2.0.0"
  )
  expect_error(
    sn_store_result(object, "trajectory", "future", future),
    "unsupported `schema_version` '2.0.0'"
  )

  partial <- future
  partial$schema_version <- "1.0.0"
  partial$tables_backup <- partial$tables
  partial$tables <- NULL
  expect_false(sn_validate_result(partial, error = FALSE)$valid)
  expect_error(
    sn_store_result(object, "trajectory", "partial", partial),
    "tables\\$primary"
  )
})

test_that("generic storage materializes canonical QC compatibility aliases", {
  object <- make_analysis_result_test_object()
  by_sample <- tibble::tibble(
    sample = c("sample_1", "sample_2"),
    retained_cells = c(10L, 12L)
  )
  overall <- tibble::tibble(metric = "retained_cells", value = 22L)
  canonical <- list(
    schema_version = "1.0.0",
    analysis_type = "qc_assessment",
    name = "canonical_qc",
    method = "qc",
    backend = "qc",
    input = list(),
    parameters = list(),
    tables = list(primary = by_sample, overall = overall),
    embeddings = list(),
    graphs = list(),
    models = list(),
    diagnostics = list(),
    warnings = character(),
    provenance = list(
      package_versions = list(Shennong = "test"),
      random_seed = NA_integer_,
      timestamp = "test"
    )
  )

  object <- sn_store_result(object, "qc_assessment", "canonical_qc", canonical)
  stored <- sn_get_result(object, "qc_assessment", "canonical_qc")

  expect_identical(stored$tables$primary, by_sample)
  expect_identical(stored$tables$overall, overall)
  expect_identical(stored$by_sample, by_sample)
  expect_identical(stored$overall, overall)
  expect_true(sn_validate_result(stored, error = FALSE)$valid)
})
