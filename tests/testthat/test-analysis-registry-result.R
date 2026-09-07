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

result_writer_names <- function() {
  c(
    "sn_annotate_de_features", "sn_assess_bulk_qc", "sn_assess_qc",
    "sn_compare_integrations", "sn_discover_programs", "sn_find_bulk_de",
    "sn_find_de", "sn_find_spatial_domains", "sn_find_spatial_features",
    "sn_integrate_spatial", "sn_interpret_annotation", "sn_interpret_de",
    "sn_interpret_enrichment", "sn_prioritize_states", "sn_run_annotation",
    "sn_run_bulk_deconvolution", "sn_run_cell_communication",
    "sn_run_clinical_association", "sn_run_cnv", "sn_run_enrichment",
    "sn_run_fate", "sn_run_grn", "sn_run_metabolism", "sn_run_milo",
    "sn_run_regulatory_activity", "sn_run_scissor",
    "sn_run_spatial_communication", "sn_run_spatial_neighborhood",
    "sn_run_survival", "sn_run_trajectory", "sn_run_velocity", "sn_run_wgcna",
    "sn_score_bulk_pathways", "sn_score_programs",
    "sn_store_cell_communication", "sn_store_deconvolution",
    "sn_store_enrichment", "sn_store_milo", "sn_store_regulatory_activity",
    "sn_store_result", "sn_test_abundance", "sn_test_programs",
    "sn_transfer_labels", "sn_write_figure_legend",
    "sn_write_presentation_summary", "sn_write_results"
  )
}

test_that("public result writers use result_id and no storage-name aliases", {
  writers <- result_writer_names()
  expect_setequal(writers, intersect(writers, getNamespaceExports("Shennong")))
  for (writer in writers) {
    parameters <- names(formals(getExportedValue("Shennong", writer)))
    expect_true("result_id" %in% parameters, info = writer)
    expect_false(any(c("name", "store_name") %in% parameters), info = writer)
  }
})

test_that("method registry remains unique and queryable", {
  methods <- sn_list_methods()
  expect_s3_class(methods, "data.frame")
  expect_equal(anyDuplicated(paste(methods$task, methods$name, sep = "::")), 0L)
  expect_true(any(methods$task == "trajectory" & methods$name == "slingshot"))
  expect_error(sn_get_method_status("not-a-method"), "Unknown Shennong method")
})

test_that("analysis-result v2 has one canonical identifier", {
  result <- list(
    schema_version = "2.0.0", analysis_type = "trajectory", result_id = "demo",
    method = "slingshot", backend = "slingshot",
    input = list(reduction = "pca", cells = 4L), parameters = list(start = "0"),
    tables = list(primary = tibble::tibble(
      cell = paste0("cell", 1:4), primary_pseudotime = 1:4
    )),
    embeddings = list(), graphs = list(), models = list(), diagnostics = list(),
    warnings = character(),
    provenance = list(
      package_versions = list(Shennong = "0.3.0.9000"), random_seed = 1L,
      timestamp = "2026-09-03 UTC", result_id = "demo",
      analysis_type = "trajectory"
    )
  )

  expect_true(sn_validate_result(result, error = FALSE)$valid)
  named <- result
  named$result_id <- NULL
  named$name <- "demo"
  expect_false(sn_validate_result(named, error = FALSE)$valid)
  old_schema <- result
  old_schema$schema_version <- "1.0.0"
  expect_false(sn_validate_result(old_schema, error = FALSE)$valid)

  malformed_de <- result
  malformed_de$analysis_type <- "de"
  malformed_de$provenance$analysis_type <- "de"
  malformed_de$tables$primary <- tibble::tibble(foo = 1)
  expect_false(sn_validate_result(malformed_de, error = FALSE)$valid)
})

test_that("schema-v2 random seeds fail closed without validator crashes", {
  result <- Shennong:::.sn_new_analysis_result(
    "trajectory", "seed_contract", "test", "test",
    tables = list(primary = tibble::tibble(
      cell = "cell1", primary_pseudotime = 1
    )),
    random_seed = 1L
  )
  invalid_seeds <- list(
    above_integer_range = as.double(.Machine$integer.max) + 1,
    enormous = 1e100,
    positive_infinity = Inf,
    negative_infinity = -Inf,
    not_a_number = NaN,
    negative = -1,
    fractional = 1.5,
    vector = c(1, 2),
    matrix = matrix(1, nrow = 1),
    character = "1",
    malformed_named_list = list(stage = as.double(.Machine$integer.max) + 1)
  )

  for (seed_name in names(invalid_seeds)) {
    malformed <- result
    malformed$provenance$random_seed <- invalid_seeds[[seed_name]]
    report <- NULL
    expect_no_error(
      report <- sn_validate_result(malformed, error = FALSE)
    )
    expect_false(report$valid, info = seed_name)
    expect_match(
      paste(report$errors, collapse = " "),
      "provenance\\$random_seed",
      info = seed_name
    )
  }

  for (seed in list(
    0L,
    as.double(.Machine$integer.max),
    NA_integer_,
    NA_real_,
    list(export = 1L, model = NA_integer_)
  )) {
    candidate <- result
    candidate$provenance$random_seed <- seed
    expect_true(sn_validate_result(candidate, error = FALSE)$valid)
  }
})

test_that("result audits continue after a malformed schema-v2 random seed", {
  object <- make_analysis_result_test_object()
  good <- Shennong:::.sn_new_analysis_result(
    "trajectory", "good", "test", "test",
    tables = list(primary = tibble::tibble(
      cell = colnames(object),
      primary_pseudotime = seq_len(ncol(object))
    )),
    random_seed = 1L
  )
  object <- sn_store_result(object, "trajectory", "good", good)
  bad <- good
  bad$result_id <- "bad"
  bad$provenance$result_id <- "bad"
  bad$provenance$random_seed <- as.double(.Machine$integer.max) + 1
  object@misc$shennong$results$trajectory <- list(bad = bad, good = good)

  audit <- NULL
  expect_no_error(
    audit <- sn_audit_results(object, include_artifacts = FALSE)
  )
  expect_setequal(audit$result_id, c("bad", "good"))
  expect_false(audit$valid[audit$result_id == "bad"])
  expect_identical(audit$status[audit$result_id == "bad"], "invalid")
  expect_match(
    audit$errors[audit$result_id == "bad"],
    "provenance\\$random_seed"
  )
  expect_true(audit$valid[audit$result_id == "good"])
  expect_identical(audit$status[audit$result_id == "good"], "valid")
})

test_that("store, list, get, audit, and delete share (type, result_id)", {
  object <- make_analysis_result_test_object()
  compact <- list(
    method = "slingshot", backend = "slingshot", input = list(cells = ncol(object)),
    parameters = list(start = "0"),
    tables = list(primary = tibble::tibble(
      cell = colnames(object), primary_pseudotime = seq_len(ncol(object))
    )),
    diagnostics = list(converged = TRUE), warnings = character(),
    provenance = list(random_seed = 1L)
  )

  object <- sn_store_result(
    object, type = "trajectory", result_id = "cd8_trajectory", result = compact
  )
  expect_identical(names(object@misc$shennong$results$trajectory), "cd8_trajectory")
  stored <- sn_get_result(object, "trajectory", "cd8_trajectory")
  expect_identical(stored$result_id, "cd8_trajectory")
  expect_null(stored$name)
  expect_identical(stored$provenance$result_id, "cd8_trajectory")
  expect_identical(stored$provenance$analysis_type, "trajectory")

  listing <- sn_list_results(object)
  expect_named(listing, c(
    "collection", "type", "result_id", "analysis", "method", "created_at",
    "n_rows", "source"
  ))
  expect_identical(listing$result_id, "cd8_trajectory")
  expect_false("name" %in% names(listing))

  audit <- sn_audit_results(object, include_artifacts = FALSE)
  expect_identical(audit$result_id, "cd8_trajectory")
  expect_identical(audit$status, "valid")
  expect_true(audit$unified)

  object <- sn_delete_result(object, "trajectory", "cd8_trajectory")
  expect_error(sn_get_result(object, "trajectory", "cd8_trajectory"), "result_id")
  expect_null(object@misc$shennong)
})

test_that("result identifiers are validated consistently", {
  object <- make_analysis_result_test_object()
  compact <- list(
    method = "custom", backend = "custom",
    tables = list(primary = tibble::tibble(value = 1))
  )
  expect_error(sn_store_result(object, "custom", "", compact), "non-empty")
  expect_error(sn_store_result(object, "custom", " padded ", compact), "whitespace")
  expect_error(sn_get_result(object, "custom", ""), "result_id")
})

test_that("future result schemas and mismatched stored identities fail closed", {
  object <- make_analysis_result_test_object()
  result <- Shennong:::.sn_new_analysis_result(
    "trajectory", "inner", "test", "test",
    tables = list(primary = tibble::tibble(
      cell = colnames(object),
      primary_pseudotime = seq_len(ncol(object))
    ))
  )

  future <- result
  future$schema_version <- "99.0.0"
  expect_error(
    sn_store_result(object, "trajectory", "future", future),
    "future `schema_version`"
  )

  object <- sn_store_result(object, "trajectory", "outer", result)
  object@misc$shennong$results$trajectory$outer$result_id <- "inner"
  object@misc$shennong$results$trajectory$outer$provenance$result_id <- "inner"
  expect_error(
    sn_get_result(object, "trajectory", "outer"),
    "does not match its storage key"
  )
  audit <- sn_audit_results(object, include_artifacts = FALSE)
  expect_false(audit$valid)
  expect_identical(audit$status, "repairable")
  expect_match(audit$errors, "storage key")

  repaired <- sn_upgrade_results(object)
  expect_identical(
    sn_get_result(repaired, "trajectory", "outer")$result_id,
    "outer"
  )
})

test_that("legacy two-part result schemas remain safely upgradeable", {
  legacy <- Shennong:::.sn_new_analysis_result(
    "de", "legacy", "test", "test",
    tables = list(primary = tibble::tibble(gene = "A"))
  )
  legacy$schema_version <- "1.0"

  upgraded <- Shennong:::.sn_upgrade_analysis_result(
    legacy,
    analysis_type = "de",
    result_id = "legacy"
  )

  expect_identical(upgraded$schema_version, "2.0.0")
  expect_identical(upgraded$provenance$migrated_from_schema_version, "1.0")
})

test_that("result upgrades accept only known legacy schema spellings", {
  legacy <- Shennong:::.sn_new_analysis_result(
    "de", "legacy", "test", "test",
    tables = list(primary = tibble::tibble(gene = "A", score = 1))
  )
  for (version in c("1", "1.0", "1.0.0")) {
    candidate <- legacy
    candidate$schema_version <- version
    upgraded <- Shennong:::.sn_upgrade_analysis_result(
      candidate, analysis_type = "de", result_id = "legacy"
    )
    expect_identical(upgraded$schema_version, "2.0.0")
  }
  for (version in c("0.1.0", "1.9.9", "2.0.0-beta")) {
    candidate <- legacy
    candidate$schema_version <- version
    expect_error(
      Shennong:::.sn_upgrade_analysis_result(
        candidate, analysis_type = "de", result_id = "legacy"
      ),
      "unsupported `schema_version`"
    )
  }
})

test_that("DE upgrades canonicalize unambiguous feature identifiers", {
  legacy <- Shennong:::.sn_new_analysis_result(
    "de", "legacy_feature", "test", "test",
    tables = list(primary = tibble::tibble(gene = c("A", "B")))
  )
  legacy$tables$primary <- tibble::tibble(feature = c("A", "B"), score = c(2, 1))

  upgraded <- Shennong:::.sn_upgrade_analysis_result(
    legacy, analysis_type = "de", result_id = "legacy_feature"
  )

  expect_identical(upgraded$tables$primary$gene, c("A", "B"))
  expect_identical(upgraded$provenance$migrated_primary_gene_from, "feature")
  expect_true(sn_validate_result(upgraded, error = FALSE)$valid)
})

test_that("DE upgrades explicit row names only when identifier sources agree", {
  legacy <- Shennong:::.sn_new_analysis_result(
    "de", "legacy_rows", "test", "test",
    tables = list(primary = tibble::tibble(gene = c("A", "B")))
  )
  row_named <- data.frame(score = c(2, 1), row.names = c("A", "B"))
  legacy$tables$primary <- row_named
  upgraded <- Shennong:::.sn_upgrade_analysis_result(
    legacy, analysis_type = "de", result_id = "legacy_rows"
  )
  expect_identical(upgraded$tables$primary$gene, c("A", "B"))
  expect_identical(upgraded$provenance$migrated_primary_gene_from, "rownames")

  legacy$tables$primary <- data.frame(
    feature = c("A", "B"), score = c(2, 1), row.names = c("X", "Y")
  )
  expect_error(
    Shennong:::.sn_upgrade_analysis_result(
      legacy, analysis_type = "de", result_id = "legacy_rows"
    ),
    "feature.*row names disagree"
  )
})

test_that("analysis results and runtime artifacts remain separate", {
  object <- make_analysis_result_test_object()
  object@misc$input_source <- list(path = "/data/example", format = "10x")
  object@misc$integration_comparison <- list(
    method = "grid", grid = tibble::tibble(method = c("harmony", "unintegrated"))
  )
  object@misc$user_payload <- list(note = "belongs to the user")

  expect_equal(nrow(sn_list_results(object)), 0L)
  with_artifacts <- sn_list_results(object, include_artifacts = TRUE)
  expect_true(all(c("input_source", "integration_comparison") %in% with_artifacts$collection))
  audit <- sn_audit_results(object)
  expect_true(any(audit$contract_scope == "artifact"))
  expect_true(any(audit$collection == "user_payload" & audit$contract_scope == "unregistered"))
})

test_that("whole legacy artifact collections require explicit confirmation", {
  object <- make_analysis_result_test_object()
  object@misc$integration <- list(user_payload = list(note = "keep me"))

  expect_error(
    sn_delete_artifact(object, "integration"),
    "confirm = TRUE"
  )
  expect_true("integration" %in% names(object@misc))

  deleted <- sn_delete_artifact(object, "integration", confirm = TRUE)
  expect_false("integration" %in% names(deleted@misc))
})

test_that("label transfer stores a canonical traceable annotation result", {
  object <- make_analysis_result_test_object()
  object$transfer_label <- c("T cell", "T cell", "B cell", "B cell")
  object$transfer_score <- c(0.9, 0.8, 0.7, 0.6)
  object <- Shennong:::.sn_store_label_transfer_result(
    object = object, prediction_prefix = "transfer",
    result_id = "pbmc_reference_transfer", method = "seurat",
    label_by = "cell_type",
    prediction_columns = c("transfer_label", "transfer_score")
  )
  stored <- sn_get_result(object, "annotation", "pbmc_reference_transfer")
  expect_identical(stored$result_id, "pbmc_reference_transfer")
  expect_identical(stored$parameters$prediction_prefix, "transfer")
  expect_equal(nrow(stored$tables$primary), ncol(object))
})

test_that("registered table stores return v2 envelopes", {
  object <- make_analysis_result_test_object()
  returned <- list(
    sn_store_deconvolution(
      object,
      tibble::tibble(sample = "sample1", cell_type = "T", fraction = 1),
      return_object = FALSE
    ),
    sn_store_regulatory_activity(
      object,
      tibble::tibble(source = "STAT1", condition = "cell1", score = 1),
      return_object = FALSE
    ),
    sn_store_milo(
      object,
      tibble::tibble(Nhood = 1, logFC = 1, SpatialFDR = 0.05),
      sample_by = "sample", group_by = "condition",
      return_object = FALSE
    ),
    sn_store_enrichment(
      object,
      tibble::tibble(ID = "GO:1", Description = "term", p.adjust = 0.05),
      return_object = FALSE
    )
  )
  expect_true(all(vapply(returned, function(result) {
    sn_validate_result(result, error = FALSE)$valid &&
      identical(result$result_id, "default")
  }, logical(1))))
})

test_that("generic storage rejects artifact-reserved namespaces", {
  object <- make_analysis_result_test_object()
  result <- list(
    method = "synthetic", backend = "synthetic",
    tables = list(primary = tibble::tibble(value = 1))
  )
  expect_error(
    sn_store_result(object, "coralysis_artifact", "bad", result),
    "reserved for a registered workflow artifact"
  )
  expect_error(
    sn_store_result(object, "input_source", "bad", result),
    "reserved for a registered workflow artifact"
  )
})
