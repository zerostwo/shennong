# Conformance evidence for the annotation::popv backend.
#
# C1: dispatch/contract on the committed tiny fixture (check-safe, skips
#     cleanly when the popv pixi environment is not installed).
# C2: direct parity against a fresh-process PopV oracle on the committed
#     integration fixture. The oracle calls the upstream public API without
#     any Shennong wrapper code.

.conformance_popv_fixture_path <- function(name) {
  file.path(.conformance_project_root(), "tests", "conformance", "fixtures", name)
}

.conformance_popv_load_fixture <- function(name) {
  readRDS(.conformance_popv_fixture_path(name))
}

# Reduced algorithm set: keeps C2 runtime bounded while exercising one
# classifier from each supported family (linear kernel SVM, boosted trees,
# CellTypist's logistic regression). scVI-backed algorithms are covered by the
# scheduled C3 tier, not this PR gate.
.popv_test_methods <- c("Support_Vector", "XGboost", "CELLTYPIST")

.popv_test_control <- function(output_dir) {
  list(
    output_dir = output_dir,
    methods = .popv_test_methods,
    hvg = NULL,
    n_samples_per_label = 50,
    quiet = TRUE
  )
}

.conformance_popv_export_inputs <- function(query, reference, input_dir) {
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  query_info <- Shennong:::.sn_write_python_object_input(
    object = query,
    input_dir = file.path(input_dir, "query"),
    assay = "RNA",
    layer = "counts"
  )
  reference_info <- Shennong:::.sn_write_python_object_input(
    object = reference,
    input_dir = file.path(input_dir, "reference"),
    assay = "RNA",
    layer = "counts"
  )
  list(input_dir = input_dir, query = query_info, reference = reference_info)
}

test_that("annotation::popv contract is admitted and structurally complete", {
  contract <- .conformance_contract("annotation::popv")
  expect_identical(contract$status, "admitted")
  expect_identical(contract$shennong$registry_key, "annotation::popv")
  expect_identical(contract$evidence$verdict, "pass")
})

test_that("annotation::popv dispatches to the pixi popv backend (C1)", {
  skip_if_not_installed("Seurat")
  .conformance_require_pixi_environment("popv")
  fixture <- .conformance_popv_load_fixture("popv-pbmc3k-tiny-v1.rds")
  run_dir <- file.path(tempdir(), paste0("sn_popv_c1_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  result <- sn_run_annotation(
    fixture$query,
    group_by = "cell_type",
    method = "popv",
    reference = fixture$reference,
    reference_label_by = "cell_type",
    species = "human",
    ontology = FALSE,
    result_id = "popv_c1",
    return_object = FALSE,
    backend_control = list(popv = .popv_test_control(file.path(run_dir, "c1")))
  )
  expect_identical(result$backend, "popv")
  expect_identical(result$method, "popv")
  expect_equal(nrow(result$tables$backend_predictions), ncol(fixture$query))
  expect_setequal(
    unique(result$tables$backend_predictions$prediction),
    c("B cells", "CD4+ T cells", "Monocytes")
  )
  stored_labels <- result$tables$cells$prediction
  expect_false(any(is.na(stored_labels)))
})

test_that("annotation::popv matches the direct PopV oracle (C2)", {
  skip_if_not_installed("Seurat")
  .conformance_require_pixi_environment("popv")
  fixture <- .conformance_popv_load_fixture("popv-pbmc3k-integration-v1.rds")
  seed <- 717L
  run_dir <- file.path(tempdir(), paste0("sn_popv_c2_", format(Sys.time(), "%Y%m%d_%H%M%S")))

  # Candidate: Shennong public entry point.
  candidate_result <- sn_run_annotation(
    fixture$query,
    group_by = "cell_type",
    method = "popv",
    reference = fixture$reference,
    reference_label_by = "cell_type",
    species = "human",
    ontology = FALSE,
    result_id = "popv_c2",
    return_object = FALSE,
    backend_control = list(popv = .popv_test_control(file.path(run_dir, "candidate")))
  )

  # Oracle: direct upstream call in the same pixi environment, same inputs,
  # same seed, no Shennong analysis code involved.
  exported <- .conformance_popv_export_inputs(
    fixture$query,
    fixture$reference,
    file.path(run_dir, "oracle", "input")
  )
  oracle_config <- list(
    label_key = "cell_type",
    prediction_mode = "retrain",
    unknown_celltype_label = "unknown",
    n_samples_per_label = 50,
    cl_obo_folder = FALSE,
    min_shared_genes = 2,
    methods = .popv_test_methods,
    hvg = NULL,
    seed = seed
  )
  config_path <- file.path(run_dir, "oracle", "config.json")
  Shennong:::.sn_write_json_file(oracle_config, config_path)
  Shennong::sn_call_pixi_environment(
    environment = "popv",
    command = "python",
    args = c(
      shQuote(file.path(
        .conformance_project_root(),
        "tests", "conformance", "oracles", "popv_oracle.py"
      )),
      "--input-dir", shQuote(exported$input_dir),
      "--output-dir", shQuote(file.path(run_dir, "oracle", "output")),
      "--config", shQuote(config_path)
    ),
    install_pixi = FALSE,
    quiet = TRUE
  )
  oracle_predictions <- utils::read.csv(
    file.path(run_dir, "oracle", "output", "predictions.csv"),
    row.names = 1,
    check.names = FALSE
  )
  oracle_labels <- setNames(
    as.character(oracle_predictions[["popv_prediction"]]),
    rownames(oracle_predictions)
  )
  oracle_agreement <- setNames(
    as.numeric(oracle_predictions[["popv_majority_vote_score"]]),
    rownames(oracle_predictions)
  )

  candidate_raw <- candidate_result$tables$backend_predictions
  candidate_labels <- setNames(candidate_raw$prediction, candidate_raw$cell)
  shared <- intersect(names(candidate_labels), names(oracle_labels))
  expect_gt(length(shared), 0.9 * ncol(fixture$query))

  compared_labels <- candidate_labels[shared] == oracle_labels[shared]
  expect_true(all(compared_labels), info = paste(
    "label mismatches:",
    paste(shared[!compared_labels], collapse = ", ")
  ))
  expect_identical(
    unname(candidate_raw$agreement[match(shared, candidate_raw$cell)]),
    unname(oracle_agreement[shared])
  )

  # Wrapper-owned additions: standardized result envelope stays intact.
  expect_identical(candidate_result$schema_version, "2.0.0")
  expect_identical(candidate_result$analysis_type, "annotation")
  expect_true("popv" %in% candidate_result$tables$evidence$method)
})
