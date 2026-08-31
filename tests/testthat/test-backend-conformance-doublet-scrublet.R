# Conformance evidence for the doublets::scrublet backend (pilot).
#
# C1: dispatch/structure on the committed tiny fixture (skips cleanly when the
#     scrublet pixi environment is not installed).
# C2: exact parity against a fresh-process scanpy oracle on the committed
#     integration fixture. The oracle calls sc.pp.scrublet directly with the
#     same seed and inputs, without any Shennong runner code.

.conformance_scrublet_fixture_path <- function(name) {
  file.path(.conformance_project_root(), "tests", "conformance", "fixtures", name)
}

.conformance_scrublet_load_fixture <- function(name) {
  readRDS(.conformance_scrublet_fixture_path(name))
}

.conformance_scrublet_export <- function(object, input_dir) {
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  Shennong:::.sn_write_python_object_input(
    object = object,
    input_dir = input_dir,
    assay = "RNA",
    layer = "counts"
  )
}

.scrublet_test_control <- function(output_dir) {
  list(
    output_dir = output_dir,
    seed = 717L,
    quiet = TRUE
  )
}

test_that("doublets::scrublet pilot contract is structurally complete", {
  contract <- .conformance_contract("doublets::scrublet")
  expect_identical(contract$status, "pilot")
  expect_true(is.null(contract$shennong$registry_key))
  expect_match(contract$ci$required_tiers[[length(contract$ci$required_tiers)]], "^C2")
})

test_that("sn_find_doublets(method='scrublet') matches the direct scanpy oracle (C2)", {
  skip_if_not_installed("Seurat")
  .conformance_require_pixi_environment("scrublet")

  fixture <- .conformance_scrublet_load_fixture("scrublet-pbmc3k-integration-v1.rds")
  run_dir <- file.path(tempdir(), paste0("sn_scrublet_c2_", format(Sys.time(), "%Y%m%d_%H%M%S")))

  # Candidate: Shennong public entry point.
  candidate <- sn_find_doublets(
    fixture,
    method = "scrublet",
    min_features = 1,
    backend_control = .scrublet_test_control(file.path(run_dir, "candidate"))
  )
  candidate_scores <- setNames(
    as.numeric(candidate$scrublet.score),
    colnames(candidate)
  )
  candidate_calls <- setNames(
    as.character(candidate$scrublet.class),
    colnames(candidate)
  )

  # Oracle: direct upstream call in the same pixi environment.
  exported <- .conformance_scrublet_export(fixture, file.path(run_dir, "oracle", "input"))
  config_path <- Shennong:::.sn_write_json_file(
    list(seed = 717),
    file.path(run_dir, "oracle", "config.json")
  )
  Shennong::sn_call_pixi_environment(
    environment = "scrublet",
    command = "python",
    args = c(
      shQuote(file.path(
        .conformance_project_root(),
        "tests", "conformance", "oracles", "scrublet_oracle.py"
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
  oracle_scores <- setNames(
    as.numeric(oracle_predictions$doublet_score),
    rownames(oracle_predictions)
  )
  oracle_calls <- ifelse(
    as.logical(oracle_predictions$is_doublet),
    "doublet",
    "singlet"
  )
  names(oracle_calls) <- rownames(oracle_predictions)

  shared <- intersect(names(candidate_scores), names(oracle_scores))
  expect_gt(length(shared), 0.9 * ncol(fixture))
  expect_identical(
    unname(candidate_scores[shared]),
    unname(oracle_scores[shared])
  )
  expect_identical(
    unname(candidate_calls[shared]),
    unname(oracle_calls[shared])
  )
})

test_that("doublets::scrublet dispatch keeps skipped-cell contract (C1)", {
  skip_if_not_installed("Seurat")
  .conformance_require_pixi_environment("scrublet")

  fixture <- .conformance_scrublet_load_fixture("scrublet-pbmc3k-tiny-v1.rds")
  counts <- SeuratObject::LayerData(fixture, layer = "counts")
  counts[, 1] <- 0
  SeuratObject::LayerData(fixture, layer = "counts") <- counts

  updated <- sn_find_doublets(
    fixture,
    method = "scrublet",
    min_features = 1,
    backend_control = .scrublet_test_control(file.path(
      tempdir(),
      paste0("sn_scrublet_c1_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    ))
  )
  expect_s4_class(updated, "Seurat")
  expect_identical(as.character(updated$scrublet.class[[1]]), "unresolved")
  expect_true(is.na(updated$scrublet.score[[1]]))
  expect_length(unique(as.character(updated$scrublet.class[-1])), 2L)
})
