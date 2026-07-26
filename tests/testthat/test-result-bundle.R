library(testthat)

make_result_bundle_test_result <- function() {
  list(
    schema_version = "1.0.0",
    analysis_type = "bulk_de",
    name = "treated_vs_control",
    method = "limma",
    backend = "limma",
    input = list(features = 2L, samples = 4L),
    parameters = list(design = ~ condition),
    tables = list(primary = tibble::tibble(
      feature = c("gene1", "gene2"),
      estimate = c(1.5, -0.5)
    )),
    embeddings = list(),
    graphs = list(),
    models = list(),
    diagnostics = list(converged = TRUE),
    warnings = character(),
    provenance = list(
      package_versions = list(Shennong = "0.2.0.9000", R = "4.6.0"),
      random_seed = 717L,
      timestamp = "2026-07-26 UTC"
    )
  )
}

test_that("result bundle v1 carries canonical result and immutable references", {
  result <- make_result_bundle_test_result()
  input_digest <- paste(rep("a", 64L), collapse = "")
  artifact_digest <- paste(rep("b", 64L), collapse = "")

  bundle <- sn_build_result_bundle(
    result,
    inputs = list(list(
      role = "expression",
      resource_id = "toil-rna",
      revision = "revision-7",
      digest = list(algorithm = "sha256", value = input_digest)
    )),
    execution = list(
      job_id = "job-42",
      runtime = "shennong-runtime",
      image_digest = paste0("sha256:", paste(rep("c", 64L), collapse = "")),
      token_count = 2048L
    ),
    artifacts = list(list(
      role = "primary_table",
      path = "tables/primary.csv",
      media_type = "text/csv",
      size_bytes = 128,
      digest = list(algorithm = "sha256", value = artifact_digest)
    )),
    created_at = "2026-07-26T00:00:00Z"
  )

  expect_s3_class(bundle, "sn_result_bundle")
  expect_identical(bundle$schema, "shennong.dev/analysis-result-bundle/v1")
  expect_identical(bundle$result$schema_version, "1.0.0")
  expect_identical(bundle$result$analysis_type, "bulk_de")
  expect_identical(bundle$result$parameters$design, "~condition")
  expect_true(bundle$validation$valid)
  expect_identical(bundle$inputs[[1L]]$revision, "revision-7")
  expect_identical(bundle$inputs[[1L]]$digest$value, input_digest)
  expect_identical(bundle$provenance$execution$job_id, "job-42")
  expect_identical(bundle$artifacts[[1L]]$role, "primary_table")
  expect_true(sn_validate_result_bundle(bundle, error = FALSE)$valid)
})

test_that("result bundle validation rejects mutable or sensitive references", {
  result <- make_result_bundle_test_result()
  bad_digest <- list(algorithm = "sha256", value = "not-a-sha256")

  expect_error(
    sn_build_result_bundle(
      result,
      inputs = list(list(
        role = "expression",
        resource_id = "toil-rna",
        revision = "revision-7",
        digest = bad_digest
      ))
    ),
    "64 hexadecimal"
  )
  expect_error(
    sn_build_result_bundle(
      result,
      inputs = list(list(
        role = "expression",
        revision = "revision-7",
        digest = list(
          algorithm = "sha256",
          value = paste(rep("a", 64L), collapse = "")
        )
      ))
    ),
    "resource_id.*artifact_id"
  )
  expect_error(
    sn_build_result_bundle(
      result,
      execution = list(Authorization = "Bearer must-not-be-bundled")
    ),
    "credentials"
  )
  expect_error(
    sn_build_result_bundle(result, created_at = "2026-07-26 UTC"),
    "RFC 3339"
  )

  binary_result <- result
  binary_result$models$sparse <- Matrix::Matrix(diag(2), sparse = TRUE)
  expect_error(
    sn_build_result_bundle(binary_result),
    "S4 object"
  )
})

test_that("result bundles export atomically with a SHA-256 handoff digest", {
  result <- make_result_bundle_test_result()
  bundle <- sn_build_result_bundle(result, created_at = "2026-07-26T00:00:00Z")
  output <- tempfile(fileext = ".json")

  exported <- sn_export_result_bundle(bundle, output)
  expect_true(file.exists(output))
  expect_identical(exported$path, normalizePath(output, winslash = "/", mustWork = TRUE))
  expect_identical(exported$digest$algorithm, "sha256")
  expect_match(exported$digest$value, "^[0-9a-f]{64}$")
  expect_gt(exported$size_bytes, 0)

  parsed <- jsonlite::read_json(output, simplifyVector = FALSE)
  expect_identical(parsed$schema, "shennong.dev/analysis-result-bundle/v1")
  expect_identical(parsed$result$schema_version, "1.0.0")
  expect_identical(parsed$validation$valid, TRUE)

  expect_error(
    sn_export_result_bundle(bundle, output),
    "already exists"
  )
  overwritten <- sn_export_result_bundle(bundle, output, overwrite = TRUE)
  expect_identical(overwritten$digest$algorithm, "sha256")
})

test_that("failed cross-platform replacement restores the previous bundle", {
  directory <- tempfile("result-bundle-publish-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE, force = TRUE), add = TRUE)
  target <- file.path(directory, "bundle.json")
  temporary <- file.path(directory, "candidate.json")
  writeLines("previous", target)
  writeLines("candidate", temporary)

  rename_with_failed_candidate <- function(from, to) {
    if (identical(from, temporary) && identical(to, target)) return(FALSE)
    file.rename(from, to)
  }
  expect_error(
    Shennong:::.sn_result_bundle_publish(
      temporary,
      target,
      overwrite = TRUE,
      rename_file = rename_with_failed_candidate
    ),
    "restored"
  )
  expect_identical(readLines(target), "previous")
  expect_identical(readLines(temporary), "candidate")
})
