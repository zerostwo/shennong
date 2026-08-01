test_that("single-cell AutoZyme benchmark evidence is complete and exact", {
  path <- system.file(
    "benchmarks",
    "single-cell-autozyme-benchmark.json",
    package = "Shennong"
  )
  expect_true(nzchar(path))
  expect_true(file.exists(path))

  report <- jsonlite::read_json(path, simplifyVector = TRUE)
  expect_identical(report$schema_version, "1.0.0")
  expect_identical(report$methodology$repetitions_per_condition, 3L)
  expect_true(report$methodology$fresh_process_per_run)
  expect_identical(report$dataset$cells, 2000L)
  expect_identical(report$dataset$features, 32738L)

  summary <- as.data.frame(report$summary, stringsAsFactors = FALSE)
  expect_setequal(
    summary$operation,
    c("merge", "joinlayers", "scdblfinder")
  )
  expect_setequal(
    summary$patch,
    c("seurat_merge", "seurat_joinlayers", "scdblfinder")
  )
  expect_true(all(summary$baseline_seconds > 0))
  expect_true(all(summary$accelerated_seconds > 0))
  expect_true(all(summary$speedup > 1))
  expect_true(all(summary$baseline_peak_rss_mb > 0))
  expect_true(all(summary$accelerated_peak_rss_mb > 0))
  expect_true(all(summary$output_identical))
  expect_true(all(summary$patch_inactive_after))

  runs <- as.data.frame(report$runs, stringsAsFactors = FALSE)
  expect_equal(nrow(runs), 18L)
  expect_setequal(runs$condition, c("baseline", "accelerated"))
  expect_true(all(runs$valid))
  expect_false(any(runs$active_after))
  fingerprints <- split(runs$output_fingerprint, runs$operation)
  expect_true(all(vapply(
    fingerprints,
    function(value) length(unique(value)) == 1L,
    logical(1)
  )))
})
