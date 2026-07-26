runtime_coverage_paths <- function() {
  root <- test_path("..", "..")
  list(
    root = root,
    runner = file.path(root, "scripts", "real-data", "run-runtime-coverage.R"),
    validator = file.path(root, "scripts", "real-data", "validate-real-data.R"),
    coverage = file.path(root, "scripts", "real-data", "coverage.csv"),
    exclusions = file.path(
      root, "scripts", "real-data", "coverage-exclusions.csv"
    )
  )
}

run_runtime_coverage_cli <- function(script, arguments) {
  output <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c(shQuote(script), vapply(arguments, shQuote, character(1))),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(output, "status")
  list(
    status = if (is.null(status)) 0L else as.integer(status),
    output = paste(output, collapse = "\n")
  )
}

test_that("runtime coverage CLI documents its local-only profiles", {
  paths <- runtime_coverage_paths()
  skip_if_not(file.exists(paths$runner), "Repository-only helper is excluded from source packages.")

  help <- run_runtime_coverage_cli(paths$runner, "--help")
  expect_identical(help$status, 0L)
  expect_match(help$output, "--data-root PATH", fixed = TRUE)
  expect_match(help$output, "--extended", fixed = TRUE)
  expect_match(help$output, "never downloads data", fixed = TRUE)
  expect_match(help$output, "current Shennong checkout namespace", fixed = TRUE)
})

test_that("runtime coverage rejects a bad fixture root before rendering", {
  paths <- runtime_coverage_paths()
  skip_if_not(
    file.exists(paths$runner) && file.exists(paths$validator),
    "Repository-only real-data helpers are excluded from source packages."
  )

  missing_root <- tempfile("shennong-runtime-coverage-missing-")
  unlink(missing_root, recursive = TRUE, force = TRUE)
  result <- run_runtime_coverage_cli(
    paths$runner,
    c("--data-root", missing_root)
  )
  expect_gt(result$status, 0L)
  expect_match(
    result$output,
    "Validating all four local real-data bundles before runtime tracing.",
    fixed = TRUE
  )
  expect_match(
    result$output,
    "runtime coverage was not rendered and no data were downloaded",
    fixed = TRUE
  )
  expect_false(dir.exists(missing_root))
})

test_that("coverage article mappings are static and source-backed", {
  paths <- runtime_coverage_paths()
  skip_if_not(file.exists(paths$coverage), "Repository-only coverage matrix is excluded from source packages.")

  coverage <- utils::read.csv(
    paths$coverage,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = character()
  )
  expect_true("article" %in% names(coverage))
  mapped <- nzchar(trimws(coverage[["article"]]))
  expect_gt(sum(mapped), 0L)

  for (index in which(mapped)) {
    article_path <- file.path(paths$root, coverage[["article"]][[index]])
    expect_true(
      file.exists(article_path),
      info = paste("missing mapped article for", coverage[["function"]][[index]])
    )
    if (!file.exists(article_path)) next
    article_source <- readLines(article_path, warn = FALSE)
    pattern <- paste0(
      "\\b", coverage[["function"]][[index]], "\\s*\\("
    )
    expect_true(
      any(grepl(pattern, article_source, perl = TRUE)),
      info = paste(
        coverage[["function"]][[index]],
        "does not appear in",
        coverage[["article"]][[index]]
      )
    )
  }
})

test_that("every exported sn function is covered or explicitly out of scope", {
  paths <- runtime_coverage_paths()
  skip_if_not(
    file.exists(paths$coverage) && file.exists(paths$exclusions),
    "Repository-only coverage classifications are excluded from source packages."
  )

  coverage <- utils::read.csv(
    paths$coverage, stringsAsFactors = FALSE, check.names = FALSE
  )
  exclusions <- utils::read.csv(
    paths$exclusions, stringsAsFactors = FALSE, check.names = FALSE
  )
  namespace <- readLines(file.path(paths$root, "NAMESPACE"), warn = FALSE)
  exports <- sub(
    "^export\\(([^)]+)\\)$", "\\1",
    grep("^export\\(sn_", namespace, value = TRUE)
  )

  expect_identical(anyDuplicated(exclusions[["function"]]), 0L)
  expect_true(all(nzchar(trimws(exclusions[["reason"]]))))
  expect_length(
    intersect(coverage[["function"]], exclusions[["function"]]), 0L
  )
  expect_setequal(
    c(coverage[["function"]], exclusions[["function"]]), exports
  )
})

test_that("runtime coverage traces and reports without enabling external downloads", {
  paths <- runtime_coverage_paths()
  skip_if_not(file.exists(paths$runner), "Repository-only helper is excluded from source packages.")
  source <- paste(readLines(paths$runner, warn = FALSE), collapse = "\n")

  expect_match(source, "validate-real-data.R", fixed = TRUE)
  expect_match(source, "pkgload::load_all", fixed = TRUE)
  expect_match(source, "detach(\"package:Shennong\"", fixed = TRUE)
  expect_match(source, "get(\"trace\", envir = asNamespace(\"methods\"))", fixed = TRUE)
  expect_match(source, "rmarkdown::render", fixed = TRUE)
  expect_match(source, "runtime-coverage.csv", fixed = TRUE)
  expect_match(source, "runtime-coverage-summary.json", fixed = TRUE)
  expect_match(source, "runtime-coverage-\", profile_label", fixed = TRUE)
  expect_match(source, "runtime-articles", fixed = TRUE)
  expect_match(source, "SHENNONG_RUN_EXTERNAL_BACKENDS = \"false\"", fixed = TRUE)
  expect_match(source, "runtime$level == \"core\"", fixed = TRUE)
  expect_match(source, "\"missing\", \"not_run\"", fixed = TRUE)
})
