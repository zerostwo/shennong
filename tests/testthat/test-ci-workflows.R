test_that("routine package checks avoid duplicate full tests", {
  path <- test_path("..", "..", ".github", "workflows", "R-CMD-check.yaml")
  skip_if_not(file.exists(path), "Repository-only workflow is excluded from source packages.")
  workflow <- readLines(path, warn = FALSE)
  text <- paste(workflow, collapse = "\n")

  expect_match(text, "workflow_dispatch:", fixed = TRUE)
  expect_match(text, "--no-tests", fixed = TRUE)
  expect_match(text, "cancel-in-progress: true", fixed = TRUE)
})

test_that("pkgdown supports incremental and explicit full builds", {
  workflow_path <- test_path("..", "..", ".github", "workflows", "pkgdown.yaml")
  helper_path <- test_path("..", "..", "scripts", "build-pkgdown.R")
  skip_if_not(file.exists(workflow_path) && file.exists(helper_path), "Repository-only CI files are excluded from source packages.")
  workflow <- readLines(workflow_path, warn = FALSE)
  helper <- readLines(helper_path, warn = FALSE)
  workflow_text <- paste(workflow, collapse = "\n")
  helper_text <- paste(helper, collapse = "\n")

  expect_match(workflow_text, "lazy = !full", fixed = TRUE)
  expect_match(workflow_text, "clean = full", fixed = TRUE)
  expect_match(workflow_text, "examples = full", fixed = TRUE)
  expect_match(helper_text, "lazy = !full", fixed = TRUE)
  expect_match(helper_text, "--full", fixed = TRUE)
})
