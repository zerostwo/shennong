library(testthat)

test_that("sn_check_file reports missing files without stopping when requested", {
  existing <- tempfile("existing-")
  file.create(existing)
  missing <- tempfile("missing-")

  expect_equal(
    sn_check_file(c(existing, missing), stop = FALSE),
    missing
  )
})

test_that("sn_check_file errors for missing files by default", {
  expect_error(
    sn_check_file(tempfile("missing-")),
    "does not exist"
  )
})

test_that("sn_set_path creates a directory and returns its path", {
  dir_path <- tempfile("shennong-dir-")

  expect_equal(sn_set_path(dir_path), dir_path)
  expect_true(dir.exists(dir_path))
})

test_that("sn_get_species returns the explicit species when provided", {
  expect_equal(sn_get_species(object = NULL, species = "human"), "human")
})
