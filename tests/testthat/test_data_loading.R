library(testthat)

test_that("sn_load_data returns cached file path when return_object is FALSE", {
  save_dir <- tempfile("shennong-cache-")
  dir.create(save_dir)
  local_h5 <- file.path(save_dir, "pbmc1k_filtered_feature_bc_matrix.h5")
  file.create(local_h5)

  expect_equal(
    sn_load_data(
      dataset = "pbmc1k",
      save_dir = save_dir,
      return_object = FALSE
    ),
    local_h5
  )
})

test_that("sn_load_pbmc delegates to sn_load_data with a deprecation warning", {
  save_dir <- tempfile("shennong-cache-")
  dir.create(save_dir)
  local_h5 <- file.path(save_dir, "pbmc3k_filtered_feature_bc_matrix.h5")
  file.create(local_h5)

  expect_warning(
    path <- sn_load_pbmc(
      dataset = "pbmc3k",
      save_dir = save_dir,
      return_object = FALSE
    ),
    "sn_load_data"
  )

  expect_equal(path, local_h5)
})
