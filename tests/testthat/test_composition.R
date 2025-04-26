library(testthat)
library(dplyr)

test_that("sn_calculate_composition works with data frames", {
  meta_df <- data.frame(
    sample = rep(c("A", "B", "C"), each = 50),
    ctype = sample(c("Tcell", "Bcell", "Myeloid"), 150, replace = TRUE),
    condition = rep(c("Control", "Treated", "Control"), each = 50),
    stringsAsFactors = FALSE
  )

  result <- sn_calculate_composition(
    x = meta_df,
    group_by = "sample",
    variable = "ctype",
    min_cells = 5,
    additional_cols = "condition"
  )

  expect_s3_class(result, "data.frame")
  expect_true(all(c("sample", "ctype", "proportion", "condition") %in% colnames(result)))
  expect_true(all(result$proportion >= 0 & result$proportion <= 100))
})

test_that("sn_calculate_composition filters groups with fewer cells than min_cells", {
  small_meta_df <- data.frame(
    sample = c("A", "A", "B"),
    ctype = c("Tcell", "Bcell", "Tcell"),
    stringsAsFactors = FALSE
  )

  result <- sn_calculate_composition(
    x = small_meta_df,
    group_by = "sample",
    variable = "ctype",
    min_cells = 2
  )

  expect_equal(nrow(result), 2)
})

test_that("sn_calculate_composition handles missing columns gracefully", {
  meta_df <- data.frame(
    sample = rep(c("A", "B", "C"), each = 50),
    ctype = sample(c("Tcell", "Bcell", "Myeloid"), 150, replace = TRUE),
    stringsAsFactors = FALSE
  )

  expect_error(
    sn_calculate_composition(
      x = meta_df,
      group_by = "sample",
      variable = "missing_column"
    ),
    "Missing required columns"
  )
})

test_that("sn_calculate_composition works with lightweight data frames", {
  meta_df <- data.frame(
    sample_id = rep(c("A", "B"), each = 10),
    cell_type = sample(c("Tcell", "Bcell"), 20, replace = TRUE),
    stringsAsFactors = FALSE
  )

  result <- sn_calculate_composition(
    x = meta_df,
    group_by = "sample_id",
    variable = "cell_type",
    min_cells = 5
  )

  expect_s3_class(result, "data.frame")
  expect_true(all(c("sample_id", "cell_type", "proportion") %in% colnames(result)))
})
