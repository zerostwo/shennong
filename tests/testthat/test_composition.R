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

test_that("sn_calculate_composition warns when additional columns vary within groups", {
  meta_df <- data.frame(
    sample = c("A", "A", "A", "B", "B"),
    ctype = c("Tcell", "Bcell", "Tcell", "Tcell", "Bcell"),
    condition = c("Control", "Treated", "Control", "Control", "Control"),
    stringsAsFactors = FALSE
  )

  expect_warning(
    result <- sn_calculate_composition(
      x = meta_df,
      group_by = "sample",
      variable = "ctype",
      min_cells = 1,
      additional_cols = "condition"
    ),
    "Additional columns are not constant"
  )

  expect_equal(unique(result$condition[result$sample == "A"]), "Control")
})

test_that("sn_calculate_composition errors when min_cells removes all groups", {
  meta_df <- data.frame(
    sample = c("A", "A", "B"),
    ctype = c("Tcell", "Bcell", "Tcell"),
    stringsAsFactors = FALSE
  )

  expect_error(
    sn_calculate_composition(
      x = meta_df,
      group_by = "sample",
      variable = "ctype",
      min_cells = 3
    ),
    "No groups remaining after filtering by `min_cells`"
  )
})

test_that("sn_calculate_composition supports Seurat-object metadata directly", {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(matrix(rpois(60, lambda = 2), nrow = 10, ncol = 6), sparse = TRUE)
  rownames(counts) <- paste0("gene", seq_len(10))
  colnames(counts) <- paste0("cell", seq_len(6))

  object <- sn_initialize_seurat_object(x = counts, project = "composition-seurat")
  object$sample <- rep(c("A", "B"), each = 3)
  object$cell_type <- c("Tcell", "Tcell", "Bcell", "Bcell", "Myeloid", "Bcell")
  object$condition <- rep(c("ctrl", "tx"), each = 3)

  result <- sn_calculate_composition(
    x = object,
    group_by = "sample",
    variable = "cell_type",
    min_cells = 1,
    additional_cols = "condition"
  )

  expect_s3_class(result, "data.frame")
  expect_true(all(c("sample", "cell_type", "proportion", "condition") %in% colnames(result)))
  expect_equal(
    as.numeric(tapply(result$proportion, result$sample, sum)),
    c(100, 100)
  )
})

test_that("sn_calculate_composition validates Seurat/data-frame inputs and additional columns", {
  expect_error(
    sn_calculate_composition(
      x = list(sample = "A"),
      group_by = "sample",
      variable = "cell_type"
    ),
    "must be a Seurat object or a data frame"
  )

  meta_df <- data.frame(
    sample = c("A", "A", "B", "B"),
    cell_type = c("Tcell", "Bcell", "Tcell", "Bcell"),
    stringsAsFactors = FALSE
  )

  expect_error(
    sn_calculate_composition(
      x = meta_df,
      group_by = "sample",
      variable = "cell_type",
      additional_cols = "missing_col"
    ),
    "Missing additional columns"
  )
})
