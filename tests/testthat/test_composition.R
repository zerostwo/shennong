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

test_that("sn_calculate_composition supports multi-column grouping", {
  meta_df <- data.frame(
    sample = c("A", "A", "A", "A", "B", "B"),
    cell_type = c("Tcell", "Tcell", "Bcell", "Bcell", "Tcell", "Tcell"),
    Phase = c("G1", "S", "G1", "G1", "S", "S"),
    stringsAsFactors = FALSE
  )

  result <- sn_calculate_composition(
    x = meta_df,
    group_by = c("sample", "cell_type"),
    variable = "Phase",
    min_cells = 1
  )

  expect_true(all(c("sample", "cell_type", "Phase", "proportion") %in% colnames(result)))
  expect_equal(
    result |>
      dplyr::group_by(sample, cell_type) |>
      dplyr::summarise(total = sum(proportion), .groups = "drop") |>
      dplyr::pull(total),
    c(100, 100, 100)
  )
  expect_equal(
    result$proportion[result$sample == "A" & result$cell_type == "Tcell" & result$Phase == "G1"],
    50
  )
  expect_equal(
    result$proportion[result$sample == "A" & result$cell_type == "Tcell" & result$Phase == "S"],
    50
  )
})

test_that("sn_calculate_composition can return counts or both counts and proportions", {
  meta_df <- data.frame(
    sample = c("A", "A", "A", "A"),
    cell_type = c("Tcell", "Tcell", "Tcell", "Bcell"),
    Phase = c("G1", "S", "S", "G1"),
    stringsAsFactors = FALSE
  )

  count_result <- sn_calculate_composition(
    x = meta_df,
    group_by = c("sample", "cell_type"),
    variable = "Phase",
    min_cells = 1,
    measure = "count"
  )

  both_result <- sn_calculate_composition(
    x = meta_df,
    group_by = c("sample", "cell_type"),
    variable = "Phase",
    min_cells = 1,
    measure = "both"
  )

  expect_true(all(c("sample", "cell_type", "Phase", "count") %in% colnames(count_result)))
  expect_false("proportion" %in% colnames(count_result))
  expect_equal(
    count_result$count[count_result$sample == "A" & count_result$cell_type == "Tcell" & count_result$Phase == "S"],
    2
  )

  expect_true(all(c("count", "group_total", "proportion") %in% colnames(both_result)))
  expect_equal(
    both_result$group_total[both_result$sample == "A" & both_result$cell_type == "Tcell" & both_result$Phase == "G1"],
    3
  )
  expect_equal(
    both_result$proportion[both_result$sample == "A" & both_result$cell_type == "Tcell" & both_result$Phase == "S"],
    2 / 3 * 100
  )
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

test_that("sn_calculate_composition does not warn when varying columns are grouped explicitly", {
  meta_df <- data.frame(
    sample = c("A", "A", "A", "A"),
    cell_type = c("Tcell", "Tcell", "Bcell", "Bcell"),
    Phase = c("G1", "S", "G1", "G1"),
    group = c("ctrl", "ctrl", "ctrl", "ctrl"),
    stringsAsFactors = FALSE
  )

  expect_no_warning(
    result <- sn_calculate_composition(
      x = meta_df,
      group_by = c("sample", "cell_type"),
      variable = "Phase",
      min_cells = 1,
      additional_cols = "group"
    )
  )

  expect_true(all(c("sample", "cell_type", "Phase", "group") %in% colnames(result)))
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

test_that("sn_calculate_composition preserves factor columns and can sort groups", {
  meta_df <- data.frame(
    cell_type = factor(c("cDC1", "cDC1", "Mono", "Mono", "B", "B"), levels = c("Mono", "cDC1", "B")),
    Mutation = factor(c("WT", "PPP2R1A", "WT", "WT", "PPP2R1A", "PPP2R1A"), levels = c("WT", "PPP2R1A")),
    stringsAsFactors = TRUE
  )

  result <- sn_calculate_composition(
    x = meta_df,
    group_by = "cell_type",
    variable = "Mutation",
    min_cells = 1,
    sort_by = "proportion",
    sort_value = "WT"
  )

  expect_true(is.factor(result$cell_type))
  expect_true(is.factor(result$Mutation))
  expect_identical(levels(result$Mutation), c("WT", "PPP2R1A"))
  expect_identical(levels(result$cell_type), c("B", "cDC1", "Mono"))
})

test_that("sn_calculate_composition rejects sorting for multi-column grouping", {
  meta_df <- data.frame(
    sample = c("A", "A", "B", "B"),
    cell_type = c("Tcell", "Bcell", "Tcell", "Bcell"),
    Phase = c("G1", "S", "G1", "S"),
    stringsAsFactors = FALSE
  )

  expect_error(
    sn_calculate_composition(
      x = meta_df,
      group_by = c("sample", "cell_type"),
      variable = "Phase",
      min_cells = 1,
      sort_by = "proportion"
    ),
    "`sort_by` is only supported"
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
