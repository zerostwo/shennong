library(testthat)

make_program_test_object <- function() {
  skip_if_not_installed("SeuratObject")
  set.seed(11)
  counts <- matrix(rpois(8 * 24, lambda = 2), nrow = 8)
  rownames(counts) <- paste0("G", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  counts[c("G1", "G2", "G3"), 1:12] <- counts[c("G1", "G2", "G3"), 1:12] + 8
  counts[c("G4", "G5", "G6"), 13:24] <- counts[c("G4", "G5", "G6"), 13:24] + 8
  object <- SeuratObject::CreateSeuratObject(counts = Matrix::Matrix(counts, sparse = TRUE), project = "program-test")
  object$sample <- rep(paste0("S", 1:4), each = 6)
  object$condition <- rep(c("control", "control", "treated", "treated"), each = 6)
  object$cell_type <- rep(c("type_a", "type_b"), 12)
  Seurat::NormalizeData(object, verbose = FALSE)
}

program_test_signatures <- function() {
  list(program_a = c("G1", "G2", "G3", "missing"), program_b = c("G4", "G5", "G6"))
}

test_that("mean scoring is sparse-aware, stored, and added to metadata", {
  object <- make_program_test_object()
  object <- sn_score_programs(
    object,
    signatures = program_test_signatures(),
    method = "mean",
    name = "programs",
    min_genes = 2
  )

  expect_true(all(c("programs_program_a", "programs_program_b") %in% colnames(object[[]])))
  result <- sn_get_result(object, "program_scoring", "programs")
  expect_equal(result$method, "mean")
  expect_equal(nrow(result$tables$scores), 2 * ncol(object))
  expect_equal(unique(result$tables$scores$level), "cell")
  expect_equal(result$tables$coverage$n_matched, c(3L, 3L))
  expect_match(result$tables$coverage$missing_genes[[1]], "missing")
  expect_gt(mean(object$programs_program_a[object$condition == "control"]), mean(object$programs_program_a[object$condition == "treated"]))
})

test_that("signature data frames and feature coverage are validated", {
  object <- make_program_test_object()
  signatures <- data.frame(
    program = c("a", "a", "drop"),
    gene = c("G1", "G2", "not_present")
  )
  result <- sn_score_programs(
    object,
    signatures = signatures,
    method = "mean",
    min_genes = 2,
    return_object = FALSE
  )
  expect_equal(result$diagnostics$retained_programs, "a")
  expect_equal(result$diagnostics$dropped_programs, "drop")
  expect_error(
    sn_score_programs(object, signatures = list(x = "not_present"), method = "mean"),
    "No signature retained"
  )
})

test_that("UCell and AUCell backends return cell-level score contracts", {
  object <- make_program_test_object()
  if (requireNamespace("UCell", quietly = TRUE)) {
    ucell <- sn_score_programs(
      object,
      signatures = program_test_signatures(),
      method = "ucell",
      name = "ucell",
      return_object = FALSE
    )
    expect_equal(nrow(ucell$tables$scores), 2 * ncol(object))
    expect_true(all(is.finite(ucell$tables$scores$score)))
  }
  if (requireNamespace("AUCell", quietly = TRUE)) {
    aucell <- sn_score_programs(
      object,
      signatures = program_test_signatures(),
      method = "aucell",
      name = "aucell",
      return_object = FALSE
    )
    expect_equal(nrow(aucell$tables$scores), 2 * ncol(object))
    expect_true(all(is.finite(aucell$tables$scores$score)))
  }
})

test_that("GSVA and ssGSEA score aggregated sample expression", {
  skip_if_not_installed("GSVA")
  object <- make_program_test_object()
  for (method in c("gsva", "ssgsea")) {
    result <- suppressWarnings(sn_score_programs(
      object,
      signatures = program_test_signatures(),
      method = method,
      group_by = "sample",
      name = paste0("sample_", method),
      return_object = FALSE,
      backend_control = stats::setNames(list(list(parameters = list(kcdf = "Gaussian"))), method)
    ))
    expect_equal(unique(result$tables$scores$level), "group")
    expect_equal(length(unique(result$tables$scores$entity)), 4L)
    expect_equal(nrow(result$tables$scores), 8L)
  }
})

test_that("program tests preserve sample as the inferential unit", {
  object <- make_program_test_object()
  object <- sn_score_programs(object, program_test_signatures(), method = "mean", name = "programs")
  comparison <- sn_test_programs(
    object,
    score_name = "programs",
    condition_by = "condition",
    sample_by = "sample",
    group_by = "cell_type",
    contrast = c("control", "treated"),
    method = "wilcox",
    store_name = "program_condition",
    return_object = FALSE
  )

  expect_equal(comparison$analysis_type, "program_comparison")
  expect_equal(comparison$diagnostics$inferential_unit, "sample")
  expect_true(all(c("estimate", "p_value", "adjusted_p_value", "n_1", "n_2") %in% colnames(comparison$tables$primary)))
  expect_true(all(comparison$tables$primary$n_1 == 2L))
  expect_true(all(comparison$tables$primary$n_2 == 2L))
})

test_that("program activity and heatmap plots render", {
  object <- make_program_test_object()
  object <- sn_score_programs(object, program_test_signatures(), method = "mean", name = "programs")
  activity <- sn_plot_program_activity(object, "programs", group_by = "condition")
  heatmap <- sn_plot_program_heatmap(object, "programs", group_by = "condition")
  expect_s3_class(activity, "ggplot")
  expect_s3_class(heatmap, "ggplot")
  expect_silent(ggplot2::ggplotGrob(activity))
  expect_silent(ggplot2::ggplotGrob(heatmap))
})
