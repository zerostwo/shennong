library(testthat)

make_qc_assessment_object <- function() {
  skip_if_not_installed("Seurat")

  counts <- Matrix::Matrix(
    matrix(rpois(40 * 8, lambda = 3), nrow = 40, ncol = 8),
    sparse = TRUE
  )
  rownames(counts) <- paste0("gene", seq_len(40))
  colnames(counts) <- paste0("cell", seq_len(8))

  object <- SeuratObject::CreateSeuratObject(
    counts = counts,
    project = "qc-assessment"
  )
  object$sample <- rep(c("s1", "s2"), each = 4)
  object$percent.mt <- c(4, 5, 18, 20, 3, 4, 8, 9)
  object$nCount_RNA_qc <- c("Passed", "Passed", "Failed", "Failed", "Passed", "Passed", "Passed", "Passed")
  object$nFeature_RNA_qc <- c("Passed", "Passed", "Failed", "Failed", "Passed", "Passed", "Passed", "Passed")
  object$scDblFinder.class <- c("singlet", "singlet", "doublet", "doublet", "singlet", "singlet", "singlet", "doublet")
  object$decontaminated_counts_zero_count <- c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE)
  object
}

test_that("sn_assess_qc summarizes current QC status by sample", {
  object <- make_qc_assessment_object()

  report <- sn_assess_qc(object, sample_by = "sample", verbose = FALSE)

  expect_type(report, "list")
  expect_true(all(c("overall", "by_sample", "messages") %in% names(report)))
  expect_equal(nrow(report$by_sample), 2)
  expect_true(all(c(
    "sample", "n_cells", "median_nCount", "median_nFeature",
    "failed_qc_fraction", "doublet_fraction", "zero_count_fraction",
    "qc_score", "qc_label"
  ) %in% colnames(report$by_sample)))
  expect_true(is.character(report$messages))
  expect_equal(report$sample_col, "sample")
})

test_that("sn_assess_qc compares filtered objects to a reference and stores reports", {
  reference <- make_qc_assessment_object()
  current <- reference[, colnames(reference)[-c(3, 4, 8)]]

  report <- sn_assess_qc(
    object = current,
    reference = reference,
    sample_by = "sample",
    verbose = FALSE
  )

  expect_true("comparison" %in% names(report))
  expect_true(all(c(
    "retention_fraction",
    "low_quality_removed_fraction",
    "doublet_removed_fraction",
    "clean_retained_fraction"
  ) %in% colnames(report$by_sample)))
  expect_true(any(report$by_sample$low_quality_removed_fraction > 0))
  expect_true(any(report$by_sample$doublet_removed_fraction > 0))

  stored <- sn_assess_qc(
    object = current,
    reference = reference,
    sample_by = "sample",
    store_name = "post_filter",
    return_object = TRUE,
    verbose = FALSE
  )

  expect_true("qc_assessments" %in% names(stored@misc))
  expect_true("post_filter" %in% names(stored@misc$qc_assessments))
})
