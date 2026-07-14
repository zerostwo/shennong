library(testthat)

make_abundance_test_object <- function() {
  skip_if_not_installed("SeuratObject")
  set.seed(41)
  samples <- paste0("S", seq_len(8L))
  conditions <- rep(c("control", "treated"), each = 4L)
  cells_per_sample <- 40L
  sample <- rep(samples, each = cells_per_sample)
  condition <- rep(conditions, each = cells_per_sample)
  batch <- rep(rep(c("B1", "B2"), 4L), each = cells_per_sample)
  cell_type <- unlist(lapply(conditions, function(group) {
    probabilities <- if (group == "treated") c(T = 0.65, B = 0.2, Myeloid = 0.15) else c(T = 0.25, B = 0.45, Myeloid = 0.3)
    sample(names(probabilities), cells_per_sample, replace = TRUE, prob = probabilities)
  }), use.names = FALSE)
  counts <- matrix(stats::rpois(20L * length(sample), lambda = 3), nrow = 20L)
  counts[1:3, cell_type == "T"] <- counts[1:3, cell_type == "T"] + 5
  counts[4:6, cell_type == "B"] <- counts[4:6, cell_type == "B"] + 5
  rownames(counts) <- paste0("G", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  object$sample <- sample
  object$condition <- condition
  object$batch <- batch
  object$cell_type <- cell_type
  embedding <- cbind(
    PC_1 = seq_len(ncol(object)) / ncol(object) + stats::rnorm(ncol(object), 0, 0.05),
    PC_2 = as.numeric(factor(cell_type)) + stats::rnorm(ncol(object), 0, 0.1),
    PC_3 = stats::rnorm(ncol(object))
  )
  rownames(embedding) <- colnames(object)
  object[["pca"]] <- SeuratObject::CreateDimReducObject(embedding, key = "PC_", assay = "RNA")
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  Seurat::ScaleData(object, verbose = FALSE)
}

test_that("Propeller abundance uses biological samples and stores a unified result", {
  skip_if_not_installed("speckle")
  object <- make_abundance_test_object()
  object <- suppressMessages(sn_test_abundance(
    object,
    method = "propeller",
    sample_by = "sample",
    condition_by = "condition",
    cell_type_by = "cell_type",
    contrast = c("treated", "control"),
    design = "batch",
    store_name = "propeller"
  ))
  result <- sn_get_result(object, "differential_abundance", "propeller")

  expect_equal(result$method, "propeller")
  expect_equal(result$diagnostics$inferential_unit, "sample")
  expect_equal(nrow(result$tables$primary), 3L)
  expect_true(all(c("estimate", "log2_ratio", "p_value", "adjusted_p_value") %in% names(result$tables$primary)))
  expect_equal(nrow(result$tables$sample_proportions), 8L * 3L)
  expect_equal(sum(result$tables$sample_contributions$contribution[result$tables$sample_contributions$cell_type == "T"]),
               result$tables$primary$estimate[result$tables$primary$feature == "T"])
  expect_true(sn_validate_result(result, error = FALSE)$valid)
})

test_that("sample-label permutation retains its null and uncertainty", {
  object <- make_abundance_test_object()
  result <- sn_test_abundance(
    object,
    method = "permutation",
    sample_by = "sample",
    condition_by = "condition",
    cell_type_by = "cell_type",
    contrast = c("treated", "control"),
    permutations = 199L,
    seed = 9L,
    return_object = FALSE
  )

  expect_equal(nrow(result$tables$permutation_null), 199L * 3L)
  expect_true(all(result$tables$primary$p_value > 0 & result$tables$primary$p_value <= 1))
  expect_true(all(c("null_sd", "adjusted_p_value") %in% names(result$tables$primary)))
  expect_equal(result$provenance$random_seed, 9L)
})

test_that("abundance inputs reject cell-level pseudoreplication hazards", {
  object <- make_abundance_test_object()
  object$bad_condition <- object$condition
  object$bad_condition[1] <- "other"
  expect_error(
    sn_test_abundance(
      object, method = "permutation", sample_by = "sample",
      condition_by = "bad_condition", cell_type_by = "cell_type",
      contrast = c("treated", "control"), permutations = 99L
    ),
    "not constant within samples"
  )
  object$condition <- "control"
  expect_error(
    sn_test_abundance(
      object, method = "permutation", sample_by = "sample",
      condition_by = "condition", cell_type_by = "cell_type", permutations = 99L
    ),
    "exactly two levels"
  )
})

test_that("abundance effect plots render", {
  object <- make_abundance_test_object()
  result <- sn_test_abundance(
    object,
    method = "permutation",
    sample_by = "sample",
    condition_by = "condition",
    cell_type_by = "cell_type",
    contrast = c("treated", "control"),
    permutations = 99L,
    return_object = FALSE
  )
  plot <- sn_plot_abundance(result)
  expect_s3_class(plot, "ggplot")
  expect_silent(ggplot2::ggplotGrob(plot))
})

test_that("Milo is available through the stable abundance entry point", {
  skip_if_not_installed("miloR")
  object <- make_abundance_test_object()
  result <- suppressMessages(sn_test_abundance(
    object,
    method = "milo",
    sample_by = "sample",
    condition_by = "condition",
    cell_type_by = "cell_type",
    contrast = c("treated", "control"),
    backend_control = list(milo = list(k = 10L, d = 3L, prop = 0.2, verbose = FALSE)),
    return_object = FALSE
  ))
  expect_equal(result$backend, "miloR")
  expect_gt(nrow(result$tables$primary), 0L)
  expect_true(all(c("feature", "estimate", "p_value", "adjusted_p_value") %in% names(result$tables$primary)))
})

test_that("sample-aware Augur prioritization returns AUC, nulls, and contributions", {
  object <- make_abundance_test_object()
  result <- sn_prioritize_states(
    object,
    method = "augur",
    phenotype = "condition",
    sample_by = "sample",
    state_by = "cell_type",
    contrast = c("treated", "control"),
    features = paste0("G", 1:10),
    max_cells_per_state = 120L,
    permutations = 9L,
    seed = 19L,
    return_object = FALSE
  )
  expect_true(all(c("state", "auc", "priority_score", "uncertainty", "p_value", "adjusted_p_value") %in% names(result$tables$primary)))
  expect_equal(nrow(result$tables$permutation_null), 9L * nrow(result$tables$primary))
  expect_gt(nrow(result$tables$sample_contributions), 0L)
  expect_true(all(result$tables$primary$auc >= 0 & result$tables$primary$auc <= 1))
  plot <- sn_plot_state_priority(result)
  expect_s3_class(plot, "ggplot")
  expect_silent(ggplot2::ggplotGrob(plot))
})

test_that("RareQ discovers states before sample-level phenotype association", {
  skip_if_not_installed("RareQ")
  object <- make_abundance_test_object()
  result <- suppressMessages(sn_prioritize_states(
    object,
    method = "rareq",
    phenotype = "condition",
    sample_by = "sample",
    contrast = c("treated", "control"),
    reduction = "pca",
    dims = 1:3,
    backend_control = list(rareq = list(
      permutations = 99L,
      find_rare = list(k = 4L, Q_cut = 0.4, ratio = 0.3)
    )),
    return_object = FALSE
  ))
  expect_equal(result$backend, "RareQ")
  expect_gt(nrow(result$tables$cells), 0L)
  expect_true(all(c("state_fraction", "phenotype_association", "priority_score") %in% names(result$tables$primary)))
  expect_gt(nrow(result$tables$permutation_null), 0L)
})

test_that("Scissor enforces its bulk-expression phenotype contract", {
  object <- make_abundance_test_object()
  expect_error(
    sn_prioritize_states(
      object, method = "scissor", phenotype = "bulk_response",
      sample_by = "sample", state_by = "cell_type", return_object = FALSE
    ),
    "bulk_expression"
  )
})

test_that("Scissor runs against explicit bulk expression and phenotype", {
  skip_if_not_installed("Scissor")
  set.seed(8)
  genes <- paste0("G", 1:20)
  cells <- paste0("c", 1:60)
  bulk_samples <- paste0("b", 1:12)
  single_counts <- matrix(stats::rpois(20L * 60L, 3), nrow = 20L, dimnames = list(genes, cells))
  bulk <- matrix(stats::rpois(20L * 12L, 20), nrow = 20L, dimnames = list(genes, bulk_samples))
  bulk_phenotype <- rep(c(0, 1), each = 6L)
  bulk[1:4, bulk_phenotype == 1] <- bulk[1:4, bulk_phenotype == 1] + 20
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(single_counts, sparse = TRUE))
  object$state <- rep(c("A", "B"), each = 30L)
  object$sample <- rep(paste0("S", 1:6), each = 10L)
  object <- Seurat::NormalizeData(object, verbose = FALSE)

  result <- suppressMessages(sn_prioritize_states(
    object,
    method = "scissor",
    phenotype = "bulk_response",
    sample_by = "sample",
    state_by = "state",
    bulk_expression = bulk,
    bulk_phenotype = bulk_phenotype,
    family = "binomial",
    backend_control = list(scissor = list(alpha = 0.1, cutoff = 0.5, tag = c("control", "case"))),
    return_object = FALSE
  ))
  expect_equal(result$backend, "Scissor")
  expect_equal(nrow(result$tables$cells), 60L)
  expect_true(all(result$tables$cells$selection %in% c("Scissor+", "Scissor-", "Unselected")))
  expect_equal(sort(result$tables$primary$state), c("A", "B"))
  expect_true(length(result$warnings) >= 1L)
})
