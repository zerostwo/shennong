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
  withr::local_pdf(tempfile("shennong-scissor-test-", fileext = ".pdf"))
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
    backend_control = list(scissor = list(cutoff = 0.99, tag = c("control", "case"))),
    return_object = FALSE
  ))
  expect_equal(result$backend, "Scissor")
  expect_equal(nrow(result$tables$cells), 60L)
  expect_true(all(result$tables$cells$selection %in% c("Scissor+", "Scissor-", "Unselected")))
  expect_equal(sort(result$tables$primary$state), c("A", "B"))
  expect_equal(result$tables$states, result$tables$primary)
  expect_equal(nrow(result$tables$correlations), 60L)
  expect_equal(nrow(result$tables$model), 1L)
  expect_equal(result$tables$model$phenotype_alignment, "positional")
  expect_true(is.na(result$parameters$alpha))
  expect_true(result$tables$model$alpha_search_requested)
  expect_true(is.finite(result$tables$model$alpha))
  expect_true(all(c("sample", "state", "contribution") %in% colnames(result$tables$samples)))
  expect_true(is.data.frame(result$tables$reliability))
  expect_true(length(result$warnings) >= 1L)
})

make_scissor_adapter_inputs <- function() {
  object <- make_abundance_test_object()
  bulk_samples <- paste0("bulk", seq_len(8L))
  genes <- rownames(object)
  bulk <- matrix(
    stats::rpois(length(genes) * length(bulk_samples), lambda = 20),
    nrow = length(genes),
    dimnames = list(genes, bulk_samples)
  )
  phenotype <- setNames(rep(c("control", "case"), each = 4L), bulk_samples)
  coefficients <- rep(0, ncol(object))
  coefficients[seq_len(6L)] <- c(-0.4, -0.2, 0.3, 0.5, 0.7, -0.1)
  names(coefficients) <- colnames(object)
  correlations <- matrix(
    seq(-0.4, 0.6, length.out = length(bulk_samples) * ncol(object)),
    nrow = length(bulk_samples),
    dimnames = list(bulk_samples, colnames(object))
  )
  reliability <- data.frame(
    Coefficient = coefficients[seq_len(6L)],
    `Probability zero` = seq(0.02, 0.12, length.out = 6L),
    check.names = FALSE,
    row.names = colnames(object)[seq_len(6L)]
  )
  list(
    object = object,
    bulk = bulk,
    phenotype = phenotype,
    coefficients = coefficients,
    correlations = correlations,
    reliability = reliability
  )
}

test_that("direct Scissor entry point preserves cells and exposes unified tables", {
  inputs <- make_scissor_adapter_inputs()
  shuffled <- rev(names(inputs$phenotype))
  updated <- sn_run_scissor(
    inputs$object,
    bulk_expression = inputs$bulk,
    bulk_phenotype = inputs$phenotype[shuffled],
    state_by = "cell_type",
    sample_by = "sample",
    family = "binomial",
    backend_control = list(
      result = list(
        Coefs = inputs$coefficients,
        para = list(alpha = 0.05, lambda = 0.2, family = "binomial")
      ),
      correlation_matrix = inputs$correlations,
      reliability_result = inputs$reliability,
      retain_correlation_matrix = TRUE
    )
  )
  result <- sn_get_result(updated, "scissor", "scissor")

  expect_equal(result$method, "scissor")
  expect_equal(result$analysis_type, "scissor")
  expect_true(sn_validate_result(result, error = FALSE)$valid)
  expect_equal(nrow(result$tables$cells), ncol(inputs$object))
  expect_setequal(result$tables$cells$cell, colnames(inputs$object))
  expect_equal(result$tables$cells, result$tables$primary)
  expect_equal(nrow(result$tables$states), length(unique(inputs$object$cell_type)))
  expect_equal(nrow(result$tables$correlations), ncol(inputs$object))
  expect_equal(nrow(result$tables$correlation_matrix), length(inputs$correlations))
  expect_equal(nrow(result$tables$reliability), 6L)
  expect_true(all(c("cell", "coefficient", "probability_zero", "state", "selection") %in%
    colnames(result$tables$reliability)))
  expect_equal(result$tables$model$phenotype_alignment, "named")
  expect_equal(result$input$bulk_samples, colnames(inputs$bulk))
  expect_equal(result$diagnostics$selected_cell_count, 6L)

  for (type in c("states", "cells", "samples", "correlations", "reliability")) {
    plot <- sn_plot_scissor(result, type = type, n = 100L)
    expect_s3_class(plot, "ggplot")
    expect_s3_class(sn_figure_spec(plot), "sn_figure_spec")
    expect_false(is.null(attr(plot, "shennong_figure_data", exact = TRUE)))
    expect_silent(ggplot2::ggplotGrob(plot))
  }
  expect_s3_class(sn_plot_scissor(updated, type = "states"), "ggplot")
})

test_that("direct Scissor defaults remain cell-first and record upstream-selected alpha", {
  inputs <- make_scissor_adapter_inputs()
  result <- sn_run_scissor(
    inputs$object,
    bulk_expression = inputs$bulk,
    bulk_phenotype = inputs$phenotype,
    seed = 29L,
    backend_control = list(
      result = list(
        Coefs = inputs$coefficients,
        para = list(alpha = 0.3, lambda = 0.1, family = "binomial")
      )
    ),
    return_object = FALSE
  )

  expect_equal(result$tables$states$state, "all_cells")
  expect_false(result$diagnostics$sample_aware)
  expect_equal(result$diagnostics$inferential_unit, "bulk_sample")
  expect_false(result$diagnostics$single_cell_sample_summary)
  expect_true(is.na(result$parameters$alpha))
  expect_equal(result$tables$model$alpha, 0.3)
  expect_true(result$tables$model$alpha_search_requested)
  expect_true(result$diagnostics$selection_cutoff_satisfied)
  expect_equal(
    result$provenance$random_seed,
    list(preprocessing = 29L, scissor_backend = 123L)
  )
})

test_that("Scissor reports when its selected fraction misses the configured cutoff", {
  inputs <- make_scissor_adapter_inputs()
  coefficients <- rep(0.1, ncol(inputs$object))
  names(coefficients) <- colnames(inputs$object)
  result <- sn_run_scissor(
    inputs$object,
    bulk_expression = inputs$bulk,
    bulk_phenotype = inputs$phenotype,
    backend_control = list(
      result = list(Coefs = coefficients, para = list(alpha = 0.9, lambda = 0.1)),
      cutoff = 0.2
    ),
    return_object = FALSE
  )

  expect_equal(result$diagnostics$selected_fraction, 1)
  expect_false(result$diagnostics$selection_cutoff_satisfied)
  expect_false(result$tables$model$selection_cutoff_satisfied)
  expect_match(paste(result$warnings, collapse = " "), "did not fall below")
})

test_that("Scissor strictly validates sample alignment and phenotype family contracts", {
  inputs <- make_scissor_adapter_inputs()
  adapter <- list(result = list(Coefs = inputs$coefficients, para = list(alpha = 0.05)))

  mismatched <- inputs$phenotype
  names(mismatched)[[1L]] <- "not-a-bulk-sample"
  expect_error(
    sn_run_scissor(
      inputs$object, inputs$bulk, mismatched,
      state_by = "cell_type", backend_control = adapter, return_object = FALSE
    ),
    "match `bulk_expression` sample names exactly"
  )

  partial <- unname(inputs$phenotype)
  names(partial) <- c(colnames(inputs$bulk)[-1L], "")
  expect_error(
    sn_run_scissor(
      inputs$object, inputs$bulk, partial,
      state_by = "cell_type", backend_control = adapter, return_object = FALSE
    ),
    "unique, non-empty sample names"
  )

  three_classes <- setNames(rep(c("a", "b", "c", "c"), each = 2L), colnames(inputs$bulk))
  expect_error(
    sn_run_scissor(
      inputs$object, inputs$bulk, three_classes,
      state_by = "cell_type", backend_control = adapter, return_object = FALSE
    ),
    "exactly two classes"
  )

  constant <- setNames(rep(1, ncol(inputs$bulk)), colnames(inputs$bulk))
  expect_error(
    sn_run_scissor(
      inputs$object, inputs$bulk, constant,
      state_by = "cell_type", family = "gaussian",
      backend_control = adapter, return_object = FALSE
    ),
    "must vary"
  )

  cox <- cbind(time = seq_len(ncol(inputs$bulk)), event = 1)
  rownames(cox) <- colnames(inputs$bulk)
  expect_error(
    sn_run_scissor(
      inputs$object, inputs$bulk, cox,
      state_by = "cell_type", family = "cox",
      backend_control = adapter, return_object = FALSE
    ),
    "both 0 .* and 1"
  )

  numeric_ids <- c("2", "1")
  numeric_id_phenotype <- matrix(
    c(10, 0, 20, 1),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("1", "2"), c("time", "event"))
  )
  aligned_numeric_ids <- Shennong:::.sn_scissor_align_phenotype(
    numeric_id_phenotype,
    bulk_samples = numeric_ids,
    family = "cox",
    control = list()
  )
  expect_identical(aligned_numeric_ids$mode, "named")
  expect_identical(rownames(aligned_numeric_ids$values), numeric_ids)
  expect_identical(
    aligned_numeric_ids$values[, "time"],
    stats::setNames(c(20, 10), numeric_ids)
  )

  expect_error(
    sn_run_scissor(
      inputs$object, inputs$bulk, inputs$phenotype,
      seed = NA_integer_, backend_control = adapter, return_object = FALSE
    ),
    "one non-negative integer"
  )
  expect_error(
    sn_plot_scissor(
      sn_run_scissor(
        inputs$object, inputs$bulk, inputs$phenotype,
        backend_control = adapter, return_object = FALSE
      ),
      n = Inf
    ),
    "one positive integer"
  )
})

test_that("Scissor preprocessing uses the selected expression layer exactly", {
  skip_if_not_installed("Scissor")
  set.seed(91)
  genes <- paste0("LayerGene", seq_len(12L))
  cells <- paste0("LayerCell", seq_len(10L))
  counts <- matrix(
    stats::rpois(length(genes) * length(cells), lambda = 4),
    nrow = length(genes),
    dimnames = list(genes, cells)
  )
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  custom <- log1p(counts + matrix(
    seq_len(length(counts)) / length(counts),
    nrow = nrow(counts)
  ))
  dimnames(custom) <- dimnames(counts)
  SeuratObject::LayerData(object[["RNA"]], layer = "scissor_custom") <- custom

  bulk_samples <- paste0("bulk_", seq_len(8L))
  bulk <- matrix(
    stats::rpois(length(genes) * length(bulk_samples), lambda = 20),
    nrow = length(genes),
    dimnames = list(genes, bulk_samples)
  )
  phenotype <- stats::setNames(rep(c("control", "case"), each = 4L), bulk_samples)
  captured <- new.env(parent = emptyenv())
  runner <- function(bulk_dataset, sc_dataset, phenotype, tag, alpha, cutoff,
                     family, Save_file) {
    captured$expression <- as.matrix(
      SeuratObject::LayerData(sc_dataset[["RNA"]], layer = "data")
    )
    list(
      Coefs = stats::setNames(rep(0, ncol(sc_dataset)), colnames(sc_dataset)),
      para = list(alpha = 0.1, lambda = 0.2, family = family)
    )
  }

  result <- suppressWarnings(suppressMessages(sn_run_scissor(
    object,
    bulk_expression = bulk,
    bulk_phenotype = phenotype,
    assay = "RNA",
    layer = "scissor_custom",
    backend_control = list(
      runner = runner,
      min_shared_features = 10L,
      nfeatures = 10L,
      npcs = 5L
    ),
    return_object = FALSE
  )))

  expect_equal(captured$expression, custom, tolerance = 0)
  expect_identical(result$input$assay, "RNA")
  expect_identical(result$input$layer, "scissor_custom")
  expect_error(
    sn_run_scissor(
      object,
      bulk_expression = bulk,
      bulk_phenotype = phenotype,
      layer = "missing_layer",
      backend_control = list(
        result = list(Coefs = rep(0, ncol(object)), para = list(alpha = 0.1))
      ),
      return_object = FALSE
    ),
    "Layer 'missing_layer' was not found"
  )
})
