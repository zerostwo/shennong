bulk_fixture <- function(features = 120L, samples = 12L) {
  set.seed(412)
  sample_names <- paste0("sample_", seq_len(samples))
  metadata <- data.frame(
    condition = factor(rep(c("normal", "tumor"), each = samples / 2), levels = c("normal", "tumor")),
    batch = factor(rep(seq_len(3), length.out = samples)),
    patient = factor(rep(seq_len(samples / 2), each = 2)),
    age = seq(42, length.out = samples),
    stage = factor(rep(c("I", "II", "III"), length.out = samples)),
    time = seq(12, by = 3, length.out = samples),
    event = rep(c(1, 0, 1), length.out = samples),
    row.names = sample_names
  )
  counts <- matrix(
    stats::rpois(features * samples, lambda = 30), nrow = features,
    dimnames = list(paste0("gene_", seq_len(features)), sample_names)
  )
  counts[seq_len(12), metadata$condition == "tumor"] <- counts[seq_len(12), metadata$condition == "tumor"] + 45L
  list(counts = counts, metadata = metadata)
}

survival_fixture <- function(samples = 48L) {
  set.seed(20260715)
  sample_names <- paste0("survival_", seq_len(samples))
  risk_score <- stats::rnorm(samples)
  age_years <- stats::rnorm(samples, mean = 60, sd = 8)
  event_time <- stats::rexp(samples, rate = 0.04 * exp(0.75 * risk_score))
  censor_time <- stats::rexp(samples, rate = 0.015)
  metadata <- data.frame(
    time = pmin(event_time, censor_time),
    event = as.integer(event_time <= censor_time),
    risk_score = risk_score,
    constant = 1,
    row.names = sample_names,
    check.names = FALSE
  )
  metadata[["age years"]] <- age_years
  counts <- matrix(
    stats::rpois(30L * samples, lambda = 35),
    nrow = 30L,
    dimnames = list(paste0("survival_gene_", seq_len(30L)), sample_names)
  )
  list(counts = counts, metadata = metadata)
}

test_that("bulk input and QC expose aligned sample diagnostics", {
  fixture <- bulk_fixture()
  result <- sn_assess_bulk_qc(fixture$counts, fixture$metadata, top_variable = 80)
  expect_true(sn_validate_result(result)$valid)
  expect_identical(result$analysis_type, "bulk_qc")
  expect_equal(result$tables$samples$sample, colnames(fixture$counts))
  expect_equal(dim(result$embeddings$pca)[1], ncol(fixture$counts))
  expect_equal(dim(result$tables$correlation), c(12L, 12L))
  expect_true(all(c("library_size", "detected_features", "outlier") %in% names(result$tables$samples)))
  unnamed_samples <- fixture$counts
  colnames(unnamed_samples) <- NULL
  expect_error(sn_assess_bulk_qc(unnamed_samples, fixture$metadata), "sample column names")
})

test_that("SummarizedExperiment bulk input is supported", {
  skip_if_not_installed("SummarizedExperiment")
  fixture <- bulk_fixture()
  object <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = fixture$counts),
    colData = S4Vectors::DataFrame(fixture$metadata)
  )
  result <- sn_assess_bulk_qc(object, assay = "counts")
  expect_true(sn_validate_result(result)$valid)
  expect_identical(result$input$source, "SummarizedExperiment")
})

test_that("bulk DE validates design and returns standardized edgeR results", {
  skip_if_not_installed("edgeR")
  fixture <- bulk_fixture()
  result <- sn_find_bulk_de(
    fixture$counts, fixture$metadata,
    design = ~batch + condition,
    contrast = c("condition", "tumor", "normal"), method = "auto"
  )
  expect_true(sn_validate_result(result)$valid)
  expect_identical(result$method, "edger")
  expect_true(all(c("gene", "log2_fold_change", "p_value", "adjusted_p_value") %in% names(result$tables$primary)))
  expect_gt(stats::median(result$tables$primary$log2_fold_change[result$tables$primary$gene %in% paste0("gene_", 1:12)]), 0)
  expect_error(
    sn_find_bulk_de(fixture$counts, fixture$metadata, design = ~missing, contrast = c("condition", "tumor", "normal")),
    "Design variable"
  )
  expect_error(
    sn_find_bulk_de(fixture$counts, fixture$metadata, design = ~condition, contrast = c("condition", "other", "normal")),
    "Contrast level"
  )
})

test_that("sn_find_de dispatches bulk and preserves the legacy wrapper", {
  fixture <- bulk_fixture(features = 20L, samples = 6L)
  backend <- list(result = data.frame(
    gene = rownames(fixture$counts),
    logFC = seq(-1, 1, length.out = nrow(fixture$counts)),
    PValue = seq(0.001, 0.2, length.out = nrow(fixture$counts)),
    FDR = seq(0.01, 0.5, length.out = nrow(fixture$counts))
  ))

  unified <- sn_find_de(
    fixture$counts,
    metadata = fixture$metadata,
    design = ~condition,
    contrast = c("condition", "tumor", "normal"),
    backend_control = backend
  )
  legacy <- sn_find_bulk_de(
    fixture$counts,
    metadata = fixture$metadata,
    design = ~condition,
    contrast = c("condition", "tumor", "normal"),
    backend_control = backend
  )

  expect_true(sn_validate_result(unified)$valid)
  expect_identical(unified$analysis_type, "bulk_de")
  expect_identical(unified$name, "bulk_de")
  expect_equal(unified$tables$primary, legacy$tables$primary)
  expect_error(
    sn_find_de(fixture$counts, metadata = fixture$metadata),
    "contrast.*required"
  )
  expect_error(
    sn_find_de(fixture$counts, modality = "single_cell"),
    "requires a Seurat object"
  )
})

test_that("limma supports normalized continuous bulk expression", {
  skip_if_not_installed("limma")
  fixture <- bulk_fixture()
  expression <- log2(fixture$counts + 1) + 0.001
  result <- sn_find_bulk_de(
    expression, fixture$metadata, design = ~batch + condition,
    contrast = c("condition", "tumor", "normal"), method = "auto"
  )
  expect_identical(result$method, "limma")
  expect_equal(nrow(result$tables$primary), nrow(expression))
  expect_true(sn_validate_result(result)$valid)
})

test_that("DESeq2 shrinkage and dream repeated measures run through real backends", {
  skip_if_not_installed("DESeq2")
  skip_if_not_installed("variancePartition")
  set.seed(908)
  metadata <- data.frame(
    condition = factor(rep(c("normal", "tumor"), 6), levels = c("normal", "tumor")),
    patient = factor(rep(seq_len(6), each = 2)),
    row.names = paste0("paired_", seq_len(12))
  )
  counts <- matrix(
    stats::rpois(150 * 12, 30), nrow = 150,
    dimnames = list(paste0("paired_gene_", seq_len(150)), rownames(metadata))
  )
  counts[seq_len(15), metadata$condition == "tumor"] <- counts[seq_len(15), metadata$condition == "tumor"] + 25L
  deseq <- suppressWarnings(suppressMessages(sn_find_bulk_de(
    counts, metadata, design = ~patient + condition,
    contrast = c("condition", "tumor", "normal"), method = "deseq2"
  )))
  dream <- suppressWarnings(sn_find_bulk_de(
    counts, metadata, design = ~condition + (1 | patient),
    contrast = c("condition", "tumor", "normal"), method = "dream"
  ))
  expect_true(sn_validate_result(deseq)$valid)
  expect_true(isTRUE(deseq$diagnostics$shrink_applied))
  expect_true(sn_validate_result(dream)$valid)
  expect_identical(dream$diagnostics$coefficient, "conditiontumor")
})

test_that("bulk pathway scoring records coverage and sample scores", {
  fixture <- bulk_fixture()
  signatures <- list(inflammatory = paste0("gene_", 1:10), stromal = paste0("gene_", 20:30))
  result <- sn_score_bulk_pathways(fixture$counts, signatures, min_genes = 3)
  expect_true(sn_validate_result(result)$valid)
  expect_equal(nrow(result$tables$scores), 2L * ncol(fixture$counts))
  expect_equal(sort(unique(result$tables$scores$pathway)), sort(names(signatures)))
  expect_true(all(result$tables$coverage$n_matched >= 3))
})

test_that("WGCNA adapters standardize modules and trait associations", {
  skip_if_not_installed("WGCNA")
  fixture <- bulk_fixture(features = 40)
  colors <- rep(c("blue", "brown"), each = 20)
  names(colors) <- rownames(fixture$counts)
  eigengenes <- data.frame(
    MEblue = scale(seq_len(12))[, 1],
    MEbrown = scale(rev(seq_len(12)))[, 1],
    row.names = colnames(fixture$counts)
  )
  result <- sn_run_wgcna(
    fixture$counts, fixture$metadata, traits = c("age", "condition"), power = 6,
    min_module_size = 5,
    backend_control = list(result = list(colors = colors, eigengenes = eigengenes))
  )
  expect_true(sn_validate_result(result)$valid)
  expect_equal(nrow(result$tables$modules), 40L)
  expect_gt(nrow(result$tables$trait_associations), 0L)
  expect_equal(result$parameters$power, 6)
})

test_that("WGCNA constructs a real network on variable expression", {
  skip_if_not_installed("WGCNA")
  fixture <- bulk_fixture(features = 80)
  result <- suppressMessages(sn_run_wgcna(
    fixture$counts, fixture$metadata, traits = "age", power = 6,
    min_module_size = 5, backend_control = list(blockwise = list(maxBlockSize = 100))
  ))
  expect_true(sn_validate_result(result)$valid)
  expect_equal(nrow(result$tables$modules), 80L)
  expect_true(all(c("gene", "module") %in% names(result$tables$modules)))
})

test_that("survival and clinical association use sample-level features", {
  skip_if_not_installed("survival")
  fixture <- bulk_fixture()
  survival <- suppressWarnings(sn_run_survival(
    fixture$counts, time = "time", event = "event",
    features = c("gene_1", "age"), covariates = "batch", metadata = fixture$metadata
  ))
  clinical <- sn_run_clinical_association(
    fixture$counts, features = c("gene_1", "gene_2"),
    clinical_vars = c("age", "stage"), metadata = fixture$metadata
  )
  expect_true(sn_validate_result(survival)$valid)
  expect_identical(survival$schema_version, "1.0.0")
  expect_equal(nrow(survival$tables$survival), 2L)
  expect_true(all(c("hazard_ratio", "conf_low", "conf_high", "ph_p_value", "status", "error") %in% names(survival$tables$survival)))
  expect_true(all(c(
    "curves", "risk", "cumulative_hazard", "groups", "log_rank",
    "proportional_hazards", "proportional_hazards_residuals", "performance"
  ) %in% names(survival$tables)))
  expect_true(sn_validate_result(clinical)$valid)
  expect_equal(nrow(clinical$tables$associations), 4L)
})

test_that("survival validates time and event contracts before fitting", {
  skip_if_not_installed("survival")
  fixture <- survival_fixture()

  negative <- fixture$metadata
  negative$time[[1]] <- -1
  expect_error(
    sn_run_survival(fixture$counts, "time", "event", "risk_score", metadata = negative),
    "cannot be negative"
  )

  invalid_event <- fixture$metadata
  invalid_event$event[[1]] <- 2
  expect_error(
    sn_run_survival(fixture$counts, "time", "event", "risk_score", metadata = invalid_event),
    "only 0.*1"
  )

  nonnumeric_time <- fixture$metadata
  nonnumeric_time$time <- as.character(nonnumeric_time$time)
  expect_error(
    sn_run_survival(fixture$counts, "time", "event", "risk_score", metadata = nonnumeric_time),
    "time must be numeric"
  )

  no_events <- fixture$metadata
  no_events$event <- 0L
  expect_error(
    sn_run_survival(fixture$counts, "time", "event", "risk_score", metadata = no_events),
    "observed event"
  )
  expect_error(
    sn_run_survival(
      fixture$counts, "time", "event", "risk_score", metadata = fixture$metadata,
      group_method = "fixed"
    ),
    "group_cutpoint"
  )
  expect_error(
    sn_run_survival(
      fixture$counts, "time", "event", "event", metadata = fixture$metadata
    ),
    "endpoints cannot also be used"
  )
  expect_error(
    sn_run_survival(
      fixture$counts, "time", "event", "risk_score", covariates = "time",
      metadata = fixture$metadata
    ),
    "endpoints cannot also be used"
  )
  expect_error(
    sn_run_survival(
      fixture$counts, "time", "event", "risk_score", covariates = "risk_score",
      metadata = fixture$metadata
    ),
    "features and covariates must be distinct"
  )
  expect_error(
    sn_run_survival(
      fixture$counts, "time", "event", "risk_score", metadata = fixture$metadata,
      risk_times = c(0, -1)
    ),
    "risk_times"
  )
})

test_that("survival isolates failed features and supports non-syntactic covariates", {
  skip_if_not_installed("survival")
  fixture <- survival_fixture()
  result <- sn_run_survival(
    fixture$counts,
    time = "time",
    event = "event",
    features = c("risk_score", "constant"),
    covariates = "age years",
    metadata = fixture$metadata
  )

  successful <- result$tables$survival[result$tables$survival$feature == "risk_score", , drop = FALSE]
  failed <- result$tables$survival[result$tables$survival$feature == "constant", , drop = FALSE]
  expect_identical(successful$status, "ok")
  expect_true(is.finite(successful$hazard_ratio))
  expect_identical(failed$status, "error")
  expect_match(failed$error, "constant")
  expect_equal(result$diagnostics$requested_models, 2L)
  expect_equal(result$diagnostics$failed_models, 1L)
  expect_equal(nrow(result$tables$performance), 2L)
})

test_that("survival returns deterministic grouping, curves, diagnostics, and performance", {
  skip_if_not_installed("survival")
  fixture <- survival_fixture()
  median_result <- sn_run_survival(
    fixture$counts, "time", "event", "risk_score", metadata = fixture$metadata
  )
  expected_median <- stats::median(fixture$metadata$risk_score)
  expect_equal(unique(median_result$tables$groups$cutpoint), expected_median)
  expect_setequal(unique(median_result$tables$groups$group), c("Low", "High"))
  expect_gt(nrow(median_result$tables$curves), 0L)
  risk_grid <- unique(median_result$tables$risk$time)
  expect_lte(length(risk_grid), 8L)
  expect_equal(nrow(median_result$tables$risk), 2L * length(risk_grid))
  expect_lt(nrow(median_result$tables$risk), nrow(median_result$tables$curves))
  expect_equal(nrow(median_result$tables$cumulative_hazard), nrow(median_result$tables$curves))
  expect_true(any(median_result$tables$proportional_hazards$scope == "feature"))
  expect_true(any(median_result$tables$proportional_hazards$scope == "global"))
  expect_true(all(c(
    "feature", "term", "time", "transformed_time", "scaled_schoenfeld"
  ) %in% names(median_result$tables$proportional_hazards_residuals)))
  expect_gt(nrow(median_result$tables$proportional_hazards_residuals), 0L)
  expect_true(is.finite(median_result$tables$performance$concordance))
  expect_identical(median_result$tables$log_rank$status, "ok")

  quantile_result <- sn_run_survival(
    fixture$counts, "time", "event", "risk_score", metadata = fixture$metadata,
    group_method = "quantile", group_quantile = 0.25
  )
  expect_equal(
    unique(quantile_result$tables$groups$cutpoint),
    as.numeric(stats::quantile(fixture$metadata$risk_score, 0.25, type = 7))
  )

  fixed_result <- sn_run_survival(
    fixture$counts, "time", "event", "risk_score", metadata = fixture$metadata,
    group_method = "fixed", group_cutpoint = 0
  )
  expect_equal(unique(fixed_result$tables$groups$cutpoint), 0)

  explicit_risk <- sn_run_survival(
    fixture$counts, "time", "event", "risk_score", metadata = fixture$metadata,
    risk_times = c(0, 10, 20)
  )
  expect_setequal(unique(explicit_risk$tables$risk$time), c(0, 10, 20))
})

test_that("survival grouping uses only endpoint-complete finite feature values", {
  skip_if_not_installed("survival")
  fixture <- survival_fixture()
  metadata <- fixture$metadata
  endpoint_missing <- order(metadata$risk_score, decreasing = TRUE)[seq_len(10L)]
  metadata$time[endpoint_missing] <- NA_real_
  metadata$event[endpoint_missing] <- NA_integer_
  complete <- !is.na(metadata$time) & !is.na(metadata$event)

  result <- sn_run_survival(
    fixture$counts, "time", "event", "risk_score", metadata = metadata
  )
  expect_equal(
    unique(result$tables$groups$cutpoint),
    stats::median(metadata$risk_score[complete])
  )
  expect_false(any(result$tables$groups$included[endpoint_missing]))
  expect_identical(result$tables$log_rank$status, "ok")

  nonfinite <- fixture$metadata
  nonfinite$risk_score[[1]] <- Inf
  finite_result <- sn_run_survival(
    fixture$counts, "time", "event", "risk_score", metadata = nonfinite
  )
  first <- finite_result$tables$groups$sample == rownames(nonfinite)[[1]]
  expect_false(finite_result$tables$groups$included[first])
  expect_true(is.na(finite_result$tables$groups$group[first]))
  expect_equal(finite_result$tables$survival$n, nrow(nonfinite) - 1L)
  expect_equal(finite_result$tables$survival$n_dropped, 1L)
  expect_equal(finite_result$tables$log_rank$n, nrow(nonfinite) - 1L)
})

test_that("survival plots expose forest, KM, risk, PH, and cumulative-hazard views", {
  skip_if_not_installed("survival")
  fixture <- survival_fixture()
  result <- sn_run_survival(
    fixture$counts, "time", "event", "risk_score", metadata = fixture$metadata
  )
  views <- c("forest", "km", "risk_table", "ph", "ph_test", "cumulative_hazard")
  plots <- lapply(views, function(view) sn_plot_survival(result, view = view, feature = "risk_score"))
  expect_true(all(vapply(plots, inherits, logical(1), "ggplot")))
  expect_true(all(vapply(plots, function(plot) !is.null(attr(plot, "shennong_figure_spec")), logical(1))))
  expect_silent(lapply(plots, ggplot2::ggplot_build))
  expect_gte(length(ggplot2::ggplot_build(plots[[4]])$data), 3L)
  expect_s3_class(sn_plot_survival(result), "ggplot")
  expect_error(sn_plot_survival(result, feature = "not_a_feature"), "not found")
})

test_that("Kaplan-Meier confidence bands use right-continuous step coordinates", {
  curves <- tibble::tibble(
    feature = rep("risk_score", 3L),
    group = rep("Low", 3L),
    time = c(0, 2, 5),
    conf_low = c(1, 0.7, 0.4),
    conf_high = c(1, 0.9, 0.8)
  )

  stepped <- Shennong:::.sn_survival_step_band(curves)

  expect_equal(stepped$time, c(0, 2, 2, 5, 5))
  expect_equal(stepped$conf_low, c(1, 1, 0.7, 0.7, 0.4))
  expect_equal(stepped$conf_high, c(1, 1, 0.9, 0.9, 0.8))
})

test_that("survival forest plots remain compatible with legacy result tables", {
  skip_if_not_installed("survival")
  fixture <- survival_fixture()
  result <- sn_run_survival(
    fixture$counts, "time", "event", "risk_score", metadata = fixture$metadata
  )
  legacy_table <- result$tables$survival
  legacy_table$status <- NULL
  legacy_table$error <- NULL
  legacy_table$km_status <- NULL
  legacy_table$km_error <- NULL
  legacy <- result
  legacy$schema_version <- "1.0"
  legacy$tables <- list(survival = legacy_table)
  legacy$table <- NULL

  plot <- sn_plot_survival(legacy)
  expect_s3_class(plot, "ggplot")
  expect_silent(ggplot2::ggplot_build(plot))
})

test_that("bulk dispatcher and plots expose the documented surface", {
  skip_if_not_installed("edgeR")
  fixture <- bulk_fixture()
  qc <- sn_run_bulk(fixture$counts, "qc", metadata = fixture$metadata)
  de <- sn_run_bulk(
    fixture$counts, "de", metadata = fixture$metadata,
    design = ~condition, contrast = c("condition", "tumor", "normal")
  )
  survival <- suppressWarnings(sn_run_survival(
    fixture$counts, "time", "event", "gene_1", metadata = fixture$metadata
  ))
  plots <- list(
    sn_plot_bulk_qc(qc), sn_plot_bulk_pca(qc, fixture$metadata, "condition"),
    sn_plot_sample_correlation(qc), sn_plot_bulk_de(de), sn_plot_survival(survival)
  )
  expect_true(all(vapply(plots, inherits, logical(1), "ggplot")))
})

test_that("external backends cannot leave the session directory deleted", {
  directory <- withr::local_tempdir()
  value <- Shennong:::.sn_with_preserved_directory(directory, {
    unlink(directory, recursive = TRUE)
    "finished"
  })
  expect_identical(value, "finished")
  expect_true(dir.exists(directory))

  unlink(directory, recursive = TRUE)
  dir.create(directory, recursive = TRUE)
  expect_error(
    Shennong:::.sn_with_preserved_directory(directory, {
      unlink(directory, recursive = TRUE)
      stop("backend failed")
    }),
    "backend failed"
  )
  expect_true(dir.exists(directory))
})
