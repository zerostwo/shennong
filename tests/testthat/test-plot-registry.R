test_that("plot method registry exposes canonical result views", {
  methods <- sn_list_plot_methods()
  expect_s3_class(methods, "tbl_df")
  expect_true(all(c(
    "analysis_type", "view", "default", "required_parameters", "accepted_input"
  ) %in% names(methods)))
  expect_true(all(c("de", "bulk_qc", "cell_communication", "trajectory") %in% methods$analysis_type))
  expect_equal(methods$view[methods$analysis_type == "de" & methods$default], "volcano")
  expect_true("sample_scatter" %in% sn_list_plot_methods("bulk_qc")$view)
  expect_equal(
    methods$required_parameters[methods$analysis_type == "trajectory" & methods$view == "gene_trend"],
    "features"
  )
  expect_equal(sn_list_plot_methods("communication")$analysis_type, rep("cell_communication", 7L))
  expect_error(sn_list_plot_methods("not_real"), "Unknown plot analysis type")
})

test_that("sn_plot_result dispatches compatible tables and legacy assessments", {
  de <- data.frame(
    gene = paste0("gene_", seq_len(12L)),
    log2_fold_change = seq(-2, 2, length.out = 12L),
    adjusted_p_value = seq(0.01, 0.12, length.out = 12L)
  )
  qc <- list(by_sample = data.frame(
    sample = c("S1", "S2"), qc_score = c(0.9, 0.7),
    n_cells = c(100L, 80L), retention_fraction = c(0.95, 0.8)
  ))

  volcano <- sn_plot_result(de, analysis_type = "de")
  qc_plot <- sn_plot_result(qc, analysis_type = "qc", view = "n_cells")
  expect_s3_class(volcano, "ggplot")
  expect_s3_class(qc_plot, "ggplot")
  expect_error(sn_plot_result(de), "analysis_type.*required")
  expect_error(
    sn_plot_result(de, analysis_type = "de", view = "bubble"),
    "Unknown view"
  )
})

test_that("sn_plot_result infers one stored Seurat result without an ID", {
  skip_if_not_installed("SeuratObject")
  counts <- Matrix::Matrix(
    matrix(
      stats::rpois(60L, 2), nrow = 6L,
      dimnames = list(paste0("gene", 1:6), paste0("cell", 1:10))
    ),
    sparse = TRUE
  )
  object <- SeuratObject::CreateSeuratObject(counts)
  table <- data.frame(
    gene = paste0("gene", 1:6),
    log2_fold_change = seq(-2, 2, length.out = 6L),
    adjusted_p_value = seq(0.01, 0.06, length.out = 6L)
  )
  result <- Shennong:::.sn_new_analysis_result(
    "de", "markers", "test", "test", tables = list(primary = table)
  )
  object <- sn_store_result(object, "de", "markers", result)

  plot <- sn_plot_result(object, view = "volcano")
  expect_s3_class(plot, "ggplot")
  expect_error(
    sn_plot_result(object, analysis_type = "enrichment"),
    "No stored results"
  )
})

test_that("canonical stored-result plots retain Seurat metadata context", {
  skip_if_not_installed("SeuratObject")
  counts <- Matrix::Matrix(
    matrix(
      seq_len(24L), nrow = 6L,
      dimnames = list(paste0("gene", 1:6), paste0("cell", 1:4))
    ),
    sparse = TRUE
  )
  object <- SeuratObject::CreateSeuratObject(counts)
  object$group <- c("A", "A", "B", "B")
  object$known <- c("T", "T", "B", "B")

  scores <- tidyr::expand_grid(
    entity = colnames(object),
    program = c("p1", "p2")
  ) |>
    dplyr::mutate(score = seq_len(dplyr::n()), level = "cell")
  program <- Shennong:::.sn_new_analysis_result(
    "program_scoring", "programs", "mean", "Shennong",
    tables = list(primary = scores, scores = scores)
  )
  object <- sn_store_result(object, "program_scoring", "programs", program)
  program_plot <- sn_plot_result(
    object,
    "program_scoring",
    "programs",
    group_by = "group"
  )
  expect_setequal(unique(program_plot$data$display_group), c("A", "B"))

  cells <- tibble::tibble(
    cell = colnames(object),
    prediction = c("T", "T", "B", "B"),
    prediction_score = 1,
    low_confidence = FALSE
  )
  annotation <- Shennong:::.sn_new_analysis_result(
    "annotation", "labels", "test", "test",
    tables = list(primary = cells, cells = cells)
  )
  object <- sn_store_result(object, "annotation", "labels", annotation)
  confusion <- sn_plot_result(
    object,
    "annotation",
    "labels",
    view = "confusion",
    truth = "known"
  )
  expect_s3_class(confusion, "ggplot")
  expect_equal(sum(confusion$data$Freq), ncol(object))
})

test_that("sn_plot_association supports observation and sample-level data", {
  data <- data.frame(
    sample = rep(paste0("S", 1:4), each = 3L),
    condition = rep(c("control", "control", "treated", "treated"), each = 3L),
    score_a = seq_len(12L),
    score_b = seq_len(12L) + rep(c(-0.2, 0.1), 6L)
  )
  plot <- sn_plot_association(
    data,
    x = "score_a",
    y = "score_b",
    sample_by = "sample",
    group_by = "condition",
    method = "pearson"
  )
  source <- attr(plot, "shennong_figure_data")
  expect_s3_class(plot, "ggplot")
  expect_equal(nrow(source), 4L)
  expect_match(plot$labels$subtitle, "r =")

  inconsistent <- data
  inconsistent$condition[[2L]] <- "treated"
  expect_error(
    sn_plot_association(
      inconsistent, "score_a", "score_b",
      sample_by = "sample", group_by = "condition"
    ),
    "constant within each sample"
  )
})

test_that("plot subsampling is reproducible without changing caller RNG state", {
  data <- data.frame(x = seq_len(50), y = rev(seq_len(50)))
  set.seed(812L)
  before <- .Random.seed

  first <- sn_plot_association(
    data, x = "x", y = "y", max_points = 12L, seed = 19L
  )
  after <- .Random.seed
  second <- sn_plot_association(
    data, x = "x", y = "y", max_points = 12L, seed = 19L
  )

  expect_identical(after, before)
  expect_identical(first$data, second$data)
})

test_that("sample-to-sample correlation scatter uses stored bulk expression", {
  set.seed(92)
  expression <- matrix(
    stats::rpois(120L, 20), nrow = 30L,
    dimnames = list(paste0("gene_", 1:30), paste0("sample_", 1:4))
  )
  qc <- sn_assess_bulk_qc(expression, top_variable = 20L)
  expect_equal(dim(qc$tables$expression), c(20L, 4L))

  direct <- sn_plot_sample_correlation(
    qc,
    view = "scatter",
    sample_x = "sample_1",
    sample_y = "sample_2"
  )
  dispatched <- sn_plot_result(
    qc,
    view = "sample_scatter",
    sample_x = "sample_1",
    sample_y = "sample_2"
  )
  expect_s3_class(direct, "ggplot")
  expect_s3_class(dispatched, "ggplot")
  expect_equal(nrow(attr(direct, "shennong_figure_data")), 20L)
})

test_that("sn_plot_distribution consolidates grouped and QC distributions", {
  data <- data.frame(
    sample = rep(paste0("S", 1:4), each = 4L),
    condition = rep(c("control", "control", "treated", "treated"), each = 4L),
    nFeature_RNA = seq(100, 250, length.out = 16L),
    percent.mt = seq(1, 8, length.out = 16L)
  )
  plot <- sn_plot_distribution(
    data,
    features = c("nFeature_RNA", "percent.mt"),
    group_by = "condition",
    sample_by = "sample",
    view = "box",
    thresholds = list(percent.mt = 5),
    show_points = TRUE
  )
  source <- attr(plot, "shennong_figure_data")
  expect_s3_class(plot, "ggplot")
  expect_equal(nrow(source$values), 8L)
  expect_equal(nrow(source$thresholds), 1L)

  qc_plot <- sn_plot_qc_thresholds(
    data,
    features = c("nFeature_RNA", "percent.mt"),
    thresholds = list(percent.mt = c(3, 6)),
    sample_by = "sample"
  )
  expect_s3_class(qc_plot, "ggplot")
  expect_true(inherits(sn_get_figure_spec(qc_plot), "sn_get_figure_spec"))
})
