library(testthat)

make_trajectory_test_object <- function(branching = TRUE) {
  skip_if_not_installed("SeuratObject")
  set.seed(29)
  if (branching) {
    n_each <- 30L
    root_x <- seq(0, 1, length.out = n_each)
    upper_x <- seq(1, 2, length.out = n_each)
    lower_x <- seq(1, 2, length.out = n_each)
    embedding <- rbind(
      cbind(root_x, rnorm(n_each, 0, 0.025)),
      cbind(upper_x, upper_x - 1 + rnorm(n_each, 0, 0.025)),
      cbind(lower_x, -(lower_x - 1) + rnorm(n_each, 0, 0.025))
    )
    clusters <- rep(c("root", "effector", "memory"), each = n_each)
  } else {
    n_each <- 20L
    x <- seq(0, 2, length.out = 3L * n_each)
    embedding <- cbind(x, sin(x * pi) / 10 + rnorm(length(x), 0, 0.02))
    clusters <- rep(c("early", "middle", "late"), each = n_each)
  }
  cells <- paste0("cell", seq_len(nrow(embedding)))
  rownames(embedding) <- cells
  colnames(embedding) <- c("PC_1", "PC_2")
  progress <- scales::rescale(embedding[, 1])
  means <- rbind(
    1 + 10 * progress,
    1 + 10 * (1 - progress),
    2 + 8 * abs(embedding[, 2]),
    2 + 5 * (embedding[, 2] > 0),
    2 + 5 * (embedding[, 2] < 0),
    matrix(3, nrow = 7L, ncol = length(progress))
  )
  counts <- matrix(stats::rpois(length(means), lambda = as.numeric(means)), nrow = nrow(means))
  rownames(counts) <- paste0("G", seq_len(nrow(counts)))
  colnames(counts) <- cells
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  object[["pca"]] <- SeuratObject::CreateDimReducObject(
    embeddings = embedding,
    key = "PC_",
    assay = "RNA"
  )
  object$seurat_clusters <- clusters
  SeuratObject::Idents(object) <- clusters
  object
}

test_that("Slingshot trajectory stores complete per-cell and graph contracts", {
  skip_if_not_installed("slingshot")
  object <- make_trajectory_test_object()
  object <- sn_run_trajectory(
    object,
    reduction = "pca",
    cluster_by = "seurat_clusters",
    start = "root",
    end = c("effector", "memory"),
    store_name = "branching",
    test_dynamic = FALSE
  )
  result <- sn_get_result(object, "trajectory", "branching")

  expect_equal(result$analysis_type, "trajectory")
  expect_equal(result$method, "slingshot")
  expect_gte(result$diagnostics$n_lineages, 2L)
  expect_equal(nrow(result$tables$cells), ncol(object))
  expect_true(all(c("primary_lineage", "primary_pseudotime") %in% names(result$tables$cells)))
  expect_true(any(grepl("^weight_", names(result$tables$cells))))
  expect_equal(nrow(result$tables$terminal_states), result$diagnostics$n_lineages)
  expect_gt(nrow(result$tables$curves), 0L)
  expect_equal(names(result$graphs$lineages), unique(result$tables$curves$lineage))
  expect_true(all(c("branching_pseudotime", "branching_lineage") %in% colnames(object[[]])))
  expect_true(sn_validate_result(result, error = FALSE)$valid)
})

test_that("trajectory endpoints and expected lineage paths are validated", {
  skip_if_not_installed("slingshot")
  object <- make_trajectory_test_object()
  expect_error(
    sn_run_trajectory(object, reduction = "pca", start = "missing", test_dynamic = FALSE),
    "not found"
  )
  expect_error(
    sn_run_trajectory(object, reduction = "pca", lineages = list(c("root", "missing")), test_dynamic = FALSE),
    "not found"
  )
  result <- sn_run_trajectory(
    object,
    reduction = "pca",
    lineages = list(to_effector = c("root", "effector"), to_memory = c("root", "memory")),
    test_dynamic = FALSE,
    return_object = FALSE
  )
  expect_equal(length(result$parameters$requested_lineages), 2L)
})

test_that("tradeSeq dynamic tests retain tests, trends, and convergence", {
  skip_if_not_installed("slingshot")
  skip_if_not_installed("tradeSeq")
  object <- make_trajectory_test_object(branching = FALSE)
  result <- suppressMessages(sn_run_trajectory(
    object,
    reduction = "pca",
    start = "early",
    end = "late",
    test_dynamic = TRUE,
    dynamic_features = paste0("G", 1:8),
    max_dynamic_features = 8L,
    nknots = 3L,
    trend_features = c("G1", "G2"),
    trend_points = 20L,
    return_object = FALSE
  ))

  expect_equal(result$backend, "slingshot+tradeSeq")
  expect_equal(nrow(result$tables$dynamic_genes), 8L)
  expect_true(all(c("feature", "test", "pvalue", "adjusted_p_value") %in% names(result$tables$dynamic_genes)))
  expect_equal(sort(unique(result$tables$fitted_trends$gene)), c("G1", "G2"))
  expect_equal(nrow(result$tables$fitted_trends), 2L * 20L)
  expect_equal(nrow(result$tables$convergence), 8L)
  expect_true(all(result$tables$convergence$converged))
  expect_equal(result$diagnostics$dynamic$tested_features, 8L)
})

test_that("trajectory and dynamic plots render from stored results", {
  skip_if_not_installed("slingshot")
  skip_if_not_installed("tradeSeq")
  object <- make_trajectory_test_object(branching = FALSE)
  object <- suppressMessages(sn_run_trajectory(
    object,
    reduction = "pca",
    start = "early",
    end = "late",
    dynamic_features = paste0("G", 1:8),
    max_dynamic_features = 8L,
    nknots = 3L,
    trend_features = c("G1", "G2"),
    trend_points = 20L
  ))
  result <- sn_get_result(object, "trajectory", "trajectory")
  plots <- list(
    sn_plot_trajectory(object, "trajectory"),
    sn_plot_pseudotime(result),
    sn_plot_lineage_probability(result, lineage = "Lineage1"),
    sn_plot_dynamic_heatmap(result),
    sn_plot_gene_trend(result, features = c("G1", "G2"))
  )
  expect_true(all(vapply(plots, inherits, logical(1), "ggplot")))
  for (plot in plots) expect_silent(ggplot2::ggplotGrob(plot))
  expect_error(sn_plot_branch_comparison(result), "No branch-specific")
})
