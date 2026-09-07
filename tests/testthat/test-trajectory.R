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
    result_id = "branching",
    test_dynamic = FALSE
  )
  result <- sn_get_result(object, "trajectory", "branching")

  expect_equal(result$analysis_type, "trajectory")
  expect_equal(result$method, "slingshot")
  expect_gte(result$diagnostics$n_lineages, 2L)
  expect_equal(nrow(result$tables$cells), ncol(object))
  expect_identical(result$tables$primary, result$tables$cells)
  expect_true(all(c("primary_lineage", "primary_pseudotime") %in% names(result$tables$cells)))
  expect_true(any(grepl("^weight_", names(result$tables$cells))))
  expect_equal(nrow(result$tables$terminal_states), result$diagnostics$n_lineages)
  expect_gt(nrow(result$tables$curves), 0L)
  expect_equal(names(result$graphs$lineages), unique(result$tables$curves$lineage))
  expect_true(all(c("branching_pseudotime", "branching_lineage") %in% colnames(object[[]])))
  expect_true(sn_validate_result(result, error = FALSE)$valid)
})

test_that("trajectory result provenance includes the scoped Slingshot patch", {
  object <- make_trajectory_test_object(branching = FALSE)
  testthat::local_mocked_bindings(
    .sn_trajectory_slingshot = function(embedding, clusters, ...) {
      cells <- rownames(embedding)
      pseudotime <- matrix(
        seq(0, 1, length.out = length(cells)),
        ncol = 1L,
        dimnames = list(cells, "Lineage1")
      )
      list(
        pseudotime = pseudotime,
        weights = matrix(1, nrow = length(cells), ncol = 1L,
                         dimnames = dimnames(pseudotime)),
        lineages = list(Lineage1 = unique(as.character(clusters))),
        curves = tibble::tibble(),
        backend = "slingshot"
      )
    },
    .package = "Shennong"
  )

  result <- sn_run_trajectory(
    object,
    reduction = "pca",
    cluster_by = "seurat_clusters",
    method = "slingshot",
    test_dynamic = FALSE,
    return_object = FALSE
  )

  expect_identical(
    result$provenance$acceleration$suppressed_patches,
    "slingshot"
  )
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

test_that("trajectory lineage assignment requires finite pseudotime", {
  embedding <- matrix(
    c(0, 0, 1, 1), ncol = 2L, byrow = TRUE,
    dimnames = list(c("c1", "c2"), c("x", "y"))
  )
  pseudotime <- matrix(
    c(NA, 0.2, NA, 0.8), nrow = 2L, byrow = TRUE,
    dimnames = list(c("c1", "c2"), c("L1", "L2"))
  )
  weights <- matrix(
    c(0.9, 0.1, 0.8, 0.2), nrow = 2L, byrow = TRUE,
    dimnames = dimnames(pseudotime)
  )

  cells <- Shennong:::.sn_trajectory_cell_table(
    embedding,
    clusters = stats::setNames(c("a", "b"), c("c1", "c2")),
    pseudotime = pseudotime,
    weights = weights
  )

  expect_identical(cells$primary_lineage, c("L2", "L2"))
  expect_equal(cells$primary_pseudotime, c(0.2, 0.8))

  pseudotime[,] <- NA_real_
  unassigned <- Shennong:::.sn_trajectory_cell_table(
    embedding,
    clusters = stats::setNames(c("a", "b"), c("c1", "c2")),
    pseudotime = pseudotime,
    weights = weights
  )
  expect_true(all(is.na(unassigned$primary_lineage)))
  expect_true(all(is.na(unassigned$primary_pseudotime)))
})

test_that("trajectory weights align by lineage name rather than column position", {
  object <- make_trajectory_test_object(branching = FALSE)
  embedding <- SeuratObject::Embeddings(object[["pca"]])
  cells <- rownames(embedding)
  clusters <- stats::setNames(as.character(object$seurat_clusters), cells)
  pseudotime <- cbind(
    LineageA = seq(0, 1, length.out = length(cells)),
    LineageB = seq(1, 0, length.out = length(cells))
  )
  rownames(pseudotime) <- cells
  weights <- cbind(LineageB = rep(0.2, length(cells)), LineageA = rep(0.8, length(cells)))
  rownames(weights) <- cells

  standardized <- Shennong:::.sn_standardize_trajectory_backend(
    list(
      pseudotime = pseudotime,
      weights = weights,
      lineages = list(LineageB = c("middle", "late"), LineageA = c("early", "middle"))
    ),
    embedding, clusters, requested = NULL, start = NULL, end = NULL, method = "mock"
  )

  expect_identical(colnames(standardized$weights), c("LineageA", "LineageB"))
  expect_equal(unname(standardized$weights[, "LineageA"]), rep(0.8, length(cells)))
  expect_identical(names(standardized$lineages), c("LineageA", "LineageB"))
  expect_error(
    Shennong:::.sn_standardize_trajectory_backend(
      list(pseudotime = pseudotime, weights = unname(weights)),
      embedding, clusters, requested = NULL, start = NULL, end = NULL, method = "mock"
    ),
    "lineage names must match"
  )
})

test_that("multi-lineage trajectory outputs require explicit coherent weights", {
  object <- make_trajectory_test_object(branching = FALSE)
  embedding <- SeuratObject::Embeddings(object[["pca"]])
  cells <- rownames(embedding)
  clusters <- stats::setNames(as.character(object$seurat_clusters), cells)
  pseudotime <- cbind(
    LineageA = seq(0, 1, length.out = length(cells)),
    LineageB = seq(1, 0, length.out = length(cells))
  )
  rownames(pseudotime) <- cells

  expect_error(
    Shennong:::.sn_standardize_trajectory_backend(
      list(pseudotime = pseudotime), embedding, clusters,
      requested = NULL, start = NULL, end = NULL, method = "mock"
    ),
    "requires explicit lineage `weights`"
  )

  weights <- matrix(0.5, nrow = nrow(pseudotime), ncol = 2L, dimnames = dimnames(pseudotime))
  pseudotime[1, 1] <- NA_real_
  weights[1, 1] <- 0.5
  expect_error(
    Shennong:::.sn_standardize_trajectory_backend(
      list(pseudotime = pseudotime, weights = weights), embedding, clusters,
      requested = NULL, start = NULL, end = NULL, method = "mock"
    ),
    "positive trajectory lineage weight.*finite pseudotime"
  )
})

test_that("trajectory adapters reject ambiguous identities and lineage columns", {
  object <- make_trajectory_test_object(branching = FALSE)
  embedding <- SeuratObject::Embeddings(object[["pca"]])
  cells <- rownames(embedding)
  clusters <- stats::setNames(as.character(object$seurat_clusters), cells)

  duplicated_long <- tibble::tibble(
    cell = c(cells, cells[[1L]]),
    lineage = "LineageA",
    pseudotime = seq_len(length(cells) + 1L)
  )
  expect_error(
    Shennong:::.sn_trajectory_backend_matrix(
      duplicated_long, cells, "pseudotime"
    ),
    "duplicate cell-lineage"
  )

  extra <- matrix(
    seq_len(length(cells) + 1L), ncol = 1L,
    dimnames = list(c(cells, "unknown-cell"), "LineageA")
  )
  expect_error(
    Shennong:::.sn_trajectory_backend_matrix(extra, cells, "pseudotime"),
    "match the analyzed cells exactly"
  )

  pseudotime <- cbind(
    `A-B` = seq(0, 1, length.out = length(cells)),
    A.B = seq(1, 0, length.out = length(cells))
  )
  rownames(pseudotime) <- cells
  weights <- matrix(0.5, nrow = length(cells), ncol = 2L, dimnames = dimnames(pseudotime))
  expect_error(
    Shennong:::.sn_standardize_trajectory_backend(
      list(pseudotime = pseudotime, weights = weights),
      embedding, clusters, requested = NULL, start = NULL, end = NULL,
      method = "mock"
    ),
    "remain unique after column-name sanitization"
  )
})

test_that("Monocle 3 defaults to an existing UMAP reduction", {
  object <- make_trajectory_test_object(branching = FALSE)
  umap <- SeuratObject::Embeddings(object[["pca"]])
  colnames(umap) <- c("UMAP_1", "UMAP_2")
  object[["umap"]] <- SeuratObject::CreateDimReducObject(
    embeddings = umap, key = "UMAP_", assay = "RNA"
  )
  runner <- function(object, method, embedding, clusters, ...) {
    expect_equal(colnames(embedding), c("UMAP_1", "UMAP_2"))
    cells <- rownames(embedding)
    pseudotime <- stats::setNames(seq(0, 1, length.out = length(cells)), cells)
    list(
      pseudotime = pseudotime,
      lineages = list(Lineage1 = unique(as.character(clusters))),
      backend = "mock_monocle3"
    )
  }
  result <- sn_run_trajectory(
    object, method = "monocle3", start = "early", test_dynamic = FALSE,
    backend_control = list(runner = runner), return_object = FALSE
  )
  expect_identical(result$input$reduction, "umap")
})

test_that("direct Monocle backend refuses to relabel PCA as UMAP", {
  object <- make_trajectory_test_object(branching = FALSE)
  embedding <- SeuratObject::Embeddings(object[["pca"]])
  clusters <- stats::setNames(as.character(object$seurat_clusters), rownames(embedding))
  expect_error(
    Shennong:::.sn_run_monocle3_trajectory(
      object, embedding, reduction = "pca", clusters = clusters,
      start = "early", assay = "RNA", counts_layer = "counts", backend_control = list()
    ),
    "genuine UMAP"
  )
})

test_that("direct Monocle backend fails closed for unsupported ends and partitions", {
  object <- make_trajectory_test_object(branching = FALSE)
  embedding <- SeuratObject::Embeddings(object[["pca"]])
  clusters <- stats::setNames(as.character(object$seurat_clusters), rownames(embedding))
  expect_error(
    Shennong:::.sn_run_monocle3_trajectory(
      object, embedding, reduction = "umap", clusters = clusters,
      start = "early", assay = "RNA", counts_layer = "counts",
      backend_control = list(), end = "late"
    ),
    "cannot enforce `end`"
  )
  expect_error(
    Shennong:::.sn_validate_monocle3_partitions(c("1", "2", "1")),
    "multiple disconnected partitions"
  )
  expect_silent(Shennong:::.sn_validate_monocle3_partitions(rep("1", 3)))
})

test_that("tradeSeq count validation rejects invalid count semantics", {
  valid <- Matrix::Matrix(matrix(c(0, 1, 2, 3), nrow = 2), sparse = TRUE)
  expect_identical(Shennong:::.sn_validate_tradeseq_counts(valid), valid)
  expect_error(Shennong:::.sn_validate_tradeseq_counts(matrix(c(0, -1))), "non-negative")
  expect_error(Shennong:::.sn_validate_tradeseq_counts(matrix(c(0, Inf))), "finite")
  expect_error(Shennong:::.sn_validate_tradeseq_counts(matrix(c(0, 1.5))), "integer-valued")
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
