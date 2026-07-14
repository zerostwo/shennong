library(testthat)

make_spatial_test_object <- function() {
  set.seed(43)
  counts <- matrix(rpois(8 * 36, 2), nrow = 8)
  rownames(counts) <- paste0("G", seq_len(nrow(counts)))
  colnames(counts) <- paste0("spot", seq_len(ncol(counts)))
  grid <- expand.grid(x = 1:6, y = 1:6)
  counts["G1", grid$x <= 3] <- counts["G1", grid$x <= 3] + 10
  counts["G2", grid$x > 3] <- counts["G2", grid$x > 3] + 10
  object <- SeuratObject::CreateSeuratObject(Matrix::Matrix(counts, sparse = TRUE))
  object$x <- grid$x
  object$y <- grid$y
  object$region <- ifelse(grid$x <= 3, "left", "right")
  object$sample <- rep(c("S1", "S2"), each = 18)
  Seurat::NormalizeData(object, verbose = FALSE)
}

test_that("Moran's I spatial features retain graph and permutation evidence", {
  object <- make_spatial_test_object()
  result <- sn_find_spatial_features(
    object, features = paste0("G", 1:4), return_object = FALSE,
    backend_control = list(k = 4, n_permutations = 19, seed = 9)
  )
  expect_true(sn_validate_result(result, error = FALSE)$valid)
  expect_equal(nrow(result$tables$features), 4L)
  expect_true(all(c("score", "p_value", "adjusted_p_value", "rank") %in% names(result$tables$features)))
  expect_equal(nrow(result$graphs$spatial), ncol(object) * 4L)
  expect_gt(result$tables$features$score[result$tables$features$feature == "G1"], 0)
})

test_that("nnSVG and SPARK-X feature adapters standardize results", {
  object <- make_spatial_test_object()
  output <- list(table = tibble::tibble(gene = c("G1", "G2"), LR_stat = c(8, 5), pval = c(0.01, 0.03), padj = c(0.02, 0.04)))
  for (method in c("nnsvg", "sparkx")) {
    result <- sn_find_spatial_features(object, method = method, backend_control = list(result = output), return_object = FALSE)
    expect_equal(result$method, method)
    expect_equal(result$tables$features$feature, c("G1", "G2"))
  }
})

test_that("spatial domain adapters store assignments in metadata", {
  object <- make_spatial_test_object()
  output <- list(domains = tibble::tibble(cell_id = colnames(object), cluster = rep(c("D1", "D2"), each = 18)))
  for (method in c("banksy", "stlearn", "bayesspace", "cellcharter")) {
    updated <- sn_find_spatial_domains(object, method = method, store_name = paste0("domain_", method), backend_control = list(result = output))
    result <- sn_get_result(updated, "spatial_domains", paste0("domain_", method))
    expect_equal(result$diagnostics$domains, 2L)
    expect_true(paste0("domain_", method) %in% colnames(updated[[]]))
  }
})

test_that("spatial neighborhoods preserve enrichment and co-occurrence", {
  object <- make_spatial_test_object()
  result <- sn_run_spatial_neighborhood(
    object, group_by = "region", return_object = FALSE,
    backend_control = list(k = 4, n_permutations = 19, seed = 12, distance_bins = 3)
  )
  expect_equal(nrow(result$graphs$spatial), ncol(object) * 4L)
  expect_equal(nrow(result$tables$enrichment), 4L)
  expect_true(all(c("z_score", "p_value", "adjusted_p_value") %in% names(result$tables$enrichment)))
  expect_gt(nrow(result$tables$cooccurrence), 0L)
})

make_communication_test_result <- function() {
  list(
    schema_version = "1.0", analysis_type = "communication", name = "communication",
    method = "synthetic", backend = "synthetic", input = list(), parameters = list(),
    tables = list(primary = tibble::tibble(
      source = c("left", "left", "right"), target = c("left", "right", "right"),
      ligand = c("L1", "L2", "L3"), receptor = c("R1", "R2", "R3"), score = c(1, 0.7, 0.9)
    )), embeddings = list(), graphs = list(), models = list(), diagnostics = list(), warnings = character(),
    provenance = list(package_versions = list(), random_seed = 1L, timestamp = "2026-01-01 UTC")
  )
}

test_that("spatial communication adds distance evidence and filters", {
  object <- make_spatial_test_object()
  unfiltered <- sn_run_spatial_communication(
    object, communication = make_communication_test_result(), group_by = "region", return_object = FALSE
  )
  expect_equal(nrow(unfiltered$tables$primary), 3L)
  expect_true(all(is.finite(unfiltered$tables$primary$spatial_distance)))
  filtered <- sn_run_spatial_communication(
    object, communication = make_communication_test_result(), group_by = "region",
    max_distance = min(unfiltered$tables$primary$spatial_distance), return_object = FALSE
  )
  expect_lte(nrow(filtered$tables$primary), nrow(unfiltered$tables$primary))
})

test_that("spatial integration and dispatcher standardize adapter embeddings", {
  object <- make_spatial_test_object()
  embedding <- tibble::tibble(cell = colnames(object), latent_1 = seq_len(ncol(object)), latent_2 = rev(seq_len(ncol(object))))
  result <- sn_integrate_spatial(object, backend_control = list(result = list(embedding = embedding)), return_object = FALSE)
  expect_equal(dim(result$embeddings$integrated), c(ncol(object), 2L))
  qc <- sn_run_spatial(object, task = "qc", k = 3)
  expect_equal(qc$diagnostics$locations, ncol(object))
  svg <- sn_run_spatial(
    object, task = "svg", method = "morans_i", features = c("G1", "G2"),
    backend_control = list(n_permutations = 3), return_object = FALSE
  )
  expect_equal(svg$analysis_type, "spatial_features")
})

test_that("spatial plots preserve aspect and render", {
  object <- make_spatial_test_object()
  features <- sn_find_spatial_features(object, features = c("G1", "G2"), backend_control = list(n_permutations = 3), return_object = FALSE)
  domains <- sn_find_spatial_domains(
    object, backend_control = list(result = list(domains = tibble::tibble(cell = colnames(object), domain = object$region))),
    return_object = FALSE
  )
  neighborhoods <- sn_run_spatial_neighborhood(object, group_by = "region", backend_control = list(n_permutations = 3), return_object = FALSE)
  communication <- sn_run_spatial_communication(object, communication = make_communication_test_result(), group_by = "region", return_object = FALSE)
  deconv <- transform(
    expand.grid(cell = colnames(object), cell_type = c("A", "B"), stringsAsFactors = FALSE),
    spatial_x = rep(object$x, 2), spatial_y = rep(object$y, 2), proportion = rep(c(0.7, 0.3), each = ncol(object))
  )
  plots <- list(
    sn_plot_spatial(object, "region"), sn_plot_spatial_feature(object, c("G1", "G2")),
    sn_plot_spatial_svg(features), sn_plot_spatial_domain(domains),
    sn_plot_spatial_neighborhood(neighborhoods), sn_plot_spatial_neighborhood(neighborhoods, type = "cooccurrence"),
    sn_plot_spatial_communication(communication), sn_plot_spatial_deconvolution(deconv)
  )
  for (plot in plots) {
    expect_s3_class(plot, "ggplot")
    expect_silent(ggplot2::ggplotGrob(plot))
  }
})

test_that("spatial methods are registered as implemented", {
  expect_true(all(sn_list_methods("spatial_svg")$implemented))
  expect_true(all(sn_list_methods("spatial_domain")$implemented))
  expect_true(sn_method_status("distance", "spatial_communication")$implemented)
})
