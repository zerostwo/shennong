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

test_that("Moran permutations change values and respect sample boundaries", {
  object <- make_spatial_test_object()
  coordinates <- Shennong:::.sn_spatial_coordinates(
    object, spatial_cols = c("x", "y"), sample_by = "sample"
  )$table
  graph <- Shennong:::.sn_spatial_knn_graph(coordinates, k = 4L)
  samples <- stats::setNames(coordinates$spatial_sample, coordinates$cell)
  cells <- coordinates$cell
  values <- stats::setNames(seq_along(cells), cells)

  set.seed(9)
  null <- replicate(19L, {
    permuted <- Shennong:::.sn_spatial_permute_values(values, cells, samples)
    expect_equal(
      tapply(unname(permuted), samples[names(permuted)], sort),
      tapply(unname(values), samples[names(values)], sort)
    )
    Shennong:::.sn_morans_i(permuted, graph, cells)
  })

  expect_gt(length(unique(null)), 1L)
  expect_true(all(samples[graph$source] == samples[graph$target]))
})

test_that("Moran permutation p-values are centered on the empirical null", {
  calls <- 0L
  statistics <- c(10, 8, 9, 10)
  testthat::local_mocked_bindings(
    .sn_morans_i = function(...) {
      calls <<- calls + 1L
      statistics[[calls]]
    },
    .package = "Shennong"
  )
  expression <- list(matrix = matrix(
    c(1, 2), nrow = 1L,
    dimnames = list("G1", c("C1", "C2"))
  ))
  coordinates <- tibble::tibble(
    cell = c("C1", "C2"), spatial_x = c(0, 1), spatial_y = c(0, 0),
    spatial_sample = "S1"
  )
  graph <- tibble::tibble(source = "C1", target = "C2", distance = 1, weight = 1)
  result <- Shennong:::.sn_spatial_morans_table(
    expression, coordinates, graph, n_permutations = 3L, seed = 9L
  )

  expect_equal(result$null_mean, 9)
  expect_equal(result$p_value, 0.75)
})

test_that("spatial permutation counts are validated before computation", {
  for (invalid in list(-1, 1.5, Inf, NA_real_, c(1, 2))) {
    expect_error(
      Shennong:::.sn_validate_spatial_permutation_count(invalid),
      "one non-negative integer"
    )
  }
  expect_identical(Shennong:::.sn_validate_spatial_permutation_count(0), 0L)
})

test_that("spatial permutation seeds do not change caller RNG state", {
  set.seed(101)
  before <- .Random.seed
  Shennong:::.sn_with_seed(9, runif(10))
  expect_identical(.Random.seed, before)
})

test_that("spatial enrichment supports one permutation without dropping matrix dimensions", {
  coordinates <- tibble::tibble(
    cell = paste0("C", 1:4), spatial_x = 1:4, spatial_y = 0,
    spatial_sample = "S1"
  )
  graph <- Shennong:::.sn_spatial_knn_graph(coordinates, k = 1L)
  labels <- stats::setNames(c("A", "A", "B", "B"), coordinates$cell)
  result <- Shennong:::.sn_spatial_enrichment(
    graph, labels, n_permutations = 1L, seed = 9L
  )
  expect_equal(nrow(result), 4L)
  expect_true(all(is.finite(result$p_value)))
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
    updated <- sn_find_spatial_domains(object, method = method, result_id = paste0("domain_", method), backend_control = list(result = output))
    result <- sn_get_result(updated, "spatial_domains", paste0("domain_", method))
    expect_equal(result$diagnostics$domains, 2L)
    expect_true(paste0("domain_", method) %in% colnames(updated[[]]))
  }
})

test_that("spatial adapter identities and numeric outputs fail closed", {
  object <- make_spatial_test_object()
  cells <- colnames(object)
  valid_domains <- tibble::tibble(cell = cells, domain = object$region)

  duplicated <- valid_domains
  duplicated$cell[[2L]] <- duplicated$cell[[1L]]
  expect_error(
    sn_find_spatial_domains(
      object,
      backend_control = list(result = list(domains = duplicated)),
      return_object = FALSE
    ),
    "unique"
  )
  unknown <- valid_domains
  unknown$cell[[1L]] <- "not-an-object-cell"
  expect_error(
    sn_find_spatial_domains(
      object,
      backend_control = list(result = list(domains = unknown)),
      return_object = FALSE
    ),
    "match every object cell exactly"
  )

  embedding <- tibble::tibble(
    cell = cells,
    latent_1 = seq_along(cells),
    latent_2 = rev(seq_along(cells))
  )
  bad_embedding <- embedding
  bad_embedding$cell[[1L]] <- "unknown-cell"
  expect_error(
    sn_integrate_spatial(
      object,
      backend_control = list(result = list(embedding = bad_embedding)),
      return_object = FALSE
    ),
    "uniquely match every"
  )
  bad_embedding <- embedding
  bad_embedding$latent_2[[1L]] <- Inf
  expect_error(
    sn_integrate_spatial(
      object,
      backend_control = list(result = list(embedding = bad_embedding)),
      return_object = FALSE
    ),
    "finite numeric"
  )

  bad_graph <- tibble::tibble(
    source = cells[[1L]], target = "unknown-cell", distance = 1
  )
  expect_error(
    sn_run_spatial_neighborhood(
      object,
      group_by = "region",
      backend_control = list(result = list(graph = bad_graph)),
      return_object = FALSE
    ),
    "belong to the analyzed object"
  )
  bad_graph$target <- cells[[2L]]
  bad_graph$distance <- -1
  expect_error(
    sn_run_spatial_neighborhood(
      object,
      group_by = "region",
      backend_control = list(result = list(graph = bad_graph)),
      return_object = FALSE
    ),
    "finite non-negative"
  )
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

test_that("spatial neighborhoods never connect independent sections", {
  object <- make_spatial_test_object()
  result <- sn_run_spatial_neighborhood(
    object, group_by = "region", sample_by = "sample", return_object = FALSE,
    backend_control = list(k = 4, n_permutations = 9, seed = 12)
  )
  sample_map <- stats::setNames(object$sample, colnames(object))
  graph <- result$graphs$spatial
  expect_true(all(sample_map[graph$source] == sample_map[graph$target]))
  expect_identical(result$input$sample_by, "sample")

  qc <- sn_run_spatial(object, task = "qc", sample_by = "sample", k = 4L)
  expect_true(all(sample_map[qc$graph$source] == sample_map[qc$graph$target]))
  expect_equal(qc$diagnostics$samples, 2L)
})

test_that("single-section spatial backends fail closed on pooled sections", {
  object <- make_spatial_test_object()
  expect_error(
    sn_find_spatial_features(
      object, method = "nnsvg", sample_by = "sample", return_object = FALSE
    ),
    "one tissue section"
  )
  expect_error(
    sn_find_spatial_domains(
      object, method = "banksy", sample_by = "sample", return_object = FALSE
    ),
    "one tissue section"
  )
})

make_communication_test_result <- function() {
  list(
    schema_version = "2.0.0", analysis_type = "cell_communication", result_id = "communication",
    method = "synthetic", backend = "synthetic", input = list(), parameters = list(),
    tables = list(primary = tibble::tibble(
      source = c("left", "left", "right"), target = c("left", "right", "right"),
      ligand = c("L1", "L2", "L3"), receptor = c("R1", "R2", "R3"), score = c(1, 0.7, 0.9)
    )), embeddings = list(), graphs = list(), models = list(), diagnostics = list(), warnings = character(),
    provenance = list(
      package_versions = list(), random_seed = 1L, timestamp = "2026-01-01 UTC",
      result_id = "communication", analysis_type = "cell_communication"
    )
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

test_that("same-group spatial distance excludes self and aggregates within sections", {
  object <- make_spatial_test_object()
  result <- sn_run_spatial_communication(
    object, communication = make_communication_test_result(),
    group_by = "region", sample_by = "sample", return_object = FALSE
  )
  same_group <- result$tables$group_distances[
    result$tables$group_distances$source == result$tables$group_distances$target,
  ]
  expect_true(all(same_group$spatial_distance > 0))
  expect_true(all(result$tables$group_distances_by_sample$spatial_distance > 0))
  expect_true(all(result$tables$group_distances$contributing_samples >= 1L))
})

test_that("spatial communication joins sample-resolved distances by sample", {
  object <- make_spatial_test_object()
  x <- object$x
  shifted <- object$sample == "S2" & object$region == "right"
  x[shifted] <- x[shifted] + 20
  object$x <- x
  communication <- make_communication_test_result()
  communication$tables$primary <- tibble::tibble(
    source = c("left", "left"), target = c("right", "right"),
    ligand = c("L1", "L1"), receptor = c("R1", "R1"),
    score = c(1, 1), sample = c("S1", "S2")
  )

  result <- sn_run_spatial_communication(
    object, communication = communication, group_by = "region",
    sample_by = "sample", return_object = FALSE
  )
  observed <- result$tables$all_interactions
  expected <- result$tables$group_distances_by_sample
  expected <- expected[expected$source == "left" & expected$target == "right", ]
  expected <- expected$spatial_distance[match(observed$sample, expected$spatial_sample)]

  expect_equal(observed$spatial_distance, expected)
  expect_gt(observed$spatial_distance[observed$sample == "S2"], observed$spatial_distance[observed$sample == "S1"])
  expect_identical(result$diagnostics$distance_matching, "source_target_sample")
})

test_that("spatial communication validates maximum distance", {
  object <- make_spatial_test_object()
  for (invalid in list(-1, Inf, NA_real_, c(1, 2))) {
    expect_error(
      sn_run_spatial_communication(
        object, communication = make_communication_test_result(),
        group_by = "region", max_distance = invalid, return_object = FALSE
      ),
      "finite non-negative"
    )
  }
})

test_that("spatial communication reuses stored cell-communication results", {
  object <- make_spatial_test_object()
  object <- sn_store_result(
    object,
    type = "cell_communication",
    result_id = "communication",
    result = make_communication_test_result()
  )
  result <- sn_run_spatial_communication(
    object,
    source_result_id = "communication",
    group_by = "region",
    return_object = FALSE
  )
  expect_equal(nrow(result$tables$all_interactions), 3L)
  expect_equal(result$models$source_result$result_id, "communication")
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

test_that("spatial registry distinguishes implemented methods from roadmap entries", {
  expect_true(all(sn_list_methods("spatial_svg")$implemented))
  domains <- sn_list_methods("spatial_domain")
  expect_true(all(domains$implemented[domains$name != "stlearn"]))
  expect_false(domains$implemented[domains$name == "stlearn"])
  expect_false(sn_get_method_status("stlearn", "spatial_domain")$runnable)
  expect_true(sn_get_method_status("distance", "spatial_communication")$implemented)
})

test_that("legacy spatial aliases warn and forward to canonical workflows", {
  object <- make_spatial_test_object()
  local_mocked_bindings(
    sn_run_cell2location = function(object, marker = NULL, ...) {
      list(object = object, marker = marker)
    },
    sn_run_tangram = function(object, reference_object = NULL, marker = NULL, ...) {
      list(object = object, reference_object = reference_object, marker = marker)
    },
    .package = "Shennong"
  )

  expect_warning(
    deconvolution <- sn_run_spatial_deconvolution(object, marker = "deconv"),
    "deprecated"
  )
  expect_identical(deconvolution$marker, "deconv")

  reference <- make_spatial_test_object()
  expect_warning(
    mapping <- sn_run_spatial_mapping(
      object,
      reference_object = reference,
      marker = "mapping"
    ),
    "deprecated"
  )
  expect_identical(mapping$reference_object, reference)
  expect_identical(mapping$marker, "mapping")
})
