library(testthat)

make_integration_benchmark_object <- function() {
  skip_if_not_installed("Seurat")
  set.seed(41)
  counts <- Matrix::Matrix(matrix(rpois(20 * 24, 3), nrow = 20), sparse = TRUE)
  rownames(counts) <- paste0("gene", seq_len(nrow(counts)))
  colnames(counts) <- paste0("cell", seq_len(ncol(counts)))
  object <- SeuratObject::CreateSeuratObject(counts)
  object$sample <- rep(c("pbmc1k", "pbmc3k", "pbmc4k"), each = 8)
  object$cell_type <- rep(c("T", "B"), 12)
  object <- Seurat::NormalizeData(object, verbose = FALSE)
  SeuratObject::VariableFeatures(object[["RNA"]]) <- rownames(object)[1:12]
  embeddings <- matrix(rnorm(24 * 5), nrow = 24, dimnames = list(colnames(object), paste0("D", 1:5)))
  for (reduction in c("pca", "harmony", "scanvi")) {
    object[[reduction]] <- Seurat::CreateDimReducObject(
      embeddings = embeddings + match(reduction, c("pca", "harmony", "scanvi")),
      key = paste0(toupper(reduction), "_"), assay = "RNA"
    )
  }
  result <- function(method, reduction, label_by = NULL) list(
    method = method,
    integration_reduction = reduction,
    graph_names = paste0(method, "_snn"),
    integration_control = if (is.null(label_by)) list() else list(label_by = label_by),
    performance = list(
      workflow = list(elapsed_seconds = 12, peak_memory_mb = 256),
      integration = list(
        elapsed_seconds = 7, peak_memory_mb = 192,
        peak_r_memory_mb = 128, backend_peak_rss_mb = 192
      ),
      reused_embedding = FALSE
    )
  )
  object@misc$integration_comparison <- list(
    methods = c("unintegrated", "harmony", "scanvi", "bbknn"),
    batch_by = "sample", assay = "RNA", layer = "counts",
    results = list(
      unintegrated = result("unintegrated", "pca"),
      harmony = result("harmony", "harmony"),
      scanvi = result("scanvi", "scanvi", "cell_type"),
      bbknn = result("bbknn", "pca")
    )
  )
  object
}

fake_scib_metrics_runner <- function(input_dir, output_dir, config, config_path) {
  expect_true(file.exists(file.path(input_dir, "normalized.mtx")))
  expect_equal(config$embedding_methods, c("unintegrated", "harmony", "scanvi"))
  expect_equal(config$baseline_method, "unintegrated")
  n_cells <- nrow(utils::read.csv(file.path(input_dir, "cells.csv")))
  n_features <- nrow(utils::read.csv(file.path(input_dir, "features.csv")))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  summary <- data.frame(
    method = config$embedding_methods,
    `Bio conservation` = c(0.7, 0.8, 0.9),
    `Batch correction` = c(0.2, 0.9, 0.8),
    Total = c(0.5, 0.84, 0.86),
    check.names = FALSE
  )
  metrics <- data.frame(
    method = rep(config$embedding_methods, each = 2),
    metric = rep(c("Silhouette label", "iLISI"), 3),
    score = seq(0.2, 0.7, length.out = 6),
    category = rep(c("Bio conservation", "Batch correction"), 3)
  )
  ranking <- data.frame(rank = 1:3, method = c("scanvi", "harmony", "unintegrated"), score = c(0.86, 0.84, 0.5))
  utils::write.csv(summary, file.path(output_dir, "summary.csv"), row.names = FALSE)
  utils::write.csv(metrics, file.path(output_dir, "metrics.csv"), row.names = FALSE)
  utils::write.csv(ranking, file.path(output_dir, "ranking.csv"), row.names = FALSE)
  jsonlite::write_json(
    list(
      scib_metrics_version = "0.5.8", jax_version = "0.test",
      jax_devices = "TFRT_CPU_0", n_cells = n_cells, n_features = n_features
    ),
    file.path(output_dir, "manifest.json"), auto_unbox = TRUE
  )
}

test_that("sn_compare_integrations stores a discoverable benchmark result", {
  object <- make_integration_benchmark_object()
  run_dir <- tempfile("shennong-scib-test-")
  updated <- sn_compare_integrations(
    object,
    label_by = "cell_type",
    accelerator = "cpu",
    backend_control = list(run_dir = run_dir, runner = fake_scib_metrics_runner),
    verbose = FALSE
  )
  result <- sn_get_result(updated, "integration_benchmark", "integration_benchmark")
  expect_true(sn_validate_result(result, error = FALSE)$valid)
  expect_equal(result$tables$ranking$method[[1]], "scanvi")
  expect_true(result$tables$summary$supervised_label_leakage[result$tables$summary$method == "scanvi"])
  expect_false(result$tables$summary$comparable[result$tables$summary$method == "bbknn"])
  expect_match(result$tables$summary$reason[result$tables$summary$method == "bbknn"], "graph_only")
  expect_equal(result$input$cells, 24L)
  expect_equal(result$input$features, 12L)
  expect_true(all(result$tables$summary$elapsed_seconds == 12))
  expect_true(all(result$tables$ranking$integration_peak_memory_mb == 192))
  listing <- sn_list_results(updated, type = "integration_benchmark")
  expect_equal(listing$name, "integration_benchmark")
})

test_that("sn_compare_integrations returns results and uses shared stratified cells", {
  object <- make_integration_benchmark_object()
  result <- sn_compare_integrations(
    object,
    label_by = "cell_type",
    max_cells = 12,
    accelerator = "cpu",
    return_object = FALSE,
    backend_control = list(run_dir = tempfile("shennong-scib-test-"), runner = fake_scib_metrics_runner),
    verbose = FALSE
  )
  expect_equal(length(result$parameters$selected_cells), 12L)
  selected <- object[[]][result$parameters$selected_cells, , drop = FALSE]
  expect_equal(length(unique(selected$sample)), 3L)
  expect_equal(length(unique(selected$cell_type)), 2L)
})

test_that("sn_compare_integrations records a selected GPU environment that fell back to CPU", {
  object <- make_integration_benchmark_object()
  result <- sn_compare_integrations(
    object,
    label_by = "cell_type",
    accelerator = "cpu",
    return_object = FALSE,
    backend_control = list(
      environment = "gpu",
      run_dir = tempfile("shennong-scib-test-"),
      runner = fake_scib_metrics_runner
    ),
    verbose = FALSE
  )
  expect_true(any(grepl("did not report a JAX GPU device", result$warnings, fixed = TRUE)))
})

test_that("sn_compare_integrations benchmarks unique embeddings within preprocessing groups", {
  object <- make_integration_benchmark_object()
  base_embedding <- Seurat::Embeddings(object, "pca")
  for (reduction in c("pca.p1", "harmony.p1", "pca.p2", "harmony.p2")) {
    object[[reduction]] <- Seurat::CreateDimReducObject(
      embeddings = base_embedding + match(reduction, c("pca.p1", "harmony.p1", "pca.p2", "harmony.p2")),
      key = paste0(gsub("[^A-Za-z0-9]", "", reduction), "_"),
      assay = "RNA"
    )
  }
  result_info <- function(run_id, embedding_id, preprocess_id, method, reduction, features) {
    list(
      run_id = run_id,
      embedding_id = embedding_id,
      preprocess_id = preprocess_id,
      method = method,
      integration_reduction = reduction,
      input_features = features,
      normalized_layer = "data",
      integration_control = list()
    )
  }
  object@misc$integration_comparison <- list(
    schema_version = "2.0.0",
    methods = c("unintegrated", "harmony"),
    batch_by = "sample", assay = "RNA", layer = "counts",
    results = list(
      u1 = result_info("u1", "uemb1", "p1", "unintegrated", "pca.p1", rownames(object)[1:8]),
      h1a = result_info("h1a", "hemb1", "p1", "harmony", "harmony.p1", rownames(object)[1:8]),
      h1b = result_info("h1b", "hemb1", "p1", "harmony", "harmony.p1", rownames(object)[1:8]),
      u2 = result_info("u2", "uemb2", "p2", "unintegrated", "pca.p2", rownames(object)[1:10]),
      h2 = result_info("h2", "hemb2", "p2", "harmony", "harmony.p2", rownames(object)[1:10])
    )
  )

  calls <- list()
  grid_runner <- function(input_dir, output_dir, config, config_path) {
    calls[[length(calls) + 1L]] <<- list(
      methods = config$embedding_methods,
      baseline = config$baseline_method,
      cells = utils::read.csv(file.path(input_dir, "cells.csv"))$cell_id
    )
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    methods <- config$embedding_methods
    scores <- seq(0.5, 0.8, length.out = length(methods))
    utils::write.csv(
      data.frame(method = methods, `Bio conservation` = scores, `Batch correction` = rev(scores), Total = scores, check.names = FALSE),
      file.path(output_dir, "summary.csv"), row.names = FALSE
    )
    utils::write.csv(
      data.frame(method = methods, metric = "Silhouette label", score = scores, category = "Bio conservation"),
      file.path(output_dir, "metrics.csv"), row.names = FALSE
    )
    utils::write.csv(
      data.frame(rank = seq_along(methods), method = methods, score = scores),
      file.path(output_dir, "ranking.csv"), row.names = FALSE
    )
    jsonlite::write_json(
      list(scib_metrics_version = "0.test", jax_version = "0.test", jax_devices = "TFRT_CPU_0", n_cells = 24L, n_features = 8L),
      file.path(output_dir, "manifest.json"), auto_unbox = TRUE
    )
  }

  result <- sn_compare_integrations(
    object,
    label_by = "cell_type",
    accelerator = "cpu",
    return_object = FALSE,
    backend_control = list(run_dir = tempfile("shennong-scib-grid-"), runner = grid_runner),
    verbose = FALSE
  )
  expect_length(calls, 2L)
  expect_equal(calls[[1]]$cells, calls[[2]]$cells)
  expect_true(all(vapply(calls, function(x) length(x$methods) == 2L, logical(1))))
  expect_true(all(vapply(calls, function(x) x$baseline %in% x$methods, logical(1))))
  expect_setequal(result$tables$summary$run_id, c("u1", "h1a", "h1b", "u2", "h2"))
  expect_equal(result$tables$summary$Total[result$tables$summary$run_id == "h1a"], result$tables$summary$Total[result$tables$summary$run_id == "h1b"])
  expect_equal(result$tables$summary$rank_within_preprocess[result$tables$summary$run_id == "h1a"], result$tables$summary$rank_within_preprocess[result$tables$summary$run_id == "h1b"])
  expect_setequal(result$tables$summary$rank_within_preprocess, c(1, 1, 1, 2, 2))
  expect_true(all(c("rank_within_preprocess", "rank_overall", "baseline_run_id") %in% colnames(result$tables$summary)))
  expect_equal(length(result$parameters$selected_features$p1), 8L)
  expect_equal(length(result$parameters$selected_features$p2), 10L)
})

test_that("sn_compare_integrations validates its comparison contract", {
  object <- make_integration_benchmark_object()
  object@misc$integration_comparison <- NULL
  expect_error(sn_compare_integrations(object, label_by = "cell_type"), "No multi-method")

  object <- make_integration_benchmark_object()
  expect_error(
    sn_compare_integrations(object, label_by = "missing", backend_control = list(runner = fake_scib_metrics_runner)),
    "Missing metadata"
  )
  expect_error(
    sn_compare_integrations(
      object, label_by = "cell_type", methods = c("harmony", "scanvi"),
      backend_control = list(runner = fake_scib_metrics_runner)
    ),
    "unintegrated"
  )
})

test_that("scib-metrics pixi environment is discoverable and GPU-aware", {
  expect_true("scib-metrics" %in% sn_list_pixi_environments())
  expect_match(sn_pixi_config_path("scib_metrics"), "scib-metrics/pixi.toml$")
  expect_true(Shennong:::.sn_pixi_gpu_aware_environment("scib-metrics"))
})
