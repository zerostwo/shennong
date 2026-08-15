.sn_scib_metrics_script_path <- function(script = NULL) {
  if (!is.null(script) && nzchar(script)) {
    script <- path.expand(script)
    if (!file.exists(script)) stop("`backend_control$script` does not exist: ", script, call. = FALSE)
    return(normalizePath(script, winslash = "/", mustWork = TRUE))
  }
  installed <- system.file("pixi", "scib-metrics", "scripts", "scib_metrics_benchmark.py", package = "Shennong")
  if (nzchar(installed) && file.exists(installed)) {
    return(normalizePath(installed, winslash = "/", mustWork = TRUE))
  }
  source_path <- file.path(getwd(), "inst", "pixi", "scib-metrics", "scripts", "scib_metrics_benchmark.py")
  if (file.exists(source_path)) return(normalizePath(source_path, winslash = "/", mustWork = TRUE))
  stop("Could not locate Shennong's scib-metrics runner.", call. = FALSE)
}

.sn_scib_comparison_manifest <- function(object) {
  manifest <- object@misc$integration_comparison %||% NULL
  if (is.null(manifest) || !is.list(manifest$results) || length(manifest$results) < 1L) {
    stop(
      "No multi-method integration manifest was found. Run `sn_run_cluster()` with multiple ",
      "`integration_method` values first.",
      call. = FALSE
    )
  }
  manifest
}

.sn_scib_method_table <- function(object, manifest, methods = NULL) {
  run_ids <- names(manifest$results)
  if (!is.null(methods)) {
    requested <- unique(as.character(methods))
    known_methods <- unique(vapply(
      manifest$results,
      function(info) info$method %||% NA_character_,
      character(1)
    ))
    unknown <- setdiff(requested, c(run_ids, known_methods))
    if (length(unknown)) {
      stop("Unknown comparison run or method(s): ", paste(unknown, collapse = ", "), ".", call. = FALSE)
    }
    run_ids <- run_ids[vapply(
      run_ids,
      function(run_id) run_id %in% requested || (manifest$results[[run_id]]$method %||% run_id) %in% requested,
      logical(1)
    )]
  }
  rows <- lapply(run_ids, function(run_id) {
    info <- manifest$results[[run_id]]
    method <- info$method %||% run_id
    reduction <- info$integration_reduction %||% if (identical(method, "unintegrated")) "pca" else method
    graph_only <- identical(method, "bbknn") || is.null(reduction) || !reduction %in% names(object@reductions)
    label_by <- info$integration_control$label_by %||% info$integration$label_by %||% NA_character_
    data.frame(
      run_id = run_id,
      embedding_id = info$embedding_id %||% run_id,
      preprocess_id = info$preprocess_id %||% "shared",
      method = method,
      reduction = if (is.null(reduction)) NA_character_ else reduction,
      graph_only = graph_only,
      training_label_by = as.character(label_by),
      stringsAsFactors = FALSE
    )
  })
  table <- do.call(rbind, rows)
  table$benchmark_id <- table$embedding_id
  table$benchmark_representative <- !duplicated(table$embedding_id)
  table$baseline_run_id <- NA_character_
  for (preprocess_id in unique(table$preprocess_id)) {
    in_group <- table$preprocess_id == preprocess_id
    baselines <- table$run_id[in_group & table$method == "unintegrated" & !table$graph_only]
    if (length(baselines) > 0L) {
      table$baseline_run_id[in_group] <- baselines[[1L]]
    }
  }
  table
}

.sn_write_scib_metrics_input <- function(object,
                                         input_dir,
                                         method_table,
                                         batch_by,
                                         label_by,
                                         assay,
                                         normalized_layer,
                                         features,
                                         max_cells,
                                         seed,
                                         cells = NULL) {
  metadata <- object[[]]
  missing_columns <- setdiff(c(batch_by, label_by), colnames(metadata))
  if (length(missing_columns)) {
    stop("Missing metadata column(s): ", paste(missing_columns, collapse = ", "), ".", call. = FALSE)
  }
  if (is.null(cells)) {
    keep <- stats::complete.cases(metadata[, c(batch_by, label_by), drop = FALSE])
    if (!all(keep)) warning("Dropping ", sum(!keep), " cells with missing batch/label values.", call. = FALSE)
    cells <- rownames(metadata)[keep]
    if (!is.null(max_cells) && length(cells) > max_cells) {
      strata <- interaction(metadata[cells, batch_by], metadata[cells, label_by], drop = TRUE, lex.order = TRUE)
      sampling_metadata <- data.frame(.sn_scib_stratum = strata, row.names = cells)
      cells <- .sn_subsample_metric_cells(
        cells = cells,
        metadata = sampling_metadata,
        max_cells = as.integer(max_cells),
        stratify_by = ".sn_scib_stratum",
        seed = seed
      )
    }
  } else {
    cells <- intersect(as.character(cells), rownames(metadata))
    if (length(cells) < 2L) stop("The shared scIB cell set is no longer present in the object.", call. = FALSE)
  }
  if (length(unique(metadata[cells, batch_by])) < 2L) stop("At least two batches are required.", call. = FALSE)
  if (length(unique(metadata[cells, label_by])) < 2L) stop("At least two biological labels are required.", call. = FALSE)

  matrix <- .sn_get_seurat_layer_data(object = object, assay = assay, layer = normalized_layer)
  features <- features %||% SeuratObject::VariableFeatures(object[[assay]])
  features <- unique(intersect(as.character(features), rownames(matrix)))
  if (length(features) < 2L) stop("At least two shared normalized features are required.", call. = FALSE)
  matrix <- .sn_as_sparse_matrix(matrix[features, cells, drop = FALSE])

  dir.create(file.path(input_dir, "embeddings"), recursive = TRUE, showWarnings = FALSE)
  Matrix::writeMM(matrix, file.path(input_dir, "normalized.mtx"))
  utils::write.csv(data.frame(cell_id = cells), file.path(input_dir, "cells.csv"), row.names = FALSE)
  utils::write.csv(data.frame(feature_id = features), file.path(input_dir, "features.csv"), row.names = FALSE)
  obs <- data.frame(
    cell_id = cells,
    batch = as.character(metadata[cells, batch_by]),
    label = as.character(metadata[cells, label_by]),
    stringsAsFactors = FALSE
  )
  utils::write.csv(obs, file.path(input_dir, "obs.csv"), row.names = FALSE)

  benchmark_table <- method_table[!method_table$graph_only & method_table$benchmark_representative, , drop = FALSE]
  embedding_methods <- benchmark_table$benchmark_id
  for (row_index in seq_len(nrow(benchmark_table))) {
    method <- benchmark_table$benchmark_id[[row_index]]
    reduction <- benchmark_table$reduction[[row_index]]
    embedding <- Seurat::Embeddings(object, reduction = reduction)[cells, , drop = FALSE]
    utils::write.csv(embedding, file.path(input_dir, "embeddings", paste0(method, ".csv")))
  }
  list(cells = cells, features = features, embedding_methods = embedding_methods)
}

.sn_scib_embedding_rank <- function(table, score, group = NULL) {
  key_columns <- c("embedding_id", group)
  representative <- !duplicated(table[, key_columns, drop = FALSE])
  unique_table <- table[representative, , drop = FALSE]
  unique_scores <- score[representative]
  if (is.null(group)) {
    unique_ranks <- rank(-unique_scores, ties.method = "min", na.last = "keep")
    return(unique_ranks[match(table$embedding_id, unique_table$embedding_id)])
  }
  unique_ranks <- stats::ave(
    -unique_scores,
    unique_table[[group]],
    FUN = function(x) rank(x, ties.method = "min", na.last = "keep")
  )
  unique_keys <- paste(unique_table[[group]], unique_table$embedding_id, sep = "\r")
  keys <- paste(table[[group]], table$embedding_id, sep = "\r")
  unique_ranks[match(keys, unique_keys)]
}

.sn_import_scib_metrics_result <- function(output_dir, method_table, label_by, name, parameters) {
  required <- file.path(output_dir, c("summary.csv", "metrics.csv", "ranking.csv", "manifest.json"))
  if (any(!file.exists(required))) {
    stop("scib-metrics output is incomplete: ", paste(basename(required[!file.exists(required)]), collapse = ", "), ".", call. = FALSE)
  }
  summary <- utils::read.csv(required[[1]], check.names = FALSE)
  metrics <- utils::read.csv(required[[2]], check.names = FALSE)
  ranking <- utils::read.csv(required[[3]], check.names = FALSE)
  backend_manifest <- jsonlite::read_json(required[[4]], simplifyVector = TRUE)
  names(summary)[names(summary) == "method"] <- "benchmark_id"
  comparable_table <- method_table[!method_table$graph_only, , drop = FALSE]
  summary <- merge(summary, comparable_table, by = "benchmark_id", all.x = TRUE, sort = FALSE)
  summary$supervised_label_leakage <- !is.na(summary$training_label_by) & summary$training_label_by == label_by
  graph_only <- method_table[method_table$graph_only, , drop = FALSE]
  if (nrow(graph_only)) {
    graph_only$comparable <- FALSE
    graph_only$reason <- "graph_only_method"
    graph_only$supervised_label_leakage <- !is.na(graph_only$training_label_by) &
      graph_only$training_label_by == label_by
    summary$comparable <- TRUE
    summary$reason <- NA_character_
    missing_graph_columns <- setdiff(colnames(summary), colnames(graph_only))
    for (column in missing_graph_columns) graph_only[[column]] <- NA
    graph_only <- graph_only[, colnames(summary), drop = FALSE]
    summary <- rbind(summary, graph_only)
  } else {
    summary$comparable <- TRUE
    summary$reason <- NA_character_
  }
  warnings <- character()
  if (any(summary$supervised_label_leakage %||% FALSE, na.rm = TRUE)) {
    warnings <- c(warnings, "One or more methods used the evaluation label during supervised training; compare them separately from unsupervised methods.")
  }
  if (nrow(graph_only)) warnings <- c(warnings, "Graph-only methods were excluded from the embedding benchmark.")
  jax_devices <- tolower(paste(backend_manifest$jax_devices %||% character(), collapse = ","))
  gpu_used <- grepl("gpu|cuda", jax_devices) || identical(backend_manifest$jax_default_backend %||% NULL, "gpu")
  if (identical(parameters$pixi_environment, "gpu") && !isTRUE(gpu_used)) {
    warnings <- c(
      warnings,
      "The GPU environment was selected, but the scib-metrics manifest did not report a JAX GPU device."
    )
  }
  names(metrics)[names(metrics) == "method"] <- "benchmark_id"
  metrics <- merge(metrics, comparable_table, by = "benchmark_id", all.x = TRUE, sort = FALSE)
  names(ranking)[names(ranking) == "method"] <- "benchmark_id"
  ranking <- merge(ranking, comparable_table, by = "benchmark_id", all.x = TRUE, sort = FALSE)
  summary$rank_within_preprocess <- .sn_scib_embedding_rank(
    summary, summary$Total, group = "preprocess_id"
  )
  ranking$rank_within_preprocess <- .sn_scib_embedding_rank(
    ranking, ranking$score, group = "preprocess_id"
  )
  ranking <- ranking[order(ranking$rank_within_preprocess, ranking$run_id), , drop = FALSE]
  .sn_new_analysis_result(
    analysis_type = "integration_benchmark",
    name = name,
    method = "scib_metrics",
    backend = "scib-metrics",
    input = list(cells = backend_manifest$n_cells, features = backend_manifest$n_features),
    parameters = parameters,
    tables = list(primary = summary, summary = summary, metrics = metrics, ranking = ranking),
    diagnostics = list(backend_manifest = backend_manifest, method_table = method_table),
    warnings = warnings,
    provenance = .sn_analysis_provenance(
      result = list(provenance = list(package_versions = list(
        Shennong = tryCatch(as.character(utils::packageVersion("Shennong")), error = function(e) NA_character_),
        scib_metrics = backend_manifest$scib_metrics_version %||% NA_character_,
        jax = backend_manifest$jax_version %||% NA_character_
      ))),
      random_seed = parameters$seed,
      capture_acceleration = FALSE
    )
  )
}

#' Compare integration methods with scib-metrics
#'
#' Benchmarks the integration reductions registered by a multi-method
#' [sn_run_cluster()] call. UMAP and t-SNE coordinates are never used as metric
#' inputs. The normalized expression matrix is exported sparsely and every
#' method is evaluated on the same cells. Parameter-grid objects are grouped by
#' preprocessing identity, each group uses its matching unintegrated baseline,
#' and each unique native embedding is scored once before its score is mapped
#' back to all cluster-resolution runs that share it.
#'
#' @param object A Seurat object containing `misc$integration_comparison`.
#' @param batch_by Metadata column containing batch identities.
#' @param label_by Metadata column containing biological labels.
#' @param methods Optional method names or run IDs to evaluate. A method name
#'   selects all of its registered parameter-grid runs.
#' @param assay RNA assay containing normalized expression.
#' @param normalized_layer Normalized, unintegrated expression layer used by
#'   scib-metrics for its baseline preparation.
#' @param features Optional common features. Defaults to variable features.
#' @param max_cells Optional shared stratified cell cap.
#' @param accelerator One of `"auto"`, `"cpu"`, or `"gpu"`.
#' @param n_jobs Neighbor-search workers.
#' @param name Stored result name.
#' @param return_object Return the updated Seurat object when `TRUE`, otherwise
#'   return the analysis-result object.
#' @param backend_control Advanced pixi, runner, metric, and filesystem controls.
#' @param seed Sampling seed.
#' @param verbose Show backend progress.
#'
#' @return A Seurat object with an `integration_benchmark` result, or the result.
#' @export
sn_compare_integrations <- function(object,
                                    batch_by = NULL,
                                    label_by,
                                    methods = NULL,
                                    assay = "RNA",
                                    normalized_layer = "data",
                                    features = NULL,
                                    max_cells = 100000L,
                                    accelerator = c("auto", "cpu", "gpu"),
                                    n_jobs = 1L,
                                    name = "integration_benchmark",
                                    return_object = TRUE,
                                    backend_control = list(),
                                    seed = 717L,
                                    verbose = TRUE) {
  check_installed(c("Seurat", "Matrix", "jsonlite"))
  if (!inherits(object, "Seurat")) stop("`object` must be a Seurat object.", call. = FALSE)
  if (!is.list(backend_control)) stop("`backend_control` must be a list.", call. = FALSE)
  accelerator <- match.arg(accelerator)
  manifest <- .sn_scib_comparison_manifest(object)
  batch_by <- batch_by %||% manifest$batch_by
  if (is.null(batch_by) || !nzchar(batch_by)) stop("`batch_by` must be supplied.", call. = FALSE)
  if (missing(label_by) || is.null(label_by) || !nzchar(label_by)) stop("`label_by` must be supplied.", call. = FALSE)
  method_table <- .sn_scib_method_table(object, manifest, methods)
  preprocess_ids <- unique(method_table$preprocess_id)
  for (preprocess_id in preprocess_ids) {
    group <- method_table[method_table$preprocess_id == preprocess_id, , drop = FALSE]
    if (!any(group$method == "unintegrated" & !group$graph_only)) {
      stop(
        "Each preprocessing group must include an `unintegrated` PCA baseline; missing for `",
        preprocess_id, "`.",
        call. = FALSE
      )
    }
    comparable_embeddings <- unique(group$embedding_id[!group$graph_only])
    if (length(comparable_embeddings) < 2L) {
      stop("At least two embedding-based methods are required in preprocessing group `", preprocess_id, "`.", call. = FALSE)
    }
  }

  runtime_dir <- .sn_shennong_runtime_dir(backend_control$runtime_dir %||% NULL)
  run_dir <- backend_control$run_dir %||% .sn_default_scvi_run_dir("scib_metrics", runtime_dir)
  selected <- .sn_resolve_scvi_accelerator(accelerator)
  pixi_environment <- backend_control$environment %||% selected$environment
  runner <- backend_control$runner %||% NULL
  shared_cells <- NULL
  selected_features <- list()
  group_results <- list()
  for (group_index in seq_along(preprocess_ids)) {
    preprocess_id <- preprocess_ids[[group_index]]
    group_table <- method_table[method_table$preprocess_id == preprocess_id, , drop = FALSE]
    group_dir <- if (length(preprocess_ids) == 1L) run_dir else file.path(run_dir, preprocess_id)
    input_dir <- file.path(group_dir, "input")
    output_dir <- file.path(group_dir, "output")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    baseline_row <- which(group_table$method == "unintegrated" & !group_table$graph_only)[[1L]]
    baseline_run_id <- group_table$run_id[[baseline_row]]
    group_features <- features %||% manifest$results[[baseline_run_id]]$input_features %||% NULL
    group_layer <- manifest$results[[baseline_run_id]]$normalized_layer %||% normalized_layer
    input <- .sn_write_scib_metrics_input(
      object, input_dir, group_table, batch_by, label_by, assay, group_layer,
      group_features, max_cells, seed, cells = shared_cells
    )
    shared_cells <- input$cells
    selected_features[[preprocess_id]] <- input$features
    baseline_method <- group_table$benchmark_id[[baseline_row]]
    config <- list(
      batch_key = "batch", label_key = "label",
      embedding_methods = input$embedding_methods,
      baseline_method = baseline_method,
      n_jobs = as.integer(n_jobs), progress_bar = isTRUE(verbose),
      min_max_scale = backend_control$min_max_scale %||% FALSE,
      solver = backend_control$solver %||% "arpack",
      bio_conservation_metrics = backend_control$bio_conservation_metrics %||% list(),
      batch_correction_metrics = backend_control$batch_correction_metrics %||% list(),
      accelerator = selected$requested
    )
    config_path <- .sn_write_json_file(config, file.path(group_dir, "config.json"))
    if (is.function(runner)) {
      runner(input_dir = input_dir, output_dir = output_dir, config = config, config_path = config_path)
    } else {
      paths <- sn_pixi_paths("scib-metrics", runtime_dir = runtime_dir)
      manifest_path <- .sn_prepare_scvi_pixi_project(
        project_dir = backend_control$pixi_project %||% paths$project_dir,
        environment = "scib-metrics",
        manifest_path = backend_control$manifest_path %||% NULL,
        manifest_lines = backend_control$manifest_lines %||% NULL,
        overwrite = isTRUE(backend_control$overwrite_manifest),
        cuda_version = backend_control$cuda_version %||% selected$cuda_version,
        platforms = backend_control$platforms %||% NULL
      )
      .sn_execute_scvi_pixi(
        pixi = backend_control$pixi %||% NULL,
        manifest_path = manifest_path,
        script = .sn_scib_metrics_script_path(backend_control$script %||% NULL),
        input_dir = input_dir, output_dir = output_dir, config_path = config_path,
        environment = pixi_environment, pixi_home = backend_control$pixi_home %||% paths$pixi_home,
        install_pixi = backend_control$install_pixi %||% TRUE,
        pixi_version = backend_control$pixi_version %||% "latest",
        pixi_download_url = backend_control$pixi_download_url %||% NULL,
        verbose = verbose, backend_label = "scib-metrics"
      )
    }
    group_parameters <- list(
      batch_by = batch_by, label_by = label_by, methods = group_table$method,
      assay = assay, normalized_layer = group_layer, max_cells = max_cells,
      accelerator = selected$requested, pixi_environment = pixi_environment,
      n_jobs = as.integer(n_jobs), seed = as.integer(seed), run_dir = group_dir,
      preprocess_id = preprocess_id, baseline_run_id = baseline_run_id,
      selected_cells = input$cells, selected_features = input$features
    )
    group_results[[preprocess_id]] <- .sn_import_scib_metrics_result(
      output_dir, group_table, label_by, name, group_parameters
    )
  }
  result <- group_results[[1L]]
  if (length(group_results) > 1L) {
    result$tables$summary <- do.call(rbind, lapply(group_results, function(x) x$tables$summary))
    result$tables$primary <- result$tables$summary
    result$tables$metrics <- do.call(rbind, lapply(group_results, function(x) x$tables$metrics))
    result$tables$ranking <- do.call(rbind, lapply(group_results, function(x) x$tables$ranking))
    result$warnings <- unique(unlist(lapply(group_results, function(x) x$warnings), use.names = FALSE))
    result$diagnostics$backend_manifest <- lapply(group_results, function(x) x$diagnostics$backend_manifest)
    result$diagnostics$method_table <- method_table
  }
  comparable_summary <- !result$tables$summary$graph_only
  result$tables$summary$rank_overall <- NA_real_
  result$tables$summary$rank_overall[comparable_summary] <- .sn_scib_embedding_rank(
    result$tables$summary[comparable_summary, , drop = FALSE],
    result$tables$summary$Total[comparable_summary]
  )
  result$tables$primary <- result$tables$summary
  result$tables$ranking$rank_overall <- .sn_scib_embedding_rank(
    result$tables$ranking,
    result$tables$ranking$score
  )
  result$tables$ranking <- result$tables$ranking[
    order(result$tables$ranking$rank_overall, result$tables$ranking$run_id),
    , drop = FALSE
  ]
  result$parameters <- list(
    batch_by = batch_by, label_by = label_by, methods = method_table$method,
    run_ids = method_table$run_id, preprocess_ids = preprocess_ids,
    assay = assay, normalized_layer = normalized_layer, max_cells = max_cells,
    accelerator = selected$requested, pixi_environment = pixi_environment,
    n_jobs = as.integer(n_jobs), seed = as.integer(seed), run_dir = run_dir,
    selected_cells = shared_cells,
    selected_features = if (length(selected_features) == 1L) selected_features[[1L]] else selected_features
  )
  result$input$cells <- length(shared_cells)
  result$input$features <- if (length(selected_features) == 1L) {
    length(selected_features[[1L]])
  } else {
    vapply(selected_features, length, integer(1))
  }
  if (!isTRUE(return_object)) return(result)
  sn_store_result(object, "integration_benchmark", name, result)
}
