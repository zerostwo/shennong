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
  methods <- methods %||% manifest$methods %||% names(manifest$results)
  methods <- unique(as.character(methods))
  unknown <- setdiff(methods, names(manifest$results))
  if (length(unknown)) stop("Unknown comparison method(s): ", paste(unknown, collapse = ", "), ".", call. = FALSE)
  rows <- lapply(methods, function(method) {
    info <- manifest$results[[method]]
    reduction <- info$integration_reduction %||% if (identical(method, "unintegrated")) "pca" else method
    graph_only <- identical(method, "bbknn") || is.null(reduction) || !reduction %in% names(object@reductions)
    label_by <- info$integration_control$label_by %||% info$integration$label_by %||% NA_character_
    data.frame(
      method = method,
      reduction = if (is.null(reduction)) NA_character_ else reduction,
      graph_only = graph_only,
      training_label_by = as.character(label_by),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
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
                                         seed) {
  metadata <- object[[]]
  missing_columns <- setdiff(c(batch_by, label_by), colnames(metadata))
  if (length(missing_columns)) {
    stop("Missing metadata column(s): ", paste(missing_columns, collapse = ", "), ".", call. = FALSE)
  }
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

  embedding_methods <- method_table$method[!method_table$graph_only]
  for (method in embedding_methods) {
    reduction <- method_table$reduction[method_table$method == method][[1L]]
    embedding <- Seurat::Embeddings(object, reduction = reduction)[cells, , drop = FALSE]
    utils::write.csv(embedding, file.path(input_dir, "embeddings", paste0(method, ".csv")))
  }
  list(cells = cells, features = features, embedding_methods = embedding_methods)
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
  summary <- merge(summary, method_table, by = "method", all.x = TRUE, sort = FALSE)
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
#' method is evaluated on the same cells and features.
#'
#' @param object A Seurat object containing `misc$integration_comparison`.
#' @param batch_by Metadata column containing batch identities.
#' @param label_by Metadata column containing biological labels.
#' @param methods Optional methods to evaluate. Defaults to all registered methods.
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
  if (!"unintegrated" %in% method_table$method) {
    stop("The comparison must include the `unintegrated` PCA baseline.", call. = FALSE)
  }
  comparable <- method_table$method[!method_table$graph_only]
  if (length(comparable) < 2L) stop("At least two embedding-based methods are required.", call. = FALSE)

  runtime_dir <- .sn_shennong_runtime_dir(backend_control$runtime_dir %||% NULL)
  run_dir <- backend_control$run_dir %||% .sn_default_scvi_run_dir("scib_metrics", runtime_dir)
  input_dir <- file.path(run_dir, "input")
  output_dir <- file.path(run_dir, "output")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  input <- .sn_write_scib_metrics_input(
    object, input_dir, method_table, batch_by, label_by, assay, normalized_layer,
    features, max_cells, seed
  )
  selected <- .sn_resolve_scvi_accelerator(accelerator)
  pixi_environment <- backend_control$environment %||% selected$environment
  config <- list(
    batch_key = "batch", label_key = "label",
    embedding_methods = input$embedding_methods,
    baseline_method = "unintegrated",
    n_jobs = as.integer(n_jobs), progress_bar = isTRUE(verbose),
    min_max_scale = backend_control$min_max_scale %||% FALSE,
    solver = backend_control$solver %||% "arpack",
    bio_conservation_metrics = backend_control$bio_conservation_metrics %||% list(),
    batch_correction_metrics = backend_control$batch_correction_metrics %||% list(),
    accelerator = selected$requested
  )
  config_path <- .sn_write_json_file(config, file.path(run_dir, "config.json"))
  runner <- backend_control$runner %||% NULL
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
  parameters <- list(
    batch_by = batch_by, label_by = label_by, methods = method_table$method,
    assay = assay, normalized_layer = normalized_layer, max_cells = max_cells,
    accelerator = selected$requested, pixi_environment = pixi_environment,
    n_jobs = as.integer(n_jobs), seed = as.integer(seed), run_dir = run_dir,
    selected_cells = input$cells, selected_features = input$features
  )
  result <- .sn_import_scib_metrics_result(output_dir, method_table, label_by, name, parameters)
  if (!isTRUE(return_object)) return(result)
  sn_store_result(object, "integration_benchmark", name, result)
}
