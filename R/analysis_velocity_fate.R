.sn_trajectory_pixi_script <- function() {
  installed <- system.file("pixi", "trajectory", "scripts", "trajectory_run.py", package = "Shennong")
  if (nzchar(installed) && file.exists(installed)) return(installed)
  source <- file.path(getwd(), "inst", "pixi", "trajectory", "scripts", "trajectory_run.py")
  if (file.exists(source)) return(normalizePath(source, winslash = "/", mustWork = TRUE))
  stop("Could not locate the bundled trajectory backend script.", call. = FALSE)
}

.sn_velocity_embedding <- function(object, reduction = NULL, dims = 1:2) {
  available <- names(object@reductions)
  reduction <- reduction %||% if ("umap" %in% available) {
    "umap"
  } else if (length(available) > 0L) {
    available[[1]]
  } else {
    NULL
  }
  if (is_null(reduction) || !reduction %in% names(object@reductions)) {
    stop("Velocity requires a stored dimensional reduction.", call. = FALSE)
  }
  embedding <- SeuratObject::Embeddings(object[[reduction]])
  dims <- as.integer(dims)
  if (length(dims) < 2L || any(dims < 1L) || max(dims) > ncol(embedding)) {
    stop("`dims` must select at least two available reduction dimensions.", call. = FALSE)
  }
  list(matrix = embedding[, dims, drop = FALSE], reduction = reduction, dims = dims)
}

.sn_write_velocity_input <- function(object,
                                     spliced_assay,
                                     spliced_layer,
                                     unspliced_assay,
                                     unspliced_layer,
                                     embedding,
                                     input_dir) {
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  spliced <- .sn_get_seurat_layer_data(object, assay = spliced_assay, layer = spliced_layer)
  unspliced <- .sn_get_seurat_layer_data(object, assay = unspliced_assay, layer = unspliced_layer)
  features <- intersect(rownames(spliced), rownames(unspliced))
  cells <- Reduce(intersect, list(colnames(spliced), colnames(unspliced), rownames(embedding$matrix)))
  if (length(features) < 3L || length(cells) < 5L) {
    stop("Velocity requires at least three shared features and five shared cells.", call. = FALSE)
  }
  Matrix::writeMM(.sn_as_sparse_matrix(spliced[features, cells, drop = FALSE]), file.path(input_dir, "spliced.mtx"))
  Matrix::writeMM(.sn_as_sparse_matrix(unspliced[features, cells, drop = FALSE]), file.path(input_dir, "unspliced.mtx"))
  utils::write.csv(data.frame(cell = cells, object[[]][cells, , drop = FALSE], check.names = FALSE), file.path(input_dir, "obs.csv"), row.names = FALSE)
  utils::write.csv(data.frame(feature = features), file.path(input_dir, "var.csv"), row.names = FALSE)
  coordinates <- data.frame(cell = cells, embedding$matrix[cells, , drop = FALSE], check.names = FALSE)
  utils::write.csv(coordinates, file.path(input_dir, "embedding.csv"), row.names = FALSE)
  list(features = features, cells = cells)
}

.sn_run_velocity_pixi <- function(object,
                                  spliced_assay,
                                  spliced_layer,
                                  unspliced_assay,
                                  unspliced_layer,
                                  embedding,
                                  backend_control) {
  root <- backend_control$run_dir %||% tempfile("shennong-velocity-")
  input_dir <- file.path(root, "input")
  output_dir <- file.path(root, "output")
  exported <- .sn_write_velocity_input(
    object, spliced_assay, spliced_layer, unspliced_assay, unspliced_layer,
    embedding, input_dir
  )
  config <- list(
    mode = "velocity", velocity_mode = backend_control$velocity_mode %||% "stochastic",
    min_shared_counts = backend_control$min_shared_counts %||% 10L,
    n_top_genes = backend_control$n_top_genes %||% min(2000L, length(exported$features)),
    n_neighbors = backend_control$n_neighbors %||% min(30L, length(exported$cells) - 1L),
    n_pcs = backend_control$n_pcs %||% min(30L, ncol(embedding$matrix)),
    max_graph_edges = backend_control$max_graph_edges %||% 100000L,
    write_h5ad = backend_control$write_h5ad %||% TRUE,
    random_seed = backend_control$seed %||% 717L
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  config_path <- file.path(root, "config.json")
  jsonlite::write_json(config, config_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  call_control <- backend_control$pixi %||% list()
  call <- c(list(
    environment = "trajectory", command = "python",
    args = c(.sn_trajectory_pixi_script(), "--input-dir", input_dir, "--output-dir", output_dir, "--config", config_path)
  ), call_control)
  do.call(sn_call_pixi_environment, call)
  cells <- utils::read.csv(file.path(output_dir, "velocity_cells.csv"), check.names = FALSE)
  graph <- if (file.exists(file.path(output_dir, "velocity_graph.csv"))) utils::read.csv(file.path(output_dir, "velocity_graph.csv"), check.names = FALSE) else data.frame()
  manifest <- jsonlite::read_json(file.path(output_dir, "manifest.json"), simplifyVector = TRUE)
  list(cells = cells, graph = graph, artifacts = manifest, expression = list(assay = spliced_assay, layer = spliced_layer))
}

.sn_standardize_velocity <- function(output, object, embedding, method) {
  if (!is.list(output)) stop("Velocity backend output must be a list.", call. = FALSE)
  cells <- tibble::as_tibble(output$cells %||% output$velocity)
  cell_column <- intersect(c("cell", "cell_id", "entity"), names(cells))[[1]] %||% NULL
  x_column <- intersect(c("velocity_1", "velocity_x", "vx", "V1"), names(cells))[[1]] %||% NULL
  y_column <- intersect(c("velocity_2", "velocity_y", "vy", "V2"), names(cells))[[1]] %||% NULL
  if (is_null(cell_column) || is_null(x_column) || is_null(y_column)) {
    stop("Velocity cells require cell, velocity_1/x, and velocity_2/y columns.", call. = FALSE)
  }
  cell_ids <- as.character(cells[[cell_column]])
  keep <- cell_ids %in% rownames(embedding$matrix)
  cell_ids <- cell_ids[keep]
  base <- tibble::tibble(
    cell = cell_ids,
    dimension_1 = as.numeric(embedding$matrix[cell_ids, 1]),
    dimension_2 = as.numeric(embedding$matrix[cell_ids, 2]),
    velocity_1 = suppressWarnings(as.numeric(cells[[x_column]][keep])),
    velocity_2 = suppressWarnings(as.numeric(cells[[y_column]][keep]))
  )
  optional <- list(
    pseudotime = c("pseudotime", "velocity_pseudotime", "latent_time"),
    confidence = c("confidence", "velocity_confidence"),
    velocity_length = c("velocity_length", "length")
  )
  for (name in names(optional)) {
    hits <- intersect(optional[[name]], names(cells))
    column <- if (length(hits) == 0L) NULL else hits[[1]]
    base[[name]] <- if (is_null(column)) NA_real_ else suppressWarnings(as.numeric(cells[[column]][keep]))
  }
  base$method <- method
  graph <- tibble::as_tibble(output$graph %||% output$edges %||% tibble::tibble())
  if (nrow(graph) > 0L) {
    source <- intersect(c("source", "from", "cell"), names(graph))[[1]] %||% NULL
    target <- intersect(c("target", "to", "neighbor"), names(graph))[[1]] %||% NULL
    weight <- intersect(c("weight", "probability", "score"), names(graph))[[1]] %||% NULL
    if (is_null(source) || is_null(target) || is_null(weight)) stop("Velocity graph requires source, target, and weight columns.", call. = FALSE)
    graph <- tibble::tibble(
      source = as.character(graph[[source]]), target = as.character(graph[[target]]),
      weight = suppressWarnings(as.numeric(graph[[weight]]))
    )
  }
  list(cells = base, graph = graph, artifacts = output$artifacts %||% list(), warnings = output$warnings %||% character())
}

#' Run RNA velocity with the managed scVelo backend
#'
#' @param object A Seurat object containing spliced and unspliced layers.
#' @param method Velocity backend; currently scVelo.
#' @param spliced_assay,unspliced_assay Assays containing count layers.
#' @param spliced_layer,unspliced_layer Layer names.
#' @param reduction,dims Embedding and dimensions used for projected vectors.
#' @param store_name Stored result name.
#' @param backend_control scVelo/pixi controls or an explicit `runner`/`result`.
#' @param return_object Return the modified object or unified velocity result.
#' @return A Seurat object or velocity result.
#' @examples
#' \dontrun{
#' object <- sn_run_velocity(object, spliced_layer = "spliced", unspliced_layer = "unspliced")
#' velocity <- sn_get_result(object, "velocity", "velocity")
#' }
#' @export
sn_run_velocity <- function(object,
                            method = "scvelo",
                            spliced_assay = NULL,
                            spliced_layer = "spliced",
                            unspliced_assay = NULL,
                            unspliced_layer = "unspliced",
                            reduction = NULL,
                            dims = 1:2,
                            store_name = "velocity",
                            backend_control = list(),
                            return_object = TRUE) {
  .sn_validate_result_object(object)
  method <- match.arg(method, "scvelo")
  spliced_assay <- spliced_assay %||% SeuratObject::DefaultAssay(object)
  unspliced_assay <- unspliced_assay %||% spliced_assay
  embedding <- .sn_velocity_embedding(object, reduction, dims)
  output <- if (is.function(backend_control$runner)) {
    backend_control$runner(
      object = object, method = method, spliced_assay = spliced_assay,
      spliced_layer = spliced_layer, unspliced_assay = unspliced_assay,
      unspliced_layer = unspliced_layer, reduction = embedding$reduction,
      dims = embedding$dims, backend_control = backend_control
    )
  } else if (!is_null(backend_control$result)) {
    backend_control$result
  } else {
    .sn_run_velocity_pixi(
      object, spliced_assay, spliced_layer, unspliced_assay,
      unspliced_layer, embedding, backend_control
    )
  }
  standardized <- .sn_standardize_velocity(output, object, embedding, method)
  cells <- standardized$cells
  object[[paste0(store_name, "_pseudotime")]] <- stats::setNames(cells$pseudotime, cells$cell)[colnames(object)]
  object[[paste0(store_name, "_confidence")]] <- stats::setNames(cells$confidence, cells$cell)[colnames(object)]
  result <- list(
    schema_version = "1.0", analysis_type = "velocity", name = store_name,
    method = method, backend = "scvelo-pixi",
    input = list(
      cells = nrow(cells), spliced_assay = spliced_assay, spliced_layer = spliced_layer,
      unspliced_assay = unspliced_assay, unspliced_layer = unspliced_layer,
      reduction = embedding$reduction, dimensions = embedding$dims
    ),
    parameters = list(velocity_mode = backend_control$velocity_mode %||% "stochastic"),
    tables = list(primary = cells, cells = cells, transition_edges = standardized$graph),
    embeddings = list(reduction = embedding$matrix, velocity = as.matrix(cells[, c("velocity_1", "velocity_2")])),
    graphs = list(transition_edges = standardized$graph),
    models = list(artifacts = standardized$artifacts),
    diagnostics = list(
      finite_vectors = sum(is.finite(cells$velocity_1) & is.finite(cells$velocity_2)),
      median_confidence = stats::median(cells$confidence, na.rm = TRUE),
      transition_edges = nrow(standardized$graph)
    ),
    warnings = as.character(standardized$warnings),
    provenance = .sn_analysis_provenance(random_seed = backend_control$seed %||% 717L)
  )
  sn_validate_result(result)
  object <- sn_store_result(object, "velocity", store_name, result)
  object <- .sn_log_seurat_command(object, assay = spliced_assay, name = "sn_run_velocity")
  if (isTRUE(return_object)) object else sn_get_result(object, "velocity", store_name)
}

.sn_run_fate_pixi <- function(velocity_result, backend_control) {
  h5ad <- backend_control$h5ad %||% velocity_result$models$artifacts$output_h5ad %||% NULL
  if (is_null(h5ad) || !file.exists(h5ad)) {
    stop("CellRank requires the scVelo backend H5AD artifact or `backend_control$h5ad`.", call. = FALSE)
  }
  root <- backend_control$run_dir %||% tempfile("shennong-fate-")
  output_dir <- file.path(root, "output")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  config <- list(
    mode = "fate", h5ad = normalizePath(h5ad, winslash = "/", mustWork = TRUE),
    n_states = backend_control$n_states %||% NULL,
    terminal_states = backend_control$terminal_states %||% NULL,
    terminal_method = backend_control$terminal_method %||% "stability",
    terminal_n_states = backend_control$terminal_n_states %||% NULL,
    stability_threshold = backend_control$stability_threshold %||% 0.96,
    compute_drivers = backend_control$compute_drivers %||% TRUE,
    n_jobs = backend_control$n_jobs %||% 1L,
    random_seed = backend_control$seed %||% 717L
  )
  config_path <- file.path(root, "config.json")
  jsonlite::write_json(config, config_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  do.call(sn_call_pixi_environment, c(list(
    environment = "trajectory", command = "python",
    args = c(.sn_trajectory_pixi_script(), "--output-dir", output_dir, "--config", config_path)
  ), backend_control$pixi %||% list()))
  list(
    probabilities = utils::read.csv(file.path(output_dir, "fate_probabilities.csv"), check.names = FALSE),
    terminal_states = utils::read.csv(file.path(output_dir, "terminal_states.csv"), check.names = FALSE),
    drivers = if (file.exists(file.path(output_dir, "lineage_drivers.csv"))) utils::read.csv(file.path(output_dir, "lineage_drivers.csv"), check.names = FALSE) else data.frame(),
    artifacts = jsonlite::read_json(file.path(output_dir, "manifest.json"), simplifyVector = TRUE)
  )
}

.sn_standardize_fate <- function(output, object) {
  if (!is.list(output)) stop("Fate backend output must be a list.", call. = FALSE)
  probabilities <- tibble::as_tibble(output$probabilities %||% output$fate_probabilities)
  cell <- intersect(c("cell", "cell_id", "entity"), names(probabilities))[[1]] %||% NULL
  state <- intersect(c("state", "lineage", "terminal_state"), names(probabilities))[[1]] %||% NULL
  probability <- intersect(c("probability", "fate_probability", "score"), names(probabilities))[[1]] %||% NULL
  if (is_null(cell) || is_null(state) || is_null(probability)) stop("Fate probabilities require cell, state, and probability columns.", call. = FALSE)
  probabilities <- tibble::tibble(
    cell = as.character(probabilities[[cell]]), state = as.character(probabilities[[state]]),
    probability = suppressWarnings(as.numeric(probabilities[[probability]]))
  )
  probabilities <- probabilities[probabilities$cell %in% colnames(object) & is.finite(probabilities$probability), , drop = FALSE]
  if (nrow(probabilities) == 0L) stop("No finite fate probabilities matched cells in `object`.", call. = FALSE)
  terminals <- tibble::as_tibble(output$terminal_states %||% tibble::tibble())
  drivers <- tibble::as_tibble(output$drivers %||% output$lineage_drivers %||% tibble::tibble())
  list(probabilities = probabilities, terminal_states = terminals, drivers = drivers, artifacts = output$artifacts %||% list(), warnings = output$warnings %||% character())
}

#' Infer terminal states and fate probabilities with CellRank
#'
#' @param object A Seurat object.
#' @param method Fate backend; currently CellRank.
#' @param velocity_name Stored velocity result used by the default pixi backend.
#' @param reduction,dims Embedding and dimensions used for plots.
#' @param store_name Stored fate result name.
#' @param backend_control CellRank/pixi controls or an explicit `runner`/`result`.
#' @param return_object Return the modified object or unified fate result.
#' @return A Seurat object or fate result.
#' @examples
#' \dontrun{
#' object <- sn_run_fate(object, velocity_name = "velocity")
#' fate <- sn_get_result(object, "fate", "fate")
#' }
#' @export
sn_run_fate <- function(object,
                        method = "cellrank",
                        velocity_name = "velocity",
                        reduction = NULL,
                        dims = 1:2,
                        store_name = "fate",
                        backend_control = list(),
                        return_object = TRUE) {
  .sn_validate_result_object(object)
  method <- match.arg(method, "cellrank")
  embedding <- .sn_velocity_embedding(object, reduction, dims)
  velocity <- tryCatch(sn_get_result(object, "velocity", velocity_name), error = function(e) NULL)
  output <- if (is.function(backend_control$runner)) {
    backend_control$runner(
      object = object, method = method, velocity_result = velocity,
      reduction = embedding$reduction, dims = embedding$dims,
      backend_control = backend_control
    )
  } else if (!is_null(backend_control$result)) {
    backend_control$result
  } else {
    if (is_null(velocity)) stop("Run `sn_run_velocity()` first or supply a CellRank runner/result.", call. = FALSE)
    .sn_run_fate_pixi(velocity, backend_control)
  }
  standardized <- .sn_standardize_fate(output, object)
  probabilities <- standardized$probabilities
  prefix <- gsub("[^[:alnum:]_]+", "_", store_name)
  metadata <- data.frame(row.names = colnames(object))
  for (state in unique(probabilities$state)) {
    values <- stats::setNames(probabilities$probability[probabilities$state == state], probabilities$cell[probabilities$state == state])
    metadata[[paste(prefix, "fate", gsub("[^[:alnum:]_]+", "_", state), sep = "_")]] <- as.numeric(values[colnames(object)])
  }
  object <- SeuratObject::AddMetaData(object, metadata = metadata)
  result <- list(
    schema_version = "1.0", analysis_type = "fate", name = store_name,
    method = method, backend = "cellrank-pixi",
    input = list(cells = length(unique(probabilities$cell)), velocity_name = velocity_name, reduction = embedding$reduction, dimensions = embedding$dims),
    parameters = list(
      n_states = backend_control$n_states %||% NULL,
      terminal_states = backend_control$terminal_states %||% NULL,
      terminal_method = backend_control$terminal_method %||% "stability",
      stability_threshold = backend_control$stability_threshold %||% 0.96
    ),
    tables = list(primary = probabilities, probabilities = probabilities, terminal_states = standardized$terminal_states, lineage_drivers = standardized$drivers),
    embeddings = list(reduction = embedding$matrix), graphs = list(),
    models = list(artifacts = standardized$artifacts),
    diagnostics = list(states = length(unique(probabilities$state)), probability_sum_range = range(tapply(probabilities$probability, probabilities$cell, sum), na.rm = TRUE)),
    warnings = as.character(standardized$warnings),
    provenance = .sn_analysis_provenance(random_seed = backend_control$seed %||% 717L)
  )
  sn_validate_result(result)
  object <- sn_store_result(object, "fate", store_name, result)
  object <- .sn_log_seurat_command(object, assay = SeuratObject::DefaultAssay(object), name = "sn_run_fate")
  if (isTRUE(return_object)) object else sn_get_result(object, "fate", store_name)
}
