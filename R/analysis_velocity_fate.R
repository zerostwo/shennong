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

.sn_write_regvelo_prior_grn <- function(prior_grn, path) {
  if (is_null(prior_grn)) {
    stop(
      "RegVelo requires `backend_control$prior_grn` as an edge table, named matrix, or CSV path.",
      call. = FALSE
    )
  }

  if (is.character(prior_grn) && length(prior_grn) == 1L) {
    if (!file.exists(prior_grn)) {
      stop("The RegVelo prior-GRN file does not exist: ", prior_grn, call. = FALSE)
    }
    if (!file.copy(prior_grn, path, overwrite = TRUE)) {
      stop("Could not stage the RegVelo prior-GRN file.", call. = FALSE)
    }
    return(invisible(path))
  }

  if (is.matrix(prior_grn) || inherits(prior_grn, "Matrix")) {
    if (is_null(rownames(prior_grn)) || is_null(colnames(prior_grn))) {
      stop("A RegVelo prior-GRN matrix requires target row names and regulator column names.", call. = FALSE)
    }
    sparse <- .sn_as_sparse_matrix(prior_grn)
    entries <- Matrix::summary(sparse)
    prior_grn <- data.frame(
      regulator = colnames(sparse)[entries$j],
      target = rownames(sparse)[entries$i],
      weight = entries$x,
      stringsAsFactors = FALSE
    )
  } else if (is.data.frame(prior_grn)) {
    regulator_hits <- intersect(c("regulator", "tf", "source", "from"), tolower(names(prior_grn)))
    target_hits <- intersect(c("target", "gene", "to"), tolower(names(prior_grn)))
    regulator_col <- if (length(regulator_hits) == 0L) NULL else regulator_hits[[1]]
    target_col <- if (length(target_hits) == 0L) NULL else target_hits[[1]]
    if (is_null(regulator_col) || is_null(target_col)) {
      stop("A RegVelo prior-GRN table requires `regulator` and `target` columns.", call. = FALSE)
    }
    original_names <- names(prior_grn)
    regulator_col <- original_names[match(regulator_col, tolower(original_names))]
    target_col <- original_names[match(target_col, tolower(original_names))]
    weight_hits <- intersect(c("weight", "score", "importance"), tolower(original_names))
    weight_name <- if (length(weight_hits) == 0L) NULL else weight_hits[[1]]
    weight_col <- if (is_null(weight_name)) NULL else original_names[match(weight_name, tolower(original_names))]
    prior_grn <- data.frame(
      regulator = as.character(prior_grn[[regulator_col]]),
      target = as.character(prior_grn[[target_col]]),
      weight = if (is_null(weight_col)) 1 else suppressWarnings(as.numeric(prior_grn[[weight_col]])),
      stringsAsFactors = FALSE
    )
  } else {
    stop("Unsupported RegVelo prior-GRN input.", call. = FALSE)
  }

  keep <- nzchar(prior_grn$regulator) & nzchar(prior_grn$target) &
    is.finite(prior_grn$weight) & prior_grn$weight != 0
  prior_grn <- unique(prior_grn[keep, , drop = FALSE])
  if (nrow(prior_grn) == 0L) {
    stop("The RegVelo prior GRN contains no finite non-zero edges.", call. = FALSE)
  }
  utils::write.csv(prior_grn, path, row.names = FALSE)
  invisible(path)
}

.sn_velocity_inference_n_pcs <- function(exported, requested = NULL) {
  maximum <- min(length(exported$features) - 1L, length(exported$cells) - 1L)
  if (!is.finite(maximum) || maximum < 2L) {
    stop(
      "Velocity inference requires at least three exported features and three exported cells to use two principal components.",
      call. = FALSE
    )
  }
  requested <- requested %||% min(30L, maximum)
  if (!is.numeric(requested) || length(requested) != 1L || is.na(requested) ||
      !is.finite(requested) || requested != as.integer(requested) || requested < 2L) {
    stop("`backend_control$n_pcs` must be one integer of at least 2.", call. = FALSE)
  }
  resolved <- as.integer(min(requested, maximum))
  if (resolved < 2L) {
    stop("The resolved velocity `n_pcs` must be at least 2 after input-size clamping.", call. = FALSE)
  }
  resolved
}

.sn_run_velocity_pixi <- function(object,
                                  method,
                                  spliced_assay,
                                  spliced_layer,
                                  unspliced_assay,
                                  unspliced_layer,
                                  embedding,
                                  backend_control) {
  run_dir_supplied <- !is_null(backend_control$run_dir)
  keep_run_dir <- backend_control$keep_run_dir %||% TRUE
  if (!is.logical(keep_run_dir) || length(keep_run_dir) != 1L || is.na(keep_run_dir)) {
    stop("`backend_control$keep_run_dir` must be TRUE or FALSE.", call. = FALSE)
  }
  runtime_dir <- .sn_shennong_runtime_dir(backend_control$runtime_dir %||% NULL)
  root <- if (isTRUE(keep_run_dir) && !run_dir_supplied) {
    .sn_create_owned_run_dir(
      parent = file.path(runtime_dir, "runs"),
      prefix = paste0("velocity-", method, "_")
    )
  } else {
    .sn_resolve_python_run_directory(
      path = backend_control$run_dir,
      method = paste0("velocity-", method),
      runtime_dir = runtime_dir,
      keep_run_dir = keep_run_dir,
      supplied = run_dir_supplied
    )
  }
  run_complete <- FALSE
  on.exit({
    if (!run_complete && .sn_is_owned_run_dir(root) && dir.exists(root)) {
      .sn_sanitize_failed_python_run(root, method = method, stage = "velocity")
    }
  }, add = TRUE)
  input_dir <- file.path(root, "input")
  output_dir <- file.path(root, "output")
  exported <- .sn_write_velocity_input(
    object, spliced_assay, spliced_layer, unspliced_assay, unspliced_layer,
    embedding, input_dir
  )
  prior_grn_path <- NULL
  if (identical(method, "regvelo")) {
    prior_grn_path <- file.path(input_dir, "prior_grn.csv")
    .sn_write_regvelo_prior_grn(backend_control$prior_grn, prior_grn_path)
  }
  config <- list(
    mode = method, velocity_mode = backend_control$velocity_mode %||% "stochastic",
    min_shared_counts = backend_control$min_shared_counts %||% 10L,
    enforce_normalization = backend_control$enforce_normalization %||% TRUE,
    log1p_transform = backend_control$log1p_transform %||% TRUE,
    n_top_genes = backend_control$n_top_genes %||% min(2000L, length(exported$features)),
    n_neighbors = backend_control$n_neighbors %||% min(30L, length(exported$cells) - 1L),
    n_pcs = .sn_velocity_inference_n_pcs(exported, backend_control$n_pcs),
    max_graph_edges = backend_control$max_graph_edges %||% 100000L,
    write_h5ad = backend_control$write_h5ad %||% TRUE,
    prior_grn = prior_grn_path,
    soft_constraint = backend_control$soft_constraint %||% TRUE,
    lam = backend_control$lam %||% 1,
    lam2 = backend_control$lam2 %||% 0,
    max_epochs = backend_control$max_epochs %||% 1500L,
    learning_rate = backend_control$learning_rate %||% 0.01,
    train_size = backend_control$train_size %||% 0.9,
    batch_size = backend_control$batch_size %||% NULL,
    early_stopping = backend_control$early_stopping %||% TRUE,
    min_max_scale = backend_control$min_max_scale %||% TRUE,
    filter_on_r2 = backend_control$filter_on_r2 %||% TRUE,
    posterior_samples = backend_control$posterior_samples %||% 30L,
    save_model = backend_control$save_model %||% TRUE,
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
  manifest <- .sn_read_integration_backend_manifest(
    output_dir,
    method = method,
    n_cells = length(exported$cells)
  )
  manifest$output_h5ad <- .sn_integration_manifest_artifact(
    manifest,
    field = "output_h5ad",
    output_dir = output_dir
  )
  if (isTRUE(keep_run_dir)) {
    manifest$run_dir <- normalizePath(root, winslash = "/", mustWork = TRUE)
    manifest$output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
    manifest$run_dir_retained <- TRUE
    manifest$retention_reason <- "Retained because CellRank fate inference consumes the velocity H5AD artifact."
    unlink(input_dir, recursive = TRUE, force = TRUE)
    unlink(config_path, force = TRUE)
  } else {
    manifest <- .sn_sanitize_python_backend_manifest(manifest, retain_paths = FALSE)
    .sn_remove_python_run_directory(root, label = paste0(method, " velocity run directory"))
  }
  run_complete <- TRUE
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
  if (length(cell_ids) == 0L || anyNA(cell_ids) || any(!nzchar(cell_ids))) {
    stop("Velocity cell identifiers cannot be missing or empty.", call. = FALSE)
  }
  if (anyDuplicated(cell_ids)) {
    stop("Velocity output must contain one unique row per cell.", call. = FALSE)
  }
  unknown_cells <- setdiff(cell_ids, rownames(embedding$matrix))
  if (length(unknown_cells) > 0L) {
    stop(
      "Velocity output contains cell(s) absent from the selected embedding: ",
      paste(utils::head(unknown_cells, 5L), collapse = ", "), ".",
      call. = FALSE
    )
  }
  keep <- rep(TRUE, length(cell_ids))
  base <- tibble::tibble(
    cell = cell_ids,
    dimension_1 = as.numeric(embedding$matrix[cell_ids, 1]),
    dimension_2 = as.numeric(embedding$matrix[cell_ids, 2]),
    velocity_1 = suppressWarnings(as.numeric(cells[[x_column]][keep])),
    velocity_2 = suppressWarnings(as.numeric(cells[[y_column]][keep]))
  )
  finite_vectors <- is.finite(base$velocity_1) & is.finite(base$velocity_2)
  partial_vectors <- xor(is.finite(base$velocity_1), is.finite(base$velocity_2))
  if (any(partial_vectors)) {
    stop("Velocity vectors must provide both finite components or neither component.", call. = FALSE)
  }
  if (!any(finite_vectors)) {
    stop("Velocity output contains no finite two-dimensional vectors.", call. = FALSE)
  }
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
    if (anyNA(graph$source) || any(!nzchar(graph$source)) ||
        anyNA(graph$target) || any(!nzchar(graph$target))) {
      stop("Velocity graph endpoints cannot be missing or empty.", call. = FALSE)
    }
    unknown_endpoints <- setdiff(unique(c(graph$source, graph$target)), cell_ids)
    if (length(unknown_endpoints) > 0L) {
      stop(
        "Velocity graph contains endpoint(s) absent from velocity cells: ",
        paste(utils::head(unknown_endpoints, 5L), collapse = ", "), ".",
        call. = FALSE
      )
    }
    if (any(!is.finite(graph$weight)) || any(graph$weight < 0)) {
      stop("Velocity transition weights must be finite and non-negative.", call. = FALSE)
    }
    if (anyDuplicated(paste(graph$source, graph$target, sep = "\r"))) {
      stop("Velocity graph must contain one unique transition per source-target pair.", call. = FALSE)
    }
  }
  list(
    cells = base, graph = graph, artifacts = output$artifacts %||% list(),
    warnings = output$warnings %||% character(), backend = output$backend %||% NULL
  )
}

#' Run RNA velocity with managed scVelo or RegVelo backends
#'
#' @param object A Seurat object containing spliced and unspliced layers.
#' @param method Velocity backend: \code{"scvelo"} or \code{"regvelo"}.
#' @param spliced_assay,unspliced_assay Assays containing count layers.
#' @param spliced_layer,unspliced_layer Layer names.
#' @param reduction,dims Embedding and dimensions used for projected vectors.
#' @param result_id Stored result name.
#' @param backend_control Backend/pixi controls or an explicit `runner`/`result`.
#'   RegVelo requires \code{prior_grn}, supplied as a regulator-target edge
#'   table, a target-by-regulator named matrix, or a CSV path. Shared scVelo
#'   preprocessing defaults to \code{enforce_normalization = TRUE} so
#'   non-integer source splicing estimates are normalized before HVG selection;
#'   \code{log1p_transform = TRUE} prepares the expression matrix for Scanpy's
#'   Seurat-flavor HVG calculation. Managed runs retain their output directory
#'   by default because \code{sn_run_fate()} consumes the generated H5AD; raw
#'   export files are removed after successful import. Set
#'   \code{keep_run_dir = FALSE} when CellRank chaining is not needed, or supply
#'   an empty \code{run_dir} to choose the retained location explicitly.
#' @param return_object Return the modified object or unified velocity result.
#' @param seed Top-level reproducibility seed. Precedence: \code{seed} >
#'   \code{backend_control$seed} > task default; the resolved value is stamped
#'   into result provenance.
#' @param verbose Top-level progress switch forwarded through
#'   \code{backend_control$verbose} when explicitly supplied.
#' @return A Seurat object or velocity result.
#' @references RegVelo documentation: \url{https://regvelo.readthedocs.io/}.
#'   Wang et al. (2026), Cell, \doi{10.1016/j.cell.2026.04.022}.
#' @examples
#' \dontrun{
#' object <- sn_run_velocity(object, spliced_layer = "spliced", unspliced_layer = "unspliced")
#' velocity <- sn_get_result(object, "velocity", "velocity")
#' }
#' @export
sn_run_velocity <- function(object,
                            method = c("scvelo", "regvelo"),
                            spliced_assay = NULL,
                            spliced_layer = "spliced",
                            unspliced_assay = NULL,
                            unspliced_layer = "unspliced",
                            reduction = NULL,
                            dims = 1:2,
                            result_id = "velocity",
                            backend_control = list(),
                            return_object = TRUE,
                            seed = NULL,
                            verbose = TRUE) {
  result_id <- .sn_validate_result_id(result_id)
  .sn_validate_seurat_object(object)
  method <- match.arg(method)
  backend_control$seed <- seed %||% backend_control$seed
  if (!missing(verbose)) {
    backend_control$verbose <- isTRUE(verbose)
  }
  spliced_assay <- spliced_assay %||% SeuratObject::DefaultAssay(object)
  unspliced_assay <- unspliced_assay %||% spliced_assay
  embedding <- .sn_velocity_embedding(object, reduction, dims)
  custom_runner <- is.function(backend_control$runner)
  supplied_result <- !is_null(backend_control$result)
  backend_label <- if (custom_runner) {
    paste0(method, "-custom-runner")
  } else if (supplied_result) {
    paste0(method, "-provided-result")
  } else {
    paste0(method, "-pixi")
  }
  output <- if (custom_runner) {
    backend_control$runner(
      object = object, method = method, spliced_assay = spliced_assay,
      spliced_layer = spliced_layer, unspliced_assay = unspliced_assay,
      unspliced_layer = unspliced_layer, reduction = embedding$reduction,
      dims = embedding$dims, backend_control = backend_control
    )
  } else if (supplied_result) {
    backend_control$result
  } else {
    .sn_run_velocity_pixi(
      object, method, spliced_assay, spliced_layer, unspliced_assay,
      unspliced_layer, embedding, backend_control
    )
  }
  standardized <- .sn_standardize_velocity(output, object, embedding, method)
  cells <- standardized$cells
  object[[paste0(result_id, "_pseudotime")]] <- stats::setNames(cells$pseudotime, cells$cell)[colnames(object)]
  object[[paste0(result_id, "_confidence")]] <- stats::setNames(cells$confidence, cells$cell)[colnames(object)]
  result <- list(
    schema_version = .sn_analysis_result_schema_version(), analysis_type = "velocity", result_id = result_id,
    method = method, backend = standardized$backend %||% backend_label,
    input = list(
      cells = nrow(cells), spliced_assay = spliced_assay, spliced_layer = spliced_layer,
      unspliced_assay = unspliced_assay, unspliced_layer = unspliced_layer,
      reduction = embedding$reduction, dimensions = embedding$dims
    ),
    parameters = if (identical(method, "regvelo")) {
      list(
        soft_constraint = backend_control$soft_constraint %||% TRUE,
        lam = backend_control$lam %||% 1,
        lam2 = backend_control$lam2 %||% 0,
        max_epochs = backend_control$max_epochs %||% 1500L,
        filter_on_r2 = backend_control$filter_on_r2 %||% TRUE,
        enforce_normalization = backend_control$enforce_normalization %||% TRUE,
        log1p_transform = backend_control$log1p_transform %||% TRUE
      )
    } else {
      list(
        velocity_mode = backend_control$velocity_mode %||% "stochastic",
        enforce_normalization = backend_control$enforce_normalization %||% TRUE,
        log1p_transform = backend_control$log1p_transform %||% TRUE
      )
    },
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
  object <- sn_store_result(object, "velocity", result_id, result)
  object <- .sn_log_seurat_command(object, assay = spliced_assay, name = "sn_run_velocity")
  if (isTRUE(return_object)) object else sn_get_result(object, "velocity", result_id)
}

.sn_run_fate_pixi <- function(velocity_result, backend_control) {
  h5ad <- backend_control$h5ad %||% velocity_result$models$artifacts$output_h5ad %||% NULL
  if (is_null(h5ad) || !file.exists(h5ad)) {
    stop("CellRank requires the scVelo backend H5AD artifact or `backend_control$h5ad`.", call. = FALSE)
  }
  run_dir_supplied <- !is_null(backend_control$run_dir)
  keep_run_dir <- backend_control$keep_run_dir %||% run_dir_supplied
  if (!is.logical(keep_run_dir) || length(keep_run_dir) != 1L || is.na(keep_run_dir)) {
    stop("`backend_control$keep_run_dir` must be TRUE or FALSE.", call. = FALSE)
  }
  run_parent <- backend_control$run_dir %||% tempdir()
  if (!dir.exists(run_parent)) {
    created <- dir.create(run_parent, recursive = TRUE, showWarnings = FALSE)
    if (!isTRUE(created) && !dir.exists(run_parent)) {
      stop("Could not create the CellRank run directory: ", run_parent, call. = FALSE)
    }
  }
  if (!isTRUE(keep_run_dir) || !run_dir_supplied) {
    # An explicit `run_dir` is a user-owned parent when cleanup is requested.
    # A retained implicit run also needs its own directory rather than using
    # the process-wide temporary directory itself.
    if (!isTRUE(keep_run_dir)) {
      root <- .sn_create_owned_run_dir(run_parent, "shennong-fate-run-")
    } else {
      root <- tempfile("shennong-fate-run-", tmpdir = run_parent)
      if (!dir.create(root, recursive = FALSE, showWarnings = FALSE)) {
        stop("Could not create a package-owned CellRank run directory.", call. = FALSE)
      }
    }
  } else {
    root <- run_parent
  }
  if (!isTRUE(keep_run_dir)) {
    on.exit(.sn_cleanup_owned_run_dir(root), add = TRUE)
  }
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
  artifacts <- jsonlite::read_json(file.path(output_dir, "manifest.json"), simplifyVector = TRUE)
  artifacts$run_dir_retained <- isTRUE(keep_run_dir)
  if (isTRUE(keep_run_dir)) {
    artifacts$output_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  }
  list(
    probabilities = utils::read.csv(file.path(output_dir, "fate_probabilities.csv"), check.names = FALSE),
    terminal_states = utils::read.csv(file.path(output_dir, "terminal_states.csv"), check.names = FALSE),
    drivers = if (file.exists(file.path(output_dir, "lineage_drivers.csv"))) utils::read.csv(file.path(output_dir, "lineage_drivers.csv"), check.names = FALSE) else data.frame(),
    artifacts = artifacts,
    backend = "cellrank-pixi"
  )
}

.sn_fate_state_metadata <- function(object, result_id, probabilities) {
  states <- unique(probabilities$state)
  prefix <- gsub("[^[:alnum:]_]+", "_", result_id)
  safe_states <- gsub("[^[:alnum:]_]+", "_", states)
  safe_states[!nzchar(safe_states)] <- "state"
  candidates <- make.unique(paste(prefix, "fate", safe_states, sep = "_"), sep = "_")

  stored_results <- .sn_result_store(object)[["fate"]] %||% list()
  previous <- stored_results[[result_id]]
  if (is_null(previous)) {
    conflicts <- intersect(candidates, colnames(object[[]]))
    if (length(conflicts) > 0L) {
      stop(
        "Fate metadata column(s) already exist but are not owned by this result: ",
        paste(conflicts, collapse = ", "),
        ". Choose a different `result_id` or remove/rename the conflicting user metadata.",
        call. = FALSE
      )
    }
    return(tibble::tibble(state = states, metadata_column = candidates))
  }

  previous <- sn_get_result(object, "fate", result_id)
  mapping <- previous$tables$state_metadata
  if (!is.data.frame(mapping) ||
      !all(c("state", "metadata_column") %in% colnames(mapping))) {
    stop(
      "The existing fate result does not contain a valid `state_metadata` ownership map; use a new `result_id`.",
      call. = FALSE
    )
  }
  mapping <- tibble::as_tibble(mapping[, c("state", "metadata_column"), drop = FALSE])
  mapping$state <- as.character(mapping$state)
  mapping$metadata_column <- as.character(mapping$metadata_column)
  if (anyNA(mapping$state) || any(!nzchar(mapping$state)) ||
      anyNA(mapping$metadata_column) || any(!nzchar(mapping$metadata_column)) ||
      anyDuplicated(mapping$state) || anyDuplicated(mapping$metadata_column)) {
    stop(
      "The existing fate result has an ambiguous `state_metadata` ownership map; use a new `result_id`.",
      call. = FALSE
    )
  }
  if (!setequal(mapping$state, states)) {
    stop(
      "Re-running a fate `result_id` with a different state set is ambiguous; use a new `result_id`.",
      call. = FALSE
    )
  }
  mapping <- mapping[match(states, mapping$state), , drop = FALSE]

  previous_probabilities <- previous$tables$probabilities %||% previous$tables$primary
  if (!is.data.frame(previous_probabilities) ||
      !all(c("cell", "state", "probability") %in% colnames(previous_probabilities))) {
    stop(
      "The existing fate result cannot prove ownership of its metadata values; use a new `result_id`.",
      call. = FALSE
    )
  }
  for (index in seq_len(nrow(mapping))) {
    column <- mapping$metadata_column[[index]]
    if (!column %in% colnames(object[[]])) next
    state <- mapping$state[[index]]
    prior <- previous_probabilities[as.character(previous_probabilities$state) == state, , drop = FALSE]
    expected <- stats::setNames(as.numeric(prior$probability), as.character(prior$cell))
    expected <- as.numeric(expected[colnames(object)])
    observed <- as.numeric(object[[column, drop = TRUE]])
    if (!isTRUE(all.equal(observed, expected, check.attributes = FALSE))) {
      stop(
        "Fate metadata column `", column,
        "` no longer matches the stored result and may be user-modified; refusing to overwrite it.",
        call. = FALSE
      )
    }
  }
  mapping
}

.sn_trajectory_cellrank <- function(velocity_result, backend_control = list()) {
  .sn_run_fate_pixi(velocity_result, backend_control)
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
  if (anyNA(probabilities$cell) || any(!nzchar(probabilities$cell)) ||
      anyNA(probabilities$state) || any(!nzchar(probabilities$state))) {
    stop("Fate probability cell and state identifiers cannot be missing or empty.", call. = FALSE)
  }
  unknown_cells <- setdiff(unique(probabilities$cell), colnames(object))
  if (length(unknown_cells) > 0L) {
    stop("Fate probabilities contain cell(s) absent from `object`: ", paste(utils::head(unknown_cells, 5L), collapse = ", "), ".", call. = FALSE)
  }
  if (anyDuplicated(paste(probabilities$cell, probabilities$state, sep = "\r"))) {
    stop("Fate probabilities must contain one unique row per cell and state.", call. = FALSE)
  }
  if (any(!is.finite(probabilities$probability)) ||
      any(probabilities$probability < 0 | probabilities$probability > 1)) {
    stop("Fate probabilities must be finite values between 0 and 1.", call. = FALSE)
  }
  if (nrow(probabilities) == 0L) stop("No fate probabilities were supplied.", call. = FALSE)
  probability_sums <- tapply(probabilities$probability, probabilities$cell, sum)
  if (any(abs(probability_sums - 1) > 1e-6)) {
    stop("Fate probabilities must sum to 1 across states for every cell.", call. = FALSE)
  }
  terminals <- tibble::as_tibble(output$terminal_states %||% tibble::tibble())
  if (nrow(terminals) > 0L) {
    terminal_cell <- intersect(c("cell", "cell_id", "entity"), names(terminals))
    terminal_state <- intersect(c("state", "lineage", "terminal_state"), names(terminals))
    if (length(terminal_state) == 0L) {
      stop("Fate terminal states require a state/lineage column.", call. = FALSE)
    }
    states <- as.character(terminals[[terminal_state[[1]]]])
    keep_terminal <- !is.na(states) & nzchar(states)
    terminals <- terminals[keep_terminal, , drop = FALSE]
    states <- states[keep_terminal]
    unknown_states <- setdiff(unique(states), unique(probabilities$state))
    if (length(unknown_states) > 0L) {
      stop(
        "Fate terminal states are absent from probability states: ",
        paste(utils::head(unknown_states, 5L), collapse = ", "), ".",
        call. = FALSE
      )
    }
    if (length(terminal_cell) > 0L) {
      terminal_cells <- as.character(terminals[[terminal_cell[[1]]]])
      unknown_terminal_cells <- setdiff(unique(terminal_cells), colnames(object))
      if (anyNA(terminal_cells) || any(!nzchar(terminal_cells)) || length(unknown_terminal_cells) > 0L) {
        stop("Fate terminal-state cells must be non-empty cells present in `object`.", call. = FALSE)
      }
    }
  }
  drivers <- tibble::as_tibble(output$drivers %||% output$lineage_drivers %||% tibble::tibble())
  list(
    probabilities = probabilities, terminal_states = terminals, drivers = drivers,
    artifacts = output$artifacts %||% list(), warnings = output$warnings %||% character(),
    backend = output$backend %||% NULL
  )
}

#' Infer terminal states and fate probabilities with CellRank
#'
#' @param object A Seurat object.
#' @param method Fate backend; currently CellRank.
#' @param source_result_id Stored velocity result used by the default pixi backend.
#' @param reduction,dims Embedding and dimensions used for plots.
#' @param result_id Stored fate result name.
#' @param backend_control CellRank/pixi controls or an explicit `runner`/`result`.
#' @param return_object Return the modified object or unified fate result.
#' @param seed Top-level reproducibility seed. Precedence: \code{seed} >
#'   \code{backend_control$seed} > task default.
#' @param verbose Top-level progress switch forwarded through
#'   \code{backend_control$verbose} when explicitly supplied.
#' @return A Seurat object or fate result.
#' @examples
#' \dontrun{
#' object <- sn_run_fate(object, source_result_id = "velocity")
#' fate <- sn_get_result(object, "fate", "fate")
#' }
#' @export
sn_run_fate <- function(object,
                        method = c("cellrank"),
                        source_result_id = "velocity",
                        reduction = NULL,
                        dims = 1:2,
                        result_id = "fate",
                        backend_control = list(),
                        return_object = TRUE,
                        seed = NULL,
                        verbose = TRUE) {
  result_id <- .sn_validate_result_id(result_id)
  .sn_validate_seurat_object(object)
  method <- match.arg(method)
  backend_control$seed <- seed %||% backend_control$seed
  if (!missing(verbose)) {
    backend_control$verbose <- isTRUE(verbose)
  }
  embedding <- .sn_velocity_embedding(object, reduction, dims)
  velocity <- tryCatch(sn_get_result(object, "velocity", source_result_id), error = function(e) NULL)
  custom_runner <- is.function(backend_control$runner)
  supplied_result <- !is_null(backend_control$result)
  backend_label <- if (custom_runner) {
    "cellrank-custom-runner"
  } else if (supplied_result) {
    "cellrank-provided-result"
  } else {
    "cellrank-pixi"
  }
  output <- if (custom_runner) {
    backend_control$runner(
      object = object, method = method, velocity_result = velocity,
      reduction = embedding$reduction, dims = embedding$dims,
      backend_control = backend_control
    )
  } else if (supplied_result) {
    backend_control$result
  } else {
    if (is_null(velocity)) stop("Run `sn_run_velocity()` first or supply a CellRank runner/result.", call. = FALSE)
    .sn_trajectory_cellrank(velocity, backend_control)
  }
  standardized <- .sn_standardize_fate(output, object)
  probabilities <- standardized$probabilities
  states <- unique(probabilities$state)
  state_metadata <- .sn_fate_state_metadata(object, result_id, probabilities)
  metadata <- data.frame(row.names = colnames(object))
  for (state in states) {
    values <- stats::setNames(probabilities$probability[probabilities$state == state], probabilities$cell[probabilities$state == state])
    column <- state_metadata$metadata_column[state_metadata$state == state][[1]]
    metadata[[column]] <- as.numeric(values[colnames(object)])
  }
  object <- SeuratObject::AddMetaData(object, metadata = metadata)
  result <- list(
    schema_version = .sn_analysis_result_schema_version(), analysis_type = "fate", result_id = result_id,
    method = method, backend = standardized$backend %||% backend_label,
    input = list(cells = length(unique(probabilities$cell)), source_result_id = source_result_id, reduction = embedding$reduction, dimensions = embedding$dims),
    parameters = list(
      n_states = backend_control$n_states %||% NULL,
      terminal_states = backend_control$terminal_states %||% NULL,
      terminal_method = backend_control$terminal_method %||% "stability",
      stability_threshold = backend_control$stability_threshold %||% 0.96
    ),
    tables = list(primary = probabilities, probabilities = probabilities, terminal_states = standardized$terminal_states, lineage_drivers = standardized$drivers, state_metadata = state_metadata),
    embeddings = list(reduction = embedding$matrix), graphs = list(),
    models = list(artifacts = standardized$artifacts),
    diagnostics = list(states = length(unique(probabilities$state)), probability_sum_range = range(tapply(probabilities$probability, probabilities$cell, sum), na.rm = TRUE)),
    warnings = as.character(standardized$warnings),
    provenance = .sn_analysis_provenance(random_seed = backend_control$seed %||% 717L)
  )
  sn_validate_result(result)
  object <- sn_store_result(object, "fate", result_id, result)
  object <- .sn_log_seurat_command(object, assay = SeuratObject::DefaultAssay(object), name = "sn_run_fate")
  if (isTRUE(return_object)) object else sn_get_result(object, "fate", result_id)
}
