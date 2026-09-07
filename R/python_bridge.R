# Seurat-to-Python bridge.
#
# Extracted from package_tools.R: helpers that serialize Seurat objects for pixi Python backends, execute runner scripts, and import results (including infercnvpy).
.sn_run_python_object_method <- function(object,
                                         environment,
                                         script_name,
                                         method,
                                         assay = NULL,
                                         layer = NULL,
                                         reference_object = NULL,
                                         reference_assay = NULL,
                                         reference_layer = NULL,
                                         spatial_cols = NULL,
                                         output_dir = NULL,
                                         runtime_dir = NULL,
                                         metadata_prefix = paste0(method, "_"),
                                         result_name = method,
                                         return_object = TRUE,
                                         config = list(),
                                         metadata_columns = NULL,
                                         reference_metadata_columns = NULL,
                                         keep_run_dir = NULL,
                                         max_artifact_import_gb = 0.5,
                                         ...) {
  check_installed(pkg = "Seurat", reason = glue::glue("to run {method} on a Seurat object."))
  if (method %in% c("scarches", "stlearn")) {
    stop(
      "The packaged `", method, "` runner is disabled because Shennong does not yet ",
      "ship a faithful upstream ", method, " workflow. Use a validated upstream workflow ",
      "until an admitted backend contract is available.",
      call. = FALSE
    )
  }
  assay <- assay %||% SeuratObject::DefaultAssay(object = object)
  layer <- layer %||% .sn_select_python_object_layer(object = object, assay = assay)
  raw_count_backends <- c(cell2location = "cell2location", scpoli = "scPoli")
  requires_raw_counts <- method %in% names(raw_count_backends)
  if (requires_raw_counts && !.sn_name_declares_count_scale(layer)) {
    stop(
      unname(raw_count_backends[[method]]),
      " requires a raw/count-like input layer; `", layer,
      "` is not declared count-scale.",
      call. = FALSE
    )
  }
  output_supplied <- !is.null(output_dir)
  if (is.null(keep_run_dir)) keep_run_dir <- output_supplied
  if (!is.logical(keep_run_dir) || length(keep_run_dir) != 1L || is.na(keep_run_dir)) {
    stop("`keep_run_dir` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!isTRUE(keep_run_dir) && identical(method, "spatialdata")) {
    stop(
      "SpatialData produces a disk-backed Zarr artifact. Supply `output_dir` or ",
      "set `keep_run_dir = TRUE` to retain it explicitly.",
      call. = FALSE
    )
  }
  output_dir <- .sn_resolve_python_run_directory(
    path = output_dir,
    method = method,
    runtime_dir = runtime_dir,
    keep_run_dir = keep_run_dir,
    supplied = output_supplied
  )
  run_complete <- FALSE
  failure_stage <- "prepare"
  if (!isTRUE(keep_run_dir)) {
    on.exit({
      if (!run_complete && dir.exists(output_dir)) {
        .sn_sanitize_failed_python_run(output_dir, method = method, stage = failure_stage)
      }
    }, add = TRUE)
  }
  input_dir <- file.path(output_dir, "input")
  result_dir <- file.path(output_dir, "output")
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  spatial_cols <- .sn_resolve_spatial_cols(object = object, spatial_cols = spatial_cols, required = environment %in% c("cell2location", "tangram", "squidpy", "spatialdata", "stlearn"))
  metadata_columns <- .sn_python_required_metadata_columns(
    method = method,
    config = config,
    explicit = metadata_columns,
    reference = FALSE
  )
  query <- .sn_write_python_object_input(
    object = object,
    input_dir = file.path(input_dir, "query"),
    assay = assay,
    layer = layer,
    spatial_cols = spatial_cols,
    metadata_columns = metadata_columns,
    require_raw_counts = requires_raw_counts
  )

  reference <- NULL
  if (!is.null(reference_object)) {
    reference_assay <- reference_assay %||% SeuratObject::DefaultAssay(object = reference_object)
    reference_layer <- reference_layer %||% .sn_select_python_object_layer(object = reference_object, assay = reference_assay)
    reference_metadata_columns <- .sn_python_required_metadata_columns(
      method = method,
      config = config,
      explicit = reference_metadata_columns,
      reference = TRUE
    )
    reference <- .sn_write_python_object_input(
      object = reference_object,
      input_dir = file.path(input_dir, "reference"),
      assay = reference_assay,
      layer = reference_layer,
      spatial_cols = NULL,
      metadata_columns = reference_metadata_columns
    )
  }

  config <- .sn_prepare_python_object_config(config = config, input_dir = input_dir)
  config <- c(
    list(
      method = method,
      assay = assay,
      layer = layer,
      reference_assay = reference_assay,
      reference_layer = reference_layer
    ),
    config
  )
  config_path <- .sn_write_json_file(config, file.path(output_dir, paste0(method, "_config.json")))
  script <- .sn_pixi_script_path(environment = environment, script_name = script_name)
  failure_stage <- "execute"
  execution_error <- tryCatch({
    .sn_execute_python_object_pixi(
      environment = environment,
      script = script,
      input_dir = input_dir,
      output_dir = result_dir,
      config_path = config_path,
      ...
    )
    NULL
  }, error = identity)
  if (!is.null(execution_error)) {
    sanitization_complete <- NA
    if (!isTRUE(keep_run_dir)) {
      sanitization_complete <- .sn_sanitize_failed_python_run(
        output_dir,
        method = method,
        stage = "execute"
      )
    }
    stop(
      conditionMessage(execution_error),
      if (!isTRUE(keep_run_dir)) .sn_python_failure_suffix(output_dir, sanitization_complete),
      call. = FALSE
    )
  }
  failure_stage <- "import"
  imported <- tryCatch(
    .sn_import_python_object_results(
      object = object,
      method = method,
      result_name = result_name,
      output_dir = result_dir,
      run_dir = output_dir,
      assay = assay,
      query = query,
      reference = reference,
      config = config,
      metadata_prefix = metadata_prefix,
      return_object = return_object,
      retain_run_dir = isTRUE(keep_run_dir),
      max_artifact_import_gb = max_artifact_import_gb
    ),
    error = identity
  )
  if (inherits(imported, "error")) {
    sanitization_complete <- NA
    if (!isTRUE(keep_run_dir)) {
      sanitization_complete <- .sn_sanitize_failed_python_run(
        output_dir,
        method = method,
        stage = "import"
      )
    }
    stop(
      conditionMessage(imported),
      if (!isTRUE(keep_run_dir)) .sn_python_failure_suffix(output_dir, sanitization_complete),
      call. = FALSE
    )
  }
  if (!isTRUE(keep_run_dir)) {
    .sn_remove_python_run_directory(output_dir, label = paste0(method, " temporary run directory"))
  }
  run_complete <- TRUE
  imported
}

.sn_python_required_metadata_columns <- function(method,
                                                 config,
                                                 explicit = NULL,
                                                 reference = FALSE) {
  columns <- explicit
  if (is.null(columns)) {
    columns <- if (isTRUE(reference)) {
      if (identical(method, "tangram")) config$cell_type_key else NULL
    } else {
      switch(
        method,
        scpoli = c(config$batch_key, config$labels_key),
        cellphonedb = config$groupby,
        squidpy = config$cluster_key,
        NULL
      )
    }
  }
  columns <- unique(as.character(unlist(columns, use.names = FALSE)))
  columns[!is.na(columns) & nzchar(columns)]
}

.sn_prepare_python_run_directory <- function(path) {
  path <- path.expand(as.character(path)[[1]])
  if (!nzchar(path)) {
    stop("Python backend `output_dir` must be a non-empty path.", call. = FALSE)
  }
  if (file.exists(path) && !dir.exists(path)) {
    stop("Python backend `output_dir` exists and is not a directory: ", path, call. = FALSE)
  }
  if (dir.exists(path)) {
    existing <- list.files(path, all.files = TRUE, no.. = TRUE)
    if (length(existing) > 0L) {
      stop(
        "Python backend `output_dir` is not empty; refusing to mix a new run with stale outputs: ",
        normalizePath(path, winslash = "/", mustWork = TRUE),
        call. = FALSE
      )
    }
  } else if (!dir.create(path, recursive = TRUE, showWarnings = FALSE)) {
    stop("Could not create Python backend `output_dir`: ", path, call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.sn_resolve_python_run_directory <- function(path = NULL,
                                             method,
                                             runtime_dir = NULL,
                                             keep_run_dir,
                                             supplied = !is.null(path)) {
  if (isTRUE(keep_run_dir)) {
    retained_path <- path %||% .sn_default_python_run_dir(
      method = method,
      runtime_dir = runtime_dir,
      temporary = FALSE
    )
    return(.sn_prepare_python_run_directory(retained_path))
  }

  parent <- if (isTRUE(supplied)) {
    path
  } else {
    dirname(.sn_default_python_run_dir(
      method = method,
      runtime_dir = runtime_dir,
      temporary = TRUE
    ))
  }
  safe_method <- gsub("[^A-Za-z0-9_.-]+", "_", as.character(method)[[1L]])
  if (!nzchar(safe_method) || safe_method %in% c(".", "..")) {
    stop("`method` must identify a safe Python backend run directory.", call. = FALSE)
  }
  .sn_create_owned_run_dir(
    parent = parent,
    prefix = paste0(safe_method, "_")
  )
}

.sn_remove_python_run_directory <- function(path,
                                            label = "temporary Python run directory",
                                            unlink_fn = unlink) {
  .sn_cleanup_owned_run_dir(
    run_dir = path,
    unlink_fn = unlink_fn,
    label = label
  )
}

.sn_sanitize_failed_python_run <- function(run_dir, method, stage) {
  if (!.sn_is_owned_run_dir(run_dir)) {
    stop(
      "Refusing to sanitize a Python run directory without a Shennong ownership marker: ",
      normalizePath(run_dir, winslash = "/", mustWork = FALSE), ".",
      call. = FALSE
    )
  }
  input_dir <- file.path(run_dir, "input")
  unlink(input_dir, recursive = TRUE, force = TRUE)
  config_files <- list.files(
    run_dir,
    pattern = "(^config\\.json$|_config\\.json$)",
    full.names = TRUE
  )
  unlink(config_files, force = TRUE)
  output_dir <- file.path(run_dir, "output")
  backend_summary <- list()
  if (dir.exists(output_dir)) {
    manifest_path <- file.path(output_dir, "manifest.json")
    if (file.exists(manifest_path)) {
      manifest <- tryCatch(jsonlite::read_json(manifest_path, simplifyVector = TRUE), error = function(...) NULL)
      if (is.list(manifest)) {
        safe_fields <- names(manifest)[
          names(manifest) %in% c("method", "status") |
            grepl("^n_[A-Za-z0-9_]+$|_version$", names(manifest))
        ]
        backend_summary <- manifest[safe_fields]
        backend_summary <- backend_summary[vapply(
          backend_summary,
          function(value) length(value) <= 1L && (is.atomic(value) || is.null(value)),
          logical(1)
        )]
      }
    }
    # Backend outputs can contain expression-derived, cell-level, or model
    # artifacts under arbitrary filenames. Retain only bounded process metrics;
    # the structural manifest summary is copied into failure.json below.
    entries <- list.files(output_dir, all.files = TRUE, no.. = TRUE, full.names = TRUE)
    keep <- basename(entries) %in% "resource-usage.txt"
    unlink(entries[!keep], recursive = TRUE, force = TRUE)
  }
  unsafe_remaining <- c(
    if (dir.exists(input_dir)) input_dir else character(),
    config_files[file.exists(config_files)],
    if (dir.exists(output_dir)) {
      entries <- list.files(output_dir, all.files = TRUE, no.. = TRUE, full.names = TRUE)
      entries[!basename(entries) %in% "resource-usage.txt"]
    } else {
      character()
    }
  )
  sanitization_complete <- length(unsafe_remaining) == 0L
  jsonlite::write_json(
    list(
      method = method,
      stage = stage,
      status = "failed",
      sanitization_complete = sanitization_complete,
      timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
      backend = backend_summary,
      retained_files = if (dir.exists(output_dir)) list.files(output_dir) else character(),
      unsafe_remaining = unsafe_remaining
    ),
    file.path(run_dir, "failure.json"),
    auto_unbox = TRUE,
    pretty = TRUE
  )
  invisible(sanitization_complete)
}

.sn_python_failure_suffix <- function(run_dir, sanitization_complete) {
  if (isTRUE(sanitization_complete)) {
    paste0(" Sanitized diagnostics remain at: ", run_dir)
  } else {
    paste0(
      " Temporary run data could not be fully sanitized; inspect and remove it manually: ",
      run_dir
    )
  }
}

.sn_with_integration_python_run <- function(method,
                                            integration_control,
                                            code) {
  if (!is.list(integration_control)) {
    stop("`integration_control` must be a list.", call. = FALSE)
  }
  supplied_run_dir <- !is.null(integration_control$run_dir)
  keep_run_dir <- integration_control$keep_run_dir
  if (is.null(keep_run_dir)) keep_run_dir <- supplied_run_dir
  if (!is.logical(keep_run_dir) || length(keep_run_dir) != 1L || is.na(keep_run_dir)) {
    stop("`integration_control$keep_run_dir` must be TRUE or FALSE.", call. = FALSE)
  }

  runtime_dir <- .sn_shennong_runtime_dir(integration_control$runtime_dir %||% NULL)
  run_dir <- .sn_resolve_python_run_directory(
    path = integration_control$run_dir,
    method = method,
    runtime_dir = runtime_dir,
    keep_run_dir = keep_run_dir,
    supplied = supplied_run_dir
  )
  integration_control$run_dir <- run_dir
  integration_control$keep_run_dir <- keep_run_dir

  complete <- FALSE
  on.exit({
    if (!complete && !isTRUE(keep_run_dir) && dir.exists(run_dir)) {
      .sn_sanitize_failed_python_run(run_dir, method = method, stage = "integration")
    }
  }, add = TRUE)

  result <- tryCatch(code(integration_control), error = identity)
  if (inherits(result, "error")) {
    sanitization_complete <- NA
    if (!isTRUE(keep_run_dir)) {
      sanitization_complete <- .sn_sanitize_failed_python_run(
        run_dir,
        method = method,
        stage = "integration"
      )
    }
    complete <- TRUE
    stop(
      conditionMessage(result),
      if (!isTRUE(keep_run_dir)) .sn_python_failure_suffix(run_dir, sanitization_complete),
      call. = FALSE
    )
  }

  if (!is.null(result$object) && inherits(result$object, "Seurat")) {
    result$object@misc$integration$run_dir_retained <- isTRUE(keep_run_dir)
    result$object@misc$integration$output_retained <- isTRUE(keep_run_dir)
    if (!isTRUE(keep_run_dir)) {
      result$object@misc$integration["run_dir"] <- list(NULL)
      result$object@misc$integration["output_h5ad"] <- list(NULL)
      if (!is.null(result$object@misc$integration$backend_performance)) {
        result$object@misc$integration$backend_performance$resource_path <- NULL
      }
    }
  }
  if (!isTRUE(keep_run_dir)) {
    .sn_remove_python_run_directory(run_dir, label = paste0(method, " integration run directory"))
  }
  complete <- TRUE
  result
}

.sn_read_integration_backend_manifest <- function(output_dir,
                                                  method,
                                                  n_cells) {
  manifest_path <- file.path(output_dir, "manifest.json")
  if (!file.exists(manifest_path)) {
    stop(method, " output is missing required `manifest.json`: ", manifest_path, call. = FALSE)
  }
  manifest <- tryCatch(
    jsonlite::read_json(manifest_path, simplifyVector = TRUE),
    error = function(error) {
      stop("Could not parse ", method, " `manifest.json`: ", conditionMessage(error), call. = FALSE)
    }
  )
  if (!is.list(manifest)) {
    stop(method, " `manifest.json` must contain a JSON object.", call. = FALSE)
  }
  manifest_method <- as.character(manifest$method %||% character())
  if (length(manifest_method) != 1L || !identical(manifest_method, method)) {
    stop(
      method, " manifest method mismatch: expected `", method,
      "`, received `", paste(manifest_method, collapse = ", "), "`.",
      call. = FALSE
    )
  }
  manifest_cells <- suppressWarnings(as.integer(manifest$n_cells %||% NA_integer_))
  if (length(manifest_cells) != 1L || is.na(manifest_cells) || manifest_cells != n_cells) {
    stop(
      method, " manifest must report `n_cells = ", n_cells,
      "` for the current input.",
      call. = FALSE
    )
  }
  manifest
}

.sn_integration_manifest_artifact <- function(manifest,
                                              field,
                                              output_dir) {
  value <- manifest[[field]]
  if (is.null(value) || length(value) == 0L) {
    return(NULL)
  }
  if (!is.atomic(value) || is.factor(value) || length(value) != 1L ||
      is.na(value) || !nzchar(as.character(value))) {
    stop("Python integration manifest `", field, "` must be one scalar path.", call. = FALSE)
  }
  value <- path.expand(as.character(value))
  root <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  candidate <- if (file.exists(value) || dir.exists(value)) value else file.path(root, value)
  if (!file.exists(candidate) && !dir.exists(candidate)) {
    stop("Python integration manifest declares missing `", field, "`: ", value, call. = FALSE)
  }
  candidate <- normalizePath(candidate, winslash = "/", mustWork = TRUE)
  if (!identical(candidate, root) && !startsWith(candidate, paste0(root, "/"))) {
    stop("Python integration manifest artifact `", field, "` is outside its output directory.", call. = FALSE)
  }
  candidate
}

.sn_select_python_object_layer <- function(object, assay) {
  layers <- SeuratObject::Layers(object[[assay]])
  if ("data" %in% layers || any(grepl("^data\\.", layers))) {
    return("data")
  }
  "counts"
}

.sn_pixi_script_path <- function(environment, script_name) {
  installed <- system.file("pixi", environment, "scripts", script_name, package = "Shennong")
  if (nzchar(installed) && file.exists(installed)) {
    return(normalizePath(installed, winslash = "/", mustWork = TRUE))
  }
  source_path <- file.path(getwd(), "inst", "pixi", environment, "scripts", script_name)
  if (file.exists(source_path)) {
    return(normalizePath(source_path, winslash = "/", mustWork = TRUE))
  }
  stop("Could not locate Shennong pixi runner script: ", file.path(environment, "scripts", script_name), call. = FALSE)
}

.sn_write_python_object_input <- function(object,
                                          input_dir,
                                          assay,
                                          layer,
                                          spatial_cols = NULL,
                                          metadata_columns = NULL,
                                          require_raw_counts = FALSE) {
  if (!is.logical(require_raw_counts) || length(require_raw_counts) != 1L ||
      is.na(require_raw_counts)) {
    stop("`require_raw_counts` must be TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(require_raw_counts) && !.sn_name_declares_count_scale(layer)) {
    stop(
      "Python raw-count input requires a raw/count-like layer; `", layer,
      "` is not declared count-scale.",
      call. = FALSE
    )
  }
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  expr <- .sn_get_seurat_layer_data(object = object, assay = assay, layer = layer)
  expr <- .sn_as_sparse_matrix(expr)
  expr <- expr[, colnames(object), drop = FALSE]
  if (is.null(rownames(expr)) || is.null(colnames(expr)) ||
      anyNA(rownames(expr)) || anyNA(colnames(expr)) ||
      any(!nzchar(rownames(expr))) || any(!nzchar(colnames(expr))) ||
      anyDuplicated(rownames(expr)) || anyDuplicated(colnames(expr))) {
    stop("Python backend input requires unique, non-empty feature and cell identifiers.", call. = FALSE)
  }
  if (isTRUE(require_raw_counts)) {
    expr <- .sn_validate_python_raw_counts(expr, label = "Python raw-count input")
  } else {
    stored_values <- if (inherits(expr, "sparseMatrix")) expr@x else as.numeric(expr)
    if (anyNA(stored_values) || any(!is.finite(stored_values))) {
      stop("Python backend input expression must contain only finite values.", call. = FALSE)
    }
  }

  matrix_path <- file.path(input_dir, "matrix.mtx")
  obs_path <- file.path(input_dir, "obs.csv")
  var_path <- file.path(input_dir, "var.csv")
  Matrix::writeMM(obj = expr, file = matrix_path)

  full_metadata <- object[[]][colnames(expr), , drop = FALSE]
  missing_metadata <- setdiff(metadata_columns, colnames(full_metadata))
  if (length(missing_metadata) > 0L) {
    stop(
      "Required Python backend metadata column(s) were not found: ",
      paste(missing_metadata, collapse = ", "),
      call. = FALSE
    )
  }
  obs <- full_metadata[, metadata_columns, drop = FALSE]
  obs <- data.frame(cell_id = rownames(full_metadata), obs, check.names = FALSE)
  utils::write.csv(obs, file = obs_path, row.names = FALSE, quote = TRUE)
  var <- data.frame(feature_id = rownames(expr), stringsAsFactors = FALSE)
  utils::write.csv(var, file = var_path, row.names = FALSE, quote = TRUE)

  spatial_path <- NULL
  if (!is.null(spatial_cols)) {
    spatial <- full_metadata[, spatial_cols, drop = FALSE]
    rownames(spatial) <- rownames(full_metadata)
    spatial_path <- file.path(input_dir, "spatial.csv")
    utils::write.csv(spatial, file = spatial_path, row.names = TRUE, quote = TRUE)
  }

  list(
    input_dir = normalizePath(input_dir, winslash = "/", mustWork = TRUE),
    matrix_path = normalizePath(matrix_path, winslash = "/", mustWork = TRUE),
    obs_path = normalizePath(obs_path, winslash = "/", mustWork = TRUE),
    var_path = normalizePath(var_path, winslash = "/", mustWork = TRUE),
    spatial_path = if (!is.null(spatial_path)) normalizePath(spatial_path, winslash = "/", mustWork = TRUE) else NULL,
    assay = assay,
    layer = layer,
    features = rownames(expr),
    cells = colnames(expr),
    metadata_columns = metadata_columns,
    n_features = nrow(expr),
    n_cells = ncol(expr)
  )
}

.sn_validate_python_raw_counts <- function(counts, label = "Python backend input") {
  counts <- .sn_as_sparse_matrix(counts)
  if (is.null(rownames(counts)) || is.null(colnames(counts)) ||
      anyNA(rownames(counts)) || anyNA(colnames(counts)) ||
      any(!nzchar(rownames(counts))) || any(!nzchar(colnames(counts))) ||
      anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts))) {
    stop(label, " must have unique, non-empty feature and cell identifiers.", call. = FALSE)
  }
  values <- if (inherits(counts, "sparseMatrix")) counts@x else as.numeric(counts)
  if (anyNA(values) || any(!is.finite(values)) || any(values < 0)) {
    stop(label, " must contain finite, non-negative raw counts.", call. = FALSE)
  }
  if (any(abs(values - round(values)) > 1e-8)) {
    stop(label, " must contain integer-like raw counts; values are never rounded silently.", call. = FALSE)
  }
  counts
}

.sn_resolve_spatial_cols <- function(object, spatial_cols = NULL, required = FALSE) {
  metadata_cols <- colnames(object[[]])
  if (!is.null(spatial_cols)) {
    if (length(spatial_cols) != 2L || !all(spatial_cols %in% metadata_cols)) {
      stop("`spatial_cols` must contain two metadata columns present in `object`.", call. = FALSE)
    }
    return(spatial_cols)
  }
  candidates <- list(
    c("x", "y"),
    c("spatial_x", "spatial_y"),
    c("imagecol", "imagerow"),
    c("pxl_col_in_fullres", "pxl_row_in_fullres"),
    c("array_col", "array_row")
  )
  for (candidate in candidates) {
    if (all(candidate %in% metadata_cols)) {
      return(candidate)
    }
  }
  if (isTRUE(required)) {
    stop("Spatial object workflow requires `spatial_cols` or recognized coordinate metadata columns.", call. = FALSE)
  }
  NULL
}

.sn_prepare_python_object_config <- function(config, input_dir) {
  if ("reference_signatures" %in% names(config)) {
    signatures <- .sn_validate_cell2location_reference_signatures(
      config$reference_signatures
    )
    signatures_path <- file.path(input_dir, "reference_signatures.csv")
    utils::write.csv(signatures, file = signatures_path, row.names = TRUE, quote = TRUE)
    config$reference_signatures <- normalizePath(signatures_path, winslash = "/", mustWork = TRUE)
  }
  config
}

.sn_validate_cell2location_reference_signatures <- function(signatures) {
  if (is.character(signatures) && length(signatures) == 1L &&
      !is.na(signatures) && nzchar(signatures)) {
    path <- path.expand(signatures)
    if (!file.exists(path) || dir.exists(path)) {
      stop("`reference_signatures` CSV file does not exist: ", path, call. = FALSE)
    }
    raw <- tryCatch(
      utils::read.csv(
        path,
        check.names = FALSE,
        stringsAsFactors = FALSE,
        na.strings = c("NA", "NaN", "Inf", "-Inf")
      ),
      error = function(error) {
        stop("Could not read `reference_signatures` CSV: ", conditionMessage(error), call. = FALSE)
      }
    )
    if (ncol(raw) < 2L) {
      stop("`reference_signatures` CSV must contain feature IDs and at least one cell-state column.", call. = FALSE)
    }
    feature_ids <- as.character(raw[[1L]])
    if (anyNA(feature_ids) || any(!nzchar(feature_ids)) || anyDuplicated(feature_ids)) {
      stop("`reference_signatures` feature identifiers must be unique and non-empty.", call. = FALSE)
    }
    signatures <- raw[-1L]
    rownames(signatures) <- feature_ids
  } else if (is.data.frame(signatures)) {
    signatures <- as.data.frame(signatures, check.names = FALSE)
  } else if (is.matrix(signatures)) {
    if (!is.numeric(signatures)) {
      stop("`reference_signatures` matrix must be numeric.", call. = FALSE)
    }
    if (is.null(rownames(signatures)) || is.null(colnames(signatures))) {
      stop("`reference_signatures` matrix must have feature and cell-state identifiers.", call. = FALSE)
    }
    if (anyNA(rownames(signatures)) || any(!nzchar(rownames(signatures))) ||
        anyDuplicated(rownames(signatures)) || anyNA(colnames(signatures)) ||
        any(!nzchar(colnames(signatures))) || anyDuplicated(colnames(signatures))) {
      stop(
        "`reference_signatures` must have unique, non-empty feature and cell-state identifiers.",
        call. = FALSE
      )
    }
    signatures <- as.data.frame(signatures, check.names = FALSE)
  } else {
    stop(
      "`reference_signatures` must be a CSV path or a numeric data frame/matrix.",
      call. = FALSE
    )
  }

  feature_ids <- rownames(signatures)
  state_ids <- colnames(signatures)
  if (nrow(signatures) < 1L || ncol(signatures) < 1L ||
      is.null(feature_ids) || anyNA(feature_ids) || any(!nzchar(feature_ids)) ||
      anyDuplicated(feature_ids) || is.null(state_ids) || anyNA(state_ids) ||
      any(!nzchar(state_ids)) || anyDuplicated(state_ids)) {
    stop(
      "`reference_signatures` must have at least one row and column with unique, non-empty feature and cell-state identifiers.",
      call. = FALSE
    )
  }
  numeric_columns <- vapply(signatures, is.numeric, logical(1))
  if (!all(numeric_columns)) {
    stop("`reference_signatures` values must be numeric.", call. = FALSE)
  }
  values <- as.matrix(signatures)
  if (anyNA(values) || any(!is.finite(values)) || any(values < 0)) {
    stop("`reference_signatures` values must be finite and non-negative.", call. = FALSE)
  }
  signatures
}

.sn_execute_python_object_pixi <- function(environment,
                                          script,
                                          input_dir,
                                          output_dir,
                                          config_path,
                                          ...) {
  sn_call_pixi_environment(
    environment = environment,
    command = "python",
    args = c(
      shQuote(script),
      "--input-dir", shQuote(input_dir),
      "--output-dir", shQuote(output_dir),
      "--config", shQuote(config_path)
    ),
    ...
  )
}

.sn_import_python_object_results <- function(object,
                                            method,
                                            result_name,
                                            output_dir,
                                            run_dir,
                                            assay,
                                            query,
                                            reference = NULL,
                                            config = list(),
                                            metadata_prefix = paste0(method, "_"),
                                            return_object = TRUE,
                                            retain_run_dir = TRUE,
                                            max_artifact_import_gb = 0.5) {
  manifest_path <- file.path(output_dir, "manifest.json")
  if (!file.exists(manifest_path)) {
    stop("Python backend output is missing required `manifest.json`: ", manifest_path, call. = FALSE)
  }
  backend_manifest <- tryCatch(
    jsonlite::read_json(manifest_path, simplifyVector = TRUE),
    error = function(e) {
      stop("Could not parse Python backend manifest: ", conditionMessage(e), call. = FALSE)
    }
  )
  if (!is.list(backend_manifest)) {
    stop("Python backend `manifest.json` must contain a JSON object.", call. = FALSE)
  }
  if (is.null(names(backend_manifest)) || anyNA(names(backend_manifest)) ||
      any(!nzchar(names(backend_manifest))) || anyDuplicated(names(backend_manifest))) {
    stop("Python backend manifest fields must have unique, non-empty names.", call. = FALSE)
  }
  manifest_method <- as.character(backend_manifest$method %||% "")
  if (length(manifest_method) != 1L || !nzchar(manifest_method)) {
    stop("Python backend manifest must declare a non-empty `method`.", call. = FALSE)
  }
  if (!identical(tolower(manifest_method), tolower(method))) {
    stop(
      "Python backend manifest method mismatch: expected '", method,
      "', received '", manifest_method, "'.",
      call. = FALSE
    )
  }
  if (!is.null(backend_manifest$n_cells)) {
    manifest_cells <- suppressWarnings(as.integer(backend_manifest$n_cells))
    if (length(manifest_cells) != 1L || is.na(manifest_cells) ||
        manifest_cells != as.integer(query$n_cells)) {
      stop("Python backend manifest reports a different number of cells than the exported input.", call. = FALSE)
    }
  }
  if (!is.null(backend_manifest$n_features)) {
    manifest_features <- suppressWarnings(as.integer(backend_manifest$n_features))
    if (length(manifest_features) != 1L || is.na(manifest_features) ||
        manifest_features != as.integer(query$n_features)) {
      stop("Python backend manifest reports a different number of features than the exported input.", call. = FALSE)
    }
  }
  artifact_paths <- .sn_collect_python_artifacts(
    manifest = backend_manifest,
    output_dir = output_dir
  )

  obs_path <- file.path(output_dir, "obs.csv")
  expected_cells <- query$cells %||% colnames(object)
  imported_metadata <- character(0)
  if (file.exists(obs_path)) {
    .sn_assert_python_artifact_budget(obs_path, max_artifact_import_gb, "Python metadata output")
    metadata <- utils::read.csv(obs_path, row.names = 1, check.names = FALSE)
    .sn_validate_exact_python_ids(
      ids = rownames(metadata),
      expected = expected_cells,
      label = "Python metadata output cell"
    )
    if (ncol(metadata) > 0L) {
      if (anyNA(colnames(metadata)) || any(!nzchar(colnames(metadata))) ||
          anyDuplicated(colnames(metadata))) {
        stop("Python metadata output columns must be unique and non-empty.", call. = FALSE)
      }
      metadata <- metadata[expected_cells, , drop = FALSE]
      if (method %in% c("cell2location", "tangram")) {
        metadata[] <- lapply(names(metadata), function(column) {
          raw <- metadata[[column]]
          values <- suppressWarnings(as.numeric(as.character(raw)))
          if (any(!is.na(raw) & is.na(values)) || any(!is.finite(values)) || any(values < 0)) {
            stop(
              "Python ", method, " metadata column '", column,
              "' must contain finite non-negative numeric values.",
              call. = FALSE
            )
          }
          values
        })
      }
      colnames(metadata) <- .sn_prefix_metadata_columns(colnames(metadata), metadata_prefix)
      if (anyDuplicated(colnames(metadata))) {
        stop("Python metadata column names collide after applying `metadata_prefix`.", call. = FALSE)
      }
      imported_metadata <- colnames(metadata)
      object <- Seurat::AddMetaData(object = object, metadata = metadata)
    }
  }

  reductions <- character(0)
  embedding_files <- list.files(output_dir, pattern = "\\.csv$", full.names = TRUE)
  embedding_files <- embedding_files[grepl("(latent|pca|umap|embedding)", basename(embedding_files), ignore.case = TRUE)]
  for (embedding_file in embedding_files) {
    .sn_assert_python_artifact_budget(
      embedding_file,
      max_artifact_import_gb,
      paste0("Python embedding `", basename(embedding_file), "`")
    )
    embedding <- tryCatch(
      .sn_read_embedding_csv(embedding_file, cells = expected_cells),
      error = function(e) {
        stop(
          "Invalid Python embedding output `", basename(embedding_file), "`: ",
          conditionMessage(e),
          call. = FALSE
        )
      }
    )
    reduction_name <- paste0(metadata_prefix, tools::file_path_sans_ext(basename(embedding_file)))
    reduction_key <- paste0(gsub("[^A-Za-z0-9]", "", toupper(reduction_name)), "_")
    colnames(embedding) <- paste0(reduction_key, seq_len(ncol(embedding)))
    object[[reduction_name]] <- Seurat::CreateDimReducObject(
      embeddings = embedding,
      key = reduction_key,
      assay = assay
    )
    reductions <- c(reductions, reduction_name)
  }
  .sn_validate_python_result_contract(
    method = method,
    imported_metadata = imported_metadata,
    imported_reductions = reductions,
    artifact_paths = artifact_paths
  )
  artifact_tables <- .sn_import_python_artifact_tables(
    method = method,
    artifact_paths = artifact_paths,
    query = query,
    reference = reference,
    import = !isTRUE(retain_run_dir),
    max_import_gb = max_artifact_import_gb
  )
  safe_backend_manifest <- .sn_sanitize_python_backend_manifest(
    backend_manifest,
    retain_paths = isTRUE(retain_run_dir)
  )
  run_manifest <- list(
    method = method,
    result_name = result_name,
    assay = assay,
    source_layer = query$layer,
    run_dir = if (isTRUE(retain_run_dir)) normalizePath(run_dir, winslash = "/", mustWork = TRUE) else NULL,
    output_dir = if (isTRUE(retain_run_dir)) normalizePath(output_dir, winslash = "/", mustWork = TRUE) else NULL,
    run_dir_retained = isTRUE(retain_run_dir),
    input_features = query$features,
    n_features = query$n_features,
    n_cells = query$n_cells,
    reference = .sn_python_input_summary(reference),
    imported_metadata = imported_metadata,
    imported_reductions = reductions,
    imported_artifacts = if (isTRUE(retain_run_dir)) artifact_paths else character(),
    imported_tables = artifact_tables,
    config = .sn_safe_python_config(config)
  )
  safe_backend_manifest <- safe_backend_manifest[
    !names(safe_backend_manifest) %in% names(run_manifest)
  ]
  run_manifest <- c(run_manifest, safe_backend_manifest)
  if (anyDuplicated(names(run_manifest))) {
    stop("Internal error: Python run manifest contains duplicate field names.", call. = FALSE)
  }
  object@misc[[method]] <- object@misc[[method]] %||% list()
  object@misc[[method]][[result_name]] <- run_manifest
  object <- .sn_log_seurat_command(object = object, name = paste0("sn_run_", method))
  if (isTRUE(return_object)) {
    return(object)
  }
  run_manifest
}

.sn_safe_python_config <- function(config) {
  if (!is.list(config) || length(config) == 0L) return(list())
  stats::setNames(lapply(names(config), function(name) {
    .sn_usage_expression_summary(config[[name]], name = name)
  }), names(config))
}

.sn_python_input_summary <- function(input) {
  if (is.null(input)) return(NULL)
  input[intersect(
    c("assay", "layer", "features", "n_features", "n_cells", "metadata_columns"),
    names(input)
  )]
}

.sn_sanitize_python_backend_manifest <- function(manifest, retain_paths = TRUE) {
  if (!is.list(manifest)) return(list(output_retained = isTRUE(retain_paths)))
  explicit_fields <- c(
    "method", "status", "assay", "source_layer", "batch_key", "labels_key",
    "groupby", "cluster_key", "key_added", "graph_name", "data_key",
    "unlabeled_category", "velocity_mode", "mode", "input_n_features"
  )
  allowed <- names(manifest)[
    names(manifest) %in% explicit_fields |
      grepl("^n_[A-Za-z0-9_]+$|_version$", names(manifest))
  ]
  safe <- manifest[allowed]
  safe <- safe[vapply(
    safe,
    function(value) length(value) == 1L && (is.atomic(value) || is.null(value)),
    logical(1)
  )]
  safe$output_retained <- isTRUE(retain_paths)
  safe
}

.sn_validate_exact_python_ids <- function(ids, expected, label, ordered = FALSE) {
  ids <- as.character(ids)
  expected <- as.character(expected)
  if (length(ids) != length(expected) || anyNA(ids) || any(!nzchar(ids)) ||
      anyDuplicated(ids)) {
    stop(label, " identifiers must be unique, non-empty, and one-per-input.", call. = FALSE)
  }
  matches <- if (isTRUE(ordered)) identical(ids, expected) else setequal(ids, expected)
  if (!matches) {
    stop(label, " identifiers do not exactly match the exported input.", call. = FALSE)
  }
  invisible(TRUE)
}

.sn_assert_python_artifact_budget <- function(paths, max_import_gb, label) {
  if (!is.numeric(max_import_gb) || length(max_import_gb) != 1L ||
      !is.finite(max_import_gb) || max_import_gb <= 0) {
    stop("`max_artifact_import_gb` must be one positive finite number.", call. = FALSE)
  }
  paths <- as.character(paths)
  missing <- paths[!file.exists(paths) | dir.exists(paths)]
  if (length(missing) > 0L) {
    stop("Cannot import missing or non-file ", label, ": ", missing[[1L]], call. = FALSE)
  }
  estimated_gb <- sum(file.info(paths)$size, na.rm = TRUE) * 4 / 1024^3
  if (!is.finite(estimated_gb) || estimated_gb > max_import_gb) {
    stop(
      "Importing ", label, " is estimated to require ",
      format(round(estimated_gb, 2), nsmall = 2),
      " GiB, exceeding `max_artifact_import_gb = ", max_import_gb,
      "`. Reduce the output or increase the budget after reviewing its expected dimensions. ",
      "Persistent retention avoids materializing supported large artifact tables but does not ",
      "bypass validation of metadata or embeddings.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.sn_import_python_artifact_tables <- function(method,
                                              artifact_paths,
                                              query,
                                              reference,
                                              import = FALSE,
                                              max_import_gb = 0.5) {
  tables <- list()
  if (identical(method, "tangram")) {
    path <- unname(artifact_paths["mapping_path"])
    if (is.null(reference)) stop("Tangram mapping import requires reference input identity.", call. = FALSE)
    .sn_validate_tangram_mapping_csv(
      path,
      query_cells = query$cells,
      reference_cells = reference$cells
    )
    if (isTRUE(import)) {
      .sn_assert_python_artifact_budget(path, max_import_gb, "Tangram mapping probabilities")
      mapping <- utils::read.csv(path, row.names = 1, check.names = FALSE)
      mapping <- as.matrix(mapping[reference$cells, query$cells, drop = FALSE])
      suppressWarnings(storage.mode(mapping) <- "numeric")
      tables$mapping <- mapping
    }
  } else if (identical(method, "squidpy")) {
    graph_path <- unname(artifact_paths["spatial_graph_path"])
    enrichment_path <- unname(artifact_paths["neighborhood_enrichment_path"])
    .sn_validate_squidpy_graph_csv(graph_path, cells = query$cells)
    if (length(enrichment_path) == 1L && !is.na(enrichment_path) && nzchar(enrichment_path)) {
      .sn_validate_squidpy_enrichment_csv(enrichment_path)
    }
    if (isTRUE(import)) {
      import_paths <- c(graph_path, enrichment_path)
      import_paths <- import_paths[!is.na(import_paths) & nzchar(import_paths)]
      .sn_assert_python_artifact_budget(import_paths, max_import_gb, "Squidpy graph results")
      graph <- utils::read.csv(graph_path, check.names = FALSE)
      graph$weight <- suppressWarnings(as.numeric(as.character(graph$weight)))
      tables$spatial_graph <- graph
      if (length(enrichment_path) == 1L && !is.na(enrichment_path) && nzchar(enrichment_path)) {
        neighborhood <- utils::read.csv(enrichment_path, check.names = FALSE)
        for (column in intersect(c("zscore", "count"), names(neighborhood))) {
          values <- suppressWarnings(as.numeric(as.character(neighborhood[[column]])))
          if (any(!is.finite(values)) || (identical(column, "count") && any(values < 0))) {
            stop("Squidpy neighborhood enrichment contains invalid `", column, "` values.", call. = FALSE)
          }
          neighborhood[[column]] <- values
        }
        tables$neighborhood_enrichment <- neighborhood
      }
    }
  } else if (identical(method, "cellphonedb")) {
    paths <- unname(artifact_paths[grepl("^result_files", names(artifact_paths))])
    files <- unique(unlist(lapply(paths, function(path) {
      if (dir.exists(path)) {
        list.files(path, pattern = "\\.(csv|txt|tsv)$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
      } else {
        path
      }
    }), use.names = FALSE))
    allowed <- .sn_cellphonedb_standard_result_names()
    files <- files[basename(files) %in% allowed]
    if (length(files) == 0L) {
      stop("CellPhoneDB completed without standard result tables.", call. = FALSE)
    }
    if (isTRUE(import)) {
      .sn_assert_python_artifact_budget(files, max_import_gb, "CellPhoneDB result tables")
    }
    for (path in files) {
      .sn_validate_cellphonedb_table_header(path)
      separator <- if (grepl("\\.csv$", path, ignore.case = TRUE)) "," else "\t"
      table_name <- make.unique(c(names(tables), tools::file_path_sans_ext(basename(path))))[[length(tables) + 1L]]
      if (isTRUE(import)) {
        tables[[table_name]] <- utils::read.table(
          path,
          header = TRUE,
          sep = separator,
          quote = "\"",
          comment.char = "",
          check.names = FALSE
        )
      }
    }
    if (isTRUE(import) && length(tables) == 0L) {
      stop("CellPhoneDB completed without importable result tables.", call. = FALSE)
    }
  }
  if (isTRUE(import)) tables else list()
}

.sn_cellphonedb_standard_result_names <- function() {
  c(
    "deconvoluted.txt", "deconvoluted_percents.txt", "interaction_scores.txt",
    "means.txt", "pvalues.txt", "relevant_interactions.txt", "significant_means.txt"
  )
}

.sn_collect_python_artifacts <- function(manifest, output_dir) {
  artifact_fields <- c(
    "mapping_path", "spatial_graph_path", "neighborhood_enrichment_path",
    "co_occurrence_path", "zarr_path", "result_files", "output_h5ad", "model_dir"
  )
  root <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  artifacts <- character(0)
  for (field in intersect(artifact_fields, names(manifest))) {
    raw_values <- manifest[[field]]
    if (!is.atomic(raw_values) || is.factor(raw_values)) {
      stop("Python backend manifest path field `", field, "` must be a character value.", call. = FALSE)
    }
    values <- as.character(raw_values)
    values <- values[!is.na(values) & nzchar(values)]
    if (!identical(field, "result_files") && length(values) > 1L) {
      stop("Python backend manifest path field `", field, "` must be one scalar path.", call. = FALSE)
    }
    for (index in seq_along(values)) {
      candidate <- path.expand(values[[index]])
      if (!file.exists(candidate) && !dir.exists(candidate)) {
        candidate <- file.path(root, candidate)
      }
      if (!file.exists(candidate) && !dir.exists(candidate)) {
        stop(
          "Python backend manifest declares missing artifact `", field, "`: ",
          values[[index]],
          call. = FALSE
        )
      }
      candidate <- normalizePath(candidate, winslash = "/", mustWork = TRUE)
      if (!identical(candidate, root) && !startsWith(candidate, paste0(root, "/"))) {
        stop(
          "Python backend artifact `", field, "` is outside the run output directory: ",
          candidate,
          call. = FALSE
        )
      }
      if (identical(tolower(as.character(manifest$method %||% "")), "cellphonedb") &&
          identical(field, "result_files") &&
          (dir.exists(candidate) ||
            !basename(candidate) %in% .sn_cellphonedb_standard_result_names())) {
        stop(
          "CellPhoneDB manifest declares a non-standard result artifact: ",
          basename(candidate),
          call. = FALSE
        )
      }
      artifact_name <- if (length(values) == 1L) field else paste0(field, "_", index)
      artifacts[[artifact_name]] <- candidate
    }
  }
  artifacts
}

.sn_validate_python_result_contract <- function(method,
                                                imported_metadata,
                                                imported_reductions,
                                                artifact_paths) {
  has_metadata <- length(imported_metadata) > 0L
  has_reduction <- length(imported_reductions) > 0L
  artifact_fields <- sub("_[0-9]+$", "", names(artifact_paths))
  require_artifact <- function(field, label = field) {
    if (!field %in% artifact_fields) {
      stop(
        "Python backend '", method, "' completed without required ", label,
        " output declared in its manifest.",
        call. = FALSE
      )
    }
  }

  if (identical(method, "scpoli")) {
    if (!has_reduction) {
      stop("scPoli completed without an importable latent embedding.", call. = FALSE)
    }
  } else if (identical(method, "cell2location")) {
    if (!has_metadata) {
      stop("cell2location completed without importable abundance estimates.", call. = FALSE)
    }
  } else if (identical(method, "tangram")) {
    require_artifact("mapping_path", "mapping probability")
  } else if (identical(method, "squidpy")) {
    require_artifact("spatial_graph_path", "spatial graph")
  } else if (identical(method, "spatialdata")) {
    require_artifact("zarr_path", "SpatialData Zarr")
  } else if (identical(method, "cellphonedb")) {
    require_artifact("result_files", "CellPhoneDB result file")
  } else if (!has_metadata && !has_reduction && length(artifact_paths) == 0L) {
    stop(
      "Python backend '", method,
      "' completed without any importable metadata, embedding, or declared artifact.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.sn_run_infercnvpy_object <- function(object,
                                      assay = NULL,
                                      layer = NULL,
                                      species = NULL,
                                      reference_key = NULL,
                                      reference_cat = NULL,
                                      gene_order = NULL,
                                      gtf_file = NULL,
                                      gtf_gene_id = "gene_name",
                                      adata_gene_id = NULL,
                                      output_dir = NULL,
                                      runtime_dir = NULL,
                                      key_added = "cnv",
                                      window_size = 100,
                                      step = 10,
                                      dynamic_threshold = 1.5,
                                      exclude_chromosomes = c("chrX", "chrY"),
                                      chunksize = 5000,
                                      n_jobs = NULL,
                                      calculate_gene_values = FALSE,
                                      lfc_clip = 3,
                                      run_pca = TRUE,
                                      run_neighbors = TRUE,
                                      run_leiden = TRUE,
                                      run_umap = FALSE,
                                      score = TRUE,
                                      leiden_resolution = 1,
                                      cnv_score_group_by = NULL,
                                      metadata_prefix = "infercnvpy_",
                                      result_name = "infercnvpy",
                                      return_object = TRUE,
                                      metadata_columns = NULL,
                                      keep_run_dir = NULL,
                                      max_artifact_import_gb = 0.5,
                                      ...) {
  check_installed(pkg = "Seurat", reason = "to run infercnvpy on a Seurat object.")
  assay <- assay %||% SeuratObject::DefaultAssay(object = object)
  layer <- layer %||% .sn_select_infercnvpy_layer(object = object, assay = assay)
  .sn_validate_log_normalized_layer(
    object = object,
    assay = assay,
    layer = layer,
    backend = "infercnvpy"
  )
  if (!is.null(reference_key) && (!nzchar(reference_key) || !reference_key %in% colnames(object[[]]))) {
    stop("`reference_key` must name a metadata column in `object`.", call. = FALSE)
  }
  .sn_validate_infercnvpy_controls(
    object = object,
    reference_key = reference_key,
    reference_cat = reference_cat,
    run_pca = run_pca,
    run_neighbors = run_neighbors,
    run_leiden = run_leiden,
    run_umap = run_umap,
    score = score,
    cnv_score_group_by = cnv_score_group_by
  )

  output_supplied <- !is.null(output_dir)
  if (is.null(keep_run_dir)) keep_run_dir <- output_supplied
  if (!is.logical(keep_run_dir) || length(keep_run_dir) != 1L || is.na(keep_run_dir)) {
    stop("`keep_run_dir` must be TRUE or FALSE.", call. = FALSE)
  }
  output_dir <- .sn_resolve_python_run_directory(
    path = output_dir,
    method = "infercnvpy",
    runtime_dir = runtime_dir,
    keep_run_dir = keep_run_dir,
    supplied = output_supplied
  )
  run_complete <- FALSE
  failure_stage <- "prepare"
  if (!isTRUE(keep_run_dir)) {
    on.exit({
      if (!run_complete && dir.exists(output_dir)) {
        .sn_sanitize_failed_python_run(output_dir, method = "infercnvpy", stage = failure_stage)
      }
    }, add = TRUE)
  }
  input_dir <- file.path(output_dir, "input")
  result_dir <- file.path(output_dir, "output")
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  input <- .sn_write_infercnvpy_input(
    object = object,
    input_dir = input_dir,
    assay = assay,
    layer = layer,
    species = species,
    gene_order = gene_order,
    gtf_file = gtf_file,
    metadata_columns = unique(c(metadata_columns, reference_key, cnv_score_group_by))
  )
  script <- .sn_infercnvpy_script_path()
  config <- list(
    method = "infercnvpy",
    assay = assay,
    layer = NULL,
    source_layer = layer,
    reference_key = reference_key,
    reference_cat = reference_cat,
    gtf_file = if (!is.null(gtf_file)) normalizePath(path.expand(gtf_file), winslash = "/", mustWork = TRUE) else NULL,
    gtf_gene_id = gtf_gene_id,
    adata_gene_id = adata_gene_id,
    key_added = key_added,
    window_size = window_size,
    step = step,
    dynamic_threshold = dynamic_threshold,
    exclude_chromosomes = exclude_chromosomes,
    chunksize = chunksize,
    n_jobs = n_jobs,
    calculate_gene_values = isTRUE(calculate_gene_values),
    lfc_clip = lfc_clip,
    run_pca = isTRUE(run_pca),
    run_neighbors = isTRUE(run_neighbors),
    run_leiden = isTRUE(run_leiden),
    run_umap = isTRUE(run_umap),
    score = isTRUE(score),
    leiden_resolution = leiden_resolution,
    cnv_score_groupby = cnv_score_group_by,
    write_h5ad = FALSE
  )
  config_path <- .sn_write_json_file(config, file.path(output_dir, "infercnvpy_config.json"))

  failure_stage <- "execute"
  execution_error <- tryCatch({
    .sn_execute_infercnvpy_pixi(
      script = script,
      input_dir = input$input_dir,
      output_dir = result_dir,
      config_path = config_path,
      ...
    )
    NULL
  }, error = identity)
  if (!is.null(execution_error)) {
    sanitization_complete <- NA
    if (!isTRUE(keep_run_dir)) {
      sanitization_complete <- .sn_sanitize_failed_python_run(
        output_dir,
        method = "infercnvpy",
        stage = "execute"
      )
    }
    stop(
      conditionMessage(execution_error),
      if (!isTRUE(keep_run_dir)) .sn_python_failure_suffix(output_dir, sanitization_complete),
      call. = FALSE
    )
  }

  failure_stage <- "import"
  manifest <- tryCatch(
    .sn_import_infercnvpy_results(
      object = object,
      output_dir = result_dir,
      run_dir = output_dir,
      assay = assay,
      input = input,
      config = config,
      metadata_prefix = metadata_prefix,
      result_name = result_name,
      return_object = return_object,
      retain_run_dir = isTRUE(keep_run_dir),
      max_artifact_import_gb = max_artifact_import_gb
    ),
    error = identity
  )
  if (inherits(manifest, "error")) {
    sanitization_complete <- NA
    if (!isTRUE(keep_run_dir)) {
      sanitization_complete <- .sn_sanitize_failed_python_run(
        output_dir,
        method = "infercnvpy",
        stage = "import"
      )
    }
    stop(
      conditionMessage(manifest),
      if (!isTRUE(keep_run_dir)) .sn_python_failure_suffix(output_dir, sanitization_complete),
      call. = FALSE
    )
  }
  if (!isTRUE(keep_run_dir)) {
    .sn_remove_python_run_directory(output_dir, label = "infercnvpy temporary run directory")
  }
  run_complete <- TRUE
  manifest
}


.sn_prefix_metadata_columns <- function(columns, prefix) {
  if (is.null(prefix) || !nzchar(prefix)) {
    return(columns)
  }
  ifelse(startsWith(columns, prefix), columns, paste0(prefix, columns))
}

.sn_read_embedding_csv <- function(path, cells) {
  embedding <- utils::read.csv(path, row.names = 1, check.names = FALSE)
  .sn_validate_exact_python_ids(
    rownames(embedding), cells,
    paste0("Embedding `", basename(path), "` cell")
  )
  if (ncol(embedding) < 1L || anyDuplicated(colnames(embedding)) ||
      anyNA(colnames(embedding)) || any(!nzchar(colnames(embedding)))) {
    stop("Embedding output must have at least one uniquely named dimension: ", path, call. = FALSE)
  }
  embedding <- as.matrix(embedding[cells, , drop = FALSE])
  storage.mode(embedding) <- "numeric"
  if (any(!is.finite(embedding))) {
    stop("Embedding output must contain only finite numeric values: ", path, call. = FALSE)
  }
  embedding
}
