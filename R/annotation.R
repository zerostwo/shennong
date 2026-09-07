.sn_annotation_label_hierarchy <- function(labels) {
  utils::data("marker_genes", package = "Shennong", envir = environment())
  marker_database <- get("marker_genes", envir = environment())
  required <- c("high_hierarchy_cell_type", "low_hierarchy_cell_type")
  if (!all(required %in% colnames(marker_database))) {
    return(stats::setNames(as.character(labels), as.character(labels)))
  }
  hierarchy <- unique(marker_database[required])
  parent <- stats::setNames(
    as.character(hierarchy$high_hierarchy_cell_type),
    as.character(hierarchy$low_hierarchy_cell_type)
  )
  labels <- as.character(labels)
  level_1 <- unname(parent[labels])
  level_1[is.na(level_1) | !nzchar(level_1)] <- labels[is.na(level_1) | !nzchar(level_1)]
  stats::setNames(level_1, labels)
}

.sn_annotation_expression <- function(object, assay = NULL, layer = "data") {
  assay <- assay %||% Seurat::DefaultAssay(object)
  if (!.sn_has_seurat_layer(object, assay = assay, layer = layer)) {
    fallback <- if (.sn_has_seurat_layer(object, assay = assay, layer = "counts")) "counts" else NULL
    if (is_null(fallback)) {
      stop("Neither layer '", layer, "' nor a counts fallback was found in assay '", assay, "'.", call. = FALSE)
    }
    .sn_log_warn("Layer '{layer}' was not found; annotation is using assay '{assay}' layer 'counts'.")
    layer <- fallback
  }
  list(
    matrix = .sn_get_seurat_layer_data(object, assay = assay, layer = layer),
    assay = assay,
    layer = layer
  )
}

.sn_annotation_reference_input <- function(reference, reference_label_by, assay = NULL, layer = "data") {
  if (inherits(reference, "Seurat")) {
    if (!reference_label_by %in% colnames(reference[[]])) {
      stop("`reference_label_by` column '", reference_label_by, "' was not found in reference metadata.", call. = FALSE)
    }
    expression <- .sn_annotation_expression(reference, assay = assay, layer = layer)
    labels <- as.character(reference[[reference_label_by, drop = TRUE]])
    names(labels) <- colnames(reference)
    return(list(matrix = expression$matrix, labels = labels))
  }
  if (inherits(reference, "SummarizedExperiment")) {
    check_installed("SummarizedExperiment", reason = "to use a SummarizedExperiment annotation reference.")
    metadata <- as.data.frame(SummarizedExperiment::colData(reference))
    if (!reference_label_by %in% colnames(metadata)) {
      stop("`reference_label_by` column '", reference_label_by, "' was not found in reference colData.", call. = FALSE)
    }
    matrix <- SummarizedExperiment::assay(reference)
    labels <- as.character(metadata[[reference_label_by]])
    names(labels) <- colnames(reference)
    return(list(matrix = matrix, labels = labels))
  }
  stop("`reference` must be a Seurat or SummarizedExperiment object for SingleR.", call. = FALSE)
}

.sn_annotation_singleR <- function(object,
                                   reference,
                                   reference_label_by,
                                   assay = NULL,
                                   layer = "data",
                                   reference_assay = NULL,
                                   reference_layer = "data",
                                   backend_control = list()) {
  check_installed("SingleR", reason = "to run `method = \"singleR\"` annotation.")
  if (is_null(reference) || is_null(reference_label_by)) {
    stop("`reference` and `reference_label_by` are required for SingleR annotation.", call. = FALSE)
  }
  query <- .sn_annotation_expression(object, assay = assay, layer = layer)
  ref <- .sn_annotation_reference_input(reference, reference_label_by, assay = reference_assay, layer = reference_layer)
  common <- intersect(rownames(query$matrix), rownames(ref$matrix))
  if (length(common) < 2L) {
    stop("SingleR annotation requires at least two shared query/reference features.", call. = FALSE)
  }
  defaults <- list(
    test = query$matrix[common, , drop = FALSE],
    ref = ref$matrix[common, , drop = FALSE],
    labels = ref$labels
  )
  prediction <- do.call(SingleR::SingleR, utils::modifyList(defaults, backend_control, keep.null = TRUE))
  labels <- as.character(prediction$pruned.labels)
  missing_pruned <- is.na(labels) | !nzchar(labels)
  labels[missing_pruned] <- as.character(prediction$labels[missing_pruned])
  scores <- as.matrix(prediction$scores)
  if (is.null(rownames(scores))) rownames(scores) <- colnames(object)
  scores <- scores[colnames(object), , drop = FALSE]
  best_scores <- apply(scores, 1, function(row) {
    finite <- row[is.finite(row)]
    if (length(finite)) max(finite) else NA_real_
  })
  evidence <- tibble::tibble(
    entity = colnames(object),
    label = labels,
    score = unname(best_scores),
    method = "singleR",
    reference_coverage = length(common) / length(unique(c(rownames(query$matrix), rownames(ref$matrix))))
  )
  delta_next <- as.numeric(prediction$delta.next %||% rep(NA_real_, ncol(object)))
  list(
    object = object,
    evidence = evidence,
    raw_predictions = tibble::tibble(
      cell = colnames(object),
      prediction = labels,
      pruned = !missing_pruned,
      delta_next = delta_next
    ),
    input = list(
      assay = query$assay,
      layer = query$layer,
      shared_features = length(common),
      reference_cells = ncol(ref$matrix)
    )
  )
}

.sn_annotation_transfer <- function(object,
                                    reference,
                                    reference_label_by,
                                    method,
                                    assay = NULL,
                                    layer = "data",
                                    backend_control = list()) {
  if (is_null(reference) || is_null(reference_label_by)) {
    stop("`reference` and `reference_label_by` are required for ", method, " annotation.", call. = FALSE)
  }
  prefix <- paste0("sn_annotation_", method)
  defaults <- list(
    object = object,
    reference = reference,
    label_by = reference_label_by,
    method = method,
    prediction_prefix = prefix,
    query_assay = assay,
    query_layer = layer,
    store_prediction_scores = TRUE,
    verbose = FALSE
  )
  transferred <- do.call(sn_transfer_labels, utils::modifyList(defaults, backend_control, keep.null = TRUE))
  label_col <- paste0(prefix, "_label")
  score_col <- paste0(prefix, "_score")
  labels <- as.character(transferred[[label_col, drop = TRUE]])
  scores <- if (score_col %in% colnames(transferred[[]])) {
    as.numeric(transferred[[score_col, drop = TRUE]])
  } else {
    rep(NA_real_, ncol(transferred))
  }
  evidence <- tibble::tibble(
    entity = colnames(transferred),
    label = labels,
    score = scores,
    method = method,
    reference_coverage = NA_real_
  )
  raw_columns <- grep(paste0("^", prefix, "_"), colnames(transferred[[]]), value = TRUE)
  raw_predictions <- tibble::as_tibble(transferred[[]][, raw_columns, drop = FALSE], rownames = "cell")
  list(
    object = transferred,
    evidence = evidence,
    raw_predictions = raw_predictions,
    input = list(assay = assay, layer = layer)
  )
}

.sn_annotation_celltypist <- function(object, assay = NULL, layer = "counts", backend_control = list()) {
  before <- colnames(object[[]])
  effective_assay <- backend_control$assay %||% assay %||% Seurat::DefaultAssay(object)
  effective_layer <- backend_control$layer %||% layer
  defaults <- list(x = object, assay = effective_assay, layer = effective_layer)
  annotated <- do.call(sn_run_celltypist, utils::modifyList(defaults, backend_control, keep.null = TRUE))
  added <- setdiff(colnames(annotated[[]]), before)
  candidates <- c(grep("majority_voting$", added, value = TRUE), grep("predicted_labels$", added, value = TRUE))
  label_col <- candidates[[1]] %||% NULL
  if (is_null(label_col)) {
    stop("CellTypist annotation did not add a recognized label metadata column.", call. = FALSE)
  }
  labels <- as.character(annotated[[label_col, drop = TRUE]])
  confidence_col <- grep("_confidence$", added, value = TRUE)
  scores <- if (length(confidence_col)) {
    suppressWarnings(as.numeric(annotated[[confidence_col[[1L]], drop = TRUE]]))
  } else {
    rep(NA_real_, ncol(annotated))
  }
  evidence <- tibble::tibble(
    entity = colnames(annotated), label = labels, score = scores,
    method = "celltypist", reference_coverage = NA_real_
  )
  list(
    object = annotated,
    evidence = evidence,
    raw_predictions = tibble::as_tibble(annotated[[]][, added, drop = FALSE], rownames = "cell"),
    input = list(assay = effective_assay, layer = effective_layer)
  )
}

.sn_annotation_symphony <- function(object,
                                    reference,
                                    reference_label_by,
                                    assay = NULL,
                                    layer = "data",
                                    backend_control = list()) {
  if (is_null(reference) || is_null(reference_label_by)) {
    stop("`reference` and `reference_label_by` are required for Symphony annotation.", call. = FALSE)
  }
  check_installed("symphony", reason = "to run `method = \"symphony\"` annotation.")
  query <- .sn_annotation_expression(object, assay = assay, layer = layer)
  build_control <- backend_control$build %||% list()
  if (inherits(reference, "Seurat")) {
    reference_input <- .sn_annotation_reference_input(
      reference,
      reference_label_by,
      assay = backend_control$reference_assay %||% NULL,
      layer = backend_control$reference_layer %||% "data"
    )
    reference_metadata <- reference[[]]
    defaults <- list(
      exp_ref = reference_input$matrix,
      metadata_ref = reference_metadata,
      vars = backend_control$vars %||% NULL,
      K = min(100L, max(2L, ncol(reference) - 1L)),
      verbose = FALSE,
      do_umap = FALSE,
      do_normalize = identical(backend_control$reference_layer %||% "data", "counts")
    )
    reference_object <- do.call(symphony::buildReference, utils::modifyList(defaults, build_control, keep.null = TRUE))
    train_labels <- reference_input$labels
  } else if (is.list(reference) && !is_null(reference$Z_corr) && !is_null(reference$meta_data)) {
    reference_object <- reference
    if (!reference_label_by %in% colnames(reference$meta_data)) {
      stop("`reference_label_by` column '", reference_label_by, "' was not found in the Symphony reference metadata.", call. = FALSE)
    }
    train_labels <- as.character(reference$meta_data[[reference_label_by]])
  } else {
    stop("Symphony `reference` must be a Seurat object or a built Symphony reference.", call. = FALSE)
  }

  map_defaults <- list(
    exp_query = query$matrix,
    metadata_query = object[[]],
    ref_obj = reference_object,
    vars = backend_control$vars %||% NULL,
    verbose = FALSE,
    do_normalize = identical(query$layer, "counts"),
    do_umap = FALSE
  )
  query_object <- do.call(
    symphony::mapQuery,
    utils::modifyList(map_defaults, backend_control$map %||% list(), keep.null = TRUE)
  )
  prefix <- "sn_annotation_symphony"
  knn_defaults <- list(
    query_obj = query_object,
    ref_obj = reference_object,
    train_labels = train_labels,
    k = min(5L, length(train_labels)),
    save_as = prefix,
    confidence = TRUE,
    seed = backend_control$seed %||% 0L
  )
  query_object <- do.call(
    symphony::knnPredict,
    utils::modifyList(knn_defaults, backend_control$knn %||% list(), keep.null = TRUE)
  )
  labels <- as.character(query_object$meta_data[[prefix]])
  scores <- as.numeric(query_object$meta_data[[paste0(prefix, "_prob")]] %||% rep(NA_real_, length(labels)))
  metadata <- data.frame(row.names = colnames(object))
  metadata[[paste0(prefix, "_label")]] <- labels
  metadata[[paste0(prefix, "_score")]] <- scores
  object <- SeuratObject::AddMetaData(object, metadata = metadata)
  list(
    object = object,
    evidence = tibble::tibble(
      entity = colnames(object), label = labels, score = scores,
      method = "symphony", reference_coverage = length(intersect(rownames(query$matrix), reference_object$vargenes$symbol)) /
        length(reference_object$vargenes$symbol)
    ),
    raw_predictions = tibble::tibble(cell = colnames(object), prediction = labels, prediction_score = scores),
    embeddings = list(symphony = t(query_object$Z)),
    models = list(reference_summary = list(reference_cells = ncol(reference_object$Z_corr), dimensions = nrow(reference_object$Z_corr))),
    input = list(assay = query$assay, layer = query$layer, reference_cells = length(train_labels))
  )
}

.sn_matrix_slice <- function(x, axis) {
  if (is.matrix(x) || is.array(x)) {
    if (axis == 2) x[, 1] else x[1, ]
  } else {
    x
  }
}

.sn_annotation_scmap_sce <- function(matrix, labels = NULL, label_column = "cell_type") {  check_installed(c("SingleCellExperiment", "SummarizedExperiment", "S4Vectors"), reason = "to prepare scmap inputs.")
  col_data <- if (is_null(labels)) {
    S4Vectors::DataFrame(row.names = colnames(matrix))
  } else {
    frame <- data.frame(labels = labels, row.names = colnames(matrix), stringsAsFactors = FALSE)
    colnames(frame) <- label_column
    S4Vectors::DataFrame(frame)
  }
  object <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = matrix, normcounts = matrix, logcounts = matrix),
    colData = col_data
  )
  SummarizedExperiment::rowData(object)$feature_symbol <- rownames(object)
  object[!duplicated(rownames(object)), ]
}

.sn_annotation_scmap <- function(object,
                                 reference,
                                 reference_label_by,
                                 assay = NULL,
                                 layer = "data",
                                 backend_control = list()) {
  check_installed("scmap", reason = "to run `method = \"scmap\"` annotation.")
  if (is_null(reference) || is_null(reference_label_by)) {
    stop("`reference` and `reference_label_by` are required for scmap annotation.", call. = FALSE)
  }
  query <- .sn_annotation_expression(object, assay = assay, layer = layer)
  reference_input <- .sn_annotation_reference_input(
    reference,
    reference_label_by,
    assay = backend_control$reference_assay %||% NULL,
    layer = backend_control$reference_layer %||% "data"
  )
  common <- intersect(rownames(query$matrix), rownames(reference_input$matrix))
  if (length(common) < 2L) {
    stop("scmap annotation requires at least two shared query/reference features.", call. = FALSE)
  }
  reference_sce <- .sn_annotation_scmap_sce(
    reference_input$matrix[common, , drop = FALSE],
    labels = reference_input$labels,
    label_column = reference_label_by
  )
  query_sce <- .sn_annotation_scmap_sce(query$matrix[common, , drop = FALSE])
  if (is_null(backend_control$features)) {
    reference_sce <- scmap::selectFeatures(reference_sce, suppress_plot = TRUE)
  } else {
    reference_sce <- scmap::setFeatures(reference_sce, intersect(backend_control$features, common))
  }
  selected <- SummarizedExperiment::rowData(reference_sce)$scmap_features
  if (!is.null(selected) && any(selected)) {
    selected_features <- rownames(reference_sce)[selected]
    reference_sce <- reference_sce[selected_features, ]
    query_sce <- query_sce[selected_features, ]
  }
  reference_sce <- scmap::indexCluster(reference_sce, cluster_col = reference_label_by)
  index <- S4Vectors::metadata(reference_sce)$scmap_cluster_index
  projected <- scmap::scmapCluster(
    projection = query_sce,
    index_list = list(reference = index),
    threshold = backend_control$threshold %||% 0.7
  )

  labels <- rep(NA_character_, ncol(object))
  scores <- rep(NA_real_, ncol(object))
  if (is.list(projected) && !is_null(projected$scmap_cluster_labs)) {
    labels_matrix <- projected$scmap_cluster_labs
    scores_matrix <- projected$scmap_cluster_siml %||% matrix(NA_real_, nrow = nrow(labels_matrix), ncol = ncol(labels_matrix))
    # Current scmap returns a cells x references matrix while older builds
    # returned references x cells; pick the axis that matches the query cells.
    cells_axis <- if (length(dim(labels_matrix)) == 2L && nrow(labels_matrix) == ncol(object)) 2L else 1L
    extracted_labels <- as.character(.sn_matrix_slice(labels_matrix, cells_axis))
    extracted_scores <- as.numeric(.sn_matrix_slice(scores_matrix, cells_axis))
    if (length(extracted_labels) == ncol(object)) {
      labels <- extracted_labels
      if (length(extracted_scores) == ncol(object)) scores <- extracted_scores
    } else {
      cli::cli_warn("scmap did not return one label per query cell; labels are set to NA.")
    }
  } else if (inherits(projected, "SingleCellExperiment")) {
    projected_metadata <- as.data.frame(SummarizedExperiment::colData(projected))
    label_col <- grep("scmap.*lab", colnames(projected_metadata), value = TRUE)[1] %||% "scmap_labels"
    score_col <- grep("scmap.*(sim|prob|score)", colnames(projected_metadata), value = TRUE)[1] %||% NULL
    extracted_labels <- as.character(projected_metadata[[label_col]])
    if (length(extracted_labels) == ncol(object)) {
      labels <- extracted_labels
      extracted_scores <- if (is_null(score_col)) rep(NA_real_, length(labels)) else as.numeric(projected_metadata[[score_col]])
      if (length(extracted_scores) == ncol(object)) scores <- extracted_scores
    } else {
      cli::cli_warn("scmap did not return one label per query cell; labels are set to NA.")
    }
  } else {
    cli::cli_warn("scmap returned an unsupported projection result; labels are set to NA.")
  }
  # scmap reports rejected / unavailable assignments as NA; surface them as an
  # explicit low-confidence label so downstream consensus stays computable.
  unassigned <- is.na(labels)
  labels[unassigned] <- "unassigned"
  scores[unassigned] <- 0
  names(labels) <- colnames(object)
  metadata <- data.frame(row.names = colnames(object))
  metadata$sn_annotation_scmap_label <- labels
  metadata$sn_annotation_scmap_score <- scores
  object <- SeuratObject::AddMetaData(object, metadata = metadata)
  list(
    object = object,
    evidence = tibble::tibble(
      entity = colnames(object), label = unname(labels), score = scores,
      method = "scmap", reference_coverage = length(common) /
        length(unique(c(rownames(query$matrix), rownames(reference_input$matrix))))
    ),
    raw_predictions = tibble::tibble(cell = colnames(object), prediction = unname(labels), prediction_score = scores),
    input = list(assay = query$assay, layer = query$layer, shared_features = length(common), reference_cells = ncol(reference_input$matrix))
  )
}

.sn_annotation_popv <- function(object,
                                reference,
                                reference_label_by,
                                assay = NULL,
                                layer = "counts",
                                backend_control = list()) {
  if (is_null(reference) || is_null(reference_label_by)) {
    stop("`reference` and `reference_label_by` are required for PopV annotation.", call. = FALSE)
  }
  effective_assay <- backend_control$assay %||% assay %||% Seurat::DefaultAssay(object)
  # PopV normalizes internally and rejects non-count matrices, so it always
  # exports raw counts regardless of the wrapper-level `layer` default.
  effective_layer <- backend_control$layer %||% "counts"
  reference_assay <- backend_control$reference_assay %||% Seurat::DefaultAssay(reference)
  reference_layer <- backend_control$reference_layer %||% "counts"
  seed <- as.integer(backend_control$seed %||% 0L)

  output_supplied <- !is.null(backend_control$output_dir)
  keep_run_dir <- backend_control$keep_run_dir
  if (is.null(keep_run_dir)) keep_run_dir <- output_supplied
  if (!is.logical(keep_run_dir) || length(keep_run_dir) != 1L || is.na(keep_run_dir)) {
    stop("`backend_control$keep_run_dir` must be TRUE or FALSE.", call. = FALSE)
  }
  run_dir <- .sn_resolve_python_run_directory(
    path = backend_control$output_dir,
    method = "popv",
    runtime_dir = .sn_shennong_runtime_dir(backend_control$runtime_dir %||% NULL),
    keep_run_dir = keep_run_dir,
    supplied = output_supplied
  )
  run_complete <- FALSE
  on.exit({
    if (!run_complete && !isTRUE(keep_run_dir) && dir.exists(run_dir)) {
      .sn_sanitize_failed_python_run(run_dir, method = "popv", stage = "annotation")
    }
  }, add = TRUE)
  input_dir <- file.path(run_dir, "input")
  result_dir <- file.path(run_dir, "output")
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  query_input <- .sn_write_python_object_input(
    object = object,
    input_dir = file.path(input_dir, "query"),
    assay = effective_assay,
    layer = effective_layer,
    metadata_columns = backend_control$query_batch_key %||% NULL,
    require_raw_counts = TRUE
  )
  reference_input <- .sn_write_python_object_input(
    object = reference,
    input_dir = file.path(input_dir, "reference"),
    assay = reference_assay,
    layer = reference_layer,
    metadata_columns = c(
      reference_label_by,
      backend_control$ref_batch_key %||% character()
    ),
    require_raw_counts = TRUE
  )
  shared_features <- length(intersect(query_input$features, reference_input$features))
  union_features <- length(unique(c(query_input$features, reference_input$features)))

  config <- list(
    method = "popv",
    label_key = reference_label_by,
    query_batch_key = backend_control$query_batch_key %||% NULL,
    ref_batch_key = backend_control$ref_batch_key %||% NULL,
    prediction_mode = backend_control$prediction_mode %||% "retrain",
    unknown_celltype_label = backend_control$unknown_celltype_label %||% "unknown",
    n_samples_per_label = backend_control$n_samples_per_label %||% 300,
    hvg = backend_control$hvg,
    methods = backend_control$methods %||% NULL,
    cl_obo_folder = backend_control$cl_obo_folder %||% FALSE,
    min_shared_genes = backend_control$min_shared_genes %||% 2L,
    seed = seed
  )
  config_path <- .sn_write_json_file(config, file.path(run_dir, "popv_config.json"))
  script <- .sn_pixi_script_path(environment = "popv", script_name = "popv_run.py")
  .sn_execute_python_object_pixi(
    environment = "popv",
    script = script,
    input_dir = input_dir,
    output_dir = result_dir,
    config_path = config_path,
    quiet = isTRUE(backend_control$quiet)
  )

  predictions <- utils::read.csv(
    file.path(result_dir, "predictions.csv"),
    row.names = 1,
    check.names = FALSE
  )
  .sn_validate_exact_python_ids(
    rownames(predictions),
    colnames(object),
    "PopV prediction cell"
  )
  required_predictions <- c("popv_prediction", "popv_majority_vote_score")
  missing_predictions <- setdiff(required_predictions, colnames(predictions))
  if (length(missing_predictions) > 0L) {
    stop(
      "PopV predictions are missing required column(s): ",
      paste(missing_predictions, collapse = ", "), ".",
      call. = FALSE
    )
  }
  aligned <- predictions[colnames(object), , drop = FALSE]
  labels <- as.character(aligned[["popv_prediction"]])
  scores <- suppressWarnings(as.numeric(aligned[["popv_majority_vote_score"]]))
  if (anyNA(labels) || any(!nzchar(labels)) || any(!is.finite(scores)) || any(scores < 0)) {
    stop("PopV predictions must contain non-empty labels and finite non-negative vote scores.", call. = FALSE)
  }
  n_methods_run <- length(setdiff(colnames(predictions), c("popv_prediction", "popv_prediction_score", "popv_majority_vote_prediction", "popv_majority_vote_score", "popv_parent")))
  normalized_scores <- scores / pmax(1, n_methods_run)
  if (any(normalized_scores > 1 + sqrt(.Machine$double.eps))) {
    stop("PopV majority-vote scores exceed the number of reported prediction methods.", call. = FALSE)
  }

  metadata <- data.frame(row.names = colnames(object))
  metadata$sn_annotation_popv_label <- labels
  metadata$sn_annotation_popv_score <- normalized_scores
  object <- SeuratObject::AddMetaData(object, metadata = metadata)

  evidence <- tibble::tibble(
    entity = colnames(object),
    label = ifelse(is.na(labels), "unassigned", labels),
    score = ifelse(is.na(normalized_scores), 0, normalized_scores),
    method = "popv",
    reference_coverage = if (union_features == 0) NA_real_ else shared_features / union_features
  )
  backend_result <- list(
    object = object,
    evidence = evidence,
    raw_predictions = tibble::tibble(
      cell = colnames(object),
      prediction = labels,
      prediction_score = normalized_scores,
      agreement = scores
    ),
    models = list(manifest = tryCatch(
      jsonlite::fromJSON(file.path(result_dir, "manifest.json"), simplifyVector = TRUE),
      error = function(e) list()
    )),
    input = list(
      assay = effective_assay,
      layer = effective_layer,
      prediction_mode = config$prediction_mode,
      dropped_cells = 0L,
      shared_features = shared_features,
      run_dir = if (isTRUE(keep_run_dir)) normalizePath(run_dir, winslash = "/", mustWork = TRUE) else NULL,
      run_dir_retained = isTRUE(keep_run_dir)
    )
  )
  if (!isTRUE(keep_run_dir)) {
    .sn_remove_python_run_directory(run_dir, label = "PopV annotation run directory")
  }
  run_complete <- TRUE
  backend_result
}

.sn_annotation_backend <- function(object,
                                   method,
                                   reference = NULL,
                                   reference_label_by = NULL,
                                   assay = NULL,
                                   layer = "data",
                                   backend_control = list()) {
  switch(
    method,
    singleR = .sn_annotation_singleR(
      object, reference, reference_label_by, assay = assay, layer = layer,
      reference_assay = backend_control$reference_assay %||% NULL,
      reference_layer = backend_control$reference_layer %||% "data",
      backend_control = backend_control$singleR %||% list()
    ),
    seurat = .sn_annotation_transfer(
      object, reference, reference_label_by, method = "seurat", assay = assay,
      layer = layer, backend_control = backend_control$seurat %||% list()
    ),
    scanvi = .sn_annotation_transfer(
      object, reference, reference_label_by, method = "scanvi", assay = assay,
      layer = layer, backend_control = backend_control$scanvi %||% list()
    ),
    celltypist = .sn_annotation_celltypist(
      object, assay = assay,
      layer = backend_control$celltypist$layer %||% "counts",
      backend_control = backend_control$celltypist %||% list()
    ),
    symphony = .sn_annotation_symphony(
      object, reference, reference_label_by, assay = assay, layer = layer,
      backend_control = backend_control$symphony %||% list()
    ),
    scmap = .sn_annotation_scmap(
      object, reference, reference_label_by, assay = assay, layer = layer,
      backend_control = backend_control$scmap %||% list()
    ),
    popv = .sn_annotation_popv(
      object, reference, reference_label_by, assay = assay, layer = layer,
      backend_control = backend_control$popv %||% list()
    ),
    stop("Annotation backend '", method, "' is not implemented.", call. = FALSE)
  )
}

.sn_annotation_cell_table <- function(cell_evidence, cells) {
  if (!is.data.frame(cell_evidence) ||
      !all(c("entity", "label", "score", "method") %in% colnames(cell_evidence))) {
    stop("The annotation backend returned an invalid per-cell evidence table.", call. = FALSE)
  }
  evidence <- tibble::tibble(cell = as.character(cell_evidence$entity)) |>
    dplyr::mutate(
      prediction = as.character(cell_evidence$label),
      prediction_score = suppressWarnings(as.numeric(cell_evidence$score)),
      method = as.character(cell_evidence$method),
      reference_coverage = if ("reference_coverage" %in% colnames(cell_evidence)) {
        suppressWarnings(as.numeric(cell_evidence$reference_coverage))
      } else {
        NA_real_
      }
    )
  duplicated <- duplicated(evidence[["cell"]]) | duplicated(rev(evidence[["cell"]]))
  if (any(duplicated)) {
    stop("The annotation backend returned more than one prediction for a cell.", call. = FALSE)
  }
  missing <- setdiff(cells, evidence[["cell"]])
  if (length(missing) > 0L) {
    filler <- tibble::tibble(
      cell = missing,
      prediction = NA_character_,
      prediction_score = NA_real_,
      method = unique(evidence[["method"]])[[1]] %||% NA_character_,
      reference_coverage = NA_real_
    )
    evidence <- dplyr::bind_rows(evidence, filler)
  }
  unassigned <- is.na(evidence[["prediction"]]) |
    !nzchar(trimws(evidence[["prediction"]]))
  evidence[["prediction"]][unassigned] <- "unassigned"
  evidence[match(cells, evidence[["cell"]]), , drop = FALSE]
}

.sn_annotation_cluster_table <- function(cells_table, clusters,
                                         confidence_threshold = NULL) {
  table <- dplyr::mutate(
    cells_table,
    cluster = unname(clusters[.data$cell]),
    .after = "cell"
  )
  groups <- split(table, table$cluster)
  dplyr::bind_rows(lapply(groups, function(rows) {
    counts <- table(rows$prediction, useNA = "ifany")
    modal <- names(counts)[which.max(counts)][[1]]
    members <- !is.na(rows$prediction) & rows$prediction == modal
    scores <- rows$prediction_score[members]
    agreement <- sum(members) / max(1L, nrow(rows))
    finite_scores <- scores[is.finite(scores)]
    mean_score <- if (length(finite_scores)) mean(finite_scores) else NA_real_
    low_confidence_share <- if (any(members)) {
      mean(rows$low_confidence[members], na.rm = TRUE)
    } else {
      1
    }
    score_is_low <- !is.finite(mean_score) ||
      (!is_null(confidence_threshold) && mean_score < confidence_threshold) ||
      low_confidence_share >= 0.5
    tibble::tibble(
      cluster = rows$cluster[[1]],
      prediction = if (is.na(modal)) NA_character_ else modal,
      prediction_score = mean_score,
      agreement_share = agreement,
      low_confidence_share = low_confidence_share,
      low_confidence = is.na(modal) || agreement < 0.5 || score_is_low,
      method = rows$method[[1]],
      reference_coverage = if (all(is.na(rows$reference_coverage))) {
        NA_real_
      } else {
        max(rows$reference_coverage, na.rm = TRUE)
      },
      n_cells = nrow(rows)
    )
  }))
}

#' Run reference-based cell-type annotation
#'
#' Stable annotation entry point for SingleR, CellTypist, Seurat label
#' transfer, Symphony mapping, scmap projection, scANVI transfer, and PopV
#' consensus voting. Computational labels and raw backend predictions are
#' retained; no LLM is allowed to overwrite them. Cluster summaries report the
#' modal predicted label per \code{group_by} group.
#'
#' @param object A \code{Seurat} object.
#' @param group_by Metadata column used for the cluster-level summary.
#' @param method Annotation method. One of \code{"singleR"} (default),
#'   \code{"celltypist"}, \code{"seurat"}, \code{"symphony"},
#'   \code{"scmap"}, \code{"scanvi"}, or \code{"popv"}.
#' @param reference Annotated reference object required by all methods except
#'   CellTypist, which uses a pre-trained model instead.
#' @param reference_label_by Reference label metadata/colData column.
#' @param tissue,disease Optional biological context recorded in provenance.
#' @param species \code{"human"} or \code{"mouse"}; inferred when possible.
#' @param ontology Map labels to the bundled Cell Ontology snapshot.
#' @param result_id Stored-result and metadata prefix.
#' @param confidence_threshold Optional minimum backend score. Scores are
#'   backend-specific and are not assumed to be calibrated across methods. By
#'   default, only missing/non-finite/non-positive scores and explicit
#'   unassigned labels are flagged from score evidence.
#' @param assay,layer Query expression source. Most backends read a
#'   log-normalized `data` layer; CellTypist and PopV independently default to
#'   raw `counts` because they normalize internally. Override those backends
#'   only via `backend_control = list(celltypist = list(layer = ...))` or
#'   `backend_control = list(popv = list(layer = ...))`.
#' @param backend_control Named backend-specific control lists. PopV exports
#'   only its required label/batch metadata and verified raw counts. Its
#'   package-owned temporary run is removed after successful import; supplying
#'   `popv$output_dir` retains that empty, explicit run location unless
#'   `popv$keep_run_dir = FALSE`, in which case it is treated as a parent and
#'   only Shennong's unique child is cleaned or sanitized.
#' @param return_object If \code{TRUE}, return the annotated object; otherwise
#'   return the stored result.
#'
#' @return An annotated Seurat object or a unified annotation result.
#'
#' @examples
#' \dontrun{
#' object <- sn_run_annotation(
#'   object,
#'   group_by = "seurat_clusters",
#'   method = "singleR",
#'   reference = reference,
#'   reference_label_by = "cell_type",
#'   species = "human"
#' )
#' sn_get_result(object, "annotation", "annotation")
#' }
#'
#' @export
sn_run_annotation <- function(object,
                              group_by = "seurat_clusters",
                              method = c("singleR", "celltypist", "seurat", "symphony", "scmap", "scanvi", "popv"),
                              reference = NULL,
                              reference_label_by = NULL,
                              tissue = NULL,
                              disease = NULL,
                              species = NULL,
                              ontology = TRUE,
                              result_id = "annotation",
                              confidence_threshold = NULL,
                              assay = NULL,
                              layer = "data",
                              backend_control = list(),
                              return_object = TRUE) {
  result_id <- .sn_validate_result_id(result_id)
  .sn_validate_seurat_object(object)
  method <- match.arg(method)
  if (!group_by %in% colnames(object[[]])) {
    stop("`group_by` column '", group_by, "' was not found in object metadata.", call. = FALSE)
  }
  if (!is.list(backend_control)) {
    stop("`backend_control` must be a named list.", call. = FALSE)
  }
  if (!is_null(confidence_threshold) &&
      (!is.numeric(confidence_threshold) || length(confidence_threshold) != 1L ||
       !is.finite(confidence_threshold))) {
    stop("`confidence_threshold` must be NULL or one finite numeric value.", call. = FALSE)
  }
  clusters <- as.character(object[[group_by, drop = TRUE]])
  names(clusters) <- colnames(object)
  if (anyNA(clusters) || any(!nzchar(clusters))) {
    stop("`group_by` contains missing or empty values.", call. = FALSE)
  }

  backend <- .sn_annotation_backend(
    object,
    method = method,
    reference = reference,
    reference_label_by = reference_label_by,
    assay = assay,
    layer = layer,
    backend_control = backend_control
  )
  object <- backend$object

  cells_table <- .sn_annotation_cell_table(backend$evidence, colnames(object))
  level_1 <- .sn_annotation_label_hierarchy(cells_table$prediction)
  cells_table$level_1 <- unname(level_1[cells_table$prediction])
  cells_table$level_2 <- cells_table$prediction
  cells_table$level_3 <- cells_table$prediction
  if (isTRUE(ontology)) {
    mapped <- sn_map_cell_ontology(cells_table$prediction)
    cells_table$ontology_id <- mapped$ontology_id
    cells_table$ontology_label <- mapped$ontology_label
  } else {
    cells_table$ontology_id <- NA_character_
    cells_table$ontology_label <- NA_character_
  }
  cells_table$low_confidence <- !is.finite(cells_table$prediction_score) |
    cells_table$prediction_score <= 0 |
    (!is_null(confidence_threshold) &
       cells_table$prediction_score < (confidence_threshold %||% -Inf)) |
    cells_table$prediction %in% c("unassigned", "unknown")

  clusters_table <- .sn_annotation_cluster_table(
    cells_table,
    clusters,
    confidence_threshold = confidence_threshold
  )

  safe_result_id <- gsub("[^[:alnum:]_]+", "_", result_id)
  cell_indices <- match(colnames(object), cells_table$cell)
  metadata <- data.frame(row.names = colnames(object))
  metadata[[paste0(safe_result_id, "_label")]] <- cells_table$prediction[cell_indices]
  metadata[[paste0(safe_result_id, "_level_1")]] <- cells_table$level_1[cell_indices]
  metadata[[paste0(safe_result_id, "_level_2")]] <- cells_table$level_2[cell_indices]
  metadata[[paste0(safe_result_id, "_level_3")]] <- cells_table$level_3[cell_indices]
  metadata[[paste0(safe_result_id, "_score")]] <- cells_table$prediction_score[cell_indices]
  metadata[[paste0(safe_result_id, "_low_confidence")]] <- cells_table$low_confidence[cell_indices]
  metadata[[paste0(safe_result_id, "_ontology_id")]] <- cells_table$ontology_id[cell_indices]
  object <- SeuratObject::AddMetaData(object, metadata = metadata)

  result <- list(
    schema_version = .sn_analysis_result_schema_version(),
    analysis_type = "annotation",
    result_id = result_id,
    method = method,
    backend = method,
    input = c(
      list(
        assay = backend$input$assay %||% assay,
        layer = backend$input$layer %||% layer,
        cells = ncol(object),
        features = nrow(object),
        group_by = group_by,
        species = sn_get_species(object, species = species),
        tissue = tissue,
        disease = disease,
        reference_label_by = reference_label_by
      ),
      backend$input[names(backend$input) %in% c(
        "shared_features", "reference_cells", "prediction_mode", "dropped_cells",
        "run_dir", "run_dir_retained"
      )]
    ),
    parameters = list(
      ontology = ontology,
      confidence_threshold = confidence_threshold,
      confidence_scale = "backend_specific"
    ),
    tables = list(
      primary = tibble::as_tibble(cells_table),
      cells = tibble::as_tibble(cells_table),
      clusters = tibble::as_tibble(clusters_table),
      evidence = tibble::as_tibble(cells_table),
      backend_predictions = tibble::as_tibble(backend$raw_predictions %||% tibble::tibble())
    ),
    embeddings = backend$embeddings %||% list(),
    graphs = list(),
    models = backend$models %||% list(),
    diagnostics = list(
      low_confidence_cells = sum(cells_table$low_confidence, na.rm = TRUE),
      low_confidence_clusters = sum(clusters_table$low_confidence, na.rm = TRUE),
      unmapped_ontology_labels = unique(cells_table$prediction[is.na(cells_table$ontology_id)])
    ),
    warnings = character(),
    provenance = .sn_analysis_provenance(random_seed = backend_control$seed %||% NA_integer_)
  )
  sn_validate_result(result)
  object <- sn_store_result(object, "annotation", result_id, result)
  object <- .sn_log_seurat_command(
    object = object,
    assay = backend$input$assay %||% assay,
    name = "sn_run_annotation"
  )
  if (isTRUE(return_object)) object else sn_get_result(object, "annotation", result_id)
}

#' Review stored annotation evidence and low-confidence labels
#'
#' @param x A Seurat object or unified annotation result.
#' @param result_id Annotation result name when \code{x} is a Seurat object.
#' @param low_confidence_only Return only low-confidence rows in the review
#'   tables.
#'
#' @return A list containing cell, cluster, evidence, and diagnostic tables.
#'
#' @examples
#' \dontrun{sn_review_annotation(object, "annotation")}
#'
#' @export
sn_review_annotation <- function(x, result_id = "annotation", low_confidence_only = TRUE) {
  result <- if (inherits(x, "Seurat")) sn_get_result(x, "annotation", result_id) else x
  sn_validate_result(result)
  if (!identical(result$analysis_type, "annotation")) {
    stop("`x` is not an annotation result.", call. = FALSE)
  }
  cells <- result$tables$cells %||% tibble::tibble()
  clusters <- result$tables$clusters %||% tibble::tibble()
  if (isTRUE(low_confidence_only)) {
    if ("low_confidence" %in% colnames(cells)) cells <- cells[cells$low_confidence, , drop = FALSE]
    if ("low_confidence" %in% colnames(clusters)) clusters <- clusters[clusters$low_confidence, , drop = FALSE]
  }
  list(
    cells = tibble::as_tibble(cells),
    clusters = tibble::as_tibble(clusters),
    evidence = tibble::as_tibble(result$tables$evidence %||% tibble::tibble()),
    diagnostics = result$diagnostics,
    provenance = result$provenance
  )
}
