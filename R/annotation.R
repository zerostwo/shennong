.sn_annotation_marker_database <- function(marker_database = NULL, species = "human") {
  if (is_null(marker_database)) {
    utils::data("marker_genes", package = "Shennong", envir = environment())
    marker_database <- get("marker_genes", envir = environment())
  }
  if (!is.data.frame(marker_database)) {
    stop("`marker_database` must be a data frame.", call. = FALSE)
  }
  required <- c("high_hierarchy_cell_type", "low_hierarchy_cell_type", species)
  missing <- setdiff(required, colnames(marker_database))
  if (length(missing) > 0L) {
    stop("`marker_database` is missing column(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }
  out <- tibble::tibble(
    parent_label = as.character(marker_database$high_hierarchy_cell_type),
    label = as.character(marker_database$low_hierarchy_cell_type),
    gene = as.character(marker_database[[species]]),
    direction = if ("direction" %in% colnames(marker_database)) {
      tolower(as.character(marker_database$direction))
    } else {
      "positive"
    }
  )
  out$direction[!out$direction %in% c("positive", "negative")] <- "positive"
  unique(out[!is.na(out$gene) & nzchar(out$gene) & !is.na(out$label) & nzchar(out$label), ])
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

.sn_annotation_marker_scores <- function(object,
                                         group_by,
                                         marker_database = NULL,
                                         species = NULL,
                                         assay = NULL,
                                         layer = "data") {
  species <- sn_get_species(object, species = species)
  markers <- .sn_annotation_marker_database(marker_database, species = species)
  expression <- .sn_annotation_expression(object, assay = assay, layer = layer)
  matrix <- expression$matrix
  groups <- as.character(object[[group_by, drop = TRUE]])
  names(groups) <- colnames(object)
  if (anyNA(groups) || any(!nzchar(groups))) {
    stop("`group_by` contains missing or empty values.", call. = FALSE)
  }

  feature_lookup <- stats::setNames(rownames(matrix), toupper(sub("\\.[0-9]+$", "", rownames(matrix))))
  markers$feature <- unname(feature_lookup[toupper(sub("\\.[0-9]+$", "", markers$gene))])
  marker_labels <- split(markers, markers$label)
  group_levels <- unique(groups)
  rows <- list()
  counter <- 0L
  for (group in group_levels) {
    cells <- names(groups)[groups == group]
    averages <- Matrix::rowMeans(matrix[, cells, drop = FALSE])
    for (label in names(marker_labels)) {
      current <- marker_labels[[label]]
      positive <- current[current$direction == "positive", , drop = FALSE]
      negative <- current[current$direction == "negative", , drop = FALSE]
      positive_features <- unique(stats::na.omit(positive$feature))
      negative_features <- unique(stats::na.omit(negative$feature))
      positive_values <- averages[positive_features]
      negative_values <- averages[negative_features]
      positive_score <- if (length(positive_values) == 0L) 0 else mean(positive_values)
      negative_score <- if (length(negative_values) == 0L) 0 else mean(negative_values)
      score <- positive_score - negative_score
      supporting <- names(sort(positive_values, decreasing = TRUE))[seq_len(min(5L, length(positive_values)))]
      conflicting <- names(sort(negative_values, decreasing = TRUE))[negative_values[order(negative_values, decreasing = TRUE)] > 0]
      conflicting <- utils::head(conflicting, 5L)
      counter <- counter + 1L
      rows[[counter]] <- tibble::tibble(
        entity = group,
        label = label,
        parent_label = current$parent_label[[1]],
        score = score,
        method = "markers",
        supporting_markers = paste(supporting, collapse = ";"),
        conflicting_markers = paste(conflicting, collapse = ";"),
        reference_coverage = if (nrow(positive) == 0L) 0 else length(positive_features) / nrow(positive)
      )
    }
  }
  list(
    evidence = dplyr::bind_rows(rows),
    assay = expression$assay,
    layer = expression$layer,
    species = species,
    marker_database = markers
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
  evidence <- dplyr::bind_rows(lapply(seq_len(ncol(scores)), function(index) {
    tibble::tibble(
      entity = rownames(scores),
      label = colnames(scores)[[index]],
      score = scores[, index],
      method = "singleR",
      reference_coverage = length(common) / length(unique(c(rownames(query$matrix), rownames(ref$matrix))))
    )
  }))
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
    rep(1, ncol(transferred))
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

.sn_annotation_celltypist <- function(object, assay = NULL, layer = "data", backend_control = list()) {
  before <- colnames(object[[]])
  defaults <- list(x = object, assay = assay %||% Seurat::DefaultAssay(object), layer = layer)
  annotated <- do.call(sn_run_celltypist, utils::modifyList(defaults, backend_control, keep.null = TRUE))
  added <- setdiff(colnames(annotated[[]]), before)
  candidates <- c(grep("majority_voting$", added, value = TRUE), grep("predicted_labels$", added, value = TRUE))
  label_col <- candidates[[1]] %||% NULL
  if (is_null(label_col)) {
    stop("CellTypist annotation did not add a recognized label metadata column.", call. = FALSE)
  }
  labels <- as.character(annotated[[label_col, drop = TRUE]])
  evidence <- tibble::tibble(
    entity = colnames(annotated), label = labels, score = 1,
    method = "celltypist", reference_coverage = NA_real_
  )
  list(
    object = annotated,
    evidence = evidence,
    raw_predictions = tibble::as_tibble(annotated[[]][, added, drop = FALSE], rownames = "cell"),
    input = list(assay = assay, layer = layer)
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
      entity = colnames(object), label = labels, score = ifelse(is.na(scores), 1, scores),
      method = "symphony", reference_coverage = length(intersect(rownames(query$matrix), reference_object$vargenes$symbol)) /
        length(reference_object$vargenes$symbol)
    ),
    raw_predictions = tibble::tibble(cell = colnames(object), prediction = labels, prediction_score = scores),
    embeddings = list(symphony = t(query_object$Z)),
    models = list(reference_summary = list(reference_cells = ncol(reference_object$Z_corr), dimensions = nrow(reference_object$Z_corr))),
    input = list(assay = query$assay, layer = query$layer, reference_cells = length(train_labels))
  )
}

.sn_annotation_scmap_sce <- function(matrix, labels = NULL, label_column = "cell_type") {
  check_installed(c("SingleCellExperiment", "SummarizedExperiment", "S4Vectors"), reason = "to prepare scmap inputs.")
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

  if (is.list(projected) && !is.null(projected$scmap_cluster_labs)) {
    labels_matrix <- projected$scmap_cluster_labs
    scores_matrix <- projected$scmap_cluster_siml %||% matrix(NA_real_, nrow = nrow(labels_matrix), ncol = ncol(labels_matrix))
    labels <- as.character(labels_matrix[1, ])
    scores <- as.numeric(scores_matrix[1, ])
  } else if (inherits(projected, "SingleCellExperiment")) {
    projected_metadata <- as.data.frame(SummarizedExperiment::colData(projected))
    label_col <- grep("scmap.*lab", colnames(projected_metadata), value = TRUE)[1] %||% "scmap_labels"
    score_col <- grep("scmap.*(sim|prob|score)", colnames(projected_metadata), value = TRUE)[1] %||% NULL
    labels <- as.character(projected_metadata[[label_col]])
    scores <- if (is_null(score_col)) rep(NA_real_, length(labels)) else as.numeric(projected_metadata[[score_col]])
  } else {
    stop("scmap returned an unsupported projection result.", call. = FALSE)
  }
  names(labels) <- colnames(object)
  metadata <- data.frame(row.names = colnames(object))
  metadata$sn_annotation_scmap_label <- labels
  metadata$sn_annotation_scmap_score <- scores
  object <- SeuratObject::AddMetaData(object, metadata = metadata)
  list(
    object = object,
    evidence = tibble::tibble(
      entity = colnames(object), label = unname(labels), score = ifelse(is.na(scores), 1, scores),
      method = "scmap", reference_coverage = length(common) /
        length(unique(c(rownames(query$matrix), rownames(reference_input$matrix))))
    ),
    raw_predictions = tibble::tibble(cell = colnames(object), prediction = unname(labels), prediction_score = scores),
    input = list(assay = query$assay, layer = query$layer, shared_features = length(common), reference_cells = ncol(reference_input$matrix))
  )
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
      object, assay = assay, layer = layer,
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
    stop("Annotation backend '", method, "' is not implemented.", call. = FALSE)
  )
}

.sn_annotation_cluster_evidence <- function(cell_evidence, clusters) {
  cells <- names(clusters)
  current <- cell_evidence[cell_evidence$entity %in% cells, , drop = FALSE]
  current$cluster <- unname(clusters[current$entity])
  cluster_sizes <- table(clusters)
  groups <- split(seq_len(nrow(current)), paste(current$cluster, current$label, current$method, sep = "\r"))
  dplyr::bind_rows(lapply(groups, function(indices) {
    rows <- current[indices, , drop = FALSE]
    cluster <- rows$cluster[[1]]
    tibble::tibble(
      entity = cluster,
      label = rows$label[[1]],
      score = sum(rows$score, na.rm = TRUE) / as.numeric(cluster_sizes[[cluster]]),
      method = rows$method[[1]],
      reference_coverage = if (all(is.na(rows$reference_coverage))) NA_real_ else max(rows$reference_coverage, na.rm = TRUE)
    )
  }))
}

.sn_annotation_attach_hierarchy <- function(table, marker_evidence) {
  hierarchy <- unique(marker_evidence[c("label", "parent_label")])
  parent <- stats::setNames(hierarchy$parent_label, hierarchy$label)
  missing_level <- is.na(table$level_1) | !nzchar(table$level_1)
  table$level_1[missing_level] <- unname(parent[table$prediction[missing_level]])
  table$level_1[is.na(table$level_1) | !nzchar(table$level_1)] <- table$prediction[is.na(table$level_1) | !nzchar(table$level_1)]
  table
}

#' Run traceable cell-type annotation
#'
#' Provides a stable annotation entry point for marker-only consensus,
#' SingleR, CellTypist, Seurat label transfer, and scANVI label transfer.
#' Consensus mode always evaluates canonical marker evidence and optionally
#' combines it with a reference backend. Computational labels and raw backend
#' predictions are retained; no LLM is allowed to overwrite them.
#'
#' @param object A \code{Seurat} object.
#' @param group_by Metadata column used for cluster-level annotation.
#' @param method Annotation method.
#' @param reference Optional annotated reference object.
#' @param reference_label_by Reference label metadata/colData column.
#' @param tissue,disease Optional biological context recorded in provenance.
#' @param species \code{"human"} or \code{"mouse"}; inferred when possible.
#' @param ontology Map labels to the bundled Cell Ontology snapshot.
#' @param store_name Stored-result and metadata prefix.
#' @param assay,layer Query expression source.
#' @param marker_database Optional marker data frame with high- and
#'   low-hierarchy labels plus species gene columns.
#' @param consensus_reference_method Backend used when consensus receives a
#'   reference.
#' @param backend_control Named backend-specific control lists.
#' @param low_confidence_threshold,margin_threshold Confidence thresholds.
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
#'   method = "consensus",
#'   species = "human"
#' )
#' sn_get_result(object, "annotation", "annotation")
#' }
#'
#' @export
sn_run_annotation <- function(object,
                              group_by = "seurat_clusters",
                              method = c("consensus", "singleR", "celltypist", "seurat", "symphony", "scmap", "scanvi"),
                              reference = NULL,
                              reference_label_by = NULL,
                              tissue = NULL,
                              disease = NULL,
                              species = NULL,
                              ontology = TRUE,
                              store_name = "annotation",
                              assay = NULL,
                              layer = "data",
                              marker_database = NULL,
                              consensus_reference_method = "singleR",
                              backend_control = list(),
                              low_confidence_threshold = 0.55,
                              margin_threshold = 0.1,
                              return_object = TRUE) {
  .sn_validate_result_object(object)
  method <- match.arg(method)
  if (!group_by %in% colnames(object[[]])) {
    stop("`group_by` column '", group_by, "' was not found in object metadata.", call. = FALSE)
  }
  if (!is.list(backend_control)) {
    stop("`backend_control` must be a named list.", call. = FALSE)
  }
  clusters <- as.character(object[[group_by, drop = TRUE]])
  names(clusters) <- colnames(object)
  marker <- .sn_annotation_marker_scores(
    object, group_by = group_by, marker_database = marker_database,
    species = species, assay = assay, layer = layer
  )

  backend <- NULL
  backend_method <- method
  if (method == "consensus" && !is_null(reference)) {
    backend_method <- match.arg(consensus_reference_method, c("singleR", "seurat", "symphony", "scmap", "scanvi"))
  }
  if (method != "consensus" || !is_null(reference)) {
    backend <- .sn_annotation_backend(
      object,
      method = backend_method,
      reference = reference,
      reference_label_by = reference_label_by,
      assay = assay,
      layer = layer,
      backend_control = backend_control
    )
    object <- backend$object
  }

  cluster_evidence <- marker$evidence
  cell_evidence <- tibble::tibble()
  if (!is_null(backend)) {
    cell_evidence <- backend$evidence
    reference_cluster <- .sn_annotation_cluster_evidence(cell_evidence, clusters)
    if (method == "consensus") {
      cluster_evidence <- dplyr::bind_rows(cluster_evidence, reference_cluster)
    } else {
      cluster_evidence <- reference_cluster
    }
  }
  cluster_predictions <- sn_annotation_consensus(
    cluster_evidence,
    ontology = ontology,
    low_confidence_threshold = low_confidence_threshold,
    margin_threshold = margin_threshold
  ) |>
    .sn_annotation_attach_hierarchy(marker$evidence)
  colnames(cluster_predictions)[colnames(cluster_predictions) == "entity"] <- "cluster"

  if (nrow(cell_evidence) > 0L) {
    cell_predictions <- sn_annotation_consensus(
      cell_evidence,
      ontology = ontology,
      low_confidence_threshold = low_confidence_threshold,
      margin_threshold = margin_threshold
    ) |>
      .sn_annotation_attach_hierarchy(marker$evidence)
    colnames(cell_predictions)[colnames(cell_predictions) == "entity"] <- "cell"
    cell_predictions$cluster <- unname(clusters[cell_predictions$cell])
  } else {
    indices <- match(clusters, cluster_predictions$cluster)
    cell_predictions <- cluster_predictions[indices, , drop = FALSE]
    cell_predictions$cell <- names(clusters)
    cell_predictions <- cell_predictions[c("cell", setdiff(colnames(cell_predictions), "cell"))]
  }

  metadata <- data.frame(row.names = colnames(object))
  cell_indices <- match(colnames(object), cell_predictions$cell)
  safe_store_name <- gsub("[^[:alnum:]_]+", "_", store_name)
  metadata[[paste0(safe_store_name, "_label")]] <- cell_predictions$prediction[cell_indices]
  metadata[[paste0(safe_store_name, "_level_1")]] <- cell_predictions$level_1[cell_indices]
  metadata[[paste0(safe_store_name, "_level_2")]] <- cell_predictions$level_2[cell_indices]
  metadata[[paste0(safe_store_name, "_level_3")]] <- cell_predictions$level_3[cell_indices]
  metadata[[paste0(safe_store_name, "_score")]] <- cell_predictions$prediction_score[cell_indices]
  metadata[[paste0(safe_store_name, "_margin")]] <- cell_predictions$margin[cell_indices]
  metadata[[paste0(safe_store_name, "_low_confidence")]] <- cell_predictions$low_confidence[cell_indices]
  metadata[[paste0(safe_store_name, "_ontology_id")]] <- cell_predictions$ontology_id[cell_indices]
  object <- SeuratObject::AddMetaData(object, metadata = metadata)

  result <- list(
    schema_version = "1.0",
    analysis_type = "annotation",
    name = store_name,
    method = method,
    backend = if (method == "consensus") paste(unique(cluster_evidence$method), collapse = "+") else backend_method,
    input = c(
      list(
        assay = marker$assay,
        layer = marker$layer,
        cells = ncol(object),
        features = nrow(object),
        group_by = group_by,
        species = marker$species,
        tissue = tissue,
        disease = disease,
        reference_label_by = reference_label_by
      ),
      backend$input %||% list()
    ),
    parameters = list(
      ontology = ontology,
      low_confidence_threshold = low_confidence_threshold,
      margin_threshold = margin_threshold,
      consensus_reference_method = if (method == "consensus" && !is_null(reference)) backend_method else NULL
    ),
    tables = list(
      cells = tibble::as_tibble(cell_predictions),
      clusters = tibble::as_tibble(cluster_predictions),
      evidence = tibble::as_tibble(cluster_evidence),
      backend_evidence = tibble::as_tibble(cell_evidence),
      backend_predictions = backend$raw_predictions %||% tibble::tibble(),
      marker_database = marker$marker_database
    ),
    embeddings = backend$embeddings %||% list(),
    graphs = list(),
    models = backend$models %||% list(),
    diagnostics = list(
      low_confidence_cells = sum(cell_predictions$low_confidence, na.rm = TRUE),
      low_confidence_clusters = sum(cluster_predictions$low_confidence, na.rm = TRUE),
      unmapped_ontology_labels = unique(cluster_predictions$prediction[is.na(cluster_predictions$ontology_id)])
    ),
    warnings = character(),
    provenance = .sn_analysis_provenance(random_seed = backend_control$seed %||% NA_integer_)
  )
  sn_validate_result(result)
  object <- sn_store_result(object, "annotation", store_name, result)
  object <- .sn_log_seurat_command(object = object, assay = marker$assay, name = "sn_run_annotation")
  if (isTRUE(return_object)) object else sn_get_result(object, "annotation", store_name)
}

#' Review stored annotation evidence and low-confidence labels
#'
#' @param x A Seurat object or unified annotation result.
#' @param store_name Annotation result name when \code{x} is a Seurat object.
#' @param low_confidence_only Return only low-confidence rows in the review
#'   tables.
#'
#' @return A list containing cell, cluster, evidence, and diagnostic tables.
#'
#' @examples
#' \dontrun{sn_review_annotation(object, "annotation")}
#'
#' @export
sn_review_annotation <- function(x, store_name = "annotation", low_confidence_only = TRUE) {
  result <- if (inherits(x, "Seurat")) sn_get_result(x, "annotation", store_name) else x
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
