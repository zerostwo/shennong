.sn_find_transfer_anchors_backend <- function(...) {
  Seurat::FindTransferAnchors(...)
}

.sn_transfer_data_backend <- function(...) {
  Seurat::TransferData(...)
}

.sn_prepare_label_transfer_name <- function(prefix, cells) {
  paste0(prefix, "_", .sn_metadata_suffix(cells))
}

.sn_transfer_labels_scanvi <- function(object,
                                       reference,
                                       label_by,
                                       prediction_prefix,
                                       assay = NULL,
                                       batch_by = NULL,
                                       features = NULL,
                                       transfer_control = list(),
                                       return_anchors = FALSE,
                                       verbose = TRUE,
                                       method_name = "scanvi") {
  if (!inherits(reference, "Seurat")) {
    stop("`reference` must be a Seurat object for scANVI/scArches label_by transfer.", call. = FALSE)
  }
  if (!label_by %in% colnames(reference[[]])) {
    stop(glue("`label_by` column '{label_by}' was not found in `reference` metadata."), call. = FALSE)
  }
  if (!is.list(transfer_control)) {
    stop("`transfer_control` must be a named list.", call. = FALSE)
  }

  assay <- assay %||% SeuratObject::DefaultAssay(object = object)
  reference_assay <- transfer_control$reference_assay %||% assay
  query_assay <- transfer_control$query_assay %||% assay
  if (!reference_assay %in% names(reference@assays)) {
    stop(glue("Reference assay '{reference_assay}' was not found."), call. = FALSE)
  }
  if (!query_assay %in% names(object@assays)) {
    stop(glue("Query assay '{query_assay}' was not found."), call. = FALSE)
  }

  reference_prefix <- transfer_control$reference_prefix %||% "reference"
  query_prefix <- transfer_control$query_prefix %||% "query"
  reference_cells <- .sn_prepare_label_transfer_name(reference_prefix, colnames(reference))
  query_cells <- .sn_prepare_label_transfer_name(query_prefix, colnames(object))

  reference <- Seurat::RenameCells(reference, new.names = reference_cells)
  object_for_transfer <- Seurat::RenameCells(object, new.names = query_cells)
  reference$.sn_transfer_role <- "reference"
  object_for_transfer$.sn_transfer_role <- "query"
  transfer_label_by <- transfer_control$label_by %||% transfer_control$label_by %||% ".sn_transfer_label"
  unlabeled_category <- transfer_control$unlabeled_category %||% "Unknown"
  reference[[transfer_label_by]] <- as.character(reference[[label_by, drop = TRUE]])
  object_for_transfer[[transfer_label_by]] <- unlabeled_category

  batch_by <- batch_by %||% transfer_control$batch_by %||% ".sn_transfer_batch"
  if (!batch_by %in% colnames(reference[[]])) {
    reference[[batch_by]] <- "reference"
  }
  if (!batch_by %in% colnames(object_for_transfer[[]])) {
    object_for_transfer[[batch_by]] <- "query"
  }

  common_features <- intersect(rownames(reference[[reference_assay]]), rownames(object_for_transfer[[query_assay]]))
  feature_set <- features %||% common_features
  feature_set <- intersect(feature_set, common_features)
  if (length(feature_set) < 2L) {
    stop("scANVI/scArches label_by transfer requires at least two shared features.", call. = FALSE)
  }

  old_reference_assay <- SeuratObject::DefaultAssay(reference)
  old_query_assay <- SeuratObject::DefaultAssay(object_for_transfer)
  on.exit(SeuratObject::DefaultAssay(reference) <- old_reference_assay, add = TRUE)
  on.exit(SeuratObject::DefaultAssay(object_for_transfer) <- old_query_assay, add = TRUE)
  SeuratObject::DefaultAssay(reference) <- reference_assay
  SeuratObject::DefaultAssay(object_for_transfer) <- query_assay
  combined <- .sn_with_default_acceleration(
    merge(reference, y = object_for_transfer, merge.data = FALSE),
    patches = "seurat_merge",
    strict = TRUE,
    operation = "merge"
  )
  combined_assay <- SeuratObject::DefaultAssay(combined)

  integration_control <- transfer_control
  integration_control$label_by <- transfer_label_by
  integration_control$label_by <- transfer_label_by
  integration_control$unlabeled_category <- unlabeled_category
  integration_control$reduction <- integration_control$reduction %||% method_name

  fit <- .sn_run_scvi_integration(
    object = combined,
    method = "scanvi",
    batch = batch_by,
    features = feature_set,
    assay = combined_assay,
    integration_control = integration_control,
    verbose = verbose
  )
  combined <- fit$object
  prediction_col <- transfer_control$prediction_col %||% "scanvi_prediction"
  if (!prediction_col %in% colnames(combined[[]])) {
    stop("scANVI/scArches backend did not return a prediction column.", call. = FALSE)
  }

  query_predictions <- combined[[prediction_col, drop = TRUE]][query_cells]
  metadata <- data.frame(row.names = colnames(object))
  metadata[[paste0(prediction_prefix, "_label")]] <- as.character(query_predictions)
  object <- Seurat::AddMetaData(object = object, metadata = metadata)
  object@misc$label_transfer[[prediction_prefix]] <- list(
    method = method_name,
    label_by = label_by,
    label_col = label_by,
    batch_by = batch_by,
    transfer_label_by = transfer_label_by,
    unlabeled_category = unlabeled_category,
    prediction_columns = colnames(metadata),
    run_dir = combined@misc$integration$run_dir %||% NULL,
    output_h5ad = combined@misc$integration$output_h5ad %||% NULL,
    transfer_control = transfer_control
  )
  object <- .sn_store_label_transfer_result(
    object = object,
    prediction_prefix = prediction_prefix,
    method = method_name,
    label_by = label_by,
    prediction_columns = colnames(metadata)
  )

  if (isTRUE(return_anchors)) {
    return(list(query = object, combined = combined, prediction_col = prediction_col))
  }
  object
}

.sn_coralysis_reference_mapping_backend <- function(...) {
  Coralysis::ReferenceMapping(...)
}

.sn_metadata_suffix <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  x <- gsub("^_+|_+$", "", x)
  x[!nzchar(x)] <- "value"
  x
}

.sn_store_label_transfer_result <- function(object,
                                            prediction_prefix,
                                            method,
                                            label_by,
                                            prediction_columns) {
  metadata <- object[[]]
  prediction_columns <- intersect(as.character(prediction_columns), colnames(metadata))
  label_column <- paste0(prediction_prefix, "_label")
  if (!label_column %in% prediction_columns) {
    stop(
      "Label-transfer output is missing expected metadata column `",
      label_column, "`.",
      call. = FALSE
    )
  }
  primary <- tibble::tibble(
    cell = rownames(metadata),
    prediction = as.character(metadata[[label_column]])
  )
  score_column <- paste0(prediction_prefix, "_score")
  if (score_column %in% prediction_columns) {
    primary$prediction_score <- as.numeric(metadata[[score_column]])
  }
  extra_columns <- setdiff(prediction_columns, c(label_column, score_column))
  for (column in extra_columns) {
    primary[[column]] <- metadata[[column]]
  }
  result <- .sn_new_analysis_result(
    analysis_type = "annotation",
    name = prediction_prefix,
    method = method,
    backend = switch(
      method,
      seurat = "Seurat::TransferData",
      coralysis = "Coralysis::ReferenceMapping",
      scanvi = "scvi-tools scANVI",
      scarches = "scvi-tools scArches",
      method
    ),
    input = list(cells = nrow(primary), label_by = label_by),
    parameters = list(
      prediction_prefix = prediction_prefix,
      prediction_columns = prediction_columns
    ),
    tables = list(primary = primary),
    diagnostics = list(
      labeled_cells = sum(!is.na(primary$prediction) & nzchar(primary$prediction)),
      missing_predictions = sum(is.na(primary$prediction) | !nzchar(primary$prediction)),
      label_levels = sort(unique(stats::na.omit(primary$prediction)))
    )
  )
  sn_store_result(object, "annotation", prediction_prefix, result)
}

.sn_get_seurat_logcounts_sce <- function(object,
                                         assay = NULL,
                                         layer = "data",
                                         verbose = TRUE) {
  check_installed(c("SingleCellExperiment", "SummarizedExperiment"))

  assay <- assay %||% SeuratObject::DefaultAssay(object = object)
  if (identical(layer, "data") && !"data" %in% SeuratObject::Layers(object = object[[assay]])) {
    if (verbose) .sn_log_info("Normalizing query object before Coralysis reference mapping.")
    object <- .sn_with_default_seurat_acceleration(
      Seurat::NormalizeData(object = object, assay = assay, verbose = verbose),
      object = object,
      assay = assay
    )
  }

  expr <- SeuratObject::LayerData(object = object, assay = assay, layer = layer)
  SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = expr),
    colData = object[[]]
  )
}

.sn_sync_coralysis_reference_label <- function(ref_sce,
                                              reference = NULL,
                                              label_col = NULL) {
  if (
    inherits(reference, "Seurat") &&
      is.character(label_col) &&
      length(label_col) == 1L &&
      label_col %in% colnames(reference[[]])
  ) {
    labels <- reference[[label_col, drop = TRUE]]
    names(labels) <- colnames(reference)
    SummarizedExperiment::colData(ref_sce)[[label_col]] <- labels[colnames(ref_sce)]
  }
  ref_sce
}

.sn_get_coralysis_reference_sce <- function(reference,
                                           reference_assay = NULL,
                                           reference_layer = "data",
                                           verbose = TRUE,
                                           label_col = NULL) {
  if (inherits(reference, "SingleCellExperiment")) {
    return(reference)
  }
  if (!inherits(reference, "Seurat")) {
    stop("`reference` must be a Seurat or SingleCellExperiment object for Coralysis label_by transfer.", call. = FALSE)
  }
  if (inherits(reference@misc$coralysis, "SingleCellExperiment")) {
    return(.sn_sync_coralysis_reference_label(
      ref_sce = reference@misc$coralysis,
      reference = reference,
      label_col = label_col
    ))
  }
  stop(
    "Coralysis label_by transfer requires a native Coralysis-trained reference stored under `reference@misc$coralysis`.\n",
    "Run `sn_run_cluster(reference, batch = ..., integration_method = \"coralysis\")` first, ",
    "or avoid `integration_control = list(store_sce = FALSE)` if the object should be used as a reference.",
    call. = FALSE
  )
}

.sn_minimize_coralysis_pca_model <- function(pca_model) {
  if (!is.list(pca_model)) {
    return(pca_model)
  }
  keep <- intersect(c("x", "center", "scale", "rotation"), names(pca_model))
  pca_model[keep]
}

.sn_prepare_coralysis_label_transfer_reference <- function(object,
                                                          label_by,
                                                          metadata_columns = NULL,
                                                          keep_umap_model = FALSE,
                                                          verbose = TRUE) {
  check_installed(c("SingleCellExperiment", "SummarizedExperiment", "S4Vectors"))

  ref_sce <- .sn_get_coralysis_reference_sce(
    reference = object,
    verbose = verbose,
    label_col = label_by
  )
  coldata <- SummarizedExperiment::colData(ref_sce)
  if (!label_by %in% colnames(coldata)) {
    stop(glue("`label_by` column '{label_by}' was not found in the Coralysis reference colData."), call. = FALSE)
  }

  coralysis <- S4Vectors::metadata(ref_sce)$coralysis
  required <- c("models", "pca.model", "pca.params")
  missing_required <- required[vapply(coralysis[required], is.null, logical(1))]
  if (length(missing_required) > 0L) {
    stop(
      "The Coralysis reference is missing required field(s): ",
      paste(missing_required, collapse = ", "),
      call. = FALSE
    )
  }
  if (is.null(coralysis$pca.params$select.icp.tables)) {
    stop("The Coralysis reference is missing `pca.params$select.icp.tables`.", call. = FALSE)
  }

  keep_cols <- unique(c(label_by, metadata_columns))
  keep_cols <- intersect(keep_cols, colnames(coldata))
  coldata <- coldata[, keep_cols, drop = FALSE]

  reference <- SingleCellExperiment::SingleCellExperiment(
    assays = list(),
    rowData = S4Vectors::DataFrame(row.names = rownames(ref_sce)),
    colData = coldata
  )
  reference_metadata <- S4Vectors::metadata(reference)
  reference_metadata$coralysis <- list(
    models = coralysis$models,
    pca.model = .sn_minimize_coralysis_pca_model(coralysis$pca.model),
    pca.params = list(select.icp.tables = coralysis$pca.params$select.icp.tables)
  )
  if (isTRUE(keep_umap_model) && !is.null(coralysis$umap.model)) {
    reference_metadata$coralysis$umap.model <- coralysis$umap.model
  }
  reference_metadata$shennong_reference <- list(
    method = "coralysis",
    label_by = label_by,
    n_features = nrow(reference),
    n_cells = ncol(reference),
    assays_dropped = TRUE,
    reduced_dims_dropped = TRUE,
    joint_probability_dropped = TRUE
  )
  S4Vectors::metadata(reference) <- reference_metadata
  reference
}

.sn_prepare_seurat_label_transfer_reference <- function(object,
                                                        label_by,
                                                        method,
                                                        assay = NULL,
                                                        layers = NULL,
                                                        features = NULL,
                                                        reduction = NULL,
                                                        metadata_columns = NULL) {
  check_installed("Seurat")
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object for this label-transfer method.", call. = FALSE)
  }
  if (!label_by %in% colnames(object[[]])) {
    stop(glue("`label_by` column '{label_by}' was not found in `object` metadata."), call. = FALSE)
  }

  assay <- assay %||% SeuratObject::DefaultAssay(object = object)
  if (!assay %in% names(object@assays)) {
    stop(glue("Assay '{assay}' was not found."), call. = FALSE)
  }

  available_layers <- SeuratObject::Layers(object = object[[assay]])
  if (is.null(layers)) {
    layers <- switch(
      method,
      seurat = intersect(c("counts", "data"), available_layers),
      scanvi = intersect("counts", available_layers),
      scarches = intersect("counts", available_layers)
    )
    if (length(layers) == 0L) {
      layers <- available_layers[1]
    }
  } else if (identical(layers, "all")) {
    layers <- available_layers
  }
  missing_layers <- setdiff(layers, available_layers)
  if (length(missing_layers) > 0L) {
    stop(
      "Layer(s) not found in assay '", assay, "': ",
      paste(missing_layers, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(reduction) && identical(method, "seurat") && "pca" %in% names(object@reductions)) {
    reduction <- "pca"
  }
  dimreducs <- reduction %||% character(0)

  old_assay <- SeuratObject::DefaultAssay(object = object)
  on.exit(SeuratObject::DefaultAssay(object = object) <- old_assay, add = TRUE)
  SeuratObject::DefaultAssay(object = object) <- assay
  reference <- Seurat::DietSeurat(
    object = object,
    assays = assay,
    layers = layers,
    features = features,
    dimreducs = dimreducs,
    graphs = character(0),
    misc = FALSE
  )

  keep_cols <- unique(c(label_by, metadata_columns))
  keep_cols <- intersect(keep_cols, colnames(reference[[]]))
  reference@meta.data <- reference@meta.data[, keep_cols, drop = FALSE]
  reference@misc$label_transfer_reference <- list(
    method = method,
    label_by = label_by,
    assay = assay,
    layers = layers,
    features = features,
    reduction = reduction
  )
  reference
}

.sn_write_label_transfer_reference <- function(reference,
                                               path = NULL,
                                               overwrite = FALSE,
                                               ...) {
  if (is.null(path)) {
    return(reference)
  }
  if (file.exists(path) && !isTRUE(overwrite)) {
    stop("`path` already exists. Set `overwrite = TRUE` to replace it: ", path, call. = FALSE)
  }
  sn_write(reference, path = path, ...)
  reference
}

#' Prepare a compact label-transfer reference
#'
#' \code{sn_prepare_label_transfer_reference()} converts a full analysis object
#' into a smaller reference object for \code{\link{sn_transfer_labels}}. For
#' native Coralysis, it returns a minimal \code{SingleCellExperiment} containing
#' the trained Coralysis models, PCA model, feature names, and selected
#' reference labels, while dropping the large reference assays, reduced
#' dimensions, and stored joint probabilities. For Seurat, scANVI, and scArches
#' workflows, it returns a slim Seurat object with only the selected assay
#' layers, labels, optional features, and optional reduction.
#'
#' @param object A Seurat object. For \code{method = "coralysis"}, an existing
#'   Coralysis-trained \code{SingleCellExperiment} is also accepted.
#' @param label_by Metadata column containing reference labels.
#' @param method Label-transfer backend the reference should support.
#' @param assay Assay to keep for Seurat/scANVI/scArches references.
#' @param layers Layers to keep for Seurat/scANVI/scArches references. Defaults
#'   to \code{counts} and \code{data} for Seurat transfer, and \code{counts} for
#'   scANVI/scArches. Use \code{"all"} to retain all layers in \code{assay}.
#' @param features Optional features to keep for Seurat/scANVI/scArches
#'   references.
#' @param reduction Optional dimensional reduction to keep. Defaults to
#'   \code{"pca"} for Seurat references when available.
#' @param metadata_columns Additional metadata columns to keep alongside
#'   \code{label_by}.
#' @param keep_umap_model For Coralysis references, keep the UMAP projection
#'   model if present. This is only needed when using
#'   \code{transfer_control = list(project.umap = TRUE)}.
#' @param path Optional output path. Use \code{.qs2} for serialized reference
#'   objects.
#' @param overwrite Logical; overwrite an existing \code{path}.
#' @param verbose Whether to print progress messages.
#' @param ... Additional arguments passed to \code{\link{sn_write}} when
#'   \code{path} is supplied.
#'
#' @return A compact Seurat or SingleCellExperiment reference object.
#'
#' @examples
#' \dontrun{
#' coral_ref <- sn_prepare_label_transfer_reference(
#'   reference,
#'   label_by = "cell_type",
#'   method = "coralysis",
#'   path = "data/processed/pbmc_coralysis_reference.qs2",
#'   overwrite = TRUE
#' )
#'
#' query <- sn_transfer_labels(
#'   query,
#'   reference = coral_ref,
#'   label_by = "cell_type",
#'   method = "coralysis"
#' )
#' }
#' @export
sn_prepare_label_transfer_reference <- function(object,
                                                label_by,
                                                method = c("coralysis", "seurat", "scanvi", "scarches"),
                                                assay = NULL,
                                                layers = NULL,
                                                features = NULL,
                                                reduction = NULL,
                                                metadata_columns = NULL,
                                                keep_umap_model = FALSE,
                                                path = NULL,
                                                overwrite = FALSE,
                                                verbose = TRUE,
                                                ...) {
  method <- match.arg(method)
  if (!is.character(label_by) || length(label_by) != 1L || !nzchar(label_by)) {
    stop("`label_by` must be a non-empty metadata column name.", call. = FALSE)
  }

  reference <- if (identical(method, "coralysis")) {
    .sn_prepare_coralysis_label_transfer_reference(
      object = object,
      label_by = label_by,
      metadata_columns = metadata_columns,
      keep_umap_model = keep_umap_model,
      verbose = verbose
    )
  } else {
    .sn_prepare_seurat_label_transfer_reference(
      object = object,
      label_by = label_by,
      method = method,
      assay = assay,
      layers = layers,
      features = features,
      reduction = reduction,
      metadata_columns = metadata_columns
    )
  }

  .sn_write_label_transfer_reference(
    reference = reference,
    path = path,
    overwrite = overwrite,
    ...
  )
}

.sn_transfer_labels_coralysis <- function(object,
                                          reference,
                                          label_col,
                                          prediction_prefix,
                                          reference_assay = NULL,
                                          query_assay = NULL,
                                          reference_layer = "data",
                                          query_layer = "data",
                                          transfer_control = list(),
                                          return_anchors = FALSE,
                                          verbose = TRUE) {
  check_installed("Coralysis")
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat query object for Coralysis label_by transfer.", call. = FALSE)
  }
  if (!is.list(transfer_control)) {
    stop("`transfer_control` must be a named list.", call. = FALSE)
  }

  ref_sce <- .sn_get_coralysis_reference_sce(
    reference = reference,
    reference_assay = reference_assay,
    reference_layer = reference_layer,
    verbose = verbose,
    label_col = label_col
  )
  if (!label_col %in% colnames(SummarizedExperiment::colData(ref_sce))) {
    stop(glue("`label_col` column '{label_col}' was not found in Coralysis reference colData."), call. = FALSE)
  }
  query_sce <- .sn_get_seurat_logcounts_sce(
    object = object,
    assay = query_assay,
    layer = query_layer,
    verbose = verbose
  )
  query_sce <- Coralysis::PrepareData(object = query_sce)

  args <- .sn_merge_control_args(
    defaults = list(
      ref = ref_sce,
      query = query_sce,
      ref.label = label_col
    ),
    control = transfer_control
  )
  mapped <- do.call(.sn_coralysis_reference_mapping_backend, args)
  mapped_coldata <- as.data.frame(SummarizedExperiment::colData(mapped))
  mapped_coldata <- mapped_coldata[colnames(object), , drop = FALSE]

  label_source <- if ("pruned_coral_labels" %in% colnames(mapped_coldata)) {
    "pruned_coral_labels"
  } else {
    "coral_labels"
  }
  if (!label_source %in% colnames(mapped_coldata)) {
    stop("Coralysis::ReferenceMapping() did not return `coral_labels`.", call. = FALSE)
  }

  metadata <- data.frame(row.names = colnames(object))
  metadata[[paste0(prediction_prefix, "_label")]] <- mapped_coldata[[label_source]]
  if ("coral_probability" %in% colnames(mapped_coldata)) {
    metadata[[paste0(prediction_prefix, "_score")]] <- mapped_coldata$coral_probability
  }
  if ("coral_labels" %in% colnames(mapped_coldata)) {
    metadata[[paste0(prediction_prefix, "_raw_label")]] <- mapped_coldata$coral_labels
  }

  object <- Seurat::AddMetaData(object, metadata = metadata)
  object@misc$label_transfer[[prediction_prefix]] <- list(
    method = "coralysis",
    label_col = label_col,
    prediction_columns = colnames(metadata),
    transfer_control = transfer_control
  )
  object <- .sn_store_label_transfer_result(
    object = object,
    prediction_prefix = prediction_prefix,
    method = "coralysis",
    label_by = label_col,
    prediction_columns = colnames(metadata)
  )

  if (isTRUE(return_anchors)) {
    return(list(query = object, mapping = mapped))
  }
  object
}

#' Transfer labels from a Seurat reference to a query object
#'
#' \code{sn_transfer_labels()} is a Shennong wrapper for reference mapping. It
#' keeps the common path compact: transfer one metadata label, add the predicted
#' label and confidence score back to the query, store cell-level predictions
#' as a canonical \code{annotation} result, and retain a compact compatibility
#' record in \code{query@misc$label_transfer}. The default \code{method =
#' "seurat"} wraps Seurat's \code{FindTransferAnchors()} and
#' \code{TransferData()} workflow. \code{method = "coralysis"} projects the
#' query onto a native Coralysis-trained reference with
#' \code{Coralysis::ReferenceMapping()}. \code{method = "scanvi"} and
#' \code{method = "scarches"} use the managed scVI-family pixi backend to train
#' a semi-supervised scANVI model with reference labels and query cells marked
#' as unlabeled, then import the predicted query labels.
#'
#' @param object A Seurat query object to annotate. This argument comes first
#'   so the function can be used in pipes.
#' @param reference A labeled Seurat reference object.
#' @param label_by Metadata column in \code{reference} to transfer.
#' @param method Label-transfer backend. \code{"seurat"} uses Seurat anchors;
#'   \code{"coralysis"} uses native \code{Coralysis::ReferenceMapping()} and
#'   requires a trained Coralysis reference stored under
#'   \code{reference@misc$coralysis}; \code{"scanvi"} and \code{"scarches"} use
#'   the scVI-family pixi backend.
#' @param prediction_prefix Prefix for metadata columns added to
#'   \code{query}. Defaults to \code{paste0(label_by, "_transfer")}.
#' @param normalization_method Normalization method passed to
#'   \code{Seurat::FindTransferAnchors()}.
#' @param reference_assay,query_assay Assays passed to
#'   \code{Seurat::FindTransferAnchors()} and \code{Seurat::TransferData()}.
#'   For \code{method = "coralysis"}, \code{query_assay} controls the assay
#'   converted to query \code{logcounts}; the reference assay is ignored when a
#'   stored Coralysis reference is available.
#' @param reference_layer,query_layer Layers used as log-normalized expression
#'   for \code{method = "coralysis"}. The query defaults to the Seurat
#'   \code{"data"} layer and is normalized first if that layer is absent.
#' @param reduction Dimensional reduction strategy passed to
#'   \code{Seurat::FindTransferAnchors()}.
#' @param reference_reduction Optional reference reduction passed to
#'   \code{Seurat::FindTransferAnchors()}.
#' @param features Optional features used to find transfer anchors.
#' @param dims Dimensions used for anchor scoring and label_by transfer.
#' @param npcs Number of PCs used by \code{Seurat::FindTransferAnchors()}.
#' @param k_anchor,k_filter,k_score,k_weight Seurat anchor/weighting
#'   parameters.
#' @param store_prediction_scores If \code{TRUE}, also store per-label
#'   prediction scores as query metadata columns.
#' @param return_anchors If \code{TRUE}, return a list containing the annotated
#'   query and backend artifacts. For Coralysis, the artifact is the mapped
#'   SingleCellExperiment.
#' @param transfer_control Optional backend-specific list. For
#'   \code{method = "coralysis"}, values are forwarded to
#'   \code{Coralysis::ReferenceMapping()}. For \code{method = "scanvi"} or
#'   \code{"scarches"}, common values include \code{batch_by},
#'   \code{runtime_dir}, \code{pixi_project}, \code{max_epochs},
#'   \code{scanvi_max_epochs}, \code{accelerator}, \code{mirror}, and
#'   \code{install_pixi}.
#' @param verbose Whether to print Seurat progress messages.
#' @param ... Additional arguments passed to \code{Seurat::FindTransferAnchors()}.
#'
#' @return A Seurat query object with transferred labels, or a list when
#'   \code{return_anchors = TRUE}.
#'
#' @examples
#' \dontrun{
#' query <- sn_transfer_labels(
#'   object = query,
#'   reference = reference,
#'   label_by = "cell_type",
#'   dims = 1:30
#' )
#' }
#' @export
sn_transfer_labels <- function(object = NULL,
                               reference,
                               label_by = NULL,
                               method = c("seurat", "coralysis", "scanvi", "scarches"),
                               prediction_prefix = NULL,
                               normalization_method = "LogNormalize",
                               reference_assay = NULL,
                               query_assay = NULL,
                               reference_layer = "data",
                               query_layer = "data",
                               reduction = "pcaproject",
                               reference_reduction = NULL,
                               features = NULL,
                               dims = 1:30,
                               npcs = 30,
                               k_anchor = 5,
                               k_filter = NA,
                               k_score = 30,
                               k_weight = 50,
                               store_prediction_scores = FALSE,
                               return_anchors = FALSE,
                               transfer_control = list(),
                               verbose = TRUE,
                               ...) {
  check_installed("Seurat")
  method <- match.arg(method)
  if (is.null(object)) {
    stop("`object` must be supplied as the query Seurat object.", call. = FALSE)
  }

  if (!inherits(reference, "Seurat") && !(identical(method, "coralysis") && inherits(reference, "SingleCellExperiment"))) {
    stop("`reference` must be a Seurat object.", call. = FALSE)
  }
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat query object.", call. = FALSE)
  }
  if (!is.character(label_by) || length(label_by) != 1L || !nzchar(label_by)) {
    stop("`label_by` must be a non-empty metadata column name.", call. = FALSE)
  }
  if (method == "coralysis") {
    prediction_prefix <- prediction_prefix %||% paste0(label_by, "_coralysis")
    return(.sn_transfer_labels_coralysis(
      object = object,
      reference = reference,
      label_col = label_by,
      prediction_prefix = prediction_prefix,
      reference_assay = reference_assay,
      query_assay = query_assay,
      reference_layer = reference_layer,
      query_layer = query_layer,
      transfer_control = transfer_control,
      return_anchors = return_anchors,
      verbose = verbose
    ))
  }
  if (method %in% c("scanvi", "scarches")) {
    prediction_prefix <- prediction_prefix %||% paste0(label_by, "_", method)
    return(.sn_transfer_labels_scanvi(
      object = object,
      reference = reference,
      label_by = label_by,
      prediction_prefix = prediction_prefix,
      assay = query_assay %||% reference_assay,
      batch_by = transfer_control$batch_by %||% NULL,
      features = features,
      transfer_control = transfer_control,
      return_anchors = return_anchors,
      verbose = verbose,
      method_name = method
    ))
  }

  if (!label_by %in% colnames(reference[[]])) {
    stop(glue("`label_by` column '{label_by}' was not found in `reference` metadata."), call. = FALSE)
  }

  ref_labels <- reference[[label_by, drop = TRUE]]
  names(ref_labels) <- colnames(reference)
  prediction_prefix <- prediction_prefix %||% paste0(label_by, "_transfer")

  anchors <- .sn_find_transfer_anchors_backend(
    reference = reference,
    query = object,
    normalization.method = normalization_method,
    reference.assay = reference_assay,
    query.assay = query_assay,
    reduction = reduction,
    reference.reduction = reference_reduction,
    features = features,
    npcs = npcs,
    dims = dims,
    k.anchor = k_anchor,
    k.filter = k_filter,
    k.score = k_score,
    verbose = verbose,
    ...
  )

  predictions <- .sn_transfer_data_backend(
    anchorset = anchors,
    refdata = ref_labels,
    weight.reduction = reduction,
    dims = dims,
    k.weight = k_weight,
    verbose = verbose
  )
  predictions <- as.data.frame(predictions)
  if (!all(colnames(object) %in% rownames(predictions)) && nrow(predictions) == ncol(object)) {
    rownames(predictions) <- colnames(object)
  }
  predictions <- predictions[colnames(object), , drop = FALSE]
  if (!"predicted.id" %in% colnames(predictions)) {
    stop("Seurat::TransferData() did not return a `predicted.id` column.", call. = FALSE)
  }

  metadata <- data.frame(row.names = colnames(object))
  metadata[[paste0(prediction_prefix, "_label")]] <- predictions$predicted.id
  if ("prediction.score.max" %in% colnames(predictions)) {
    metadata[[paste0(prediction_prefix, "_score")]] <- predictions$prediction.score.max
  }
  if (isTRUE(store_prediction_scores)) {
    score_cols <- grep("^prediction\\.score\\.", colnames(predictions), value = TRUE)
    score_cols <- setdiff(score_cols, "prediction.score.max")
    for (score_col in score_cols) {
      suffix <- sub("^prediction\\.score\\.", "", score_col)
      metadata[[paste0(prediction_prefix, "_score_", .sn_metadata_suffix(suffix))]] <- predictions[[score_col]]
    }
  }

  object <- Seurat::AddMetaData(object, metadata = metadata)
  object@misc$label_transfer[[prediction_prefix]] <- list(
    label_by = label_by,
    label_col = label_by,
    normalization_method = normalization_method,
    reduction = reduction,
    dims = dims,
    features = features,
    prediction_columns = colnames(metadata)
  )
  object <- .sn_store_label_transfer_result(
    object = object,
    prediction_prefix = prediction_prefix,
    method = "seurat",
    label_by = label_by,
    prediction_columns = colnames(metadata)
  )

  if (isTRUE(return_anchors)) {
    return(list(query = object, anchors = anchors, predictions = predictions))
  }
  object
}
