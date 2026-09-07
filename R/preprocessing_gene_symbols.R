.sn_standardize_assay_layer <- function(layer_data,
                                        feature_map,
                                        count_like = FALSE) {
  feature_index <- match(rownames(layer_data), names(feature_map))
  if (anyNA(feature_index)) {
    stop("An RNA assay layer contained features that were absent from the counts layer.", call. = FALSE)
  }

  mapped <- unname(feature_map[feature_index])
  keep <- !is.na(mapped) & nzchar(mapped)
  if (!any(keep)) {
    return(NULL)
  }
  layer_data <- layer_data[keep, , drop = FALSE]
  mapped <- mapped[keep]
  if (isTRUE(count_like)) {
    rownames(layer_data) <- mapped
    if (anyDuplicated(mapped)) {
      layer_data <- .sn_aggregate_rows_by_group(layer_data, groups = mapped)
    }
    return(layer_data)
  }

  duplicated_targets <- names(which(table(feature_map, useNA = "no") > 1L))
  safe <- !mapped %in% duplicated_targets
  if (!any(safe)) {
    return(NULL)
  }
  layer_data <- layer_data[safe, , drop = FALSE]
  rownames(layer_data) <- mapped[safe]
  layer_data
}

.sn_assay_derived_component <- function(component, assay) {
  if (!"assay.used" %in% methods::slotNames(component)) {
    return(FALSE)
  }
  used <- tryCatch(
    methods::slot(component, "assay.used"),
    error = function(error) character()
  )
  any(as.character(used) %in% assay)
}

.sn_invalidate_assay_derived_state <- function(object, assay) {
  object@reductions <- object@reductions[!vapply(
    object@reductions,
    .sn_assay_derived_component,
    logical(1),
    assay = assay
  )]
  object@graphs <- object@graphs[!vapply(
    object@graphs,
    .sn_assay_derived_component,
    logical(1),
    assay = assay
  )]
  object@commands <- object@commands[!vapply(
    object@commands,
    .sn_assay_derived_component,
    logical(1),
    assay = assay
  )]
  object
}

.sn_rebuild_standardized_rna_assay <- function(object,
                                               feature_map,
                                               standardized_counts) {
  old_assay <- object[["RNA"]]
  old_layers <- SeuratObject::Layers(old_assay)
  count_layers <- grep("^counts(\\.|$)", old_layers, value = TRUE)
  if (length(count_layers) == 0L) {
    stop("The RNA assay must contain a counts layer for gene-symbol standardization.", call. = FALSE)
  }

  transformed_counts <- stats::setNames(lapply(count_layers, function(current_layer) {
    .sn_standardize_assay_layer(
      layer_data = SeuratObject::LayerData(
        object = object,
        assay = "RNA",
        layer = current_layer
      ),
      feature_map = feature_map,
      count_like = TRUE
    )
  }), sub("^counts\\.?", "", count_layers))
  names(transformed_counts)[!nzchar(names(transformed_counts))] <- "counts"

  new_assay <- if (length(transformed_counts) == 1L && identical(count_layers, "counts")) {
    SeuratObject::CreateAssay5Object(counts = standardized_counts)
  } else {
    SeuratObject::CreateAssay5Object(counts = transformed_counts)
  }

  other_layers <- setdiff(old_layers, count_layers)
  omitted_layers <- character()
  for (current_layer in other_layers) {
    layer_data <- .sn_standardize_assay_layer(
      layer_data = SeuratObject::LayerData(
        object = object,
        assay = "RNA",
        layer = current_layer
      ),
      feature_map = feature_map,
      count_like = grepl("count", current_layer, ignore.case = TRUE)
    )
    if (is_null(layer_data)) {
      omitted_layers <- c(omitted_layers, current_layer)
      next
    }
    SeuratObject::LayerData(new_assay, layer = current_layer) <- layer_data
  }

  old_metadata <- old_assay[[]]
  if (nrow(old_metadata) > 0L) {
    mapped <- unname(feature_map[match(rownames(old_metadata), names(feature_map))])
    unique_targets <- names(which(table(feature_map, useNA = "no") == 1L))
    keep_metadata <- !is.na(mapped) & mapped %in% unique_targets
    if (any(keep_metadata)) {
      metadata <- old_metadata[keep_metadata, , drop = FALSE]
      rownames(metadata) <- mapped[keep_metadata]
      new_assay <- SeuratObject::AddMetaData(new_assay, metadata = metadata)
    }
  }

  old_variable_features <- SeuratObject::VariableFeatures(old_assay)
  if (length(old_variable_features) > 0L) {
    mapped_variable_features <- unname(feature_map[match(old_variable_features, names(feature_map))])
    unique_targets <- names(which(table(feature_map, useNA = "no") == 1L))
    mapped_variable_features <- unique(mapped_variable_features[
      !is.na(mapped_variable_features) & mapped_variable_features %in% unique_targets
    ])
    if (length(mapped_variable_features) > 0L) {
      SeuratObject::VariableFeatures(new_assay) <- mapped_variable_features
    }
  }
  SeuratObject::Key(new_assay) <- SeuratObject::Key(old_assay)
  new_assay@assay.orig <- old_assay@assay.orig
  new_assay@misc <- old_assay@misc

  object[["RNA"]] <- new_assay
  object <- .sn_invalidate_assay_derived_state(object, assay = "RNA")
  if (length(omitted_layers) > 0L) {
    .sn_log_warn(
      "Dropped RNA assay layer(s) whose values could not be safely combined after gene-symbol standardization: {paste(omitted_layers, collapse = ', ')}."
    )
  }
  object
}
