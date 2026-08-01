# SeuratObject::JoinLayers.Assay5 accelerator.
#
# Lifted from autozyme task `test_seurat_joinlayers`; accepted at 55da460.
# The fast path is deliberately restricted to exact plain Assay5 objects and
# explicit counts -> counts joins. Every other call delegates to the captured
# upstream method.

if (requireNamespace("SeuratObject", quietly = TRUE) &&
    requireNamespace("Matrix", quietly = TRUE)) {

getFromNamespace <- utils::getFromNamespace
slot <- methods::slot
`slot<-` <- methods::`slot<-`
validObject <- methods::validObject

original_join_layers_assay5 <- utils::getFromNamespace(
  "JoinLayers.Assay5", "SeuratObject"
)

fast_join_layers_assay5 <- local({
  original <- original_join_layers_assay5
  stitch_matrix <- getFromNamespace("StitchMatrix", "SeuratObject")
  prep_layer_data <- getFromNamespace(".PrepLayerData2", "SeuratObject")
  check_gc <- getFromNamespace("CheckGC", "SeuratObject")
  layer_data <- getFromNamespace("LayerData", "SeuratObject")
  sparse_empty_matrix <- getFromNamespace("SparseEmptyMatrix", "SeuratObject")
  variable_features <- getFromNamespace("VariableFeatures", "SeuratObject")
  set_variable_features <- getFromNamespace("VariableFeatures<-", "SeuratObject")
  layers_fn <- getFromNamespace("Layers", "SeuratObject")

  has_bpcells_marker <- function(x) {
    cls <- class(x = x)
    identical(x = attr(x = cls, which = "package"), y = "BPCells")
  }

  has_exact_dgc_marker <- function(x) {
    cls <- class(x = x)
    identical(x = cls[[1L]], y = "dgCMatrix") &&
      identical(x = attr(x = cls, which = "package"), y = "Matrix")
  }

  fast_sparse_cbind_aligned <- function(mats) {
    if (!requireNamespace("Matrix", quietly = TRUE)) {
      return(NULL)
    }
    tryCatch(
      expr = {
        joined <- do.call(what = cbind, args = mats)
        if (!has_exact_dgc_marker(joined)) {
          return(NULL)
        }
        joined
      },
      error = function(e) NULL
    )
  }

  build_logmap <- function(map, selected_layers, layer, values) {
    mat <- slot(object = map, name = ".Data")
    keep_cols <- colnames(x = mat)[!(colnames(x = mat) %in% selected_layers)]
    mat <- mat[, keep_cols, drop = FALSE]

    idx <- as.integer(x = values)
    new_col <- matrix(
      data = FALSE,
      nrow = nrow(x = map),
      ncol = 1L,
      dimnames = list(rownames(x = map), layer)
    )
    if (length(x = idx)) {
      new_col[idx, layer] <- TRUE
    }
    mat <- cbind(mat, new_col)
    slot(object = map, name = ".Data") <- mat
    map
  }

  prove_contiguous_membership <- function(object, selected_layers) {
    feature_map <- slot(object = object, name = "features")
    cell_map <- slot(object = object, name = "cells")
    rowmap <- feature_map[, selected_layers]
    colmap <- cell_map[, selected_layers]
    row_membership <- slot(object = rowmap, name = ".Data")
    col_membership <- slot(object = colmap, name = ".Data")
    if (!identical(x = colnames(x = row_membership), y = selected_layers) ||
        !identical(x = colnames(x = col_membership), y = selected_layers)) {
      return(NULL)
    }

    if (!is.matrix(x = row_membership) ||
        !is.matrix(x = col_membership) ||
        !is.logical(x = row_membership) ||
        !is.logical(x = col_membership) ||
        !all(row_membership)) {
      return(NULL)
    }

    feature_names <- rownames(x = feature_map)
    cell_names <- rownames(x = cell_map)
    n_cells <- length(x = cell_names)
    cell_counts <- as.integer(x = colSums(x = col_membership))
    if (any(cell_counts <= 0L) ||
        !identical(x = as.integer(sum(cell_counts)), y = n_cells) ||
        !identical(
          x = as.integer(rowSums(x = col_membership)),
          y = rep.int(1L, nrow(x = col_membership))
        )) {
      return(NULL)
    }
    range_ends <- cumsum(cell_counts)
    range_starts <- range_ends - cell_counts + 1L
    expected_layer_id <- rep.int(seq_along(along.with = cell_counts), times = cell_counts)
    observed_layer_id <- max.col(m = col_membership, ties.method = "first")
    if (!identical(x = observed_layer_id, y = expected_layer_id)) {
      return(NULL)
    }

    list(
      feature_map = feature_map,
      cell_map = cell_map,
      feature_names = feature_names,
      cell_names = cell_names,
      n_features = length(x = feature_names),
      n_cells = n_cells,
      cell_counts = cell_counts,
      range_starts = range_starts,
      range_ends = range_ends
    )
  }

  prove_selected_membership <- function(object, selected_layers) {
    feature_map <- slot(object = object, name = "features")
    cell_map <- slot(object = object, name = "cells")
    rowmap <- feature_map[, selected_layers]
    colmap <- cell_map[, selected_layers]
    row_membership <- slot(object = rowmap, name = ".Data")
    col_membership <- slot(object = colmap, name = ".Data")
    if (!identical(x = colnames(x = row_membership), y = selected_layers) ||
        !identical(x = colnames(x = col_membership), y = selected_layers)) {
      return(NULL)
    }
    if (!is.matrix(x = row_membership) ||
        !is.matrix(x = col_membership) ||
        !is.logical(x = row_membership) ||
        !is.logical(x = col_membership)) {
      return(NULL)
    }

    feature_names <- rownames(x = feature_map)
    cell_names <- rownames(x = cell_map)
    feature_counts <- as.integer(x = colSums(x = row_membership))
    cell_counts <- as.integer(x = colSums(x = col_membership))
    selected_feature_rows <- rowSums(x = row_membership) > 0L
    cell_row_sums <- as.integer(rowSums(x = col_membership))
    selected_cell_rows <- cell_row_sums > 0L
    if (any(feature_counts <= 0L) ||
        any(cell_counts <= 0L) ||
        any(cell_row_sums > 1L) ||
        !any(selected_feature_rows) ||
        !any(selected_cell_rows)) {
      return(NULL)
    }

    all_layers_have_selected_features <- all(
      feature_counts == sum(selected_feature_rows)
    )
    if (isTRUE(all_layers_have_selected_features) &&
        anyDuplicated(feature_names[selected_feature_rows])) {
      return(NULL)
    }
    layer_features <- if (isTRUE(all_layers_have_selected_features)) {
      NULL
    } else {
      vector(mode = "list", length = length(x = selected_layers))
    }
    layer_cells <- vector(mode = "list", length = length(x = selected_layers))
    for (i in seq_along(along.with = selected_layers)) {
      if (!isTRUE(all_layers_have_selected_features)) {
        layer_features[[i]] <- feature_names[row_membership[, i]]
        if (anyDuplicated(layer_features[[i]])) {
          return(NULL)
        }
      }
      layer_cells[[i]] <- cell_names[col_membership[, i]]
      if (anyDuplicated(layer_cells[[i]])) {
        return(NULL)
      }
    }
    if (!is.null(x = layer_features)) {
      names(layer_features) <- selected_layers
    }
    names(layer_cells) <- selected_layers

    list(
      feature_map = feature_map,
      cell_map = cell_map,
      feature_names = feature_names,
      cell_names = cell_names,
      selected_feature_idx = which(selected_feature_rows),
      selected_feature_names = feature_names[selected_feature_rows],
      selected_cell_idx = which(selected_cell_rows),
      selected_cell_names = cell_names[selected_cell_rows],
      n_features = length(x = feature_names),
      n_cells = length(x = cell_names),
      layer_features = layer_features,
      layer_cells = layer_cells,
      concat_cells = unlist(layer_cells, use.names = FALSE),
      all_layers_have_selected_features = all_layers_have_selected_features
    )
  }

  proof_layer_features <- function(proof, i) {
    if (is.null(x = proof$layer_features)) {
      proof$selected_feature_names
    } else {
      proof$layer_features[[i]]
    }
  }

  bpcells_layer_contract_ok <- function(x, expected_features, expected_cells) {
    methods::is(x, "IterableMatrix") &&
      identical(
        x = as.integer(dim(x = x)),
        y = c(length(x = expected_features), length(x = expected_cells))
      ) &&
      identical(x = rownames(x = x), y = expected_features) &&
      identical(x = colnames(x = x), y = expected_cells)
  }

  stitch_selected_layers <- function(mats, proof, selected_layers) {
    tryCatch(
      expr = stitch_matrix(
        x = mats[[1L]],
        y = mats[-1L],
        rowmap = proof$feature_map[, selected_layers],
        colmap = proof$cell_map[, selected_layers]
      ),
      error = function(e) NULL
    )
  }

  public_stitch_iterable <- function(mats, proof) {
    target_features <- proof$selected_feature_names
    if (!isTRUE(proof$all_layers_have_selected_features)) {
      for (i in seq_along(along.with = mats)) {
        missing_row <- setdiff(
          x = target_features,
          y = proof_layer_features(proof = proof, i = i)
        )
        if (length(x = missing_row) > 0L) {
          zero_i <- sparse_empty_matrix(
            nrow = length(x = missing_row),
            ncol = ncol(x = mats[[i]]),
            colnames = proof$layer_cells[[i]],
            rownames = missing_row
          )
          zero_i <- methods::as(object = zero_i, Class = "IterableMatrix")
          mats[[i]] <- rbind(mats[[i]], zero_i)[
            target_features,
            ,
            drop = FALSE
          ]
        }
      }
    }
    tryCatch(
      expr = Reduce(f = cbind, x = mats),
      error = function(e) NULL
    )
  }

  prepare_stitched_value <- function(joined, proof) {
    cell_reorder <- match(
      x = proof$selected_cell_names,
      table = proof$concat_cells
    )
    if (anyNA(cell_reorder) ||
        anyDuplicated(cell_reorder) ||
        !identical(x = sort(cell_reorder), y = seq_along(cell_reorder))) {
      return(NULL)
    }
    value <- prep_layer_data(
      x = joined,
      target = c(proof$n_features, proof$n_cells),
      dnames = list(proof$selected_feature_names, proof$concat_cells),
      fmargin = 1L
    )
    value <- value[
      seq_len(nrow(x = value)),
      cell_reorder,
      drop = FALSE
    ]
    expected_dim <- c(
      length(x = proof$selected_feature_names),
      length(x = proof$selected_cell_names)
    )
    if (!identical(x = as.integer(dim(x = value)), y = as.integer(expected_dim))) {
      return(NULL)
    }
    value
  }

  dgc_contract_ok <- function(x, expected_features, expected_cells) {
    dnames <- dimnames(x = x)
    rownames_ok <- is.null(x = dnames[[1L]]) ||
      identical(x = dnames[[1L]], y = expected_features)
    colnames_ok <- is.null(x = dnames[[2L]]) ||
      identical(x = dnames[[2L]], y = expected_cells)
    has_exact_dgc_marker(x) &&
      identical(
        x = as.integer(dim(x = x)),
        y = c(length(x = expected_features), length(x = expected_cells))
      ) &&
      rownames_ok &&
      colnames_ok
  }

  finalize_fast_join <- function(
    object,
    layer_store,
    selected_layers,
    new,
    value,
    proof,
    var.features
  ) {
    on.exit(expr = check_gc(), add = TRUE)

    unselected_layers <- names(x = layer_store)[!(names(x = layer_store) %in% selected_layers)]
    final_layer_list <- c(stats::setNames(object = list(value), nm = new), layer_store[unselected_layers])
    final_layers <- names(x = final_layer_list)
    feature_values <- proof$selected_feature_idx
    if (is.null(x = feature_values)) {
      feature_values <- seq_len(length.out = proof$n_features)
    }

    slot(object = object, name = "layers") <- final_layer_list
    slot(object = object, name = "features") <- build_logmap(
      map = proof$feature_map,
      selected_layers = selected_layers,
      layer = new,
      values = feature_values
    )
    slot(object = object, name = "cells") <- build_logmap(
      map = proof$cell_map,
      selected_layers = selected_layers,
      layer = new,
      values = if (is.null(x = proof$selected_cell_idx)) {
        seq_len(length.out = proof$n_cells)
      } else {
        proof$selected_cell_idx
      }
    )
    slot(object = object, name = "default") <- if (length(x = final_layers)) 1L else 0L
    object <- set_variable_features(object = object, value = var.features)
    validObject(object = object)
    object
  }

  function(
    object,
    layers = NULL,
    new = NULL,
    ...
  ) {
    fallback <- function() {
      original(object = object, layers = layers, new = new, ...)
    }

    object_class <- class(x = object)
    if (length(x = object_class) != 1L ||
        !identical(x = object_class[[1L]], y = "Assay5") ||
        !identical(x = attr(x = object_class, which = "package"), y = "SeuratObject")) {
      return(fallback())
    }
    if (!identical(x = layers, y = "counts") ||
        !identical(x = new, y = "counts")) {
      return(fallback())
    }
    if (...length() != 0L) {
      return(fallback())
    }

    selected_layers <- tryCatch(
      expr = suppressWarnings(layers_fn(object = object, search = layers)),
      error = function(e) character()
    )
    if (length(x = selected_layers) <= 1L) {
      return(fallback())
    }

    layer_store <- slot(object = object, name = "layers")
    if (!all(selected_layers %in% names(x = layer_store)) ||
        new %in% names(x = layer_store)) {
      return(fallback())
    }

    selected_raw_layers <- unname(obj = layer_store[selected_layers])
    all_bpcells <- all(vapply(
      X = selected_raw_layers,
      FUN = has_bpcells_marker,
      FUN.VALUE = logical(1L)
    ))
    all_dgc <- all(vapply(
      X = selected_raw_layers,
      FUN = has_exact_dgc_marker,
      FUN.VALUE = logical(1L)
    ))
    if (!isTRUE(all_bpcells) && !isTRUE(all_dgc)) {
      return(fallback())
    }

    var.features <- variable_features(object = object)

    if (isTRUE(all_dgc)) {
      contiguous_proof <- prove_contiguous_membership(
        object = object,
        selected_layers = selected_layers
      )
      if (!is.null(x = contiguous_proof)) {
        direct_ok <- TRUE
        for (i in seq_along(along.with = selected_raw_layers)) {
          cell_slice <- contiguous_proof$cell_names[
            contiguous_proof$range_starts[[i]]:contiguous_proof$range_ends[[i]]
          ]
          if (!dgc_contract_ok(
            x = selected_raw_layers[[i]],
            expected_features = contiguous_proof$feature_names,
            expected_cells = cell_slice
          )) {
            direct_ok <- FALSE
            break
          }
        }
        if (isTRUE(direct_ok)) {
          joined <- fast_sparse_cbind_aligned(mats = selected_raw_layers)
          if (!is.null(x = joined)) {
            value <- prep_layer_data(
              x = joined,
              target = dim(x = object),
              dnames = list(
                contiguous_proof$feature_names,
                contiguous_proof$cell_names
              ),
              fmargin = 1L
            )
            value <- value[
              seq_len(nrow(x = value)),
              seq_len(ncol(x = value)),
              drop = FALSE
            ]
            if (identical(
              x = as.integer(dim(x = value)),
              y = as.integer(dim(x = object))
            )) {
              return(finalize_fast_join(
                object = object,
                layer_store = layer_store,
                selected_layers = selected_layers,
                new = new,
                value = value,
                proof = contiguous_proof,
                var.features = var.features
              ))
            }
          }
        }
      }

      proof <- prove_selected_membership(
        object = object,
        selected_layers = selected_layers
      )
      if (is.null(x = proof)) {
        return(fallback())
      }
      selected_layer_data <- vector(
        mode = "list",
        length = length(x = selected_layers)
      )
      names(selected_layer_data) <- selected_layers
      for (i in seq_along(along.with = selected_layers)) {
        mat <- tryCatch(
          expr = layer_data(object = object, layer = selected_layers[[i]]),
          error = function(e) NULL
        )
        if (is.null(x = mat) ||
            !dgc_contract_ok(
              x = mat,
              expected_features = proof_layer_features(proof = proof, i = i),
              expected_cells = proof$layer_cells[[i]]
            )) {
          return(fallback())
        }
        selected_layer_data[[i]] <- mat
      }
      joined <- stitch_selected_layers(
        mats = selected_layer_data,
        proof = proof,
        selected_layers = selected_layers
      )
      value <- if (is.null(x = joined)) NULL else prepare_stitched_value(
        joined = joined,
        proof = proof
      )
      if (is.null(x = value) || !has_exact_dgc_marker(value)) {
        return(fallback())
      }
      return(finalize_fast_join(
        object = object,
        layer_store = layer_store,
        selected_layers = selected_layers,
        new = new,
        value = value,
        proof = proof,
        var.features = var.features
      ))
    }

    if (!requireNamespace("BPCells", quietly = TRUE)) {
      return(fallback())
    }
    matrix_type <- getFromNamespace("matrix_type", "BPCells")
    proof <- prove_selected_membership(object = object, selected_layers = selected_layers)
    if (is.null(x = proof)) {
      return(fallback())
    }

    selected_layer_data <- vector(mode = "list", length = length(x = selected_layers))
    names(selected_layer_data) <- selected_layers
    for (i in seq_along(along.with = selected_raw_layers)) {
      mat <- tryCatch(
        expr = layer_data(object = object, layer = selected_layers[[i]]),
        error = function(e) NULL
      )
      if (is.null(x = mat) ||
          !bpcells_layer_contract_ok(
        x = mat,
        expected_features = proof_layer_features(proof = proof, i = i),
        expected_cells = proof$layer_cells[[i]]
      )) {
        return(fallback())
      }
      selected_layer_data[[i]] <- mat
    }

    input_types <- tryCatch(
      expr = vapply(
        X = selected_layer_data,
        FUN = matrix_type,
        FUN.VALUE = character(1L)
      ),
      error = function(e) character()
    )
    if (length(x = input_types) != length(x = selected_layer_data) ||
        anyNA(input_types) ||
        any(!nzchar(input_types)) ||
        length(x = unique(x = input_types)) != 1L) {
      return(fallback())
    }
    expected_output_type <- if (!isTRUE(proof$all_layers_have_selected_features)) {
      "double"
    } else {
      input_types[[1L]]
    }
    joined <- public_stitch_iterable(
      mats = selected_layer_data,
      proof = proof
    )
    if (is.null(x = joined)) {
      return(fallback())
    }
    value <- prepare_stitched_value(
      joined = joined,
      proof = proof
    )
    output_type <- tryCatch(
      expr = if (is.null(x = value)) character() else matrix_type(value),
      error = function(e) character()
    )
    if (!identical(x = output_type, y = expected_output_type)) {
      return(fallback())
    }
    finalize_fast_join(
      object = object,
      layer_store = layer_store,
      selected_layers = selected_layers,
      new = new,
      value = value,
      proof = proof,
      var.features = var.features
    )
  }
})

.seurat_joinlayers_smoke_load <- function(task_dir, tier) {
  task <- yaml::read_yaml(file.path(task_dir, "task.yaml"))
  datasets <- Filter(function(x) identical(x$tier, tier), task$datasets)
  if (length(datasets) != 1L) {
    stop("expected one dataset for tier '", tier, "'")
  }
  path <- datasets[[1L]]$path
  if (startsWith(path, "./")) {
    path <- file.path(task_dir, sub("^\\./", "", path))
  } else if (!startsWith(path, "/")) {
    path <- file.path(task_dir, path)
  }
  readRDS(path)
}

.seurat_joinlayers_smoke_call <- function(inputs) {
  warnings <- character()
  result <- withCallingHandlers(
    SeuratObject::JoinLayers(
      object = inputs$assay,
      layers = "counts",
      new = "counts"
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(
    assay = result,
    warnings = warnings,
    input_assay_after = inputs$assay,
    fast_route = TRUE
  )
}

.seurat_joinlayers_smoke_save <- function(result, dir, ...) {
  saveRDS(result, file.path(dir, "joinlayers_output.rds"), version = 3)
}

register_patch(
  name = "seurat_joinlayers",
  upstream = "SeuratObject",
  targets = list(JoinLayers.Assay5 = fast_join_layers_assay5),
  smoke = list(
    load = .seurat_joinlayers_smoke_load,
    call = .seurat_joinlayers_smoke_call,
    save = .seurat_joinlayers_smoke_save
  ),
  tested_against = "SeuratObject 5.4.0.9001",
  tested_upstream_versions = list(SeuratObject = "5.4.0.9001")
)
}
