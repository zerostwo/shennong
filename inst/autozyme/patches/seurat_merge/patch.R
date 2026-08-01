# SeuratObject::merge.Assay5 accelerator.
#
# Lifted from autozyme task `test_seurat_merge`; final best eeba105
# (the merge.Assay5 implementation itself was established by the audited
# b1cf81a parent).  This package patch deliberately claims one namespace
# binding only: merge.Assay5.  The task also experimented with merge.StdAssay
# and merge.Seurat, but those broader surfaces are not part of this release.

if (requireNamespace("SeuratObject", quietly = TRUE)) {
  # File-scope captures keep fallback independent of later namespace rebinding.
  .seurat_merge_orig_assay5 <-
    utils::getFromNamespace("merge.Assay5", "SeuratObject")
  .seurat_merge_LogMap <- utils::getFromNamespace("LogMap", "SeuratObject")
  .seurat_merge_EmptyDF <-
    utils::getFromNamespace("EmptyDF", "SeuratObject")
  .seurat_merge_Key <- utils::getFromNamespace("Key", "SeuratObject")
  .seurat_merge_DefaultLayer <-
    utils::getFromNamespace("DefaultLayer", "SeuratObject")
  .seurat_merge_set_DefaultLayer <-
    utils::getFromNamespace("DefaultLayer<-", "SeuratObject")
  .seurat_merge_Layers <-
    utils::getFromNamespace("Layers", "SeuratObject")
  .seurat_merge_Cells <- utils::getFromNamespace("Cells", "SeuratObject")
  .seurat_merge_Features <-
    utils::getFromNamespace("Features", "SeuratObject")
  .seurat_merge_LayerData <-
    utils::getFromNamespace("LayerData", "SeuratObject")
  .seurat_merge_Assays <-
    utils::getFromNamespace("Assays", "SeuratObject")
  .seurat_merge_DefaultAssay <-
    utils::getFromNamespace("DefaultAssay", "SeuratObject")
  .seurat_merge_Idents <-
    utils::getFromNamespace("Idents", "SeuratObject")
  .seurat_merge_Project <-
    utils::getFromNamespace("Project", "SeuratObject")

  .seurat_merge_or <- function(x, y) if (is.null(x)) y else x

  # `inherits(object, "Assay5")` is intentionally too broad here: upstream
  # constructs the result with class(x), so an untested subclass can carry
  # additional validity rules or storage semantics.  Exact plain Assay5 is
  # the only released fast-path class.
  .seurat_merge_is_plain_assay5 <- function(object) {
    cls <- class(object)
    identical(as.character(cls), "Assay5") &&
      identical(attr(cls, "package"), "SeuratObject")
  }

  .seurat_merge_valid_vector <- function(x, n) {
    is.character(x) && length(x) == n
  }

  .seurat_merge_has_extra_backend_attributes <- function(layer) {
    attribute_names <- names(attributes(layer))
    if (is.null(attribute_names)) return(FALSE)
    allowed <- if (isS4(layer)) {
      c(methods::slotNames(layer), "class")
    } else {
      c("names", "dim", "dimnames", "class", "row.names", "tsp", "levels")
    }
    length(setdiff(attribute_names, allowed)) > 0L
  }

  fast_merge_Assay5 <- function(
    x,
    y,
    labels = NULL,
    add.cell.ids = NULL,
    collapse = FALSE,
    ...,
    zyme = TRUE
  ) {
    fallback <- function() {
      .seurat_merge_orig_assay5(
        x = x,
        y = y,
        labels = labels,
        add.cell.ids = add.cell.ids,
        collapse = collapse,
        ...
      )
    }

    if (!isTRUE(zyme)) {
      return(fallback())
    }

    # Read dot names without materializing their promises.  Only the private
    # `.cell.names` value used by the accepted task implementation is ever
    # forced; unrelated dots remain lazy exactly as in upstream.
    dot_names <- ...names()
    cell_names_hits <- which(dot_names == ".cell.names")
    if (length(cell_names_hits) > 1L) {
      return(fallback())
    }
    .cell.names <- if (length(cell_names_hits) == 1L) {
      ...elt(cell_names_hits[[1L]])
    } else {
      NULL
    }

    assays <- tryCatch(c(x, y), error = function(e) NULL)
    if (is.null(assays) || !length(assays) ||
        !all(vapply(
          assays,
          .seurat_merge_is_plain_assay5,
          logical(1)
        ))) {
      return(fallback())
    }

    n_assays <- length(assays)
    if (!identical(collapse, FALSE)) {
      return(fallback())
    }
    if (!is.null(labels) &&
        !.seurat_merge_valid_vector(labels, n_assays)) {
      return(fallback())
    }
    if (!is.null(add.cell.ids) &&
        !.seurat_merge_valid_vector(add.cell.ids, n_assays)) {
      return(fallback())
    }
    if (!is.null(.cell.names) &&
        (!is.list(.cell.names) || length(.cell.names) != n_assays ||
         !all(vapply(seq_along(assays), function(i) {
           is.character(.cell.names[[i]]) &&
             length(.cell.names[[i]]) == ncol(assays[[i]]) &&
             !anyNA(.cell.names[[i]])
         }, logical(1))))) {
      return(fallback())
    }

    labels_fast <- .seurat_merge_or(
      labels,
      as.character(seq_along(assays))
    )
    labels_fast[is.na(labels_fast)] <-
      as.character(which(is.na(labels_fast)))
    add_cell_ids_fast <- add.cell.ids
    if (!is.null(add_cell_ids_fast)) {
      add_cell_ids_fast[is.na(add_cell_ids_fast)] <-
        as.character(which(is.na(add_cell_ids_fast)))
    }

    # Duplicate generated layer names have overwrite behavior that is best
    # left to the canonical setter path.  Preflight before changing local
    # assay copies so every unsupported shape reaches pristine upstream code.
    assay_layers <- lapply(assays, .seurat_merge_Layers)
    has_extra_backend_attributes <- any(vapply(assays, function(assay) {
      raw_layers <- methods::slot(assay, "layers")
      any(vapply(
        raw_layers,
        .seurat_merge_has_extra_backend_attributes,
        logical(1)
      ))
    }, logical(1)))
    if (has_extra_backend_attributes) {
      return(fallback())
    }
    generated_layer_names <- unlist(Map(
      function(layer_names, label) paste(layer_names, label, sep = "."),
      assay_layers,
      labels_fast
    ), use.names = FALSE)
    if (anyDuplicated(generated_layer_names)) {
      return(fallback())
    }

    metadata <- lapply(assays, function(assay) assay[[]])
    malformed_vf <- any(vapply(metadata, function(mf) {
      vf_names <- grep("^vf_", names(mf), value = TRUE)
      any(lengths(strsplit(vf_names, "_", fixed = TRUE)) < 4L)
    }, logical(1)))
    if (malformed_vf) {
      return(fallback())
    }

    result <- tryCatch({
      for (i in seq_along(assays)) {
        if (is.null(.cell.names) && !is.null(add_cell_ids_fast)) {
          colnames(assays[[i]]) <- paste(
            colnames(assays[[i]]),
            add_cell_ids_fast[[i]],
            sep = "."
          )
        }
      }

      assay_cells <- if (is.null(.cell.names)) {
        lapply(assays, colnames)
      } else {
        .cell.names
      }
      assay_features <- lapply(assays, rownames)

      features_vec <- assay_features[[1L]]
      if (length(assay_features) > 1L &&
          !all(vapply(
            assay_features[-1L],
            identical,
            logical(1),
            y = features_vec
          ))) {
        features_vec <- Reduce(union, assay_features)
      }

      cells_vec <- unlist(assay_cells, use.names = FALSE)
      if (anyDuplicated(cells_vec)) {
        cells_vec <- Reduce(union, assay_cells)
      }

      features_all <- .seurat_merge_LogMap(features_vec)
      combined <- methods::new(
        Class = "Assay5",
        layers = list(),
        cells = .seurat_merge_LogMap(cells_vec),
        features = features_all,
        meta.data = .seurat_merge_EmptyDF(n = nrow(features_all)),
        misc = list(),
        key = .seurat_merge_or(
          .seurat_merge_Key(x),
          character(0)
        )
      )

      default <- .seurat_merge_DefaultLayer(assays[[1L]])
      layer_count <- sum(lengths(assay_layers))
      use_batched_maps <-
        layer_count * (length(features_vec) + length(cells_vec)) >= 5e5

      if (use_batched_maps) {
        build_logmap <- function(values, entries, entry_names) {
          map <- matrix(
            FALSE,
            nrow = length(values),
            ncol = length(entry_names),
            dimnames = list(values, entry_names)
          )
          for (j in seq_along(entry_names)) {
            idx <- match(entries[[j]], values)
            if (anyNA(idx)) {
              stop("layer map contains values absent from merged map",
                   call. = FALSE)
            }
            map[idx, j] <- TRUE
          }
          methods::new("LogMap", .Data = map)
        }

        layer_values <- list()
        layer_features <- list()
        layer_cells <- list()
        for (i in seq_along(assays)) {
          for (layer in assay_layers[[i]]) {
            layer_name <- paste(layer, labels_fast[[i]], sep = ".")
            layer_values[[layer_name]] <-
              methods::slot(assays[[i]], "layers")[[layer]]
            current_cells <- .seurat_merge_Cells(assays[[i]], layer = layer)
            if (!is.null(.cell.names)) {
              cell_map <- stats::setNames(
                .cell.names[[i]],
                colnames(assays[[i]])
              )
              current_cells <- unname(cell_map[current_cells])
              if (anyNA(current_cells)) {
                stop("private cell map does not cover a layer", call. = FALSE)
              }
            }
            layer_features[[layer_name]] <-
              .seurat_merge_Features(assays[[i]], layer = layer)
            layer_cells[[layer_name]] <- current_cells
          }
        }
        methods::slot(combined, "layers") <- layer_values
        methods::slot(combined, "features") <- build_logmap(
          features_vec,
          layer_features,
          generated_layer_names
        )
        methods::slot(combined, "cells") <- build_logmap(
          cells_vec,
          layer_cells,
          generated_layer_names
        )
      } else {
        for (i in seq_along(assays)) {
          for (layer in assay_layers[[i]]) {
            layer_name <- paste(layer, labels_fast[[i]], sep = ".")
            current_cells <- .seurat_merge_Cells(assays[[i]], layer = layer)
            if (!is.null(.cell.names)) {
              cell_map <- stats::setNames(
                .cell.names[[i]],
                colnames(assays[[i]])
              )
              current_cells <- unname(cell_map[current_cells])
              if (anyNA(current_cells)) {
                stop("private cell map does not cover a layer", call. = FALSE)
              }
            }
            methods::slot(combined, "layers")[[layer_name]] <-
              methods::slot(assays[[i]], "layers")[[layer]]
            methods::slot(combined, "features")[[layer_name]] <-
              .seurat_merge_Features(assays[[i]], layer = layer)
            methods::slot(combined, "cells")[[layer_name]] <- current_cells
          }
        }
      }

      for (i in seq_along(assays)) {
        mf <- metadata[[i]]
        if (!ncol(mf)) {
          next
        }
        vf_idx <- grep("^vf_", names(mf))
        if (length(vf_idx)) {
          names(mf)[vf_idx] <- vapply(
            names(mf)[vf_idx],
            function(vf) {
              parts <- strsplit(vf, "_", fixed = TRUE)[[1L]]
              paste(
                paste(parts[1:2], collapse = "_"),
                paste(
                  paste(parts[3:(length(parts) - 1L)], collapse = "_"),
                  labels_fast[[i]],
                  sep = "."
                ),
                parts[[length(parts)]],
                sep = "_"
              )
            },
            character(1)
          )
        }
        combined[[]] <- mf
      }

      default_layers <- .seurat_merge_Layers(combined, search = default)
      layer_names <- names(methods::slot(combined, "layers"))
      if (identical(default_layers, layer_names)) {
        methods::slot(combined, "default") <- length(default_layers)
      } else {
        combined <- .seurat_merge_set_DefaultLayer(
          combined,
          value = default_layers
        )
      }
      methods::validObject(combined)
      combined
    }, error = function(e) NULL)

    if (is.null(result)) fallback() else result
  }

  # Smoke recipe for test_seurat_merge.  Loading, fixture extraction, and
  # compact input hashing are outside the timed region; `call` contains only
  # the public merge() invocation plus condition capture used by evaluate.R.
  .seurat_merge_smoke_hash <- function(object) {
    path <- tempfile("autozyme-seurat-merge-")
    on.exit(unlink(path), add = TRUE)
    saveRDS(object, path, version = 3, compress = FALSE)
    unname(tools::md5sum(path))
  }

  .seurat_merge_smoke_load <- function(task_dir, tier) {
    task <- yaml::read_yaml(file.path(task_dir, "task.yaml"))
    dataset <- Filter(function(entry) identical(entry$tier, tier), task$datasets)
    if (!length(dataset)) {
      stop("no dataset for tier '", tier, "' in task.yaml")
    }
    input <- readRDS(resolve_dataset_path(task_dir, dataset[[1L]]$path))
    list(source = input, input_before = digest::digest(input))
  }

  .seurat_merge_smoke_call <- function(inputs) {
    warnings <- character()
    messages <- character()
    error <- NULL
    source <- inputs$source
    merge_dr <- .seurat_merge_or(source$merge.dr, FALSE)
    run_call <- function() {
      if (identical(source$mode, "assay5")) {
        assay_name <- .seurat_merge_or(source$assay, "RNA")
        x <- if (inherits(source$x, "Seurat")) {
          source$x[[assay_name]]
        } else {
          source$x
        }
        y <- lapply(source$y, function(object) {
          if (inherits(object, "Seurat")) object[[assay_name]] else object
        })
        if (isTRUE(source$lazy_dots)) {
          return(base::merge(
            x = x,
            y = y,
            labels = source$add.cell.ids,
            add.cell.ids = NULL,
            collapse = FALSE,
            lazy = stop("lazy dot forced")
          ))
        }
        return(base::merge(
          x = x,
          y = y,
          labels = source$add.cell.ids,
          add.cell.ids = NULL,
          collapse = FALSE
        ))
      }
      args <- list(
        x = source$x,
        y = source$y,
        add.cell.ids = source$add.cell.ids,
        collapse = FALSE,
        merge.data = TRUE,
        merge.dr = merge_dr,
        project = source$project
      )
      if (isTRUE(source$lazy_dots)) {
        args$lazy <- quote(stop("lazy dot forced"))
      }
      do.call(base::merge, args)
    }
    value <- tryCatch(
      withCallingHandlers(
        run_call(),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        },
        message = function(m) {
          messages <<- c(messages, conditionMessage(m))
          invokeRestart("muffleMessage")
        }
      ),
      error = function(e) {
        error <<- e
        NULL
      }
    )
    list(
      value = value,
      status = if (is.null(error)) "ok" else "error",
      warnings = warnings,
      messages = messages,
      error = error,
      source = source,
      input_before = inputs$input_before
    )
  }

  .seurat_merge_has_slot <- function(object, name) {
    name %in% methods::slotNames(object)
  }

  .seurat_merge_summarize_logmap <- function(map) {
    list(
      class = class(map),
      dim = as.integer(dim(map)),
      dimnames = dimnames(map),
      data = as.matrix(map)
    )
  }

  .seurat_merge_raw_layer <- function(assay, layer) {
    if (inherits(assay, "StdAssay") &&
        .seurat_merge_has_slot(assay, "layers")) {
      layers <- methods::slot(assay, "layers")
      if (layer %in% names(layers)) return(layers[[layer]])
    }
    if (.seurat_merge_has_slot(assay, layer)) {
      return(methods::slot(assay, layer))
    }
    NULL
  }

  .seurat_merge_summarize_layer <- function(view, raw = NULL) {
    raw <- .seurat_merge_or(raw, view)
    mat <- methods::as(view, "dgCMatrix")
    list(
      backend_class = class(raw),
      accessor_class = class(view),
      accessor_attributes_digest =
        digest::digest(attributes(view)),
      raw_class = class(raw),
      raw_attributes_digest = digest::digest(attributes(raw)),
      dim = as.integer(dim(mat)),
      rownames = rownames(mat),
      colnames = colnames(mat),
      x = mat@x,
      i = mat@i,
      p = mat@p
    )
  }

  .seurat_merge_summarize_assay <- function(assay) {
    layer_names <- .seurat_merge_Layers(assay)
    layers <- stats::setNames(lapply(layer_names, function(layer) {
      .seurat_merge_summarize_layer(
        .seurat_merge_LayerData(assay, layer = layer, fast = FALSE),
        .seurat_merge_raw_layer(assay, layer)
      )
    }), layer_names)
    list(
      class = class(assay),
      features = rownames(assay),
      cells = colnames(assay),
      layer_names = layer_names,
      default_layer = if (inherits(assay, "StdAssay")) {
        .seurat_merge_DefaultLayer(assay)
      } else {
        character()
      },
      default_slot = if (.seurat_merge_has_slot(assay, "default")) {
        methods::slot(assay, "default")
      } else {
        NA_integer_
      },
      key = .seurat_merge_Key(assay),
      misc = if (.seurat_merge_has_slot(assay, "misc")) {
        methods::slot(assay, "misc")
      } else {
        list()
      },
      feature_metadata = if (.seurat_merge_has_slot(assay, "meta.data")) {
        methods::slot(assay, "meta.data")
      } else {
        data.frame(row.names = rownames(assay))
      },
      cells_logmap = if (.seurat_merge_has_slot(assay, "cells")) {
        .seurat_merge_summarize_logmap(methods::slot(assay, "cells"))
      } else {
        NULL
      },
      features_logmap = if (.seurat_merge_has_slot(assay, "features")) {
        .seurat_merge_summarize_logmap(methods::slot(assay, "features"))
      } else {
        NULL
      },
      layers = layers
    )
  }

  .seurat_merge_smoke_save <- function(result, dir, tier = "small", ...) {
    input_after <- digest::digest(result$source)
    empty <- function() {
      list(
        status = result$status,
        warnings = result$warnings,
        messages = result$messages,
        error_class = if (is.null(result$error)) character() else class(result$error),
        error_message = if (is.null(result$error)) character() else
          conditionMessage(result$error),
        input_digest_before = result$input_before,
        input_digest_after = input_after,
        input_immutable = identical(result$input_before, input_after),
        dims = integer(),
        cells = character(),
        features = character(),
        assays = character(),
        default_assay = character(),
        layer_names = list(),
        assay_features = list(),
        assay_layers = list(),
        meta_data = data.frame(),
        active_ident = character(),
        active_ident_levels = character(),
        project = character(),
        version = character(),
        misc = list(),
        valid_object = identical(result$status, "error")
      )
    }

    out <- if (!identical(result$status, "ok")) {
      empty()
    } else if (inherits(result$value, "Seurat")) {
      assay_names <- .seurat_merge_Assays(result$value)
      assay_layers <- stats::setNames(lapply(assay_names, function(assay) {
        .seurat_merge_summarize_assay(result$value[[assay]])
      }), assay_names)
      identities <- .seurat_merge_Idents(result$value)
      list(
        status = "ok",
        warnings = result$warnings,
        messages = result$messages,
        error_class = character(),
        error_message = character(),
        input_digest_before = result$input_before,
        input_digest_after = input_after,
        input_immutable = identical(result$input_before, input_after),
        dims = as.integer(dim(result$value)),
        cells = colnames(result$value),
        features = rownames(result$value),
        assays = assay_names,
        default_assay = .seurat_merge_DefaultAssay(result$value),
        layer_names = stats::setNames(
          lapply(assay_layers, `[[`, "layer_names"),
          assay_names
        ),
        assay_features = stats::setNames(
          lapply(assay_layers, `[[`, "features"),
          assay_names
        ),
        assay_layers = assay_layers,
        meta_data = result$value[[]],
        active_ident = as.character(identities),
        active_ident_levels = levels(identities),
        project = .seurat_merge_Project(result$value),
        version = as.character(methods::slot(result$value, "version")),
        misc = methods::slot(result$value, "misc"),
        valid_object = isTRUE(methods::validObject(result$value, test = TRUE))
      )
    } else {
      assay_name <- "assay_result"
      assay_layers <- stats::setNames(
        list(.seurat_merge_summarize_assay(result$value)),
        assay_name
      )
      list(
        status = "ok",
        warnings = result$warnings,
        messages = result$messages,
        error_class = character(),
        error_message = character(),
        input_digest_before = result$input_before,
        input_digest_after = input_after,
        input_immutable = identical(result$input_before, input_after),
        dims = as.integer(dim(result$value)),
        cells = colnames(result$value),
        features = rownames(result$value),
        assays = assay_name,
        default_assay = character(),
        layer_names = stats::setNames(
          lapply(assay_layers, `[[`, "layer_names"),
          assay_name
        ),
        assay_features = stats::setNames(
          lapply(assay_layers, `[[`, "features"),
          assay_name
        ),
        assay_layers = assay_layers,
        meta_data = data.frame(),
        active_ident = character(),
        active_ident_levels = character(),
        project = character(),
        version = character(),
        misc = list(),
        valid_object = isTRUE(methods::validObject(result$value, test = TRUE))
      )
    }
    saveRDS(out, file.path(dir, "merge_output.rds"))
  }

  register_patch(
    name = "seurat_merge",
    upstream = "SeuratObject",
    targets = list(merge.Assay5 = fast_merge_Assay5),
    smoke = list(
      load = .seurat_merge_smoke_load,
      call = .seurat_merge_smoke_call,
      save = .seurat_merge_smoke_save
    ),
    tested_against = "SeuratObject 5.4.0.9001",
    tested_upstream_versions = list(SeuratObject = "5.4.0.9001")
  )
}
