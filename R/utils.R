#' Check if files exist
#'
#' This function takes a vector of file paths as input and checks if each file exists.
#' If a file does not exist, the function returns a message indicating which file(s) do(es) not exist.
#'
#' @param x A vector of file paths.
#' @param stop A logical value indicating whether to stop the function if a file does not exist.
#' Default is TRUE.
#'
#' @return If stop is TRUE, the function stops and returns a message indicating which file(s) do(es) not exist.
#' If stop is FALSE, the function returns a vector of file paths that do not exist.
#'
#' @examples
#' existing <- tempfile("existing-")
#' file.create(existing)
#' missing <- tempfile("missing-")
#' sn_check_file(c(existing, missing), stop = FALSE)
#'
#' @export
#' @keywords file, existence
#' @seealso \code{\link{file.exists}}, \code{\link{stop}}
sn_check_file <- function(x, stop = TRUE) {
  not_exist <- x[!file.exists(x)]
  if (length(x = not_exist) > 0) {
    if (stop) {
      stop(
        paste0(
          "Please check the file path, the file listed below does not exist:\n"
        ),
        paste0(not_exist, collapse = "\n")
      )
    } else {
      return(not_exist)
    }
  }
}

#' Create a directory path if needed
#'
#' This helper expands glue expressions in a path string, creates the directory
#' if it does not already exist, and returns the resulting path.
#'
#' @param path A character scalar giving the directory path to create.
#'
#' @return The resulting path string.
#'
#' @examples
#' tmp_dir <- tempfile("shennong-")
#' sn_set_path(tmp_dir)
#'
#' @export
sn_set_path <- function(path) {
  path <- glue(path)
  if (!file.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  return(path)
}

check_installed_github <- function(pkg, repo, reason = NULL) {
  if (rlang::is_installed(pkg)) {
    return(invisible(TRUE))
  }

  reason <- reason %||% paste0("to use functionality that depends on `", pkg, "`.")
  stop(
    glue(
      "Package '{pkg}' is required {reason}\n",
      "Install it with:\n",
      "  remotes::install_github('{repo}')"
    ),
    call. = FALSE
  )
}

.sn_log_seurat_command <- function(object, assay = NULL, name = NULL) {
  cmd <- get("LogSeuratCommand", envir = asNamespace("SeuratObject"))(
    object = object,
    return.command = TRUE
  )
  command_name <- name %||% slot(cmd, "name")
  slot(cmd, "assay.used") <- assay %||% SeuratObject::DefaultAssay(object)
  slot(cmd, "name") <- command_name
  object[[command_name]] <- cmd
  object
}

.sn_validate_seurat_assay_layer <- function(object, assay = "RNA", layer = "counts") {
  if (!inherits(object, "Seurat")) {
    stop("Input must be a Seurat object.")
  }

  if (!assay %in% names(object@assays)) {
    stop(glue("Assay '{assay}' was not found in the Seurat object."))
  }

  matched_layers <- .sn_match_seurat_layers(object = object, assay = assay, layer = layer)
  if (length(matched_layers) == 0) {
    stop(glue("Layer '{layer}' was not found in assay '{assay}'."))
  }

  invisible(TRUE)
}

.sn_match_seurat_layers <- function(object, assay = "RNA", layer = "counts") {
  assay_layers <- SeuratObject::Layers(object[[assay]])
  if (layer %in% assay_layers) {
    return(layer)
  }

  grep(pattern = paste0("^", layer, "(\\.|$)"), x = assay_layers, value = TRUE)
}

.sn_combine_seurat_layers <- function(object, assay = "RNA", layers) {
  mats <- lapply(layers, function(current_layer) {
    SeuratObject::LayerData(object = object, assay = assay, layer = current_layer)
  })

  feature_names <- rownames(mats[[1]])
  combined <- Matrix::Matrix(
    0,
    nrow = length(feature_names),
    ncol = ncol(object),
    sparse = TRUE,
    dimnames = list(feature_names, colnames(object))
  )

  for (mat in mats) {
    combined[rownames(mat), colnames(mat)] <- mat
  }

  combined
}

.sn_get_seurat_layer_data <- function(object, assay = "RNA", layer = "counts") {
  .sn_validate_seurat_assay_layer(object = object, assay = assay, layer = layer)
  matched_layers <- .sn_match_seurat_layers(object = object, assay = assay, layer = layer)

  if (length(matched_layers) == 1) {
    return(SeuratObject::LayerData(object = object, assay = assay, layer = matched_layers))
  }

  .sn_combine_seurat_layers(object = object, assay = assay, layers = matched_layers)
}

.sn_prepare_seurat_analysis_input <- function(object, assay = "RNA", layer = "counts") {
  .sn_validate_seurat_assay_layer(object = object, assay = assay, layer = layer)

  original_assay <- SeuratObject::DefaultAssay(object)
  original_counts <- NULL
  assay_layers <- SeuratObject::Layers(object[[assay]])
  has_exact_counts <- "counts" %in% assay_layers
  analysis_counts <- .sn_get_seurat_layer_data(object = object, assay = assay, layer = layer)
  needs_temp_counts <- !identical(layer, "counts") || !has_exact_counts || length(.sn_match_seurat_layers(object, assay, layer)) > 1

  if (has_exact_counts) {
    original_counts <- SeuratObject::LayerData(object = object, assay = assay, layer = "counts")
  }

  if (needs_temp_counts) {
    SeuratObject::LayerData(object = object, assay = assay, layer = "counts") <- analysis_counts
  }

  SeuratObject::DefaultAssay(object) <- assay

  list(
    object = object,
    context = list(
      original_assay = original_assay,
      analysis_assay = assay,
      original_counts = original_counts,
      had_exact_counts = has_exact_counts,
      needs_temp_counts = needs_temp_counts
    )
  )
}

.sn_restore_seurat_analysis_input <- function(object, context) {
  if (isTRUE(context$needs_temp_counts)) {
    if (isTRUE(context$had_exact_counts)) {
      SeuratObject::LayerData(
        object = object,
        assay = context$analysis_assay,
        layer = "counts"
      ) <- context$original_counts
    } else {
      SeuratObject::LayerData(
        object = object,
        assay = context$analysis_assay,
        layer = "counts"
      ) <- NULL
    }
  }

  SeuratObject::DefaultAssay(object) <- context$original_assay
  object
}
