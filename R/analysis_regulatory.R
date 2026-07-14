.sn_dorothea_network <- function(species, confidence_levels = c("A", "B", "C")) {
  check_installed("dorothea", reason = "to use the bundled DoRothEA regulons.")
  data_name <- if (identical(species, "human")) "dorothea_hs" else "dorothea_mm"
  env <- new.env(parent = emptyenv())
  utils::data(list = data_name, package = "dorothea", envir = env)
  network <- get(data_name, envir = env)
  if ("confidence" %in% colnames(network)) {
    network <- network[network$confidence %in% confidence_levels, , drop = FALSE]
  }
  network
}

.sn_progeny_network <- function(species, top = 500) {
  check_installed("progeny", reason = "to use the PROGENy pathway model.")
  model_fun <- get("getModel", envir = asNamespace("progeny"), inherits = FALSE)
  model_fun(organism = if (identical(species, "human")) "Human" else "Mouse", top = top)
}

.sn_normalize_regulatory_network <- function(network, method) {
  network <- tibble::as_tibble(network)
  if (identical(method, "dorothea")) {
    source_col <- intersect(c("tf", "source"), colnames(network))[[1]] %||% NULL
    target_col <- intersect(c("target", "gene"), colnames(network))[[1]] %||% NULL
    mor_col <- intersect(c("mor", "mor_tf", "weight"), colnames(network))[[1]] %||% NULL
  } else {
    source_col <- intersect(c("pathway", "source"), colnames(network))[[1]] %||% NULL
    target_col <- intersect(c("target", "gene"), colnames(network))[[1]] %||% NULL
    mor_col <- intersect(c("weight", "mor"), colnames(network))[[1]] %||% NULL
  }
  if (is.null(source_col) || is.null(target_col)) {
    stop("`network` must contain source and target columns.", call. = FALSE)
  }
  if (is.null(mor_col)) {
    network$mor <- 1
    mor_col <- "mor"
  }
  out <- data.frame(
    source = as.character(network[[source_col]]),
    target = as.character(network[[target_col]]),
    mor = as.numeric(network[[mor_col]]),
    stringsAsFactors = FALSE
  )
  out <- out[!is.na(out$source) & !is.na(out$target) & !is.na(out$mor), , drop = FALSE]
  unique(out)
}

#' Infer transcription-factor or pathway activity
#'
#' \code{sn_run_regulatory_activity()} runs fast footprint-style activity
#' inference using DoRothEA regulons or PROGENy pathway models through
#' \pkg{decoupleR}. DoRothEA returns transcription-factor activities; PROGENy
#' returns pathway activities.
#'
#' @param object A Seurat object.
#' @param method One of \code{"dorothea"} or \code{"progeny"}.
#' @param assay,layer Assay and layer used to retrieve expression.
#' @param group_by Optional metadata column. When supplied, expression is
#'   averaged by group before activity inference.
#' @param species Species used to choose built-in resources.
#' @param network Optional user-supplied regulatory network. If omitted,
#'   DoRothEA uses the \pkg{dorothea} package data and PROGENy uses
#'   \pkg{progeny} or \pkg{decoupleR} resources when installed.
#' @param confidence_levels DoRothEA confidence levels to keep.
#' @param progeny_top Number of target genes per PROGENy pathway.
#' @param minsize Minimum target-set size passed to \code{decoupleR::run_ulm()}.
#' @param store_name Name used under
#'   \code{object@misc$regulatory_activity_results}.
#' @param return_object If \code{TRUE}, return the updated Seurat object.
#' @param ... Additional arguments passed to \code{decoupleR::run_ulm()}.
#'
#' @return A Seurat object or stored-result list.
#' @export
sn_run_regulatory_activity <- function(object,
                                       method = c("dorothea", "progeny"),
                                       assay = NULL,
                                       layer = "data",
                                       group_by = NULL,
                                       species = NULL,
                                       network = NULL,
                                       confidence_levels = c("A", "B", "C"),
                                       progeny_top = 500,
                                       minsize = 5L,
                                       store_name = "default",
                                       return_object = TRUE,
                                       ...) {
  .sn_validate_seurat_object(object)
  check_installed("decoupleR", reason = "to run fast DoRothEA/PROGENy activity inference.")
  method <- match.arg(method)
  species <- .sn_resolve_species(object, species)
  mat <- .sn_expression_matrix(object = object, assay = assay, layer = layer)
  if (!is.null(group_by)) {
    if (!group_by %in% colnames(object[[]])) {
      stop("`group_by` must identify a metadata column in `object`.", call. = FALSE)
    }
    mat <- .sn_group_average_matrix(mat, object[[group_by, drop = TRUE]][colnames(mat)])
  }

  network <- network %||% if (identical(method, "dorothea")) {
    .sn_dorothea_network(species, confidence_levels = confidence_levels)
  } else {
    .sn_progeny_network(species, top = progeny_top)
  }
  network <- .sn_normalize_regulatory_network(network, method = method)
  source <- target <- mor <- NULL

  table <- decoupleR::run_ulm(
    mat = mat,
    network = network,
    .source = source,
    .target = target,
    .mor = mor,
    minsize = minsize,
    ...
  ) |>
    tibble::as_tibble()
  table$analysis_type <- if (identical(method, "dorothea")) "transcription_factor" else "pathway"

  sn_store_regulatory_activity(
    object = object,
    result = table,
    store_name = store_name,
    method = method,
    group_by = group_by,
    species = species,
    network = network,
    return_object = return_object
  )
}

#' Store regulatory activity results on a Seurat object
#'
#' @param object A Seurat object.
#' @param result Regulatory activity table.
#' @param store_name Name used under
#'   \code{object@misc$regulatory_activity_results}.
#' @param method Inference method.
#' @param group_by Optional grouping column.
#' @param species Optional species label.
#' @param network Optional regulatory network used for inference.
#' @param return_object If \code{TRUE}, return the updated object.
#'
#' @return A Seurat object or stored-result list.
#' @export
sn_store_regulatory_activity <- function(object,
                                         result,
                                         store_name = "default",
                                         method = "dorothea",
                                         group_by = NULL,
                                         species = NULL,
                                         network = NULL,
                                         return_object = TRUE) {
  .sn_validate_seurat_object(object)
  stored_result <- list(
    schema_version = "1.0.0",
    package_version = as.character(utils::packageVersion("Shennong")),
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    table = tibble::as_tibble(result),
    analysis = "regulatory_activity",
    method = method,
    group_by = group_by,
    species = species,
    network = network
  )
  object <- .sn_store_misc_result(
    object = object,
    collection = "regulatory_activity_results",
    store_name = store_name,
    result = stored_result
  )
  if (isTRUE(return_object)) {
    return(.sn_log_seurat_command(object = object, name = "sn_store_regulatory_activity"))
  }
  stored_result
}

#' Retrieve stored regulatory activity results
#'
#' @param object A Seurat object.
#' @param activity_name Name of the stored result.
#' @param sources Optional TF or pathway names to keep.
#' @param conditions Optional cell or group names to keep.
#' @param with_metadata If \code{TRUE}, return the full stored-result list.
#'
#' @return A tibble or stored-result list.
#' @export
sn_get_regulatory_activity_result <- function(object,
                                              activity_name = "default",
                                              sources = NULL,
                                              conditions = NULL,
                                              with_metadata = FALSE) {
  .sn_validate_seurat_object(object)
  stored <- .sn_get_misc_result(
    object = object,
    collection = "regulatory_activity_results",
    store_name = activity_name
  )
  if (isTRUE(with_metadata)) {
    return(stored)
  }
  table <- tibble::as_tibble(stored$table)
  if (!is.null(sources) && "source" %in% colnames(table)) {
    table <- dplyr::filter(table, .data$source %in% sources)
  }
  if (!is.null(conditions) && "condition" %in% colnames(table)) {
    table <- dplyr::filter(table, .data$condition %in% conditions)
  }
  table
}
