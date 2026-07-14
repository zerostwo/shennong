.sn_analysis_result_required_fields <- function() {
  c(
    "schema_version", "analysis_type", "name", "method", "backend", "input",
    "parameters", "tables", "embeddings", "graphs", "models", "diagnostics",
    "warnings", "provenance"
  )
}

.sn_analysis_provenance <- function(result = NULL, random_seed = NULL) {
  existing <- result$provenance %||% list()
  package_version <- tryCatch(
    as.character(utils::packageVersion("Shennong")),
    error = function(e) as.character(result$package_version %||% NA_character_)
  )
  existing$package_versions <- existing$package_versions %||% list(
    Shennong = package_version,
    R = paste(R.version$major, R.version$minor, sep = ".")
  )
  existing$random_seed <- existing$random_seed %||% random_seed %||%
    result$random_seed %||% NA_integer_
  existing$timestamp <- existing$timestamp %||% result$created_at %||%
    format(Sys.time(), tz = "UTC", usetz = TRUE)
  existing
}

.sn_result_tables <- function(result) {
  tables <- result$tables %||% list()
  if (!is_null(result$table) && is_null(tables$primary)) {
    tables$primary <- result$table
  }
  if (!is_null(result$overall) && is_null(tables$overall)) {
    tables$overall <- result$overall
  }
  if (!is_null(result$by_sample) && is_null(tables$by_sample)) {
    tables$by_sample <- result$by_sample
  }
  tables
}

.sn_upgrade_analysis_result <- function(result,
                                        analysis_type,
                                        name,
                                        method = NULL,
                                        backend = NULL) {
  if (!is.list(result)) {
    stop("`result` must be a list.", call. = FALSE)
  }
  analysis_type <- as.character(analysis_type %||% result$analysis_type %||% result$analysis)
  name <- as.character(name %||% result$name)
  method <- as.character(method %||% result$method %||% "unknown")
  backend <- as.character(backend %||% result$backend %||% method)

  result$schema_version <- as.character(result$schema_version %||% "1.0")
  result$analysis_type <- analysis_type
  result$name <- name
  result$method <- method
  result$backend <- backend
  result$input <- result$input %||% list()
  result$parameters <- result$parameters %||% list()
  result$tables <- .sn_result_tables(result)
  result$embeddings <- result$embeddings %||% list()
  result$graphs <- result$graphs %||% list()
  result$models <- result$models %||% list()
  result$diagnostics <- result$diagnostics %||% list()
  result$warnings <- as.character(result$warnings %||% character())
  result$provenance <- .sn_analysis_provenance(result)

  result$package_version <- result$package_version %||%
    result$provenance$package_versions$Shennong %||% NA_character_
  result$created_at <- result$created_at %||% result$provenance$timestamp
  result$analysis <- result$analysis %||% analysis_type
  if (is_null(result$table) && is.data.frame(result$tables$primary)) {
    result$table <- result$tables$primary
  }
  result
}

.sn_result_validation <- function(result) {
  errors <- character()
  warnings <- character()
  if (!is.list(result)) {
    return(list(valid = FALSE, errors = "`result` must be a list.", warnings = warnings))
  }

  required <- .sn_analysis_result_required_fields()
  missing <- setdiff(required, names(result))
  if (length(missing) > 0L) {
    errors <- c(errors, paste0("Missing required field(s): ", paste(missing, collapse = ", "), "."))
  }
  scalar_fields <- c("schema_version", "analysis_type", "name", "method", "backend")
  for (field in intersect(scalar_fields, names(result))) {
    value <- result[[field]]
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
      errors <- c(errors, paste0("`", field, "` must be a non-empty character scalar."))
    }
  }
  list_fields <- c("input", "parameters", "tables", "embeddings", "graphs", "models", "diagnostics", "provenance")
  for (field in intersect(list_fields, names(result))) {
    if (!is.list(result[[field]])) {
      errors <- c(errors, paste0("`", field, "` must be a list."))
    }
  }
  if ("warnings" %in% names(result) && !is.character(result$warnings)) {
    errors <- c(errors, "`warnings` must be a character vector.")
  }
  if (is.list(result$provenance)) {
    provenance_missing <- setdiff(c("package_versions", "random_seed", "timestamp"), names(result$provenance))
    if (length(provenance_missing) > 0L) {
      errors <- c(
        errors,
        paste0("`provenance` is missing field(s): ", paste(provenance_missing, collapse = ", "), ".")
      )
    }
  }
  if (length(result$tables %||% list()) == 0L) {
    warnings <- c(warnings, "Result contains no tables.")
  }
  list(valid = length(errors) == 0L, errors = errors, warnings = warnings)
}

#' Validate a Shennong analysis result
#'
#' @param result A result list following the Shennong analysis-result contract.
#' @param error If \code{TRUE}, stop when validation fails. If \code{FALSE},
#'   return the validation report without stopping.
#'
#' @return A validation report with \code{valid}, \code{errors}, and
#'   \code{warnings} fields. Successful reports are returned invisibly when
#'   \code{error = TRUE}.
#'
#' @examples
#' result <- list(
#'   schema_version = "1.0", analysis_type = "demo", name = "example",
#'   method = "mean", backend = "base", input = list(), parameters = list(),
#'   tables = list(primary = data.frame(value = 1)), embeddings = list(),
#'   graphs = list(), models = list(), diagnostics = list(), warnings = character(),
#'   provenance = list(package_versions = list(), random_seed = 1L, timestamp = "2026-01-01 UTC")
#' )
#' sn_validate_result(result, error = FALSE)
#'
#' @export
sn_validate_result <- function(result, error = TRUE) {
  report <- .sn_result_validation(result)
  class(report) <- c("sn_result_validation", "list")
  if (!report$valid && isTRUE(error)) {
    stop(
      "Invalid Shennong analysis result:\n- ",
      paste(report$errors, collapse = "\n- "),
      call. = FALSE
    )
  }
  if (isTRUE(error)) invisible(report) else report
}

.sn_result_collection <- function(type) {
  registry <- .sn_misc_result_registry()
  hit <- registry[registry$type == type & registry$listable, , drop = FALSE]
  if (nrow(hit) == 0L) NULL else hit$collection[[1]]
}

.sn_validate_result_object <- function(object) {
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object.", call. = FALSE)
  }
  invisible(TRUE)
}

.sn_prepare_result_for_collection <- function(result, type, name) {
  result <- .sn_upgrade_analysis_result(result, analysis_type = type, name = name)
  collection <- .sn_result_collection(type)
  if (!is_null(collection)) {
    entry <- .sn_misc_registry_entry(collection)
    required <- entry$required_fields[[1]]
    if ("table" %in% required && is_null(result$table)) {
      stop(
        "Result type '", type, "' requires `tables$primary` to be a data frame.",
        call. = FALSE
      )
    }
    if ("database" %in% required) {
      result$database <- result$database %||% result$parameters$database %||% "unknown"
    }
  }
  sn_validate_result(result)
  result
}

#' Store a Shennong analysis result on a Seurat object
#'
#' Registered legacy result types are stored in their established
#' \code{object@misc} collection. New result types use the generic
#' \code{object@misc$analysis_results} collection.
#'
#' @param object A \code{Seurat} object.
#' @param type Analysis type, for example \code{"trajectory"} or \code{"de"}.
#' @param name Stable name used to retrieve the result.
#' @param result A result list. Missing contract fields are filled when they can
#'   be inferred without changing the analytical content.
#'
#' @return The modified \code{Seurat} object.
#'
#' @examples
#' \dontrun{
#' obj <- sn_store_result(obj, "trajectory", "cd8_slingshot", result)
#' sn_get_result(obj, "trajectory", "cd8_slingshot")
#' }
#'
#' @export
sn_store_result <- function(object, type, name, result) {
  .sn_validate_result_object(object)
  if (!is.character(type) || length(type) != 1L || !nzchar(type)) {
    stop("`type` must be a non-empty character scalar.", call. = FALSE)
  }
  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
    stop("`name` must be a non-empty character scalar.", call. = FALSE)
  }
  type <- tolower(type)
  prepared <- .sn_prepare_result_for_collection(result, type = type, name = name)
  collection <- .sn_result_collection(type)
  if (!is_null(collection)) {
    return(.sn_store_misc_result(object, collection, name, prepared))
  }

  misc_data <- methods::slot(object, "misc")
  misc_data$analysis_results <- misc_data$analysis_results %||% list()
  misc_data$analysis_results[[type]] <- misc_data$analysis_results[[type]] %||% list()
  misc_data$analysis_results[[type]][[name]] <- prepared
  methods::slot(object, "misc") <- misc_data
  object
}

#' Retrieve a stored Shennong analysis result
#'
#' @param object A \code{Seurat} object.
#' @param type Analysis type.
#' @param name Stored result name.
#'
#' @return A validated Shennong analysis-result list.
#'
#' @examples
#' \dontrun{sn_get_result(obj, "trajectory", "cd8_slingshot")}
#'
#' @export
sn_get_result <- function(object, type, name) {
  .sn_validate_result_object(object)
  type <- tolower(as.character(type))
  name <- as.character(name)
  collection <- .sn_result_collection(type)
  if (!is_null(collection)) {
    result <- .sn_get_misc_result(object, collection, name)
  } else {
    misc_data <- methods::slot(object, "misc")
    results <- misc_data$analysis_results[[type]] %||% list()
    if (!name %in% names(results)) {
      stop(
        "No stored result named '", name, "' was found for analysis type '", type, "'.",
        call. = FALSE
      )
    }
    result <- results[[name]]
  }
  result <- .sn_upgrade_analysis_result(result, analysis_type = type, name = name)
  sn_validate_result(result)
  result
}

.sn_generic_result_summary <- function(object) {
  misc_data <- methods::slot(object, "misc")
  collections <- misc_data$analysis_results %||% list()
  if (length(collections) == 0L) {
    return(tibble::tibble())
  }
  dplyr::bind_rows(lapply(names(collections), function(type) {
    entries <- collections[[type]]
    if (length(entries) == 0L) return(tibble::tibble())
    names <- names(entries)
    tibble::tibble(
      collection = "analysis_results",
      type = type,
      name = names,
      analysis = vapply(names, function(name) entries[[name]]$analysis_type %||% type, character(1)),
      method = vapply(names, function(name) entries[[name]]$method %||% NA_character_, character(1)),
      created_at = vapply(names, function(name) entries[[name]]$provenance$timestamp %||% NA_character_, character(1)),
      n_rows = unname(vapply(names, function(name) {
        tables <- entries[[name]]$tables %||% list()
        sum(vapply(tables, function(table) if (is.data.frame(table)) nrow(table) else 0L, integer(1)))
      }, integer(1))),
      source = NA_character_
    )
  }))
}

#' Delete a stored Shennong analysis result
#'
#' @param object A \code{Seurat} object.
#' @param type Analysis type.
#' @param name Stored result name.
#'
#' @return The modified \code{Seurat} object.
#'
#' @examples
#' \dontrun{obj <- sn_delete_result(obj, "trajectory", "cd8_slingshot")}
#'
#' @export
sn_delete_result <- function(object, type, name) {
  .sn_validate_result_object(object)
  type <- tolower(as.character(type))
  name <- as.character(name)
  collection <- .sn_result_collection(type)
  misc_data <- methods::slot(object, "misc")
  if (!is_null(collection)) {
    entries <- misc_data[[collection]] %||% list()
    if (!name %in% names(entries)) {
      stop("No stored result named '", name, "' was found for analysis type '", type, "'.", call. = FALSE)
    }
    entries[[name]] <- NULL
    misc_data[[collection]] <- entries
  } else {
    entries <- misc_data$analysis_results[[type]] %||% list()
    if (!name %in% names(entries)) {
      stop("No stored result named '", name, "' was found for analysis type '", type, "'.", call. = FALSE)
    }
    entries[[name]] <- NULL
    misc_data$analysis_results[[type]] <- entries
  }
  methods::slot(object, "misc") <- misc_data
  object
}
