.sn_analysis_result_required_fields <- function() {
  c(
    "schema_version", "analysis_type", "result_id", "method", "backend", "input",
    "parameters", "tables", "embeddings", "graphs", "models", "diagnostics",
    "warnings", "provenance"
  )
}

.sn_analysis_result_schema_version <- function() "2.0.0"

.sn_validate_result_id <- function(result_id) {
  if (!is.character(result_id) || length(result_id) != 1L ||
      is.na(result_id) || !nzchar(result_id)) {
    stop("`result_id` must be a non-empty character scalar.", call. = FALSE)
  }
  if (!identical(result_id, trimws(result_id))) {
    stop("`result_id` must not start or end with whitespace.", call. = FALSE)
  }
  result_id
}

.sn_is_supported_analysis_result_schema_version <- function(version) {
  version <- as.character(version %||% NA_character_)
  if (length(version) != 1L || is.na(version) ||
      !grepl("^[0-9]+\\.[0-9]+\\.[0-9]+([+-][0-9A-Za-z.-]+)?$", version)) {
    return(FALSE)
  }
  identical(version, .sn_analysis_result_schema_version())
}

.sn_analysis_result_requires_primary <- function(analysis_type) {
  !analysis_type %in% c("interpretation")
}

.sn_analysis_provenance <- function(result = NULL,
                                    random_seed = NULL,
                                    capture_acceleration = TRUE) {
  result <- result %||% list()
  existing <- result[["provenance"]] %||% list()
  package_version <- tryCatch(
    as.character(utils::packageVersion("Shennong")),
    error = function(e) as.character(result[["package_version"]] %||% NA_character_)
  )
  existing[["package_versions"]] <- existing[["package_versions"]] %||% list(
    Shennong = package_version,
    R = paste(R.version$major, R.version$minor, sep = ".")
  )
  existing[["random_seed"]] <- existing[["random_seed"]] %||% random_seed %||%
    result[["random_seed"]] %||% NA_integer_
  existing[["timestamp"]] <- existing[["timestamp"]] %||% result[["created_at"]] %||%
    format(Sys.time(), tz = "UTC", usetz = TRUE)
  acceleration <- if (isTRUE(capture_acceleration)) {
    tryCatch(.sn_acceleration_provenance(), error = function(e) list())
  } else {
    list()
  }
  if (length(acceleration) > 0L) {
    existing[["acceleration"]] <- existing[["acceleration"]] %||% acceleration
  }
  existing
}

.sn_contextual_analysis_provenance <- function(result = NULL,
                                               random_seed = NULL) {
  .sn_analysis_provenance(
    result = result,
    random_seed = random_seed,
    capture_acceleration = is.environment(
      getOption("shennong.acceleration.provenance_context")
    )
  )
}

.sn_result_tables <- function(result) {
  tables <- result[["tables"]] %||% list()
  legacy_table <- result[["table"]]
  if (!is_null(legacy_table) && is_null(tables[["primary"]])) {
    tables[["primary"]] <- legacy_table
  }
  legacy_overall <- result[["overall"]]
  if (!is_null(legacy_overall) && is_null(tables[["overall"]])) {
    tables[["overall"]] <- legacy_overall
  }
  legacy_by_sample <- result[["by_sample"]]
  if (!is_null(legacy_by_sample) && is_null(tables[["by_sample"]])) {
    tables[["by_sample"]] <- legacy_by_sample
  }
  if (is_null(tables[["primary"]]) && is.data.frame(tables[["by_sample"]])) {
    tables[["primary"]] <- tables[["by_sample"]]
  }
  primary_alias <- switch(
    as.character(result[["analysis_type"]] %||% result[["analysis"]] %||% ""),
    annotation = tables[["cells"]],
    program_scoring = tables[["scores"]],
    trajectory = tables[["cells"]],
    NULL
  )
  if (is_null(tables[["primary"]]) && is.data.frame(primary_alias)) {
    tables[["primary"]] <- primary_alias
  }
  tables
}

.sn_upgrade_de_primary_gene <- function(tables) {
  primary <- tables[["primary"]]
  if (!is.data.frame(primary) || "gene" %in% colnames(primary)) {
    return(list(tables = tables, source = NULL))
  }

  feature_available <- "feature" %in% colnames(primary)
  row_ids <- rownames(primary)
  rownames_available <- nrow(primary) > 0L && length(row_ids) == nrow(primary) &&
    !identical(as.character(row_ids), as.character(seq_len(nrow(primary))))
  candidates <- list()
  if (feature_available) candidates$feature <- as.character(primary[["feature"]])
  if (rownames_available) candidates$rownames <- as.character(row_ids)
  if (length(candidates) == 0L) return(list(tables = tables, source = NULL))

  valid <- vapply(candidates, function(value) {
    length(value) == nrow(primary) && !anyNA(value) && all(nzchar(trimws(value)))
  }, logical(1))
  if (!all(valid)) {
    stop(
      "Cannot upgrade DE identifiers: `feature`/row names contain missing or empty values.",
      call. = FALSE
    )
  }
  if (length(candidates) > 1L &&
      !identical(unname(candidates[[1L]]), unname(candidates[[2L]]))) {
    stop(
      "Cannot upgrade DE identifiers because `feature` and row names disagree.",
      call. = FALSE
    )
  }

  source <- names(candidates)[[1L]]
  primary[["gene"]] <- candidates[[source]]
  tables[["primary"]] <- primary
  list(tables = tables, source = source)
}

.sn_new_analysis_result <- function(analysis_type,
                                    result_id,
                                    method,
                                    backend = method,
                                    input = list(),
                                    parameters = list(),
                                    tables = list(),
                                    embeddings = list(),
                                    graphs = list(),
                                    models = list(),
                                    diagnostics = list(),
                                    warnings = character(),
                                    provenance = NULL,
                                    random_seed = NA_integer_) {
  result_id <- .sn_validate_result_id(result_id)
  provenance <- provenance %||% .sn_analysis_provenance(random_seed = random_seed)
  provenance[["result_id"]] <- result_id
  provenance[["analysis_type"]] <- analysis_type
  result <- list(
    schema_version = .sn_analysis_result_schema_version(),
    analysis_type = analysis_type,
    result_id = result_id,
    method = method,
    backend = backend,
    input = input,
    parameters = parameters,
    tables = tables,
    embeddings = embeddings,
    graphs = graphs,
    models = models,
    diagnostics = diagnostics,
    warnings = as.character(warnings),
    provenance = provenance
  )
  sn_validate_result(result)
  result
}

.sn_upgrade_analysis_result <- function(result,
                                        analysis_type,
                                        result_id,
                                        method = NULL,
                                        backend = NULL) {
  if (!is.list(result)) {
    stop("`result` must be a list.", call. = FALSE)
  }
  analysis_type <- as.character(
    analysis_type %||% result[["analysis_type"]] %||% result[["analysis"]]
  )
  result_id <- .sn_validate_result_id(result_id)
  method <- as.character(method %||% result[["method"]] %||% "unknown")
  backend <- as.character(backend %||% result[["backend"]] %||% method)

  source_schema_version <- result[["schema_version"]]
  if (!is_null(source_schema_version)) {
    valid_source_version <- is.character(source_schema_version) &&
      length(source_schema_version) == 1L && !is.na(source_schema_version)
    if (!valid_source_version) {
      stop(
        "Cannot upgrade result with invalid `schema_version` '",
        paste(as.character(source_schema_version), collapse = ", "), "'.",
        call. = FALSE
      )
    }
    migratable_versions <- c("1", "1.0", "1.0.0", .sn_analysis_result_schema_version())
    semantic_version <- grepl(
      "^[0-9]+\\.[0-9]+\\.[0-9]+([+-][0-9A-Za-z.-]+)?$",
      source_schema_version
    )
    if (semantic_version &&
        utils::compareVersion(
          sub("[+-].*$", "", source_schema_version),
          .sn_analysis_result_schema_version()
        ) > 0L) {
      stop(
        "Cannot read future `schema_version` '", source_schema_version,
        "' with Shennong result schema '", .sn_analysis_result_schema_version(),
        "'. Upgrade Shennong before reading or storing this result.",
        call. = FALSE
      )
    }
    if (!source_schema_version %in% migratable_versions) {
      stop(
        "Cannot safely upgrade unsupported `schema_version` '",
        source_schema_version, "'. Migratable versions are: ",
        paste(migratable_versions, collapse = ", "), ".",
        call. = FALSE
      )
    }
  }

  result[["schema_version"]] <- .sn_analysis_result_schema_version()
  result[["analysis_type"]] <- analysis_type
  result[["name"]] <- NULL
  result[["result_id"]] <- result_id
  result[["method"]] <- method
  result[["backend"]] <- backend
  result[["input"]] <- result[["input"]] %||% list()
  result[["parameters"]] <- result[["parameters"]] %||% list()
  result[["tables"]] <- .sn_result_tables(result)
  de_gene_migration <- list(tables = result[["tables"]], source = NULL)
  if (identical(analysis_type, "de")) {
    de_gene_migration <- .sn_upgrade_de_primary_gene(result[["tables"]])
    result[["tables"]] <- de_gene_migration$tables
  }
  result[["embeddings"]] <- result[["embeddings"]] %||% list()
  result[["graphs"]] <- result[["graphs"]] %||% list()
  result[["models"]] <- result[["models"]] %||% list()
  result[["diagnostics"]] <- result[["diagnostics"]] %||% list()
  result[["warnings"]] <- as.character(result[["warnings"]] %||% character())
  result[["provenance"]] <- .sn_analysis_provenance(
    result,
    capture_acceleration = FALSE
  )
  result[["provenance"]][["result_id"]] <- result_id
  result[["provenance"]][["analysis_type"]] <- analysis_type
  if (!is_null(source_schema_version) &&
      !identical(source_schema_version, .sn_analysis_result_schema_version())) {
    result[["provenance"]][["migrated_from_schema_version"]] <-
      source_schema_version
  }
  if (!is_null(de_gene_migration$source)) {
    result[["provenance"]][["migrated_primary_gene_from"]] <-
      de_gene_migration$source
  }

  result[["table"]] <- NULL
  result[["overall"]] <- NULL
  result[["by_sample"]] <- NULL
  result[["package_version"]] <- NULL
  result[["created_at"]] <- NULL
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
  scalar_fields <- c("schema_version", "analysis_type", "result_id", "method", "backend")
  for (field in intersect(scalar_fields, names(result))) {
    value <- result[[field]]
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
      errors <- c(errors, paste0("`", field, "` must be a non-empty character scalar."))
    }
  }
  schema_version <- result[["schema_version"]] %||% NA_character_
  if (is.character(schema_version) && length(schema_version) == 1L &&
             !is.na(schema_version) && nzchar(schema_version) &&
             !grepl("^[0-9]+\\.[0-9]+\\.[0-9]+([+-][0-9A-Za-z.-]+)?$", schema_version)) {
    errors <- c(errors, "`schema_version` must use semantic version form such as '2.0.0'.")
  }
  if (is.character(schema_version) &&
      length(schema_version) == 1L && !is.na(schema_version) && nzchar(schema_version) &&
      grepl("^[0-9]+\\.[0-9]+\\.[0-9]+([+-][0-9A-Za-z.-]+)?$", schema_version) &&
      !.sn_is_supported_analysis_result_schema_version(schema_version)) {
    errors <- c(
      errors,
      paste0(
        "Unsupported `schema_version` '", schema_version,
        "'; expected '", .sn_analysis_result_schema_version(), "'."
      )
    )
  }
  list_fields <- c("input", "parameters", "tables", "embeddings", "graphs", "models", "diagnostics", "provenance")
  for (field in intersect(list_fields, names(result))) {
    if (!is.list(result[[field]])) {
      errors <- c(errors, paste0("`", field, "` must be a list."))
    }
  }
  if ("warnings" %in% names(result) && !is.character(result[["warnings"]])) {
    errors <- c(errors, "`warnings` must be a character vector.")
  }
  provenance <- result[["provenance"]]
  if (is.list(provenance)) {
    provenance_missing <- setdiff(
      c("package_versions", "random_seed", "timestamp"),
      names(provenance)
    )
    if (length(provenance_missing) > 0L) {
      errors <- c(
        errors,
        paste0("`provenance` is missing field(s): ", paste(provenance_missing, collapse = ", "), ".")
      )
    }
    for (field in c("result_id", "analysis_type")) {
      value <- provenance[[field]]
      if (!is_null(value) && !identical(value, result[[field]])) {
        errors <- c(errors, paste0("`provenance$", field, "` must match `", field, "`."))
      }
    }
    package_versions <- provenance[["package_versions"]]
    if (!is_null(package_versions)) {
      valid_version_names <- is.list(package_versions) && (
        length(package_versions) == 0L ||
          (!is_null(names(package_versions)) && !anyNA(names(package_versions)) &&
             all(nzchar(names(package_versions))) && !anyDuplicated(names(package_versions)))
      )
      valid_version_values <- is.list(package_versions) && all(vapply(
        package_versions,
        function(value) is.character(value) && length(value) == 1L &&
          (is.na(value) || nzchar(value)),
        logical(1)
      ))
      if (!valid_version_names || !valid_version_values) {
        errors <- c(
          errors,
          "`provenance$package_versions` must be a named list of character scalars."
        )
      }
    }
    timestamp <- provenance[["timestamp"]]
    if (!is_null(timestamp) && (!is.character(timestamp) || length(timestamp) != 1L ||
        is.na(timestamp) || !nzchar(timestamp))) {
      errors <- c(errors, "`provenance$timestamp` must be a non-empty character scalar.")
    }
    valid_seed <- function(value) {
      if (!is.numeric(value) || length(value) != 1L || is.object(value) ||
          !is.null(dim(value))) {
        return(FALSE)
      }
      if (is.na(value)) {
        return(!is.nan(value))
      }
      is.finite(value) && value >= 0 && value <= .Machine$integer.max &&
        value == trunc(value)
    }
    random_seed <- provenance[["random_seed"]]
    if (!is_null(random_seed)) {
      valid_random_seed <- if (is.list(random_seed)) {
        (length(random_seed) == 0L ||
          (!is_null(names(random_seed)) && !anyNA(names(random_seed)) &&
             all(nzchar(names(random_seed))) && !anyDuplicated(names(random_seed)))) &&
          all(vapply(random_seed, valid_seed, logical(1)))
      } else {
        valid_seed(random_seed)
      }
      if (!valid_random_seed) {
        errors <- c(
          errors,
          paste0(
            "`provenance$random_seed` must be one non-negative integer/NA or ",
            "a named list of such values."
          )
        )
      }
    }
  }
  analysis_type <- result[["analysis_type"]] %||% ""
  analysis_type <- if (is.character(analysis_type) && length(analysis_type) == 1L &&
                       !is.na(analysis_type)) analysis_type else ""
  result_tables <- result[["tables"]]
  tables <- if (is.list(result_tables)) result_tables else list()
  primary <- tables[["primary"]]
  if (length(tables) == 0L) {
    if (.sn_analysis_result_requires_primary(analysis_type)) {
      errors <- c(errors, "Table-producing analysis results require `tables$primary`.")
    } else {
      warnings <- c(warnings, "Result contains no tables.")
    }
  } else if (.sn_analysis_result_requires_primary(analysis_type) &&
             !is.data.frame(primary)) {
    errors <- c(errors, "`tables$primary` must be a data frame or tibble.")
  }
  if (is.data.frame(primary)) {
    if (ncol(primary) == 0L) {
      errors <- c(errors, "`tables$primary` must contain at least one column.")
    }
    if (anyDuplicated(names(primary))) {
      errors <- c(errors, "`tables$primary` must not contain duplicate column names.")
    }
  }
  spec <- .sn_analysis_result_type_specs()[[analysis_type]]
  if (!is_null(spec) && is.data.frame(primary)) {
    missing_primary_columns <- setdiff(spec$primary_columns, names(primary))
    if (length(missing_primary_columns) > 0L) {
      errors <- c(
        errors,
        paste0(
          "`tables$primary` for analysis type '", analysis_type,
          "' is missing column(s): ", paste(missing_primary_columns, collapse = ", "), "."
        )
      )
    }
    errors <- c(
      errors,
      .sn_result_primary_semantic_errors(primary, analysis_type, spec)
    )
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
#'   schema_version = "2.0.0", analysis_type = "demo", result_id = "example",
#'   method = "mean", backend = "base", input = list(), parameters = list(),
#'   tables = list(primary = data.frame(value = 1)), embeddings = list(),
#'   graphs = list(), models = list(), diagnostics = list(), warnings = character(),
#'   provenance = list(package_versions = list(), random_seed = 1L,
#'     timestamp = "2026-01-01 UTC", result_id = "example", analysis_type = "demo")
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

.sn_prepare_result <- function(result, type, result_id) {
  result <- .sn_upgrade_analysis_result(
    result,
    analysis_type = type,
    result_id = result_id
  )
  sn_validate_result(result)
  result
}

.sn_result_store <- function(object) {
  misc_data <- methods::slot(object, "misc")
  shennong <- misc_data[["shennong"]] %||% list()
  shennong[["results"]] %||% list()
}

.sn_stored_result_identity_errors <- function(result, type, result_id) {
  errors <- character()
  if (!is.list(result)) return("Stored result must be a list.")
  if (!identical(result[["analysis_type"]], type)) {
    errors <- c(
      errors,
      paste0(
        "Stored result analysis type '", result[["analysis_type"]] %||% "<missing>",
        "' does not match its storage type '", type, "'."
      )
    )
  }
  if (!identical(result[["result_id"]], result_id)) {
    errors <- c(
      errors,
      paste0(
        "Stored result ID '", result[["result_id"]] %||% "<missing>",
        "' does not match its storage key '", result_id, "'."
      )
    )
  }
  errors
}

#' Store a Shennong analysis result on a Seurat object
#'
#' Every result is stored at
#' \code{object@misc$shennong$results[[analysis_type]][[result_id]]}.
#'
#' @param object A \code{Seurat} object.
#' @param type Analysis type, for example \code{"trajectory"} or \code{"de"}.
#' @param result_id Stable identifier used to store and retrieve the result.
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
sn_store_result <- function(object, type, result_id, result) {
  .sn_validate_seurat_object(object)
  result_id <- .sn_validate_result_id(result_id)
  if (!is.character(type) || length(type) != 1L || !nzchar(type)) {
    stop("`type` must be a non-empty character scalar.", call. = FALSE)
  }
  type <- tolower(type)
  artifact_types <- .sn_misc_result_registry() |>
    dplyr::filter(.data$contract_scope == "artifact") |>
    dplyr::pull(.data$type)
  if (grepl("_artifact$", type) || type %in% artifact_types) {
    stop(
      "`type` is reserved for a registered workflow artifact and cannot be stored as a unified analysis result.",
      call. = FALSE
    )
  }
  prepared <- .sn_prepare_result(result, type = type, result_id = result_id)
  misc_data <- methods::slot(object, "misc")
  misc_data[["shennong"]] <- misc_data[["shennong"]] %||% list()
  misc_data[["shennong"]][["results"]] <-
    misc_data[["shennong"]][["results"]] %||% list()
  misc_data[["shennong"]][["results"]][[type]] <-
    misc_data[["shennong"]][["results"]][[type]] %||% list()
  misc_data[["shennong"]][["results"]][[type]][[result_id]] <- prepared
  methods::slot(object, "misc") <- misc_data
  object
}

#' Retrieve a stored Shennong analysis result
#'
#' @param object A \code{Seurat} object.
#' @param type Analysis type.
#' @param result_id Stored result identifier.
#'
#' @return A validated Shennong analysis-result list.
#'
#' @examples
#' \dontrun{sn_get_result(obj, "trajectory", "cd8_slingshot")}
#'
#' @export
sn_get_result <- function(object, type, result_id) {
  .sn_validate_seurat_object(object)
  type <- tolower(as.character(type))
  result_id <- .sn_validate_result_id(result_id)
  results <- .sn_result_store(object)[[type]] %||% list()
  if (!result_id %in% names(results)) {
    stop(
      "No result with `result_id = \"", result_id,
      "\"` was found for analysis type '", type, "'.",
      call. = FALSE
    )
  }
  result <- results[[result_id]]
  sn_validate_result(result)
  identity_errors <- .sn_stored_result_identity_errors(result, type, result_id)
  if (length(identity_errors) > 0L) {
    stop(
      "Invalid stored Shennong analysis result identity:\n- ",
      paste(identity_errors, collapse = "\n- "),
      call. = FALSE
    )
  }
  result
}

.sn_generic_result_summary <- function(object) {
  collections <- .sn_result_store(object)
  if (length(collections) == 0L) {
    return(tibble::tibble())
  }
  dplyr::bind_rows(lapply(names(collections), function(type) {
    entries <- collections[[type]]
    if (length(entries) == 0L) return(tibble::tibble())
    result_ids <- names(entries)
    tibble::tibble(
      collection = "shennong.results",
      type = type,
      result_id = result_ids,
      analysis = vapply(result_ids, function(result_id) {
        entries[[result_id]][["analysis_type"]] %||% type
      }, character(1)),
      method = vapply(result_ids, function(result_id) {
        entries[[result_id]][["method"]] %||% NA_character_
      }, character(1)),
      created_at = vapply(result_ids, function(result_id) {
        provenance <- entries[[result_id]][["provenance"]] %||% list()
        provenance[["timestamp"]] %||% NA_character_
      }, character(1)),
      n_rows = unname(vapply(result_ids, function(result_id) {
        tables <- entries[[result_id]][["tables"]] %||% list()
        primary <- tables[["primary"]]
        if (is.data.frame(primary)) nrow(primary) else 0L
      }, integer(1))),
      source = NA_character_
    )
  }))
}

.sn_stored_analysis_result_entries <- function(object, type = NULL,
                                               include_artifacts = FALSE) {
  misc_data <- methods::slot(object, "misc")
  registry <- .sn_misc_result_registry()
  registry <- registry[registry$contract_scope == "artifact", , drop = FALSE]
  entries <- list()
  if (isTRUE(include_artifacts)) {
    for (index in seq_len(nrow(registry))) {
      current_type <- registry$type[[index]]
      if (!is_null(type) && !current_type %in% type) next
      collection <- registry$collection[[index]]
      collection_data <- misc_data[[collection]] %||% list()
      if (length(collection_data) > 0L) {
        entries[[length(entries) + 1L]] <- list(
          collection = collection,
          type = current_type,
          result_id = collection,
          contract_scope = "artifact",
          result = collection_data
        )
      }
    }
  }
  stored_results <- .sn_result_store(object)
  for (current_type in names(stored_results)) {
    if (!is_null(type) && !current_type %in% type) next
    for (result_id in names(stored_results[[current_type]])) {
      entries[[length(entries) + 1L]] <- list(
        collection = "shennong.results",
        type = current_type,
        result_id = result_id,
        contract_scope = "analysis_result",
        result = stored_results[[current_type]][[result_id]]
      )
    }
  }
  if (isTRUE(include_artifacts)) {
    registered_collections <- c(registry$collection, "shennong")
    unregistered_collections <- setdiff(names(misc_data), registered_collections)
    for (collection in unregistered_collections) {
      collection_data <- misc_data[[collection]]
      if (length(collection_data) == 0L) next
      current_type <- paste0(collection, "_unregistered")
      if (!is_null(type) && !current_type %in% type && !collection %in% type) next
      entries[[length(entries) + 1L]] <- list(
        collection = collection,
        type = current_type,
        result_id = collection,
        contract_scope = "unregistered",
        result = collection_data
      )
    }
  }
  entries
}

#' Audit stored Shennong analysis results
#'
#' Inspect every registered analysis result and artifact, plus unknown populated
#' top-level \code{object@misc} entries, without mutating the object. The audit
#' distinguishes results that satisfy the current contract from malformed
#' results that can be normalized safely, and reports unknown payloads as
#' \code{unregistered}.
#'
#' @param object A \code{Seurat} object.
#' @param type Optional analysis type or character vector of types to inspect.
#' @param include_artifacts If \code{TRUE}, also report registered runtime,
#'   cache, and backend artifacts that intentionally remain outside the
#'   tabular analysis-result contract, as well as unknown populated top-level
#'   \code{object@misc} entries.
#'
#' @return A tibble with one row per stored result and validation, migration,
#'   schema-version, and canonical-primary-table status.
#'
#' @examples
#' \dontrun{sn_audit_results(object)}
#'
#' @export
sn_audit_results <- function(object, type = NULL, include_artifacts = TRUE) {
  .sn_validate_seurat_object(object)
  requested_types <- if (is_null(type)) NULL else tolower(as.character(type))
  entries <- .sn_stored_analysis_result_entries(
    object,
    requested_types,
    include_artifacts = include_artifacts
  )
  if (length(entries) == 0L) {
    return(tibble::tibble(
      collection = character(), type = character(), result_id = character(),
      contract_scope = character(),
      schema_version = character(), target_schema_version = character(),
      status = character(), unified = logical(), valid = logical(),
      repairable = logical(), primary_rows = integer(),
      errors = character(), upgrade_errors = character(), warnings = character()
    ))
  }
  dplyr::bind_rows(lapply(entries, function(entry) {
    raw <- entry$result
    if (identical(entry$contract_scope, "unregistered")) {
      return(tibble::tibble(
        collection = entry$collection,
        type = entry$type,
        result_id = entry$result_id,
        contract_scope = "unregistered",
        schema_version = NA_character_,
        target_schema_version = NA_character_,
        status = "unregistered",
        unified = NA,
        valid = NA,
        repairable = FALSE,
        primary_rows = 0L,
        errors = "",
        upgrade_errors = "",
        warnings = paste0(
          "Top-level `object@misc$", entry$collection,
          "` is not registered as a Shennong analysis result or artifact."
        )
      ))
    }
    if (identical(entry$contract_scope, "artifact")) {
      raw_schema <- if (is.list(raw)) raw[["schema_version"]] %||% NA_character_ else NA_character_
      return(tibble::tibble(
        collection = entry$collection,
        type = entry$type,
        result_id = entry$result_id,
        contract_scope = "artifact",
        schema_version = as.character(raw_schema),
        target_schema_version = NA_character_,
        status = "artifact",
        unified = NA,
        valid = NA,
        repairable = FALSE,
        primary_rows = 0L,
        errors = "",
        upgrade_errors = "",
        warnings = paste0(
          "Registered runtime/cache artifact; intentionally outside the ",
          "analysis-result contract."
        )
      ))
    }
    raw_report <- sn_validate_result(raw, error = FALSE)
    identity_errors <- .sn_stored_result_identity_errors(
      raw,
      entry$type,
      entry$result_id
    )
    if (length(identity_errors) > 0L) {
      raw_report$valid <- FALSE
      raw_report$errors <- c(raw_report$errors, identity_errors)
    }
    upgraded <- tryCatch(
      .sn_prepare_result(raw, type = entry$type, result_id = entry$result_id),
      error = identity
    )
    repairable <- !inherits(upgraded, "error")
    upgraded_report <- if (repairable) {
      sn_validate_result(upgraded, error = FALSE)
    } else {
      list(valid = FALSE, errors = conditionMessage(upgraded), warnings = character())
    }
    raw_schema <- if (is.list(raw)) raw[["schema_version"]] %||% NA_character_ else NA_character_
    raw_tables <- if (is.list(raw)) raw[["tables"]] else NULL
    raw_primary <- if (is.list(raw_tables)) raw_tables[["primary"]] else NULL
    canonical_version <- .sn_is_supported_analysis_result_schema_version(raw_schema)
    unified <- repairable && isTRUE(upgraded_report$valid) &&
      isTRUE(raw_report$valid) && canonical_version &&
      (!.sn_analysis_result_requires_primary(entry$type) || is.data.frame(raw_primary))
    status <- if (unified) "valid" else if (isTRUE(upgraded_report$valid)) "repairable" else "invalid"
    primary <- if (repairable) upgraded[["tables"]][["primary"]] else NULL
    tibble::tibble(
      collection = entry$collection,
      type = entry$type,
      result_id = entry$result_id,
      contract_scope = "analysis_result",
      schema_version = as.character(raw_schema),
      target_schema_version = .sn_analysis_result_schema_version(),
      status = status,
      unified = unified,
      valid = isTRUE(raw_report$valid),
      repairable = isTRUE(upgraded_report$valid),
      primary_rows = if (is.data.frame(primary)) nrow(primary) else 0L,
      errors = paste(raw_report$errors %||% character(), collapse = "; "),
      upgrade_errors = paste(upgraded_report$errors %||% character(), collapse = "; "),
      warnings = paste(unique(c(raw_report$warnings, upgraded_report$warnings)), collapse = "; ")
    )
  })) |>
    dplyr::arrange(.data$collection, .data$type, .data$result_id)
}

#' Upgrade stored Shennong analysis results
#'
#' Normalize stored analysis results to the current canonical result envelope.
#'
#' @param object A \code{Seurat} object.
#' @param type Optional analysis type or character vector of types to normalize.
#' @param strict If \code{TRUE}, stop at the first result that cannot be safely
#'   normalized. If \code{FALSE}, leave invalid entries unchanged and warn.
#'
#' @return The modified \code{Seurat} object.
#'
#' @examples
#' \dontrun{object <- sn_upgrade_results(object)}
#'
#' @export
sn_upgrade_results <- function(object, type = NULL, strict = TRUE) {
  .sn_validate_seurat_object(object)
  requested_types <- if (is_null(type)) NULL else tolower(as.character(type))
  entries <- .sn_stored_analysis_result_entries(object, requested_types)
  for (entry in entries) {
    updated <- tryCatch(
      sn_store_result(object, entry$type, entry$result_id, entry$result),
      error = identity
    )
    if (inherits(updated, "error")) {
      message <- paste0(
        "Could not upgrade stored result '", entry$result_id, "' of type '",
        entry$type, "': ", conditionMessage(updated)
      )
      if (isTRUE(strict)) stop(message, call. = FALSE)
      warning(message, call. = FALSE)
    } else {
      object <- updated
    }
  }
  object
}

#' Delete a stored Shennong analysis result
#'
#' @param object A \code{Seurat} object.
#' @param type Analysis type.
#' @param result_id Stored result identifier.
#'
#' @return The modified \code{Seurat} object.
#'
#' @examples
#' \dontrun{obj <- sn_delete_result(obj, "trajectory", "cd8_slingshot")}
#'
#' @export
sn_delete_result <- function(object, type, result_id) {
  .sn_validate_seurat_object(object)
  type <- tolower(as.character(type))
  result_id <- .sn_validate_result_id(result_id)
  misc_data <- methods::slot(object, "misc")
  entries <- misc_data[["shennong"]][["results"]][[type]] %||% list()
  if (!result_id %in% names(entries)) {
    stop(
      "No result with `result_id = \"", result_id,
      "\"` was found for analysis type '", type, "'.",
      call. = FALSE
    )
  }
  entries[[result_id]] <- NULL
  if (length(entries) == 0L) {
    misc_data[["shennong"]][["results"]][[type]] <- NULL
  } else {
    misc_data[["shennong"]][["results"]][[type]] <- entries
  }
  if (length(misc_data[["shennong"]][["results"]]) == 0L) {
    misc_data[["shennong"]][["results"]] <- NULL
  }
  if (length(misc_data[["shennong"]]) == 0L) {
    misc_data[["shennong"]] <- NULL
  }
  methods::slot(object, "misc") <- misc_data
  object
}

#' Delete a registered workflow artifact
#'
#' Registered artifact collections (for example \code{"clustering_stage_cache"},
#' \code{"integration_comparison"}, or \code{"label_transfer"}) are removed
#' either member-by-member or as a whole container. Because these historical
#' collections live at top level in \code{object@misc}, their names alone cannot
#' prove Shennong ownership; every deletion therefore requires explicit
#' confirmation. Unknown artifact types fail closed; unified analysis results
#' must be deleted through \code{\link{sn_delete_result}} instead.
#'
#' @param object A Seurat object.
#' @param artifact_type Registered artifact type, i.e. the collection name or
#'   its \code{"*_artifact"} alias (for example \code{"integration"} or
#'   \code{"integration_artifact"}).
#' @param artifact_id Optional artifact identifier. When omitted, the entire artifact
#'   collection container is removed only when \code{confirm = TRUE}.
#' @param confirm Explicit confirmation required before deleting a member or an
#'   entire top-level artifact collection. This protects user-owned payloads
#'   that use a legacy collection name also registered by Shennong.
#'
#' @return The modified Seurat object.
#'
#' @examples
#' \dontrun{
#' obj <- sn_delete_artifact(
#'   obj, "integration_comparison", "pbmc_grid", confirm = TRUE
#' )
#' obj <- sn_delete_artifact(obj, "label_transfer", confirm = TRUE)
#' }
#'
#' @export
sn_delete_artifact <- function(object, artifact_type, artifact_id = NULL,
                               confirm = FALSE) {
  .sn_validate_seurat_object(object)
  registry <- .sn_misc_result_registry()
  artifacts <- registry[registry$contract_scope == "artifact", , drop = FALSE]
  requested <- tolower(as.character(artifact_type))
  if (!requested %in% artifacts$type && !requested %in% artifacts$collection) {
    stop(
      "'",
      artifact_type,
      "' is not a registered artifact type. Registered types: ",
      paste(sort(artifacts$type), collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  type_index <- match(requested, artifacts$type)
  if (is_na(type_index)) {
    type_index <- match(requested, artifacts$collection)
  }
  collection <- artifacts$collection[[type_index]]
  misc_data <- methods::slot(object, "misc")
  if (!isTRUE(confirm)) {
    target <- if (is_null(artifact_id)) {
      paste0("the entire top-level artifact collection '", collection, "'")
    } else {
      paste0(
        "artifact '", as.character(artifact_id), "' in the legacy top-level collection '",
        collection, "'"
      )
    }
    stop(
      "Deleting ", target,
      " requires `confirm = TRUE`; a legacy `object@misc` name does not prove ",
      "that Shennong owns the payload.",
      call. = FALSE
    )
  }
  if (is_null(artifact_id)) {
    if (is_null(misc_data[[collection]])) {
      warning("No '", collection, "' artifact collection was present.", call. = FALSE)
      return(object)
    }
    misc_data[[collection]] <- NULL
    methods::slot(object, "misc") <- misc_data
    return(object)
  }
  entries <- misc_data[[collection]]
  if (!is.list(entries) || !as.character(artifact_id) %in% names(entries)) {
    stop(
      "No artifact with artifact_id '", artifact_id, "' was found in the '", collection, "' collection.",
      call. = FALSE
    )
  }
  entries[[as.character(artifact_id)]] <- NULL
  if (length(entries) == 0L) {
    misc_data[[collection]] <- NULL
  } else {
    misc_data[[collection]] <- entries
  }
  methods::slot(object, "misc") <- misc_data
  object
}
#' List stored Shennong analysis and interpretation results on a Seurat object
#'
#' @param object A \code{Seurat} object.
#' @param type Optional analysis type used to filter the result inventory.
#'
#' @return A tibble describing registered Shennong stored-result collections,
#'   including DE, enrichment, interpretation, deconvolution, Milo,
#'   communication, regulatory activity, and QC assessment entries when present.
#'   The canonical lookup key is reported in \code{result_id}.
#'
#' @examples
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'   counts <- matrix(rpois(10 * 12, lambda = 1), nrow = 10, ncol = 12)
#'   rownames(counts) <- c(
#'     "CD3D", "CD3E", "TRAC", "LTB", "MS4A1",
#'     "CD79A", "HLA-DRA", "LYZ", "ACTB", "MALAT1"
#'   )
#'   colnames(counts) <- paste0("cell", 1:12)
#'   obj <- sn_initialize_seurat_object(counts, species = "human")
#'   obj <- Seurat::NormalizeData(obj, verbose = FALSE)
#'   obj <- sn_find_de(
#'     obj,
#'     analysis = "markers",
#'     group_by = NULL,
#'     layer = "data",
#'     min_pct = 0,
#'     logfc_threshold = 0,
#'     return_object = TRUE,
#'     verbose = FALSE
#'   )
#'   sn_list_results(obj)
#' }
#' @param include_artifacts Include registered workflow artifacts that do not
#'   implement the unified analysis-result contract.
#' @export
sn_list_results <- function(object, type = NULL, include_artifacts = FALSE) {
  .sn_validate_seurat_object(object)

  result <- .sn_generic_result_summary(object)
  if (ncol(result) == 0L) {
    result <- tibble::tibble(
      collection = character(), type = character(), result_id = character(),
      analysis = character(), method = character(), created_at = character(),
      n_rows = integer(), source = character()
    )
  }
  if (isTRUE(include_artifacts)) {
    artifact_entries <- .sn_stored_analysis_result_entries(
      object,
      include_artifacts = TRUE
    )
    artifact_entries <- Filter(
      function(entry) identical(entry$contract_scope, "artifact"),
      artifact_entries
    )
    artifact_summary <- dplyr::bind_rows(lapply(artifact_entries, function(entry) {
      artifact <- entry$result
      artifact_field <- function(name) {
        if (is.list(artifact)) artifact[[name]] else NULL
      }
      tibble::tibble(
        collection = entry$collection,
        type = entry$type,
        result_id = entry$result_id,
        analysis = as.character(artifact_field("analysis") %||% NA_character_)[[1]],
        method = as.character(artifact_field("method") %||% NA_character_)[[1]],
        created_at = as.character(artifact_field("created_at") %||% NA_character_)[[1]],
        n_rows = if (is.list(artifact)) .sn_result_n_rows(artifact) else 0L,
        source = as.character(artifact_field("source_result_id") %||% NA_character_)[[1]]
      )
    }))
    result <- dplyr::bind_rows(result, artifact_summary)
  }
  if (!is_null(type)) {
    requested_types <- tolower(as.character(type))
    result <- dplyr::filter(result, .data$type %in% .env$requested_types)
  }
  result |>
    dplyr::arrange(.data$collection, .data$result_id)
}

#' Retrieve a stored DE result from a Seurat object
#'
#' @param object A \code{Seurat} object.
#' @param result_id Identifier of the stored DE result.
#' @param top_n Optional number of rows to keep. When supplied together with a
#'   ranking column, results are reduced to the top rows overall or per group.
#' @param direction One of \code{"all"}, \code{"up"}, or \code{"down"}.
#' @param groups Optional subset of group labels to keep.
#' @param with_metadata If \code{TRUE}, return the full stored result list
#'   instead of just the result table.
#'
#' @return A tibble or stored-result list.
#'
#' @examples
#' \dontrun{
#' markers <- sn_get_de_result(seurat_obj, result_id = "cluster_markers", top_n = 5)
#' }
#' @export
sn_get_de_result <- function(object,
                             result_id = "default",
                             top_n = NULL,
                             direction = c("all", "up", "down"),
                             groups = NULL,
                             with_metadata = FALSE) {
  .sn_validate_seurat_object(object)

  stored <- sn_get_result(object, type = "de", result_id = result_id)
  if (isTRUE(with_metadata)) {
    return(stored)
  }

  .sn_subset_ranked_table(
    table = stored$tables$primary,
    rank_col = stored$rank_col,
    group_col = stored$group_col,
    top_n = top_n,
    direction = direction,
    groups = groups
  )
}

#' Retrieve a stored enrichment result from a Seurat object
#'
#' @param object A \code{Seurat} object.
#' @param result_id Identifier of the stored enrichment result.
#' @param top_n Optional number of top terms to keep.
#' @param groups Optional subset of cluster/group labels when the stored table
#'   includes a \code{Cluster} column.
#' @param with_metadata If \code{TRUE}, return the full stored result list
#'   instead of just the term table.
#'
#' @return A tibble or stored-result list.
#'
#' @examples
#' \dontrun{
#' terms <- sn_get_enrichment_result(seurat_obj, result_id = "cluster_gsea", top_n = 10)
#' }
#' @export
sn_get_enrichment_result <- function(object,
                                     result_id = "default",
                                     top_n = NULL,
                                     groups = NULL,
                                     with_metadata = FALSE) {
  .sn_validate_seurat_object(object)

  stored <- sn_get_result(object, type = "enrichment", result_id = result_id)
  if (isTRUE(with_metadata)) {
    return(stored)
  }

  table <- tibble::as_tibble(stored$tables$primary)
  group_col <- c("Cluster", "cluster", ".sign")[
    c("Cluster", "cluster", ".sign") %in% colnames(table)
  ][1] %||% NULL
  rank_col <- c("NES", "Count", "GeneRatio", "p.adjust", "pvalue")[
    c("NES", "Count", "GeneRatio", "p.adjust", "pvalue") %in% colnames(table)
  ][1] %||% NULL

  if (!is_null(groups) && !is_null(group_col)) {
    table <- dplyr::filter(table, .data[[group_col]] %in% groups)
  }

  if (is_null(top_n) || is_null(rank_col)) {
    return(table)
  }

  if (rank_col %in% c("p.adjust", "pvalue")) {
    table <- table[order(table[[rank_col]], decreasing = FALSE), , drop = FALSE]
  } else {
    table <- table[order(abs(table[[rank_col]]), decreasing = TRUE), , drop = FALSE]
  }

  utils::head(table, top_n)
}

#' Retrieve a stored interpretation result from a Seurat object
#'
#' @param object A \code{Seurat} object.
#' @param result_id Identifier of the stored interpretation result.
#'
#' @return The stored interpretation-result list.
#'
#' @examples
#' \dontrun{
#' interpretation <- sn_get_interpretation_result(seurat_obj, "annotation_note")
#' }
#' @export
sn_get_interpretation_result <- function(object, result_id = "default") {
  .sn_validate_seurat_object(object)

  sn_get_result(object, type = "interpretation", result_id = result_id)
}

.sn_resolve_result_input <- function(x, type, result_id = NULL) {
  result <- if (inherits(x, "Seurat")) {
    if (is_null(result_id)) {
      stop("`result_id` is required when `x` is a Seurat object.", call. = FALSE)
    }
    sn_get_result(x, type, result_id)
  } else {
    x
  }
  sn_validate_result(result)
  if (!identical(result$analysis_type, type)) stop("Expected a ", type, " result.", call. = FALSE)
  result
}
