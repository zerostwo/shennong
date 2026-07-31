.sn_analysis_result_required_fields <- function() {
  c(
    "schema_version", "analysis_type", "name", "method", "backend", "input",
    "parameters", "tables", "embeddings", "graphs", "models", "diagnostics",
    "warnings", "provenance"
  )
}

.sn_analysis_result_schema_version <- function() "1.0.0"

.sn_normalize_analysis_result_schema_version <- function(version) {
  version <- as.character(version %||% .sn_analysis_result_schema_version())
  if (length(version) != 1L || is.na(version) || !nzchar(version)) {
    return(version)
  }
  if (grepl("^[0-9]+$", version)) return(paste0(version, ".0.0"))
  if (grepl("^[0-9]+\\.[0-9]+$", version)) return(paste0(version, ".0"))
  version
}

.sn_is_supported_analysis_result_schema_version <- function(version) {
  version <- as.character(version %||% NA_character_)
  if (length(version) != 1L || is.na(version) ||
      !grepl("^[0-9]+\\.[0-9]+\\.[0-9]+([+-][0-9A-Za-z.-]+)?$", version)) {
    return(FALSE)
  }
  identical(version, .sn_analysis_result_schema_version())
}

.sn_analysis_result_type_specs <- function() {
  list(
    annotation = list(primary_columns = c("cell", "prediction")),
    program_scoring = list(primary_columns = c("entity", "program", "score")),
    trajectory = list(primary_columns = c("cell", "primary_pseudotime")),
    state_priority = list(primary_columns = c("state", "priority_score")),
    scissor = list(primary_columns = c("cell", "coefficient", "selection")),
    bulk_survival = list(primary_columns = c(
      "feature", "hazard_ratio", "conf_low", "conf_high", "p_value"
    ))
  )
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
    tryCatch(.sn_autozyme_provenance(), error = function(e) list())
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
      getOption("shennong.autozyme.provenance_context")
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

.sn_new_analysis_result <- function(analysis_type,
                                    name,
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
  result <- list(
    schema_version = .sn_analysis_result_schema_version(),
    analysis_type = analysis_type,
    name = name,
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
    provenance = provenance %||% .sn_analysis_provenance(random_seed = random_seed)
  )
  sn_validate_result(result)
  result
}

.sn_upgrade_analysis_result <- function(result,
                                        analysis_type,
                                        name,
                                        method = NULL,
                                        backend = NULL) {
  if (!is.list(result)) {
    stop("`result` must be a list.", call. = FALSE)
  }
  analysis_type <- as.character(
    analysis_type %||% result[["analysis_type"]] %||% result[["analysis"]]
  )
  name <- as.character(name %||% result[["name"]])
  method <- as.character(method %||% result[["method"]] %||% "unknown")
  backend <- as.character(backend %||% result[["backend"]] %||% method)

  original_schema_version <- result[["schema_version"]]
  if (!is_null(original_schema_version)) {
    if (!is.character(original_schema_version) || length(original_schema_version) != 1L ||
        is.na(original_schema_version) || !nzchar(original_schema_version)) {
      stop("Cannot upgrade a result with an invalid `schema_version`.", call. = FALSE)
    }
    normalized_schema <- .sn_normalize_analysis_result_schema_version(original_schema_version)
    if (!identical(normalized_schema, .sn_analysis_result_schema_version())) {
      stop(
        "Cannot upgrade unsupported `schema_version` '", original_schema_version,
        "'; this Shennong version only upgrades compatible 1.0 results.",
        call. = FALSE
      )
    }
  }

  result[["schema_version"]] <- .sn_analysis_result_schema_version()
  result[["analysis_type"]] <- analysis_type
  result[["name"]] <- name
  result[["method"]] <- method
  result[["backend"]] <- backend
  result[["input"]] <- result[["input"]] %||% list()
  result[["parameters"]] <- result[["parameters"]] %||% list()
  result[["tables"]] <- .sn_result_tables(result)
  result[["embeddings"]] <- result[["embeddings"]] %||% list()
  result[["graphs"]] <- result[["graphs"]] %||% list()
  result[["models"]] <- result[["models"]] %||% list()
  result[["diagnostics"]] <- result[["diagnostics"]] %||% list()
  result[["warnings"]] <- as.character(result[["warnings"]] %||% character())
  result[["provenance"]] <- .sn_analysis_provenance(
    result,
    capture_acceleration = FALSE
  )
  if (!is_null(original_schema_version) &&
      !identical(original_schema_version, result[["schema_version"]])) {
    result[["provenance"]][["source_schema_version"]] <-
      result[["provenance"]][["source_schema_version"]] %||% original_schema_version
  }

  result[["package_version"]] <- result[["package_version"]] %||%
    result[["provenance"]][["package_versions"]][["Shennong"]] %||% NA_character_
  result[["created_at"]] <- result[["created_at"]] %||%
    result[["provenance"]][["timestamp"]]
  result[["analysis"]] <- result[["analysis"]] %||% analysis_type
  if (is_null(result[["table"]]) && is.data.frame(result[["tables"]][["primary"]])) {
    result[["table"]] <- result[["tables"]][["primary"]]
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
  schema_version <- result[["schema_version"]] %||% NA_character_
  normalized_schema <- .sn_normalize_analysis_result_schema_version(schema_version)
  legacy_schema <- is.character(schema_version) && length(schema_version) == 1L &&
    !is.na(schema_version) && !identical(schema_version, normalized_schema) &&
    identical(normalized_schema, .sn_analysis_result_schema_version())
  if (isTRUE(legacy_schema)) {
    warnings <- c(
      warnings,
      paste0(
        "Legacy `schema_version` '", schema_version, "' is compatible but should be upgraded to '",
        .sn_analysis_result_schema_version(), "'."
      )
    )
  } else if (is.character(schema_version) && length(schema_version) == 1L &&
             !is.na(schema_version) && nzchar(schema_version) &&
             !grepl("^[0-9]+\\.[0-9]+\\.[0-9]+([+-][0-9A-Za-z.-]+)?$", schema_version)) {
    errors <- c(errors, "`schema_version` must use semantic version form such as '1.0.0'.")
  }
  if (!isTRUE(legacy_schema) && is.character(schema_version) &&
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
    provenance_missing <- setdiff(c("package_versions", "random_seed", "timestamp"), names(provenance))
    if (length(provenance_missing) > 0L) {
      errors <- c(
        errors,
        paste0("`provenance` is missing field(s): ", paste(provenance_missing, collapse = ", "), ".")
      )
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
      is.numeric(value) && length(value) == 1L &&
        (is.na(value) || (is.finite(value) && value >= 0 && value == as.integer(value)))
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
  if (isTRUE(legacy_schema) && !is.data.frame(primary)) {
    primary <- switch(
      analysis_type,
      annotation = tables[["cells"]],
      program_scoring = tables[["scores"]],
      trajectory = tables[["cells"]],
      NULL
    )
    if (is.data.frame(primary)) {
      warnings <- c(
        warnings,
        paste0(
          "Legacy '", analysis_type,
          "' result uses a named table alias; upgrade it to populate `tables$primary`."
        )
      )
    }
  }
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
  }
  legacy_table <- result[["table"]]
  if (is.data.frame(legacy_table) && is.data.frame(tables[["primary"]]) &&
      !identical(legacy_table, tables[["primary"]])) {
    errors <- c(errors, "Legacy `table` and canonical `tables$primary` are not synchronized.")
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
#'   schema_version = "1.0.0", analysis_type = "demo", name = "example",
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
  hit <- registry[
    registry$type == type & registry$contract_scope == "analysis_result" & registry$listable,
    ,
    drop = FALSE
  ]
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
    legacy_table <- result[["table"]]
    if ("table" %in% required && is.data.frame(legacy_table)) {
      result[["tables"]][["primary"]] <- legacy_table
    }
    if ("table" %in% required && is_null(legacy_table)) {
      stop(
        "Result type '", type, "' requires `tables$primary` to be a data frame.",
        call. = FALSE
      )
    }
    if ("database" %in% required) {
      result[["database"]] <- result[["database"]] %||%
        result[["parameters"]][["database"]] %||% "unknown"
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
    prepared <- .sn_prepare_misc_result(collection, name, prepared)
    return(.sn_store_misc_result(object, collection, name, prepared))
  }

  misc_data <- methods::slot(object, "misc")
  misc_data[["analysis_results"]] <- misc_data[["analysis_results"]] %||% list()
  misc_data[["analysis_results"]][[type]] <-
    misc_data[["analysis_results"]][[type]] %||% list()
  misc_data[["analysis_results"]][[type]][[name]] <- prepared
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
    results <- misc_data[["analysis_results"]][[type]] %||% list()
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
  collections <- misc_data[["analysis_results"]] %||% list()
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
      analysis = vapply(names, function(name) {
        entries[[name]][["analysis_type"]] %||% type
      }, character(1)),
      method = vapply(names, function(name) {
        entries[[name]][["method"]] %||% NA_character_
      }, character(1)),
      created_at = vapply(names, function(name) {
        provenance <- entries[[name]][["provenance"]] %||% list()
        provenance[["timestamp"]] %||% NA_character_
      }, character(1)),
      n_rows = unname(vapply(names, function(name) {
        tables <- entries[[name]][["tables"]] %||% list()
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
  if (!isTRUE(include_artifacts)) {
    registry <- registry[registry$contract_scope == "analysis_result", , drop = FALSE]
  }
  entries <- list()
  for (index in seq_len(nrow(registry))) {
    current_type <- registry$type[[index]]
    if (!is_null(type) && !current_type %in% type) next
    collection <- registry$collection[[index]]
    collection_data <- misc_data[[collection]] %||% list()
    contract_scope <- registry$contract_scope[[index]]
    if (identical(contract_scope, "artifact")) {
      if (length(collection_data) > 0L) {
        entries[[length(entries) + 1L]] <- list(
          collection = collection,
          type = current_type,
          name = collection,
          contract_scope = contract_scope,
          result = collection_data
        )
      }
      next
    }
    for (name in names(collection_data)) {
      entries[[length(entries) + 1L]] <- list(
        collection = collection,
        type = current_type,
        name = name,
        contract_scope = "analysis_result",
        result = collection_data[[name]]
      )
    }
  }
  generic <- misc_data[["analysis_results"]] %||% list()
  for (current_type in names(generic)) {
    if (!is_null(type) && !current_type %in% type) next
    for (name in names(generic[[current_type]])) {
      entries[[length(entries) + 1L]] <- list(
        collection = "analysis_results",
        type = current_type,
        name = name,
        contract_scope = "analysis_result",
        result = generic[[current_type]][[name]]
      )
    }
  }
  if (isTRUE(include_artifacts)) {
    registered_collections <- c(registry$collection, "analysis_results")
    unregistered_collections <- setdiff(names(misc_data), registered_collections)
    for (collection in unregistered_collections) {
      collection_data <- misc_data[[collection]]
      if (length(collection_data) == 0L) next
      current_type <- paste0(collection, "_unregistered")
      if (!is_null(type) && !current_type %in% type && !collection %in% type) next
      entries[[length(entries) + 1L]] <- list(
        collection = collection,
        type = current_type,
        name = collection,
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
#' distinguishes results that already satisfy the current contract from legacy
#' results that can be upgraded safely, and reports unknown payloads as
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
  .sn_validate_result_object(object)
  requested_types <- if (is_null(type)) NULL else tolower(as.character(type))
  entries <- .sn_stored_analysis_result_entries(
    object,
    requested_types,
    include_artifacts = include_artifacts
  )
  if (length(entries) == 0L) {
    return(tibble::tibble(
      collection = character(), type = character(), name = character(),
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
        name = entry$name,
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
        name = entry$name,
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
    upgraded <- tryCatch(
      .sn_prepare_result_for_collection(raw, type = entry$type, name = entry$name),
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
    status <- if (unified) "valid" else if (isTRUE(upgraded_report$valid)) "legacy" else "invalid"
    primary <- if (repairable) upgraded[["tables"]][["primary"]] else NULL
    tibble::tibble(
      collection = entry$collection,
      type = entry$type,
      name = entry$name,
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
    dplyr::arrange(.data$collection, .data$type, .data$name)
}

#' Upgrade stored Shennong analysis results
#'
#' Upgrade listable legacy results in place while preserving their established
#' physical \code{object@misc} collections and compatibility aliases.
#'
#' @param object A \code{Seurat} object.
#' @param type Optional analysis type or character vector of types to upgrade.
#' @param strict If \code{TRUE}, stop at the first result that cannot be safely
#'   upgraded. If \code{FALSE}, leave invalid entries unchanged and warn.
#'
#' @return The modified \code{Seurat} object.
#'
#' @examples
#' \dontrun{object <- sn_upgrade_results(object)}
#'
#' @export
sn_upgrade_results <- function(object, type = NULL, strict = TRUE) {
  .sn_validate_result_object(object)
  requested_types <- if (is_null(type)) NULL else tolower(as.character(type))
  entries <- .sn_stored_analysis_result_entries(object, requested_types)
  for (entry in entries) {
    updated <- tryCatch(
      sn_store_result(object, entry$type, entry$name, entry$result),
      error = identity
    )
    if (inherits(updated, "error")) {
      message <- paste0(
        "Could not upgrade stored result '", entry$name, "' of type '",
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
    entries <- misc_data[["analysis_results"]][[type]] %||% list()
    if (!name %in% names(entries)) {
      stop("No stored result named '", name, "' was found for analysis type '", type, "'.", call. = FALSE)
    }
    entries[[name]] <- NULL
    misc_data[["analysis_results"]][[type]] <- entries
  }
  methods::slot(object, "misc") <- misc_data
  object
}
