.sn_result_bundle_schema <- function() {
  "shennong.dev/analysis-result-bundle/v1"
}

.sn_result_bundle_sensitive_fields <- function() {
  c(
    "password", "passwd", "token", "access_token", "refresh_token",
    "api_key", "admin_key", "authorization", "client_secret", "secret",
    "access_key", "secret_key"
  )
}

.sn_result_bundle_scalar <- function(value, label) {
  if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
    stop("`", label, "` must be a non-empty character scalar.", call. = FALSE)
  }
  value
}

.sn_result_bundle_timestamp <- function(value, label) {
  value <- .sn_result_bundle_scalar(value, label)
  if (!grepl(
    "^[0-9]{4}-[0-9]{2}-[0-9]{2}T[0-9]{2}:[0-9]{2}:[0-9]{2}(\\.[0-9]+)?Z$",
    value
  )) {
    stop(
      "`", label, "` must be an RFC 3339 UTC timestamp ending in `Z`.",
      call. = FALSE
    )
  }
  parsed <- suppressWarnings(as.POSIXct(
    sub("Z$", "", value),
    format = "%Y-%m-%dT%H:%M:%OS",
    tz = "UTC"
  ))
  if (is.na(parsed)) {
    stop("`", label, "` must be a valid RFC 3339 UTC timestamp.", call. = FALSE)
  }
  value
}

.sn_result_bundle_prepare_json <- function(value, path = "bundle") {
  if (is_null(value)) return(NULL)
  if (inherits(value, "formula")) {
    return(paste(deparse(value, width.cutoff = 500L), collapse = " "))
  }
  if (inherits(value, "POSIXt")) {
    return(format(value, tz = "UTC", usetz = TRUE))
  }
  if (inherits(value, "Date")) {
    return(format(value, "%Y-%m-%d"))
  }
  if (inherits(value, "difftime")) {
    return(as.numeric(value, units = "secs"))
  }
  if (is.function(value) || is.environment(value) || typeof(value) %in% c(
    "externalptr", "weakref", "symbol"
  )) {
    stop(
      "`", path, "` contains an object that cannot be represented in JSON. ",
      "Store executable or binary state as a digest-addressed artifact instead.",
      call. = FALSE
    )
  }
  if (is.language(value)) {
    return(paste(deparse(value, width.cutoff = 500L), collapse = " "))
  }
  if (isS4(value)) {
    stop(
      "`", path, "` contains an S4 object that cannot be represented in the ",
      "Result Bundle JSON boundary. Store it as a digest-addressed artifact.",
      call. = FALSE
    )
  }
  if (is.raw(value) || is.complex(value)) {
    stop(
      "`", path, "` contains binary or complex values that are outside the ",
      "Result Bundle JSON boundary. Store them as a digest-addressed artifact.",
      call. = FALSE
    )
  }
  if (is.factor(value)) return(as.character(value))
  if (inherits(value, "package_version")) return(as.character(value))
  if (is.data.frame(value)) {
    prepared <- value
    for (index in seq_along(prepared)) {
      column_name <- names(prepared)[[index]] %||% as.character(index)
      prepared[[index]] <- .sn_result_bundle_prepare_json(
        prepared[[index]],
        paste0(path, "$", column_name)
      )
    }
    return(prepared)
  }
  if (is.matrix(value) || is.array(value)) {
    if (is.numeric(value) && any(!is.finite(value) & !is.na(value))) {
      stop("`", path, "` contains non-finite numeric values.", call. = FALSE)
    }
    return(array(
      as.vector(value),
      dim = dim(value),
      dimnames = dimnames(value)
    ))
  }
  if (is.list(value) || is.pairlist(value)) {
    value <- as.list(value)
    output <- lapply(seq_along(value), function(index) {
      element_name <- names(value)[[index]] %||% as.character(index)
      .sn_result_bundle_prepare_json(
        value[[index]],
        paste0(path, "$", element_name)
      )
    })
    names(output) <- names(value)
    return(output)
  }
  if (is.numeric(value) && any(!is.finite(value) & !is.na(value))) {
    stop("`", path, "` contains non-finite numeric values.", call. = FALSE)
  }
  if (is.atomic(value)) return(value)
  stop(
    "`", path, "` contains an unsupported value of type `", typeof(value), "`.",
    call. = FALSE
  )
}

.sn_result_bundle_credential_paths <- function(value, path = "bundle") {
  if (!is.list(value) && !is.data.frame(value)) return(character())
  value_names <- names(value)
  hits <- character()
  if (!is_null(value_names)) {
    normalized <- gsub("[^a-z0-9]+", "_", tolower(value_names))
    sensitive <- normalized %in% .sn_result_bundle_sensitive_fields()
    if (any(sensitive)) {
      hits <- paste0(path, "$", value_names[sensitive])
    }
  }
  children <- unlist(lapply(seq_along(value), function(index) {
    child_name <- value_names[[index]] %||% as.character(index)
    .sn_result_bundle_credential_paths(
      value[[index]],
      paste0(path, "$", child_name)
    )
  }), use.names = FALSE)
  unique(c(hits, children))
}

.sn_result_bundle_digest <- function(digest_value, label) {
  if (!is.list(digest_value)) {
    stop("`", label, "` must be a digest record.", call. = FALSE)
  }
  unknown <- setdiff(names(digest_value), c("algorithm", "value"))
  if (length(unknown) > 0L) {
    stop(
      "`", label, "` contains unsupported field(s): ",
      paste(unknown, collapse = ", "), ".",
      call. = FALSE
    )
  }
  algorithm <- tolower(.sn_result_bundle_scalar(
    digest_value[["algorithm"]],
    paste0(label, "$algorithm")
  ))
  value <- tolower(.sn_result_bundle_scalar(
    digest_value[["value"]],
    paste0(label, "$value")
  ))
  if (!identical(algorithm, "sha256")) {
    stop("`", label, "$algorithm` must be `sha256` for Result Bundle v1.", call. = FALSE)
  }
  if (!grepl("^[0-9a-f]{64}$", value)) {
    stop("`", label, "$value` must contain exactly 64 hexadecimal characters.", call. = FALSE)
  }
  list(algorithm = algorithm, value = value)
}

.sn_result_bundle_relative_path <- function(path, label) {
  path <- .sn_result_bundle_scalar(path, label)
  normalized <- gsub("\\\\", "/", path)
  if (startsWith(normalized, "/") ||
      grepl("^[A-Za-z]:/", normalized) ||
      any(strsplit(normalized, "/", fixed = TRUE)[[1L]] == "..")) {
    stop("`", label, "` must be a safe bundle-relative path.", call. = FALSE)
  }
  normalized
}

.sn_result_bundle_size <- function(value, label) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 0 || value != floor(value)) {
    stop("`", label, "` must be one non-negative whole-byte count.", call. = FALSE)
  }
  as.numeric(value)
}

.sn_result_bundle_normalize_records <- function(records, kind) {
  if (is_null(records)) return(list())
  if (!is.list(records) || is.data.frame(records)) {
    stop("`", kind, "` must be a list of records.", call. = FALSE)
  }
  lapply(seq_along(records), function(index) {
    record <- records[[index]]
    label <- paste0(kind, "[[", index, "]]")
    if (!is.list(record) || is.data.frame(record)) {
      stop("`", label, "` must be a record list.", call. = FALSE)
    }
    if (identical(kind, "inputs")) {
      required <- c("role", "revision", "digest")
      allowed <- c(
        required, "resource_id", "artifact_id", "media_type", "size_bytes",
        "metadata"
      )
    } else {
      required <- c("role", "digest")
      allowed <- c(
        required, "artifact_id", "path", "media_type", "size_bytes", "metadata"
      )
    }
    missing <- setdiff(required, names(record))
    if (length(missing) > 0L) {
      stop(
        "`", label, "` is missing required field(s): ",
        paste(missing, collapse = ", "), ".",
        call. = FALSE
      )
    }
    if (identical(kind, "inputs") &&
        !any(c("resource_id", "artifact_id") %in% names(record))) {
      stop(
        "`", label, "` must contain at least one immutable input identifier: ",
        "`resource_id` or `artifact_id`.",
        call. = FALSE
      )
    }
    unknown <- setdiff(names(record), allowed)
    if (length(unknown) > 0L) {
      stop(
        "`", label, "` contains unsupported field(s): ",
        paste(unknown, collapse = ", "), ".",
        call. = FALSE
      )
    }
    record[["role"]] <- .sn_result_bundle_scalar(
      record[["role"]],
      paste0(label, "$role")
    )
    if (identical(kind, "inputs")) {
      record[["revision"]] <- .sn_result_bundle_scalar(
        record[["revision"]],
        paste0(label, "$revision")
      )
    }
    for (field in intersect(c("resource_id", "artifact_id", "media_type"), names(record))) {
      record[[field]] <- .sn_result_bundle_scalar(
        record[[field]],
        paste0(label, "$", field)
      )
    }
    if ("path" %in% names(record)) {
      record[["path"]] <- .sn_result_bundle_relative_path(
        record[["path"]],
        paste0(label, "$path")
      )
    }
    if ("size_bytes" %in% names(record)) {
      record[["size_bytes"]] <- .sn_result_bundle_size(
        record[["size_bytes"]],
        paste0(label, "$size_bytes")
      )
    }
    if ("metadata" %in% names(record) && !is.list(record[["metadata"]])) {
      stop("`", label, "$metadata` must be a list.", call. = FALSE)
    }
    record[["digest"]] <- .sn_result_bundle_digest(
      record[["digest"]],
      paste0(label, "$digest")
    )
    .sn_result_bundle_prepare_json(record, label)
  })
}

.sn_result_bundle_json <- function(bundle, pretty = FALSE) {
  encoded <- jsonlite::toJSON(
    unclass(bundle),
    auto_unbox = TRUE,
    dataframe = "rows",
    matrix = "rowmajor",
    Date = "ISO8601",
    POSIXt = "ISO8601",
    factor = "string",
    null = "null",
    na = "null",
    digits = NA,
    pretty = isTRUE(pretty)
  )
  if (!jsonlite::validate(encoded)) {
    stop("Result Bundle serialization did not produce valid JSON.", call. = FALSE)
  }
  encoded
}

.sn_result_bundle_validation <- function(bundle) {
  errors <- character()
  warnings <- character()
  if (!is.list(bundle)) {
    return(list(valid = FALSE, errors = "`bundle` must be a list.", warnings = warnings))
  }
  required <- c(
    "schema", "created_at", "result", "validation", "inputs", "provenance",
    "artifacts"
  )
  missing <- setdiff(required, names(bundle))
  if (length(missing) > 0L) {
    errors <- c(
      errors,
      paste0("Missing required field(s): ", paste(missing, collapse = ", "), ".")
    )
  }
  if ("schema" %in% names(bundle) &&
      !identical(bundle[["schema"]], .sn_result_bundle_schema())) {
    errors <- c(
      errors,
      paste0("`schema` must be `", .sn_result_bundle_schema(), "`.")
    )
  }
  if ("created_at" %in% names(bundle)) {
    created_at_error <- tryCatch(
      {
        .sn_result_bundle_timestamp(bundle[["created_at"]], "created_at")
        NULL
      },
      error = conditionMessage
    )
    if (!is_null(created_at_error)) errors <- c(errors, created_at_error)
  }
  if ("result" %in% names(bundle)) {
    result_report <- tryCatch(
      sn_validate_result(bundle[["result"]], error = FALSE),
      error = identity
    )
    if (inherits(result_report, "error")) {
      errors <- c(errors, conditionMessage(result_report))
    } else if (!isTRUE(result_report$valid)) {
      errors <- c(
        errors,
        paste0("Invalid bundled analysis result: ", paste(result_report$errors, collapse = "; "))
      )
    }
  }
  validation <- bundle[["validation"]]
  if (!is_null(validation)) {
    if (!is.list(validation) || !identical(validation[["valid"]], TRUE) ||
        !is.character(validation[["errors"]] %||% character()) ||
        !is.character(validation[["warnings"]] %||% character())) {
      errors <- c(
        errors,
        "`validation` must record a successful analysis-result validation report."
      )
    }
  }
  for (kind in c("inputs", "artifacts")) {
    if (kind %in% names(bundle)) {
      record_error <- tryCatch(
        {
          .sn_result_bundle_normalize_records(bundle[[kind]], kind)
          NULL
        },
        error = conditionMessage
      )
      if (!is_null(record_error)) errors <- c(errors, record_error)
    }
  }
  provenance <- bundle[["provenance"]]
  if (!is_null(provenance) && !is.list(provenance)) {
    errors <- c(errors, "`provenance` must be a list.")
  } else if (is.list(provenance)) {
    provenance_required <- c(
      "package_versions", "random_seed", "result_timestamp", "execution"
    )
    provenance_missing <- setdiff(provenance_required, names(provenance))
    if (length(provenance_missing) > 0L) {
      errors <- c(
        errors,
        paste0(
          "`provenance` is missing field(s): ",
          paste(provenance_missing, collapse = ", "), "."
        )
      )
    }
    if (!is.list(provenance[["package_versions"]])) {
      errors <- c(errors, "`provenance$package_versions` must be a list.")
    }
    if (!is.list(provenance[["execution"]])) {
      errors <- c(errors, "`provenance$execution` must be a list.")
    }
    result_timestamp <- provenance[["result_timestamp"]]
    if (!is_null(result_timestamp) &&
        (!is.character(result_timestamp) || length(result_timestamp) != 1L ||
          is.na(result_timestamp) || !nzchar(result_timestamp))) {
      errors <- c(
        errors,
        "`provenance$result_timestamp` must be a non-empty character scalar."
      )
    }
  }
  credential_paths <- .sn_result_bundle_credential_paths(bundle)
  if (length(credential_paths) > 0L) {
    errors <- c(
      errors,
      paste0(
        "Result Bundles must not contain credentials; remove: ",
        paste(credential_paths, collapse = ", "), "."
      )
    )
  }
  json_error <- tryCatch(
    {
      .sn_result_bundle_json(bundle)
      NULL
    },
    error = conditionMessage
  )
  if (!is_null(json_error)) errors <- c(errors, json_error)
  list(valid = length(errors) == 0L, errors = unique(errors), warnings = warnings)
}

#' Build, validate, and export a Shennong Result Bundle v1
#'
#' A Result Bundle is the package-owned, JSON-compatible handoff boundary
#' between a validated Shennong analysis result and an external orchestrator.
#' It contains no service calls or credentials. Input records identify exact
#' immutable revisions and SHA-256 digests; artifact records describe candidate
#' outputs that an authorized external service may verify and promote.
#'
#' @param result A canonical Shennong analysis result accepted by
#'   [sn_validate_result()].
#' @param inputs A list of immutable input-reference records. Each record
#'   requires `role`, `revision`, at least one of `resource_id` or
#'   `artifact_id`, and `digest = list(algorithm = "sha256",
#'   value = "<64 hex characters>")`. Optional fields are `media_type`,
#'   `size_bytes`, and `metadata`.
#' @param execution A JSON-compatible list of execution provenance such as job,
#'   runtime, image, environment-lock, or source revision identifiers.
#'   Credentials are rejected.
#' @param artifacts A list of candidate output-artifact records. Each record
#'   requires `role` and a SHA-256 `digest`. Optional fields are `artifact_id`,
#'   bundle-relative `path`, `media_type`, `size_bytes`, and `metadata`.
#' @param created_at RFC 3339 UTC creation timestamp ending in `Z`. By default
#'   the current UTC time is used.
#' @param bundle A Result Bundle to validate or export.
#' @param error If `TRUE`, stop for an invalid bundle. If `FALSE`, return a
#'   validation report.
#' @param path Destination JSON file.
#' @param pretty Write indented JSON.
#' @param overwrite Replace an existing destination file.
#'
#' @return `sn_build_result_bundle()` returns an `sn_result_bundle` list.
#'   `sn_validate_result_bundle()` returns a report with `valid`, `errors`, and
#'   `warnings`. `sn_export_result_bundle()` invisibly returns the normalized
#'   file path, byte size, and SHA-256 digest of the exported JSON.
#'
#' @examples
#' result <- list(
#'   schema_version = "1.0.0", analysis_type = "demo", name = "example",
#'   method = "mean", backend = "base", input = list(), parameters = list(),
#'   tables = list(primary = data.frame(feature = "gene1", value = 1)),
#'   embeddings = list(), graphs = list(), models = list(), diagnostics = list(),
#'   warnings = character(),
#'   provenance = list(
#'     package_versions = list(Shennong = "0.2.0"),
#'     random_seed = 1L,
#'     timestamp = "2026-01-01 UTC"
#'   )
#' )
#' bundle <- sn_build_result_bundle(result)
#' sn_validate_result_bundle(bundle, error = FALSE)
#' path <- tempfile(fileext = ".json")
#' sn_export_result_bundle(bundle, path)
#' unlink(path)
#'
#' @name result_bundle
NULL

#' @rdname result_bundle
#' @export
sn_build_result_bundle <- function(result,
                                   inputs = list(),
                                   execution = list(),
                                   artifacts = list(),
                                   created_at = format(
                                     Sys.time(),
                                     "%Y-%m-%dT%H:%M:%SZ",
                                     tz = "UTC"
                                   )) {
  sn_validate_result(result)
  prepared_result <- .sn_result_bundle_prepare_json(result, "result")
  sn_validate_result(prepared_result)
  inputs <- .sn_result_bundle_normalize_records(inputs, "inputs")
  artifacts <- .sn_result_bundle_normalize_records(artifacts, "artifacts")
  if (!is.list(execution)) {
    stop("`execution` must be a JSON-compatible list.", call. = FALSE)
  }
  execution <- .sn_result_bundle_prepare_json(execution, "execution")
  created_at <- .sn_result_bundle_timestamp(created_at, "created_at")
  result_report <- sn_validate_result(prepared_result, error = FALSE)
  result_provenance <- prepared_result[["provenance"]]
  bundle <- list(
    schema = .sn_result_bundle_schema(),
    created_at = created_at,
    result = prepared_result,
    validation = unclass(result_report),
    inputs = inputs,
    provenance = list(
      package_versions = result_provenance[["package_versions"]] %||% list(),
      random_seed = result_provenance[["random_seed"]] %||% NA_integer_,
      result_timestamp = result_provenance[["timestamp"]] %||% NA_character_,
      execution = execution
    ),
    artifacts = artifacts
  )
  class(bundle) <- c("sn_result_bundle", "list")
  sn_validate_result_bundle(bundle)
  bundle
}

#' @rdname result_bundle
#' @export
sn_validate_result_bundle <- function(bundle, error = TRUE) {
  report <- .sn_result_bundle_validation(bundle)
  class(report) <- c("sn_result_bundle_validation", "list")
  if (!report$valid && isTRUE(error)) {
    stop(
      "Invalid Shennong Result Bundle:\n- ",
      paste(report$errors, collapse = "\n- "),
      call. = FALSE
    )
  }
  if (isTRUE(error)) invisible(report) else report
}

.sn_result_bundle_publish <- function(temporary,
                                      target,
                                      overwrite,
                                      rename_file = file.rename,
                                      remove_file = unlink) {
  target_exists <- file.exists(target)
  if (target_exists && !isTRUE(overwrite)) {
    stop("Result Bundle destination already exists: ", target, call. = FALSE)
  }
  if (!target_exists) {
    if (!rename_file(temporary, target)) {
      stop("Could not publish Result Bundle to: ", target, call. = FALSE)
    }
    return(invisible(target))
  }

  # POSIX rename replaces the old file atomically. Filesystems that refuse to
  # replace an existing destination use a same-directory backup and rollback.
  if (rename_file(temporary, target)) return(invisible(target))
  backup <- tempfile(
    pattern = ".shennong-result-bundle-backup-",
    tmpdir = dirname(target),
    fileext = ".json"
  )
  if (!rename_file(target, backup)) {
    stop(
      "Could not replace Result Bundle; the existing destination was not changed: ",
      target,
      call. = FALSE
    )
  }
  if (rename_file(temporary, target)) {
    remove_file(backup, force = TRUE)
    return(invisible(target))
  }
  restored <- rename_file(backup, target)
  if (!restored) {
    stop(
      "Could not replace Result Bundle or restore the destination. The previous ",
      "file remains at: ", backup,
      call. = FALSE
    )
  }
  stop(
    "Could not replace Result Bundle; the existing destination was restored.",
    call. = FALSE
  )
}

#' @rdname result_bundle
#' @export
sn_export_result_bundle <- function(bundle,
                                    path,
                                    pretty = TRUE,
                                    overwrite = FALSE) {
  sn_validate_result_bundle(bundle)
  path <- .sn_result_bundle_scalar(path, "path")
  if (!is.logical(pretty) || length(pretty) != 1L || is.na(pretty)) {
    stop("`pretty` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("`overwrite` must be TRUE or FALSE.", call. = FALSE)
  }
  if (dir.exists(path)) {
    stop("`path` must name a JSON file, not a directory.", call. = FALSE)
  }
  if (file.exists(path) && !isTRUE(overwrite)) {
    stop("Result Bundle destination already exists: ", path, call. = FALSE)
  }
  directory <- dirname(path)
  if (!dir.exists(directory) && !dir.create(directory, recursive = TRUE)) {
    stop("Could not create Result Bundle directory: ", directory, call. = FALSE)
  }
  directory <- normalizePath(directory, winslash = "/", mustWork = TRUE)
  target <- file.path(directory, basename(path))
  temporary <- tempfile(
    pattern = ".shennong-result-bundle-",
    tmpdir = directory,
    fileext = ".json"
  )
  on.exit(unlink(temporary, force = TRUE), add = TRUE)
  writeLines(
    as.character(.sn_result_bundle_json(bundle, pretty = pretty)),
    temporary,
    useBytes = TRUE
  )
  .sn_result_bundle_publish(temporary, target, overwrite = overwrite)
  normalized <- normalizePath(target, winslash = "/", mustWork = TRUE)
  output <- list(
    path = normalized,
    size_bytes = as.numeric(file.info(normalized)$size),
    digest = list(
      algorithm = "sha256",
      value = digest::digest(
        file = normalized,
        algo = "sha256",
        serialize = FALSE
      )
    )
  )
  class(output) <- c("sn_result_bundle_export", "list")
  invisible(output)
}
