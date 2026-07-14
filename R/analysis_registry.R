.sn_method_registry_dir <- function() {
  installed <- system.file("methods", package = "Shennong")
  candidates <- unique(c(
    if (nzchar(installed)) installed else character(),
    file.path(getwd(), "inst", "methods")
  ))
  candidates <- candidates[dir.exists(candidates)]
  if (length(candidates) == 0L) {
    stop("Shennong method registry files were not found.", call. = FALSE)
  }
  candidates[[1]]
}

.sn_method_registry_required_fields <- function() {
  c(
    "name", "task", "runtime", "implemented", "default", "description",
    "input_requirements", "outputs", "supports", "cpu_gpu", "install", "citation"
  )
}

.sn_read_method_registry_file <- function(path) {
  entries <- tryCatch(
    jsonlite::fromJSON(path, simplifyVector = FALSE),
    error = function(e) {
      stop(
        "Could not parse method registry '", basename(path), "': ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  if (!is.list(entries) || length(entries) == 0L) {
    stop("Method registry '", basename(path), "' must contain a non-empty list.", call. = FALSE)
  }

  required <- .sn_method_registry_required_fields()
  lapply(seq_along(entries), function(index) {
    entry <- entries[[index]]
    missing <- setdiff(required, names(entry))
    if (length(missing) > 0L) {
      stop(
        "Method registry '", basename(path), "' entry ", index,
        " is missing field(s): ", paste(missing, collapse = ", "),
        call. = FALSE
      )
    }
    entry$registry_file <- basename(path)
    entry
  })
}

.sn_method_registry <- function() {
  registry_dir <- .sn_method_registry_dir()
  paths <- sort(list.files(registry_dir, pattern = "\\.ya?ml$", full.names = TRUE))
  if (length(paths) == 0L) {
    stop("No method registry YAML files were found in ", registry_dir, ".", call. = FALSE)
  }

  entries <- unlist(lapply(paths, .sn_read_method_registry_file), recursive = FALSE)
  keys <- vapply(entries, function(entry) paste(entry$task, entry$name, sep = "::"), character(1))
  if (anyDuplicated(keys)) {
    duplicated_keys <- unique(keys[duplicated(keys)])
    stop(
      "Method registry contains duplicate task/method entries: ",
      paste(duplicated_keys, collapse = ", "),
      call. = FALSE
    )
  }
  entries
}

.sn_method_package_available <- function(package) {
  package <- as.character(package %||% "")
  if (!nzchar(package) || package %in% c("base", "Shennong")) {
    return(TRUE)
  }
  suppressWarnings(requireNamespace(package, quietly = TRUE))
}

.sn_method_availability <- function(entry) {
  implemented <- isTRUE(entry$implemented)
  runtime <- tolower(as.character(entry$runtime))
  package_available <- .sn_method_package_available(entry$package %||% "")
  executable <- as.character(entry$executable %||% "")
  executable_available <- if (nzchar(executable)) nzchar(Sys.which(executable)) else TRUE
  pixi_available <- if (identical(runtime, "pixi")) {
    isTRUE(sn_check_pixi(quiet = TRUE)$installed)
  } else {
    NA
  }
  available <- implemented && package_available && executable_available &&
    (!identical(runtime, "pixi") || isTRUE(pixi_available))

  reason <- if (!implemented) {
    "Backend is registered for roadmap discovery but is not implemented yet."
  } else if (!package_available) {
    paste0("Optional R package '", entry$package, "' is not installed.")
  } else if (!executable_available) {
    paste0("Required executable '", executable, "' is not discoverable on PATH.")
  } else if (identical(runtime, "pixi") && !isTRUE(pixi_available)) {
    "pixi is not installed or discoverable."
  } else {
    "Method is available."
  }

  list(
    available = available,
    package_available = package_available,
    executable_available = executable_available,
    pixi_available = pixi_available,
    reason = reason
  )
}

.sn_method_registry_row <- function(entry, include_status = TRUE) {
  status <- if (include_status) .sn_method_availability(entry) else list(
    available = NA,
    package_available = NA,
    pixi_available = NA,
    reason = NA_character_
  )
  tibble::tibble(
    task = as.character(entry$task),
    name = as.character(entry$name),
    description = as.character(entry$description),
    runtime = as.character(entry$runtime),
    package = as.character(entry$package %||% NA_character_),
    environment = as.character(entry$environment %||% NA_character_),
    implemented = isTRUE(entry$implemented),
    default = isTRUE(entry$default),
    available = isTRUE(status$available),
    availability_reason = as.character(status$reason),
    cpu_gpu = as.character(entry$cpu_gpu),
    install_action = as.character(entry$install),
    input_requirements = list(entry$input_requirements),
    outputs = list(entry$outputs),
    supports = list(entry$supports),
    citation = as.character(entry$citation),
    registry_file = as.character(entry$registry_file)
  )
}

#' List registered Shennong analysis methods
#'
#' Reads the shipped method registry and reports whether each backend is
#' implemented and currently available. Registered roadmap methods remain
#' discoverable even when their optional dependency or pixi runtime has not
#' been installed.
#'
#' @param task Optional task name such as \code{"trajectory"},
#'   \code{"annotation"}, or \code{"bulk"}.
#' @param available Optional logical filter. Use \code{TRUE} for methods that
#'   can run in the current session and \code{FALSE} for unavailable methods.
#'
#' @return A tibble with one row per task/method pair.
#'
#' @examples
#' sn_list_methods("trajectory")
#' sn_list_methods(available = TRUE)
#'
#' @export
sn_list_methods <- function(task = NULL, available = NULL) {
  entries <- .sn_method_registry()
  table <- dplyr::bind_rows(lapply(entries, .sn_method_registry_row))
  if (!is_null(task)) {
    requested_tasks <- tolower(as.character(task))
    table <- dplyr::filter(table, .data$task %in% .env$requested_tasks)
  }
  if (!is_null(available)) {
    if (!is.logical(available) || length(available) != 1L || is.na(available)) {
      stop("`available` must be NULL, TRUE, or FALSE.", call. = FALSE)
    }
    table <- dplyr::filter(table, .data$available == available)
  }
  dplyr::arrange(table, .data$task, dplyr::desc(.data$default), .data$name)
}

#' Report the status of a registered Shennong method
#'
#' @param method A registered method name.
#' @param task Optional task used to disambiguate a method registered for more
#'   than one workflow.
#'
#' @return A named list describing availability, runtime, environment,
#'   installation action, requirements, outputs, and citation.
#'
#' @examples
#' sn_method_status("slingshot", task = "trajectory")
#'
#' @export
sn_method_status <- function(method, task = NULL) {
  if (!is.character(method) || length(method) != 1L || !nzchar(method)) {
    stop("`method` must be a non-empty character scalar.", call. = FALSE)
  }
  table <- sn_list_methods(task = task)
  hit <- table[tolower(table$name) == tolower(method), , drop = FALSE]
  if (nrow(hit) == 0L) {
    tasks <- sort(unique(table$task))
    stop(
      "Unknown Shennong method '", method, "'",
      if (!is_null(task)) paste0(" for task '", task, "'") else "",
      ". Available task(s): ", paste(tasks, collapse = ", "),
      call. = FALSE
    )
  }
  if (nrow(hit) > 1L) {
    stop(
      "Method '", method, "' is registered for multiple tasks: ",
      paste(hit$task, collapse = ", "), ". Supply `task`.",
      call. = FALSE
    )
  }
  row <- hit[1, , drop = FALSE]
  list(
    method = row$name[[1]],
    task = row$task[[1]],
    available = row$available[[1]],
    implemented = row$implemented[[1]],
    reason = row$availability_reason[[1]],
    runtime = row$runtime[[1]],
    package = row$package[[1]],
    environment = row$environment[[1]],
    default = row$default[[1]],
    install_action = row$install_action[[1]],
    input_requirements = row$input_requirements[[1]],
    outputs = row$outputs[[1]],
    supports = row$supports[[1]],
    cpu_gpu = row$cpu_gpu[[1]],
    citation = row$citation[[1]]
  )
}
