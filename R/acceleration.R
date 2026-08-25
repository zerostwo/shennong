# Optional ShennongOpt acceleration ------------------------------------------

.sn_acceleration_option_name <- "shennong.acceleration"

# Legacy wrapper patch keys -> ShennongOpt patch names. Keys without a mapping
# run upstream unaccelerated until a matching ShennongOpt patch exists.
.sn_acceleration_patch_map <- c(
  "coralysis" = "coralysis",
  "decontx" = "decontx",
  "decontx_standalone" = "decontx",
  "lisi" = "lisi",
  "rogue" = "rogue",
  "scdblfinder" = "scdblfinder",
  "scran" = "scran",
  "seurat" = "seurat",
  "ucell" = "ucell"
)

.sn_acceleration_available <- function() {
  requireNamespace("ShennongOpt", quietly = TRUE)
}

.sn_acceleration_legacy_exports <- c(
  sn_list_accelerations = "sno_list_patches",
  sn_check_acceleration = "sno_status",
  sn_enable_acceleration = "sno_activate",
  sn_disable_acceleration = "sno_deactivate",
  sn_is_acceleration_disabled = "sno_is_disabled",
  sn_with_acceleration_disabled = "sno_with_disabled"
)

.sn_acceleration_get_export <- function(fun) {
  exports <- getNamespaceExports("ShennongOpt")
  legacy <- unname(.sn_acceleration_legacy_exports[fun])
  resolved <- if (fun %in% exports) {
    fun
  } else {
    legacy
  }
  if (length(resolved) != 1L || is.na(resolved) || !resolved %in% exports) {
    stop(
      "ShennongOpt does not export the required acceleration API `", fun, "()`.",
      call. = FALSE
    )
  }
  base::getExportedValue("ShennongOpt", resolved)
}

.sn_acceleration_call <- function(fun, ...) {
  do.call(.sn_acceleration_get_export(fun), list(...))
}

.sn_acceleration_registered_patches <- function(installed = TRUE) {
  if (!.sn_acceleration_available()) {
    return(character(0))
  }
  tryCatch(
    as.character(.sn_acceleration_call("sn_list_accelerations", installed = installed)),
    error = function(e) character(0)
  )
}

.sn_map_acceleration_patches <- function(patches) {
  patches <- unique(tolower(trimws(as.character(patches))))
  patches <- patches[!is.na(patches) & nzchar(patches)]
  mapped <- unname(.sn_acceleration_patch_map[patches])
  supported <- unique(mapped[!is.na(mapped)])
  registered <- .sn_acceleration_registered_patches(installed = FALSE)
  installable <- .sn_acceleration_registered_patches(installed = TRUE)
  supported <- intersect(supported, registered)
  list(
    requested = patches,
    supported = intersect(supported, installable),
    unsupported = c(
      setdiff(patches, names(.sn_acceleration_patch_map)),
      setdiff(supported, installable)
    )
  )
}

.sn_acceleration_default_enabled <- function() {
  if (!isTRUE(getOption(.sn_acceleration_option_name, TRUE))) {
    return(FALSE)
  }
  disabled <- Sys.getenv("SHENNONG_ACCELERATION_DISABLED", unset = "")
  truthy <- tolower(trimws(disabled)) %in% c("1", "true", "t", "yes", "y", "on")
  !any(truthy)
}

#' Check the current acceleration status
#'
#' Reports which ShennongOpt patches are registered, installable, and active in
#' the current session.
#'
#' @return A named character vector mapping patch name to `"active"` or
#'   `"inactive"`, or an empty vector when ShennongOpt is not installed.
#'
#' @examples
#' \dontrun{sn_check_acceleration()}
#'
#' @export
sn_check_acceleration <- function() {
  if (!.sn_acceleration_available()) {
    return(stats::setNames(character(0), character(0)))
  }
  tryCatch(.sn_acceleration_call("sn_check_acceleration"), error = function(e) {
    stats::setNames(character(0), character(0))
  })
}

#' Activate ShennongOpt acceleration patches
#'
#' Rebinds accelerated implementations into their upstream namespaces.
#' Activation is process-local and idempotent; pass patch names to activate a
#' subset, or nothing to activate every registered patch whose upstream package
#' is installed. See `ShennongOpt::sn_list_accelerations()` for available names.
#'
#' @param name Optional patch name(s) to activate.
#' @return Invisible named logical vector of activation results.
#'
#' @examples
#' \dontrun{
#' sn_enable_acceleration()
#' sn_disable_acceleration()
#' }
#'
#' @export
sn_enable_acceleration <- function(name = NULL) {
  if (!.sn_acceleration_available()) {
    stop("ShennongOpt is not installed; acceleration is unavailable.", call. = FALSE)
  }
  invisible(.sn_acceleration_call("sn_enable_acceleration", name))
}

#' Deactivate ShennongOpt acceleration patches
#'
#' Restores upstream implementations. Idempotent.
#'
#' @param name Optional patch name(s) to deactivate; default restores all.
#' @return Invisible NULL.
#'
#' @export
sn_disable_acceleration <- function(name = NULL) {
  if (!.sn_acceleration_available()) {
    return(invisible(NULL))
  }
  invisible(.sn_acceleration_call("sn_disable_acceleration", name))
}

#' Run an expression with selected acceleration patches active
#'
#' Activates the requested ShennongOpt patches for the duration of
#' \code{expr} and restores the previous activation state afterwards. Patches
#' without a ShennongOpt counterpart run unaccelerated and are recorded as
#' suppressed in the acceleration provenance context.
#'
#' @param expr Expression to evaluate.
#' @param name Patch name(s) to activate.
#' @return The value of \code{expr}.
#'
#' @examples
#' \dontrun{sn_with_acceleration(obj <- Seurat::RunPCA(obj), name = "seurat")}
#'
#' @export
sn_with_acceleration <- function(expr, name = NULL) {
  .sn_with_acceleration_impl(expr, name)
}

.sn_with_acceleration_impl <- function(expr, patches) {
  mapped <- .sn_map_acceleration_patches(patches)
  enabled <- .sn_acceleration_default_enabled() &&
    .sn_acceleration_available() &&
    length(mapped$supported) > 0L
  previously_active <- character(0)
  activated <- character(0)
  if (enabled) {
    previously_active <- names(which(unname(
      .sn_acceleration_call("sn_check_acceleration")
    ) == "active"))
    newly <- setdiff(mapped$supported, previously_active)
    if (length(newly) > 0L) {
      result <- suppressMessages(tryCatch(
        .sn_acceleration_call("sn_enable_acceleration", newly),
        error = function(e) NULL
      ))
      result <- if (is.null(result)) logical(0) else as.logical(result)
      names(result) <- if (is.null(names(result))) newly else names(result)
      activated <- names(result)[result & !is.na(result)]
    }
  }
  used <- intersect(mapped$supported, union(previously_active, activated))
  suppressed <- c(mapped$unsupported, setdiff(mapped$supported, used))
  .sn_record_acceleration_usage(used)
  .sn_record_acceleration_suppression(suppressed)
  if (length(activated) > 0L) {
    on.exit({
      .sn_acceleration_call("sn_disable_acceleration", setdiff(activated, previously_active))
    }, add = TRUE)
  }
  force(expr)
}

.sn_with_default_acceleration <- function(expr, patches = character(), ...) {
  # `...` accepts legacy provenance hints (strict, operation) from callers.
  .sn_with_acceleration_impl(expr, patches)
}

.sn_with_default_seurat_acceleration <- function(expr, ...) {
  dots <- list(...)
  patches <- unique(c("seurat", as.character(dots$patches)))
  .sn_with_acceleration_impl(expr, patches)
}

.sn_with_explicit_acceleration_or_disabled <- function(expr, patches) {
  .sn_with_explicit_acceleration_context({
    .sn_with_acceleration_impl(expr, patches)
  }, patches)
}

.sn_with_acceleration_provenance_context <- function(expr, patches) {
  option_name <- "shennong.acceleration.provenance_context"
  current_options <- options()
  previous_present <- option_name %in% names(current_options)
  previous <- current_options[[option_name]]
  owns_context <- !is.environment(previous)
  context <- if (owns_context) new.env(parent = emptyenv()) else previous
  if (is_null(context$patches)) {
    context$patches <- character()
  }
  if (is_null(context$used_patches)) {
    context$used_patches <- character()
  }
  if (is_null(context$suppressed_patches)) {
    context$suppressed_patches <- character()
  }
  context$patches <- union(context$patches, patches)

  if (owns_context) {
    options(stats::setNames(list(context), option_name))
    on.exit({
      value <- if (previous_present) previous else NULL
      options(stats::setNames(list(value), option_name))
    }, add = TRUE)
  }

  force(expr)
}

.sn_record_acceleration_usage <- function(patches) {
  context <- getOption("shennong.acceleration.provenance_context")
  if (is.environment(context)) {
    context$used_patches <- union(
      context$used_patches %||% character(),
      intersect(as.character(patches), context$patches %||% character())
    )
  }
  if (exists(".sn_usage_record_acceleration", mode = "function")) {
    .sn_usage_record_acceleration(patches)
  }
  invisible(NULL)
}

.sn_record_acceleration_suppression <- function(patches) {
  context <- getOption("shennong.acceleration.provenance_context")
  if (is.environment(context)) {
    context$suppressed_patches <- union(
      context$suppressed_patches %||% character(),
      intersect(as.character(patches), context$patches %||% character())
    )
  }
  invisible(NULL)
}

.sn_acceleration_explicit_patches <- function() {
  patches <- getOption("shennong.acceleration.explicit_patches", character())
  unique(tolower(trimws(as.character(patches))))
}

.sn_with_explicit_acceleration_context <- function(expr, patches) {
  option_name <- "shennong.acceleration.explicit_patches"
  current_options <- options()
  previous_present <- option_name %in% names(current_options)
  previous <- current_options[[option_name]]
  patches <- unique(tolower(trimws(as.character(patches))))
  options(stats::setNames(list(union(
    .sn_acceleration_explicit_patches(),
    patches
  )), option_name))
  on.exit({
    value <- if (previous_present) previous else NULL
    options(stats::setNames(list(value), option_name))
  }, add = TRUE)
  force(expr)
}

#' (internal) Evaluate an expression on the unaccelerated upstream code path
#'
#' @keywords internal
.sn_with_acceleration_disabled <- function(expr) {
  if (!.sn_acceleration_available() || isTRUE(.sn_acceleration_call("sn_is_acceleration_disabled"))) {
    return(force(expr))
  }
  with_disabled <- .sn_acceleration_get_export("sn_with_acceleration_disabled")
  with_disabled(expr)
}

.sn_acceleration_effective_active_patches <- function() {
  status <- sn_check_acceleration()
  if (length(status) == 0L) {
    return(character(0))
  }
  disabled <- tryCatch(
    .sn_acceleration_call("sn_is_acceleration_disabled"),
    error = function(e) FALSE
  )
  if (isTRUE(disabled)) character(0) else names(status)[status == "active"]
}

.sn_acceleration_provenance <- function() {
  context <- getOption("shennong.acceleration.provenance_context")
  if (!is.environment(context)) {
    return(list())
  }
  patches <- context$patches %||% character()
  if (length(patches) == 0L) {
    return(list())
  }
  used <- intersect(context$used_patches %||% character(), patches)
  explicitly_suppressed <- intersect(context$suppressed_patches %||% character(), patches)
  suppressed <- union(explicitly_suppressed, setdiff(patches, used))
  list(
    patches = patches,
    used_patches = used,
    suppressed_patches = suppressed
  )
}
