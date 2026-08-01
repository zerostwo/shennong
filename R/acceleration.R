# Optional AutoZyme acceleration ------------------------------------------------

.sn_autozyme_expected_version <- "0.3.1"
.sn_autozyme_expected_sha <- "718541d9489596c7c1d75f52e9b3a8b2a429d1f9"
.sn_autozyme_default_patches <- c(
  "cellchat",
  "clusterprofiler",
  "fgsea",
  "nichenetr",
  "scdblfinder",
  "seurat",
  "seurat_joinlayers",
  "seurat_merge",
  "soupx",
  "tradeseq",
  "wgcna"
)
.sn_autozyme_clusterprofiler_cache <- new.env(parent = emptyenv())

# This manifest is intentionally conservative. Versions are the exact upstream
# versions against which the pinned AutoZyme revision declares its patches.
# `approximate = TRUE` means the accelerated implementation can change more than
# floating-point round-off and therefore requires a second, explicit opt-in.
.sn_autozyme_patch_manifest <- list(
  bayesspace = list(
    upstream = "BayesSpace", versions = "1.21.2",
    equivalence = "approximate", approximate = TRUE
  ),
  cellchat = list(
    upstream = "CellChat", versions = "2.2.0.9001",
    equivalence = "exact_scoped", approximate = FALSE
  ),
  clusterprofiler = list(
    upstream = "clusterProfiler", versions = "4.20.0",
    equivalence = "exact_scoped", approximate = FALSE,
    source_sha256 = "97711d3821dabfe497c3cecb9932c6e1829af01ba86e5289487996de1b544e77"
  ),
  decontx = list(
    upstream = "celda", versions = "1.24.0",
    equivalence = "approximate", approximate = TRUE
  ),
  fgsea = list(
    upstream = "fgsea", versions = "1.34.2",
    equivalence = "numeric_tolerance", approximate = FALSE
  ),
  infercnv = list(
    upstream = "infercnv", versions = "1.24.0",
    equivalence = "approximate", approximate = TRUE
  ),
  maftools = list(
    upstream = "maftools", versions = "2.24.0",
    equivalence = "exact", approximate = FALSE
  ),
  mast = list(
    upstream = "MAST", versions = "1.35.2",
    equivalence = "approximate", approximate = TRUE
  ),
  nichenetr = list(
    upstream = "nichenetr", versions = "2.2.1.1",
    equivalence = "exact_scoped", approximate = FALSE
  ),
  rctd = list(
    upstream = "spacexr", versions = "2.2.1",
    equivalence = "approximate", approximate = TRUE
  ),
  scriabin = list(
    upstream = "scriabin", versions = "0.0.0.9000",
    equivalence = "exact_scoped", approximate = FALSE
  ),
  scdblfinder = list(
    upstream = "scDblFinder", versions = "1.27.6",
    equivalence = "exact_scoped", approximate = FALSE,
    source_sha256 = "50f128d3e3f40ea3a25cdebef2b5c522a3cb6db1edf3c7c2d26646bdc9431423"
  ),
  seurat = list(
    upstream = "Seurat", versions = c("5.2.1", "5.4.0"),
    equivalence = "numeric_tolerance_scoped", approximate = FALSE
  ),
  seurat_joinlayers = list(
    upstream = "SeuratObject", versions = c("5.4.0", "5.4.0.9001"),
    equivalence = "exact_scoped", approximate = FALSE,
    source_sha256 = "b2379438658f8c39f8073b214ba446b72c3d8d2a059921d163eaa9abf9c761dc"
  ),
  seurat_merge = list(
    upstream = "SeuratObject", versions = c("5.4.0", "5.4.0.9001"),
    equivalence = "exact_scoped", approximate = FALSE,
    source_sha256 = "cbb4af06fa7ad8994af70627b7420154137c7d743bae0e06ce7c942fcae7f78d"
  ),
  slingshot = list(
    upstream = "slingshot", versions = "2.16.0",
    equivalence = "approximate", approximate = TRUE
  ),
  soupx = list(
    upstream = "SoupX", versions = "1.6.2",
    equivalence = "exact_scoped", approximate = FALSE,
    source_sha256 = "757185c34a34f5500733d3f5044023d9d22e27bf0de4c7670299225c1f70dc85"
  ),
  tradeseq = list(
    upstream = "tradeSeq", versions = "1.22.0",
    equivalence = "exact_scoped", approximate = FALSE
  ),
  vegan = list(
    upstream = "vegan", versions = "2.8.0",
    equivalence = "numeric_tolerance_scoped", approximate = FALSE
  ),
  wgcna = list(
    upstream = "WGCNA", versions = "1.74",
    equivalence = "numeric_tolerance_scoped", approximate = FALSE
  )
)

.sn_validate_autozyme_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", arg), call. = FALSE)
  }
  x
}

.sn_validate_autozyme_patches <- function(patches) {
  if (!is.character(patches) || length(patches) == 0L || anyNA(patches)) {
    stop("`patches` must be a non-empty character vector.", call. = FALSE)
  }

  patches <- tolower(trimws(patches))
  if (any(!nzchar(patches))) {
    stop("`patches` cannot contain empty names.", call. = FALSE)
  }

  unknown <- setdiff(patches, names(.sn_autozyme_patch_manifest))
  if (length(unknown) > 0L) {
    stop(
      sprintf(
        "Unsupported AutoZyme patch(es): %s. Supported patches are: %s.",
        paste(unknown, collapse = ", "),
        paste(names(.sn_autozyme_patch_manifest), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  unique(patches)
}

.sn_autozyme_is_installed <- function(package = "autozyme") {
  path <- tryCatch(
    base::find.package(package, quiet = TRUE),
    error = function(e) character(0)
  )
  length(path) == 1L && nzchar(path)
}

.sn_autozyme_installed_version <- function(package) {
  if (!.sn_autozyme_is_installed(package)) {
    return(NA_character_)
  }

  tryCatch(
    as.character(utils::packageVersion(package)),
    error = function(e) NA_character_
  )
}

.sn_autozyme_description <- function() {
  absent <- list(version = NA_character_, remote_sha = NA_character_)
  if (!.sn_autozyme_is_installed("autozyme")) {
    return(absent)
  }

  desc <- tryCatch(
    utils::packageDescription("autozyme"),
    error = function(e) NULL
  )
  if (is.null(desc)) {
    return(absent)
  }

  version <- desc[["Version"]]
  remote_sha <- desc[["RemoteSha"]]
  list(
    version = if (is.null(version) || !nzchar(version)) NA_character_ else version,
    remote_sha = if (is.null(remote_sha) || !nzchar(remote_sha)) NA_character_ else remote_sha
  )
}

.sn_autozyme_vendored_patch_path <- function(patch) {
  system.file(
    "autozyme", "patches", patch, "patch.R",
    package = "Shennong"
  )
}

.sn_autozyme_source_matches <- function(path, expected) {
  if (!nzchar(path) || !file.exists(path) ||
      is_null(expected) || !nzchar(expected)) {
    return(FALSE)
  }

  actual <- tryCatch(
    digest::digest(file = path, algo = "sha256"),
    error = function(error) NA_character_
  )
  !is.na(actual) && identical(actual, expected)
}

.sn_autozyme_patch_source_status <- function(patch, spec) {
  autozyme_path <- system.file(
    "patches", patch, "patch.R",
    package = "autozyme"
  )
  expected <- spec$source_sha256
  autozyme_source_match <- .sn_autozyme_source_matches(
    autozyme_path,
    expected
  )
  vendored_path <- .sn_autozyme_vendored_patch_path(patch)
  bundled_by_shennong <- nzchar(vendored_path) && file.exists(vendored_path)
  vendored_source_match <- .sn_autozyme_source_matches(
    vendored_path,
    expected
  )

  registered <- FALSE
  if (.sn_autozyme_is_installed("autozyme")) {
    registered <- tryCatch(
      patch %in% .sn_autozyme_call("list_patches"),
      error = function(error) FALSE
    )
  }

  provider <- if (isTRUE(autozyme_source_match)) {
    "autozyme"
  } else if (isTRUE(vendored_source_match)) {
    "shennong"
  } else if (registered) {
    "autozyme"
  } else {
    NA_character_
  }

  list(
    registered = registered,
    bundled_by_shennong = bundled_by_shennong,
    source_match = isTRUE(autozyme_source_match) ||
      isTRUE(vendored_source_match),
    vendored_source_match = vendored_source_match,
    provider = provider
  )
}

.sn_register_vendored_autozyme_patch <- function(patch) {
  patch_path <- .sn_autozyme_vendored_patch_path(patch)
  if (!nzchar(patch_path) || !file.exists(patch_path)) {
    stop(
      sprintf("Shennong does not bundle AutoZyme patch `%s`.", patch),
      call. = FALSE
    )
  }
  expected <- .sn_autozyme_patch_manifest[[patch]]$source_sha256
  if (!.sn_autozyme_source_matches(patch_path, expected)) {
    stop(
      sprintf(
        "Shennong's vendored AutoZyme patch `%s` failed source validation.",
        patch
      ),
      call. = FALSE
    )
  }

  patch_env <- new.env(parent = asNamespace("Shennong"))
  patch_env$register_patch <- base::getExportedValue(
    "autozyme",
    "register_patch"
  )
  sys.source(patch_path, envir = patch_env)

  available <- .sn_autozyme_call("list_patches")
  if (!patch %in% available) {
    stop(
      sprintf("Shennong's vendored AutoZyme patch `%s` did not register.", patch),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.sn_autozyme_call <- function(fun, ...) {
  do.call(base::getExportedValue("autozyme", fun), list(...))
}

.sn_autozyme_namespace_loaded <- function() {
  base::isNamespaceLoaded("autozyme")
}

.sn_autozyme_status_vector <- function() {
  if (!.sn_autozyme_is_installed("autozyme") ||
      !.sn_autozyme_namespace_loaded()) {
    return(character(0))
  }

  status <- .sn_autozyme_call("status")
  if (!is.character(status) || is.null(names(status))) {
    stop("AutoZyme returned an invalid patch status.", call. = FALSE)
  }
  status
}

.sn_autozyme_active_patches <- function() {
  status <- .sn_autozyme_status_vector()
  names(status)[!is.na(status) & status == "active"]
}

.sn_autozyme_effective_active_patches <- function() {
  active <- .sn_autozyme_active_patches()
  if (length(active) == 0L || !.sn_autozyme_namespace_loaded()) {
    return(active)
  }
  disabled <- isTRUE(tryCatch(
    .sn_autozyme_call("is_disabled"),
    error = function(error) FALSE
  ))
  if (disabled) character(0) else active
}

.sn_autozyme_default_enabled <- function() {
  if (!isTRUE(getOption("shennong.autozyme", TRUE))) {
    return(FALSE)
  }

  disabled <- Sys.getenv(
    c("AUTOZYME_DISABLED", "AUTOZYME_DISABLE"),
    unset = ""
  )
  truthy <- tolower(trimws(disabled)) %in% c("1", "true", "t", "yes", "y", "on")
  !any(truthy)
}

.sn_capture_autozyme_future_option <- function() {
  current <- options()
  name <- "future.globals.maxSize"
  list(
    present = name %in% names(current),
    value = current[[name]]
  )
}

.sn_restore_autozyme_future_option <- function(snapshot) {
  value <- if (isTRUE(snapshot$present)) snapshot$value else NULL
  options(stats::setNames(list(value), "future.globals.maxSize"))
  invisible(NULL)
}

.sn_autozyme_rollback_details <- function(rollback) {
  details <- character()
  if (length(rollback$deactivation_errors) > 0L) {
    details <- c(
      details,
      paste0(
        "deactivation error(s): ",
        paste(rollback$deactivation_errors, collapse = "; ")
      )
    )
  }
  if (!is_null(rollback$verification_error)) {
    details <- c(
      details,
      paste0(
        "post-rollback state could not be verified: ",
        rollback$verification_error
      )
    )
  }
  if (length(rollback$residual_active) > 0L) {
    details <- c(
      details,
      paste0(
        "residual active patch(es): ",
        paste(rollback$residual_active, collapse = ", ")
      )
    )
  }
  details
}

.sn_with_autozyme_provenance_context <- function(expr, patches) {
  option_name <- "shennong.autozyme.provenance_context"
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

.sn_record_autozyme_usage <- function(patches) {
  context <- getOption("shennong.autozyme.provenance_context")
  if (is.environment(context)) {
    context$used_patches <- union(
      context$used_patches %||% character(),
      intersect(as.character(patches), context$patches %||% character())
    )
  }
  invisible(NULL)
}

.sn_record_autozyme_suppression <- function(patches) {
  context <- getOption("shennong.autozyme.provenance_context")
  if (is.environment(context)) {
    context$suppressed_patches <- union(
      context$suppressed_patches %||% character(),
      intersect(as.character(patches), context$patches %||% character())
    )
  }
  invisible(NULL)
}

.sn_with_default_autozyme <- function(expr, patches, strict = TRUE) {
  patches <- unique(tolower(trimws(as.character(patches))))
  patches <- intersect(patches, .sn_autozyme_default_patches)
  .sn_with_autozyme_provenance_context(
    .sn_with_default_autozyme_impl(expr, patches = patches, strict = strict),
    patches = patches
  )
}

.sn_with_default_autozyme_impl <- function(expr, patches, strict = TRUE) {
  patches <- unique(tolower(trimws(as.character(patches))))
  patches <- intersect(patches, .sn_autozyme_default_patches)
  strict <- .sn_validate_autozyme_flag(strict, "strict")
  future_option <- .sn_capture_autozyme_future_option()
  future_option_pending <- TRUE
  restore_before_analysis <- function() {
    if (isTRUE(future_option_pending)) {
      .sn_restore_autozyme_future_option(future_option)
      future_option_pending <<- FALSE
    }
    invisible(NULL)
  }
  force_analysis <- function() {
    restore_before_analysis()
    force(expr)
  }
  on.exit(restore_before_analysis(), add = TRUE)

  if (length(patches) == 0L) {
    return(force_analysis())
  }
  if (!.sn_autozyme_default_enabled() ||
      !.sn_autozyme_is_installed("autozyme")) {
    .sn_record_autozyme_usage(intersect(
      patches,
      tryCatch(.sn_autozyme_effective_active_patches(), error = function(error) character())
    ))
    return(force_analysis())
  }

  check_error <- NULL
  checks <- tryCatch(
    sn_check_autozyme(
      patches = patches,
      strict = strict,
      allow_approximate = FALSE
    ),
    error = function(error) {
      check_error <<- error
      NULL
    }
  )
  if (!is_null(check_error)) {
    warning(
      sprintf(
        "Could not check automatic AutoZyme acceleration; continuing without it: %s",
        conditionMessage(check_error)
      ),
      call. = FALSE
    )
    return(force_analysis())
  }

  eligible <- checks$patch[checks$eligible & !checks$active]
  if (length(eligible) == 0L) {
    .sn_record_autozyme_usage(intersect(
      patches,
      tryCatch(.sn_autozyme_effective_active_patches(), error = function(error) character())
    ))
    return(force_analysis())
  }

  before_error <- NULL
  before <- tryCatch(
    .sn_autozyme_active_patches(),
    error = function(error) {
      before_error <<- error
      character(0)
    }
  )
  if (!is_null(before_error)) {
    warning(
      sprintf(
        "Could not inspect AutoZyme state; continuing without automatic acceleration: %s",
        conditionMessage(before_error)
      ),
      call. = FALSE
    )
    .sn_record_autozyme_usage(intersect(
      patches,
      tryCatch(.sn_autozyme_effective_active_patches(), error = function(error) character())
    ))
    return(force_analysis())
  }

  activation_error <- NULL
  tryCatch(
    suppressPackageStartupMessages(
      suppressMessages(
        sn_enable_autozyme(
          patches = eligible,
          strict = strict,
          allow_approximate = FALSE
        )
      )
    ),
    error = function(error) {
      activation_error <<- error
      NULL
    }
  )
  restore_before_analysis()
  if (!is_null(activation_error)) {
    after_error <- NULL
    after <- tryCatch(
      .sn_autozyme_active_patches(),
      error = function(error) {
        after_error <<- error
        character(0)
      }
    )
    residual <- if (is_null(after_error)) {
      setdiff(intersect(eligible, after), before)
    } else {
      character(0)
    }
    if (!is_null(after_error) || length(residual) > 0L) {
      details <- if (!is_null(after_error)) {
        sprintf("state verification failed: %s", conditionMessage(after_error))
      } else {
        sprintf("residual active patch(es): %s", paste(residual, collapse = ", "))
      }
      stop(
        sprintf(
          paste0(
            "Automatic AutoZyme activation failed and its state could not be ",
            "safely restored (%s): %s"
          ),
          details,
          conditionMessage(activation_error)
        ),
        call. = FALSE
      )
    }
    warning(
      sprintf(
        "Automatic AutoZyme activation failed; continuing without it: %s",
        conditionMessage(activation_error)
      ),
      call. = FALSE
    )
    .sn_record_autozyme_usage(intersect(
      patches,
      tryCatch(.sn_autozyme_effective_active_patches(), error = function(error) character())
    ))
    return(force_analysis())
  }

  newly_activated <- setdiff(eligible, before)
  active_for_call <- intersect(
    patches,
    .sn_autozyme_effective_active_patches()
  )
  .sn_record_autozyme_usage(active_for_call)
  if (length(active_for_call) > 0L) {
    .sn_log_info(sprintf(
      "[AutoZyme] Acceleration enabled for this call (patches: %s).",
      paste(active_for_call, collapse = ", ")
    ))
  }
  cleanup_pending <- length(newly_activated) > 0L
  on.exit({
    if (isTRUE(cleanup_pending)) {
      cleanup_pending <- FALSE
      emergency_rollback <- .sn_autozyme_rollback(newly_activated)
      if (!isTRUE(emergency_rollback$complete)) {
        stop(
          paste0(
            "Automatic AutoZyme scope exited before normal cleanup and could ",
            "not restore its prior state: ",
            paste(.sn_autozyme_rollback_details(emergency_rollback), collapse = "; "),
            "."
          ),
          call. = FALSE
        )
      } else if (length(emergency_rollback$deactivation_errors) > 0L) {
        warning(
          paste0(
            "Automatic AutoZyme scope restored state with diagnostics: ",
            paste(emergency_rollback$deactivation_errors, collapse = "; "),
            "."
          ),
          call. = FALSE
        )
      }
    }
  }, add = TRUE)
  value <- NULL
  analysis_error <- NULL
  tryCatch(
    value <- force_analysis(),
    error = function(error) {
      analysis_error <<- error
      NULL
    }
  )
  analysis_future_option <- .sn_capture_autozyme_future_option()
  on.exit(
    .sn_restore_autozyme_future_option(analysis_future_option),
    add = TRUE,
    after = TRUE
  )

  rollback <- .sn_autozyme_rollback(newly_activated)
  cleanup_pending <- FALSE
  rollback_details <- .sn_autozyme_rollback_details(rollback)
  if (!isTRUE(rollback$complete)) {
    message <- paste0(
      "Automatic AutoZyme scope could not restore its prior state: ",
      paste(rollback_details, collapse = "; "),
      "."
    )
    if (!is_null(analysis_error)) {
      message <- paste0(
        message,
        " The analytical call also failed: ",
        conditionMessage(analysis_error),
        "."
      )
    }
    stop(message, call. = FALSE)
  }
  if (length(rollback$deactivation_errors) > 0L) {
    warning(
      paste0(
        "Automatic AutoZyme state was restored with diagnostics: ",
        paste(rollback$deactivation_errors, collapse = "; "),
        "."
      ),
      call. = FALSE
    )
  }
  if (!is_null(analysis_error)) {
    stop(analysis_error)
  }

  value
}

.sn_autozyme_with_disabled <- function() {
  base::getExportedValue("autozyme", "with_disabled")
}

.sn_with_autozyme_disabled <- function(expr) {
  if (!.sn_autozyme_is_installed("autozyme") ||
      !.sn_autozyme_namespace_loaded()) {
    return(force(expr))
  }

  with_disabled <- tryCatch(
    .sn_autozyme_with_disabled(),
    error = function(error) {
      stop(
        paste0(
          "AutoZyme is loaded, but its patches could not be suspended safely ",
          "for this call: ", conditionMessage(error)
        ),
        call. = FALSE
      )
    }
  )

  with_disabled(expr)
}

.sn_seurat_uses_bpcells <- function(object, assay = NULL) {
  if (!inherits(object, "Seurat")) {
    return(TRUE)
  }

  assays <- assay %||% names(object@assays)
  if (!is.character(assays) || length(assays) == 0L || anyNA(assays) ||
      any(!assays %in% names(object@assays))) {
    return(TRUE)
  }

  for (current_assay in unique(assays)) {
    layers <- tryCatch(
      SeuratObject::Layers(object[[current_assay]]),
      error = function(error) NULL
    )
    if (is_null(layers)) {
      return(TRUE)
    }
    for (layer in layers) {
      matrix <- tryCatch(
        SeuratObject::LayerData(
          object = object,
          assay = current_assay,
          layer = layer
        ),
        error = function(error) NULL
      )
      if (is_null(matrix) || .sn_is_iterable_matrix(matrix)) {
        return(TRUE)
      }
    }
  }

  FALSE
}

.sn_with_default_seurat_autozyme <- function(expr, object, assay = NULL) {
  if (.sn_seurat_uses_bpcells(object = object, assay = assay)) {
    .sn_record_autozyme_suppression("seurat")
    active <- tryCatch(
      .sn_autozyme_effective_active_patches(),
      error = function(error) character()
    )
    if ("seurat" %in% active) {
      return(.sn_with_autozyme_disabled(expr))
    }
    return(.sn_with_default_autozyme(
      expr,
      patches = c("seurat_joinlayers", "seurat_merge")
    ))
  }

  .sn_with_default_autozyme(
    expr,
    patches = c("seurat", "seurat_joinlayers", "seurat_merge"),
    strict = FALSE
  )
}

#' Check optional AutoZyme acceleration compatibility
#'
#' Reports whether selected AutoZyme patches can be enabled. This check is
#' side-effect free: it does not
#' install packages, activate patches, prepare Python, or change thread counts.
#'
#' @param patches Character vector of AutoZyme patch names to inspect or manage.
#' @param strict If `TRUE`, require the installed upstream package version to
#'   exactly match a version validated by the pinned AutoZyme revision.
#' @param allow_approximate If `TRUE`, permit patches whose documented fast
#'   path can change analytical results beyond floating-point tolerance.
#'
#' @return A data frame with installation, version, activation, and eligibility
#'   information for each requested patch.
#'
#' @details
#' Shennong lazily activates compatible, non-approximate patches for the scope
#' of an integrated workflow call. Automatic activation covers CellChat,
#' clusterProfiler, fgsea, NicheNet, scDblFinder, Seurat, SeuratObject
#' `merge.Assay5`/`JoinLayers.Assay5`, SoupX, tradeSeq, and WGCNA. It
#' normally requires the pinned AutoZyme build, or an exact patch fingerprint
#' recorded by Shennong, and an exactly validated upstream version. Automatic
#' Seurat scopes relax the version-label check and retain runtime structure
#' guards; the scDblFinder patch instead requires exact target-function
#' fingerprints. Shennong bundles the
#' scDblFinder, SeuratObject merge/JoinLayers, and SoupX patches (plus
#' SoupX's compiled kernels), validates
#' each bundled fingerprint, and registers them through AutoZyme when the
#' official package does not provide `scdblfinder`, `seurat_joinlayers`,
#' `seurat_merge`, or `soupx`. The returned
#' `patch_provider` and `bundled_by_shennong` columns make this source explicit.
#' Set
#' `options(shennong.autozyme = FALSE)`
#' or the environment variable `AUTOZYME_DISABLED=true` (the legacy alias
#' `AUTOZYME_DISABLE=true` is also accepted) to prevent automatic activation.
#' Newly activated patches are restored even when the workflow errors.
#' Automatic workflow scopes log the patch names after activation succeeds and
#' before the analytical expression starts. This confirms that the patches are
#' enabled for the call; individual AutoZyme fast paths may still fall back to
#' upstream code when an input is outside their validated scope. These
#' controls do not disable patches that are already active and do not affect
#' explicit calls to [sn_enable_autozyme()] or [sn_with_autozyme()].
#' The older broad Seurat acceleration is bypassed for BPCells-backed objects
#' so an on-disk layer is never coerced to an in-memory sparse matrix. The
#' narrow `seurat_joinlayers` patch remains eligible because its validated
#' fast path supports public BPCells `IterableMatrix` layers.
#' @export
#'
#' @examples
#' sn_check_autozyme()
sn_check_autozyme <- function(
    patches = c(
      "cellchat", "clusterprofiler", "fgsea", "nichenetr", "scdblfinder",
      "seurat", "seurat_joinlayers", "seurat_merge", "soupx", "tradeseq",
      "wgcna"
    ),
    strict = TRUE,
    allow_approximate = FALSE) {
  patches <- .sn_validate_autozyme_patches(patches)
  strict <- .sn_validate_autozyme_flag(strict, "strict")
  allow_approximate <- .sn_validate_autozyme_flag(
    allow_approximate,
    "allow_approximate"
  )

  autozyme_installed <- .sn_autozyme_is_installed("autozyme")
  description <- .sn_autozyme_description()
  source_match <- if (is.na(description$remote_sha)) {
    NA
  } else {
    identical(description$remote_sha, .sn_autozyme_expected_sha)
  }
  build_version_match <- !is.na(description$version) &&
    identical(description$version, .sn_autozyme_expected_version)
  active <- if (autozyme_installed) .sn_autozyme_active_patches() else character(0)

  rows <- lapply(patches, function(patch) {
    spec <- .sn_autozyme_patch_manifest[[patch]]
    patch_source <- .sn_autozyme_patch_source_status(patch, spec)
    accepted_source <- isTRUE(source_match) || isTRUE(patch_source$source_match)
    build_match <- build_version_match && accepted_source
    installed_version <- .sn_autozyme_installed_version(spec$upstream)
    upstream_installed <- !is.na(installed_version)
    version_match <- upstream_installed && installed_version %in% spec$versions
    approximation_allowed <- !isTRUE(spec$approximate) || allow_approximate
    patch_available <- isTRUE(patch_source$registered) ||
      isTRUE(patch_source$vendored_source_match)
    eligible <- autozyme_installed && patch_available && upstream_installed &&
      approximation_allowed && (!strict || (build_match && version_match))

    reason <- if (!autozyme_installed) {
      "autozyme is not installed"
    } else if (!patch_available) {
      "patch is unavailable from both AutoZyme and Shennong"
    } else if (strict && !build_version_match) {
      sprintf(
        "installed AutoZyme version %s does not match pinned version %s",
        description$version,
        .sn_autozyme_expected_version
      )
    } else if (strict && !accepted_source) {
      paste0(
        "AutoZyme source does not match the pinned revision and no trusted ",
        "Shennong-vendored patch is available"
      )
    } else if (!upstream_installed) {
      sprintf("upstream package %s is not installed", spec$upstream)
    } else if (!approximation_allowed) {
      "approximate patch requires allow_approximate = TRUE"
    } else if (strict && !version_match) {
      sprintf(
        "installed upstream version %s is not exactly validated",
        installed_version
      )
    } else if (!strict && !version_match) {
      "eligible with explicitly allowed upstream version drift"
    } else {
      "eligible"
    }

    data.frame(
      patch = patch,
      upstream = spec$upstream,
      installed_version = installed_version,
      tested_versions = paste(spec$versions, collapse = ","),
      equivalence = spec$equivalence,
      approximate = isTRUE(spec$approximate),
      default = patch %in% .sn_autozyme_default_patches,
      autozyme_installed = autozyme_installed,
      autozyme_version = description$version,
      autozyme_version_match = build_version_match,
      autozyme_remote_sha = description$remote_sha,
      autozyme_source_match = source_match,
      autozyme_patch_source_match = patch_source$source_match,
      registered = patch_source$registered,
      bundled_by_shennong = patch_source$bundled_by_shennong,
      patch_provider = patch_source$provider,
      upstream_installed = upstream_installed,
      version_match = version_match,
      active = patch %in% active,
      eligible = eligible,
      reason = reason,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}

.sn_require_autozyme <- function() {
  if (.sn_autozyme_is_installed("autozyme")) {
    return(invisible(TRUE))
  }

  stop(
    paste0(
      "The optional package `autozyme` is not installed. Install the pinned ",
      "revision declared in Shennong's DESCRIPTION before enabling it; ",
      "Shennong will not install AutoZyme automatically."
    ),
    call. = FALSE
  )
}

.sn_autozyme_rollback <- function(patches) {
  future_option <- .sn_capture_autozyme_future_option()
  on.exit(.sn_restore_autozyme_future_option(future_option), add = TRUE)

  patches <- unique(as.character(patches))
  if (length(patches) == 0L) {
    return(list(
      complete = TRUE,
      residual_active = character(),
      deactivation_errors = character(),
      verification_error = NULL
    ))
  }

  deactivation_errors <- character()
  before_error <- NULL
  active_before <- tryCatch(
    .sn_autozyme_active_patches(),
    error = function(error) {
      before_error <<- conditionMessage(error)
      character()
    }
  )
  rollback_targets <- if (is_null(before_error)) {
    intersect(patches, active_before)
  } else {
    patches
  }

  for (patch in rev(rollback_targets)) {
    tryCatch(
      .sn_autozyme_call("deactivate", patch),
      error = function(error) {
        deactivation_errors <<- c(
          deactivation_errors,
          sprintf("%s (%s)", patch, conditionMessage(error))
        )
      }
    )
  }

  verification_error <- NULL
  active_after <- tryCatch(
    .sn_autozyme_active_patches(),
    error = function(error) {
      verification_error <<- conditionMessage(error)
      character()
    }
  )
  residual_active <- if (is_null(verification_error)) {
    intersect(patches, active_after)
  } else {
    character()
  }
  if (!is_null(before_error)) {
    deactivation_errors <- c(
      deactivation_errors,
      paste0("pre-rollback status check failed (", before_error, ")")
    )
  }

  list(
    complete = is_null(verification_error) && length(residual_active) == 0L,
    residual_active = residual_active,
    deactivation_errors = deactivation_errors,
    verification_error = verification_error
  )
}

#' Enable optional AutoZyme patches
#'
#' Manual activation is transactional. All selected patches are checked before
#' activation; if activation fails, only patches newly enabled by this call are
#' rolled back. Existing user-activated patches are left untouched.
#'
#' @inheritParams sn_check_autozyme
#'
#' @return Invisibly, the post-activation compatibility data frame returned by
#'   [sn_check_autozyme()].
#' @export
#'
#' @examples
#' \dontrun{
#' sn_enable_autozyme()
#' }
sn_enable_autozyme <- function(
    patches = c("cellchat", "nichenetr"),
    strict = TRUE,
    allow_approximate = FALSE) {
  patches <- .sn_validate_autozyme_patches(patches)
  strict <- .sn_validate_autozyme_flag(strict, "strict")
  allow_approximate <- .sn_validate_autozyme_flag(
    allow_approximate,
    "allow_approximate"
  )
  .sn_require_autozyme()

  checks <- sn_check_autozyme(
    patches = patches,
    strict = strict,
    allow_approximate = allow_approximate
  )
  blocked <- checks[!checks$eligible, , drop = FALSE]
  if (nrow(blocked) > 0L) {
    details <- paste0(blocked$patch, " (", blocked$reason, ")")
    stop(
      sprintf("Cannot enable AutoZyme: %s.", paste(details, collapse = "; ")),
      call. = FALSE
    )
  }

  vendored <- checks$patch[
    checks$patch_provider == "shennong" & !checks$active
  ]
  for (patch in stats::na.omit(vendored)) {
    .sn_register_vendored_autozyme_patch(patch)
  }

  available <- .sn_autozyme_call("list_patches")
  missing_patches <- setdiff(patches, available)
  if (length(missing_patches) > 0L) {
    stop(
      sprintf(
        "The installed AutoZyme build does not provide patch(es): %s.",
        paste(missing_patches, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  approximate <- checks$patch[checks$approximate]
  if (length(approximate) > 0L) {
    warning(
      sprintf(
        paste0(
          "Enabling approximate AutoZyme patch(es): %s. ",
          "These fast paths can change analytical results."
        ),
        paste(approximate, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  before <- .sn_autozyme_active_patches()
  requested_new <- setdiff(patches, before)
  post_checks <- NULL
  tryCatch(
    {
      for (patch in requested_new) {
        ok <- .sn_autozyme_call("activate", patch)
        if (!isTRUE(ok)) {
          stop(sprintf("AutoZyme did not activate patch `%s`.", patch), call. = FALSE)
        }
      }

      after <- .sn_autozyme_active_patches()
      missing_active <- setdiff(patches, after)
      if (length(missing_active) > 0L) {
        stop(
          sprintf(
            "AutoZyme did not report patch(es) as active: %s.",
            paste(missing_active, collapse = ", ")
          ),
          call. = FALSE
        )
      }

      post_checks <- sn_check_autozyme(
        patches = patches,
        strict = strict,
        allow_approximate = allow_approximate
      )
    },
    error = function(e) {
      rollback <- .sn_autozyme_rollback(requested_new)
      rollback_details <- character()
      if (length(rollback$deactivation_errors) > 0L) {
        rollback_details <- c(
          rollback_details,
          paste0(
            "deactivation error(s): ",
            paste(rollback$deactivation_errors, collapse = "; ")
          )
        )
      }
      if (!is_null(rollback$verification_error)) {
        rollback_details <- c(
          rollback_details,
          paste0(
            "post-rollback state could not be verified: ",
            rollback$verification_error
          )
        )
      }
      if (length(rollback$residual_active) > 0L) {
        rollback_details <- c(
          rollback_details,
          paste0(
            "residual active patch(es): ",
            paste(rollback$residual_active, collapse = ", ")
          )
        )
      }
      message <- conditionMessage(e)
      if (!isTRUE(rollback$complete) || length(rollback_details) > 0L) {
        message <- paste0(
          message,
          " AutoZyme rollback diagnostics: ",
          paste(rollback_details, collapse = "; "),
          "."
        )
      }
      stop(message, call. = FALSE)
    }
  )

  invisible(post_checks)
}

#' Disable optional AutoZyme patches
#'
#' This function is idempotent and never calls AutoZyme's global
#' `deactivate_all()` helper. Only the explicitly selected patches are changed.
#'
#' @param patches Character vector of supported patch names to disable.
#'
#' @return Invisibly, `TRUE` when at least one selected patch was active and
#'   disabled, otherwise `FALSE`.
#' @export
#'
#' @examples
#' \dontrun{
#' sn_disable_autozyme()
#' }
sn_disable_autozyme <- function(patches = c(
    "cellchat", "clusterprofiler", "fgsea", "nichenetr", "scdblfinder",
    "seurat", "seurat_joinlayers", "seurat_merge", "soupx", "tradeseq",
    "wgcna"
)) {
  patches <- .sn_validate_autozyme_patches(patches)
  if (!.sn_autozyme_is_installed("autozyme")) {
    return(invisible(FALSE))
  }

  active <- .sn_autozyme_active_patches()
  selected <- intersect(patches, active)
  for (patch in selected) {
    .sn_autozyme_call("deactivate", patch)
  }

  invisible(length(selected) > 0L)
}

#' Evaluate code with temporary AutoZyme acceleration
#'
#' Only patches newly activated by this call are disabled on exit, including
#' when `expr` raises an error. Patches that were already active remain active.
#'
#' @param expr Code to evaluate after activation.
#' @inheritParams sn_check_autozyme
#'
#' @return The value of `expr`.
#' @export
#'
#' @examples
#' \dontrun{
#' result <- sn_with_autozyme({
#'   sn_run_cell_communication(
#'     object,
#'     method = "cellchat",
#'     group_by = "cell_type"
#'   )
#' })
#' }
sn_with_autozyme <- function(
    expr,
    patches = c("cellchat", "nichenetr"),
    strict = TRUE,
    allow_approximate = FALSE) {
  patches <- .sn_validate_autozyme_patches(patches)
  before <- .sn_autozyme_active_patches()
  sn_enable_autozyme(
    patches = patches,
    strict = strict,
    allow_approximate = allow_approximate
  )
  newly_activated <- setdiff(patches, before)

  if (length(newly_activated) > 0L) {
    on.exit(
      tryCatch(
        sn_disable_autozyme(newly_activated),
        error = function(e) {
          warning(
            sprintf("Could not restore AutoZyme state: %s", conditionMessage(e)),
            call. = FALSE
          )
        }
      ),
      add = TRUE
    )
  }

  force(expr)
}

.sn_autozyme_provenance <- function() {
  description <- .sn_autozyme_description()
  installed <- .sn_autozyme_is_installed("autozyme")
  active <- if (installed) .sn_autozyme_active_patches() else character(0)
  disabled <- installed && .sn_autozyme_namespace_loaded() && isTRUE(tryCatch(
    .sn_autozyme_call("is_disabled"),
    error = function(error) FALSE
  ))
  if (disabled) {
    active <- character(0)
  }
  context <- getOption("shennong.autozyme.provenance_context")
  if (is.environment(context)) {
    relevant <- context$patches %||% character()
    used <- intersect(context$used_patches %||% character(), relevant)
    suppressed <- intersect(
      context$suppressed_patches %||% character(),
      relevant
    )
    active <- union(
      setdiff(intersect(active, relevant), suppressed),
      used
    )
  }
  if (!installed || length(active) == 0L) {
    return(list())
  }

  patch_details <- lapply(active, function(patch) {
    spec <- .sn_autozyme_patch_manifest[[patch]]
    if (is.null(spec)) {
      return(list(patch = patch))
    }
    installed_version <- .sn_autozyme_installed_version(spec$upstream)
    patch_source <- .sn_autozyme_patch_source_status(patch, spec)
    list(
      patch = patch,
      upstream = spec$upstream,
      installed_version = installed_version,
      tested_versions = spec$versions,
      version_match = !is.na(installed_version) && installed_version %in% spec$versions,
      equivalence = spec$equivalence,
      approximate = isTRUE(spec$approximate),
      provider = patch_source$provider,
      bundled_by_shennong = isTRUE(patch_source$bundled_by_shennong)
    )
  })
  names(patch_details) <- active

  list(
    available = installed,
    version = description$version,
    remote_sha = description$remote_sha,
    expected_version = .sn_autozyme_expected_version,
    expected_remote_sha = .sn_autozyme_expected_sha,
    active_patches = active,
    patches = patch_details
  )
}
