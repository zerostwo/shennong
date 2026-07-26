# Optional AutoZyme acceleration ------------------------------------------------

.sn_autozyme_expected_version <- "0.3.1"
.sn_autozyme_expected_sha <- "718541d9489596c7c1d75f52e9b3a8b2a429d1f9"
.sn_autozyme_default_patches <- c("cellchat", "nichenetr")

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
    upstream = "clusterProfiler", versions = "4.16.0",
    equivalence = "exact_scoped", approximate = FALSE
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
  seurat = list(
    upstream = "Seurat", versions = c("5.2.1", "5.4.0"),
    equivalence = "numeric_tolerance_scoped", approximate = FALSE
  ),
  slingshot = list(
    upstream = "slingshot", versions = "2.16.0",
    equivalence = "approximate", approximate = TRUE
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

.sn_autozyme_call <- function(fun, ...) {
  do.call(base::getExportedValue("autozyme", fun), list(...))
}

.sn_autozyme_status_vector <- function() {
  if (!.sn_autozyme_is_installed("autozyme") ||
      !base::isNamespaceLoaded("autozyme")) {
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

#' Check optional AutoZyme acceleration compatibility
#'
#' Reports whether selected AutoZyme patches can be enabled without changing
#' Shennong's default execution. This check is side-effect free: it does not
#' install packages, activate patches, prepare Python, or change thread counts.
#'
#' @param patches Character vector of AutoZyme patch names. The conservative
#'   defaults cover Shennong's CellChat and NicheNet communication paths.
#' @param strict If `TRUE`, require the installed upstream package version to
#'   exactly match a version validated by the pinned AutoZyme revision.
#' @param allow_approximate If `TRUE`, permit patches whose documented fast
#'   path can change analytical results beyond floating-point tolerance.
#'
#' @return A data frame with installation, version, activation, and eligibility
#'   information for each requested patch.
#' @export
#'
#' @examples
#' sn_check_autozyme()
sn_check_autozyme <- function(
    patches = c("cellchat", "nichenetr"),
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
  build_match <- build_version_match && isTRUE(source_match)
  active <- if (autozyme_installed) .sn_autozyme_active_patches() else character(0)

  rows <- lapply(patches, function(patch) {
    spec <- .sn_autozyme_patch_manifest[[patch]]
    installed_version <- .sn_autozyme_installed_version(spec$upstream)
    upstream_installed <- !is.na(installed_version)
    version_match <- upstream_installed && installed_version %in% spec$versions
    approximation_allowed <- !isTRUE(spec$approximate) || allow_approximate
    eligible <- autozyme_installed && upstream_installed &&
      approximation_allowed && (!strict || (build_match && version_match))

    reason <- if (!autozyme_installed) {
      "autozyme is not installed"
    } else if (strict && !build_version_match) {
      sprintf(
        "installed AutoZyme version %s does not match pinned version %s",
        description$version,
        .sn_autozyme_expected_version
      )
    } else if (strict && !isTRUE(source_match)) {
      "installed AutoZyme source does not match the pinned revision"
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
#' Activation is explicit and transactional. All selected patches are checked
#' before activation; if activation fails, only patches newly enabled by this
#' call are rolled back. Existing user-activated patches are left untouched.
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
sn_disable_autozyme <- function(patches = c("cellchat", "nichenetr")) {
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
#'   sn_run_cell_communication(object, method = "cellchat")
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
  if (!installed || length(active) == 0L) {
    return(list())
  }

  patch_details <- lapply(active, function(patch) {
    spec <- .sn_autozyme_patch_manifest[[patch]]
    if (is.null(spec)) {
      return(list(patch = patch))
    }
    installed_version <- .sn_autozyme_installed_version(spec$upstream)
    list(
      patch = patch,
      upstream = spec$upstream,
      installed_version = installed_version,
      tested_versions = spec$versions,
      version_match = !is.na(installed_version) && installed_version %in% spec$versions,
      equivalence = spec$equivalence,
      approximate = isTRUE(spec$approximate)
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
