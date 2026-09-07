# Pixi runtime.
#
# Extracted from package_tools.R: pixi binary install/check helpers, bundled environment discovery/selection/rendering, command execution, mirror configuration, and deprecated pixi aliases.
#' Check, install, and configure pixi for Shennong Python backends
#'
#' Shennong uses pixi for optional Python backends such as scVI/scANVI. The R
#' package ships runner scripts, but the Python environments themselves are
#' created under the user-level Shennong runtime directory
#' \code{~/.shennong/}, not inside the current analysis project and not inside
#' the installed R package.
#'
#' @param pixi Optional pixi executable path. When \code{NULL}, Shennong checks
#'   \code{options("shennong.pixi")}, \code{PATH}, and \code{~/.pixi/bin/pixi}.
#' @param quiet Logical; suppress status messages where possible.
#' @param install Logical used by \code{sn_ensure_pixi()}; install pixi when the
#'   requested exact version is not already available.
#' @param version Explicit Pixi release version downloaded from the official
#'   GitHub release. Defaults to the package-tested version \code{"0.69.0"};
#'   mutable \code{"latest"} downloads are rejected.
#' @param pixi_home Pixi home directory. The installer defaults to
#'   \code{"~/.pixi"}; Shennong runtime workflows set \code{PIXI_HOME} to a
#'   user-level path such as \code{"~/.shennong/pixi/home"}.
#' @param bin_dir Optional directory where the standalone pixi binary should be
#'   installed.
#' @param no_path_update Logical; when \code{TRUE}, ask the installer not to
#'   edit shell startup files.
#' @param download_url Optional custom pixi archive or binary URL. This is
#'   useful for institutional mirrors and requires an explicit \code{sha256}.
#' @param sha256 Optional expected SHA-256 digest. Official versioned downloads
#'   obtain the release sidecar digest automatically; custom URLs must supply
#'   this value.
#' @param force Logical; reinstall even if pixi is already available.
#'
#' @return \code{sn_check_pixi()} and \code{sn_ensure_pixi()} return a named
#'   list with install status, executable path, and version. \code{sn_install_pixi()}
#'   invisibly returns the refreshed check result.
#'
#' @examples
#' info <- sn_check_pixi(quiet = TRUE)
#' info$installed
#' \dontrun{
#' sn_ensure_pixi()
#' sn_install_pixi(pixi_home = "~/.pixi", no_path_update = TRUE)
#' }
#'
#' @export
sn_check_pixi <- function(pixi = NULL, quiet = FALSE) {
  pixi_path <- .sn_resolve_pixi_path(pixi = pixi)
  installed <- !is.null(pixi_path) && nzchar(pixi_path) && file.exists(pixi_path)
  version <- NULL
  output <- character(0)

  if (installed) {
    output <- tryCatch(
      suppressWarnings(.sn_system2_capture(pixi_path, "--version")),
      error = function(e) character(0)
    )
    status <- attr(output, "status") %||% 0L
    if (identical(status, 0L) && length(output) > 0L) {
      version <- sub("^pixi\\s+", "", output[[1]])
    }
  }

  result <- list(
    installed = installed,
    path = if (installed) pixi_path else NA_character_,
    version = version,
    output = output
  )

  if (!quiet) {
    if (installed) {
      .sn_log_info("pixi detected at {pixi_path}; version = {version %||% 'unknown'}.")
    } else {
      .sn_log_warn("pixi was not found. Use `sn_install_pixi()` or `sn_ensure_pixi()` to install it.")
    }
  }

  invisible(result)
}

#' @rdname sn_check_pixi
#' @export
sn_install_pixi <- function(version = "0.69.0",
                            pixi_home = "~/.pixi",
                            bin_dir = NULL,
                            no_path_update = TRUE,
                            download_url = NULL,
                            sha256 = NULL,
                            force = FALSE,
                            quiet = FALSE) {
  version <- .sn_normalize_requested_pixi_version(version)
  current <- sn_check_pixi(quiet = TRUE)
  current_matches <- isTRUE(current$installed) &&
    .sn_pixi_version_matches(current$version, version)
  if (current_matches && !isTRUE(force)) {
    if (!quiet) {
      .sn_log_info("pixi {version} is already installed at {current$path}.")
    }
    return(invisible(current))
  }
  if (isTRUE(current$installed) && !current_matches && !quiet) {
    .sn_log_warn(
      "Installed pixi version {current$version %||% 'unknown'} does not match required version {version}; installing the required release."
    )
  }

  # Install a versioned, checksummed release payload instead of executing the
  # mutable remote shell/PowerShell installer. `no_path_update` remains in the
  # public API for compatibility; direct installation never edits startup files.
  invisible(no_path_update)
  downloaded <- .sn_download_verified_pixi(
    version = version,
    download_url = download_url,
    sha256 = sha256,
    quiet = quiet
  )
  on.exit(unlink(downloaded$path, force = TRUE), add = TRUE)

  install_dir <- if (!is.null(bin_dir) && nzchar(bin_dir)) {
    path.expand(bin_dir)
  } else {
    file.path(path.expand(pixi_home), "bin")
  }
  installed_binary <- .sn_install_verified_pixi_archive(
    archive = downloaded$path,
    asset = downloaded$asset,
    install_dir = install_dir
  )
  refreshed <- sn_check_pixi(pixi = installed_binary, quiet = quiet)
  if (!isTRUE(refreshed$installed)) {
    refreshed <- sn_check_pixi(quiet = quiet)
  }
  if (!isTRUE(refreshed$installed)) {
    stop("Verified Pixi installation completed, but the executable could not be found.", call. = FALSE)
  }
  if (!.sn_pixi_version_matches(refreshed$version, version)) {
    stop(
      "Installed Pixi executable reports version '", refreshed$version %||% "unknown",
      "', but version '", version, "' was required.",
      call. = FALSE
    )
  }

  invisible(refreshed)
}

#' @rdname sn_check_pixi
#' @export
sn_ensure_pixi <- function(pixi = NULL,
                           install = TRUE,
                           version = "0.69.0",
                           pixi_home = "~/.pixi",
                           bin_dir = NULL,
                           no_path_update = TRUE,
                           download_url = NULL,
                           sha256 = NULL,
                           quiet = FALSE) {
  version <- .sn_normalize_requested_pixi_version(version)
  info <- sn_check_pixi(pixi = pixi, quiet = TRUE)
  version_matches <- isTRUE(info$installed) &&
    .sn_pixi_version_matches(info$version, version)
  if (version_matches) {
    if (!quiet) {
      .sn_log_info("pixi {version} detected at {info$path}.")
    }
    return(invisible(info))
  }

  if (!isTRUE(install)) {
    if (isTRUE(info$installed)) {
      stop(
        "Installed Pixi version '", info$version %||% "unknown",
        "' does not match required version '", version,
        "'. Set `install = TRUE` to install the required release.",
        call. = FALSE
      )
    }
    stop(
      "pixi is not installed. Install it with `sn_install_pixi()` or set ",
      "`install = TRUE` in `sn_ensure_pixi()`.",
      call. = FALSE
    )
  }

  sn_install_pixi(
    version = version,
    pixi_home = pixi_home,
    bin_dir = bin_dir,
    no_path_update = no_path_update,
    download_url = download_url,
    sha256 = sha256,
    force = isTRUE(info$installed),
    quiet = quiet
  )
}

#' Inspect Shennong pixi runtime paths
#'
#' Returns the user-level paths used for optional pixi-managed Python
#' environments. Shennong follows the same convention as downloaded example
#' data: runtime files are generated under \code{~/.shennong/} by default, not
#' under the current analysis project and not inside the installed R package.
#'
#' @param environment Python environment name. Use
#'   \code{sn_list_pixi_environments()} to see bundled configs.
#' @param runtime_dir Optional explicit Shennong runtime directory. Defaults to
#'   \code{getOption("shennong.runtime_dir")}, \code{SHENNONG_RUNTIME_DIR},
#'   \code{SHENNONG_HOME}, then \code{"~/.shennong"}.
#'
#' @return A named list of runtime paths.
#'
#' @examples
#' sn_get_pixi_paths("scvi", runtime_dir = tempfile("shennong-home-"))
#'
#' @export
sn_get_pixi_paths <- function(environment = NULL,
                          runtime_dir = NULL) {
  requested_environment <- tolower(as.character(environment %||% "scvi"))
  environment <- .sn_normalize_pixi_environment(requested_environment)
  runtime_dir <- .sn_shennong_runtime_dir(runtime_dir)
  pixi_root <- file.path(runtime_dir, "pixi")
  project_dir <- file.path(pixi_root, environment)

  list(
    environment = requested_environment,
    family = environment,
    runtime_dir = normalizePath(runtime_dir, winslash = "/", mustWork = FALSE),
    pixi_root = normalizePath(pixi_root, winslash = "/", mustWork = FALSE),
    pixi_home = normalizePath(file.path(pixi_root, "home"), winslash = "/", mustWork = FALSE),
    project_dir = normalizePath(project_dir, winslash = "/", mustWork = FALSE),
    source_config_path = sn_get_pixi_config_path(environment),
    manifest_path = normalizePath(file.path(project_dir, "pixi.toml"), winslash = "/", mustWork = FALSE),
    workspace_env_dir = normalizePath(file.path(project_dir, ".pixi", "envs"), winslash = "/", mustWork = FALSE),
    runs_dir = normalizePath(file.path(runtime_dir, "runs"), winslash = "/", mustWork = FALSE)
  )
}

#' List bundled pixi environment configs
#'
#' @return A character vector of environment names with pixi manifests bundled
#'   under \code{inst/pixi/}.
#'
#' @examples
#' sn_list_pixi_environments()
#'
#' @export
sn_list_pixi_environments <- function() {
  sort(.sn_pixi_environment_names())
}

#' Locate a bundled pixi config
#'
#' @param environment Python environment name.
#'
#' @return A path to the package-bundled \code{pixi.toml} template.
#'
#' @examples
#' sn_get_pixi_config_path("scvi")
#'
#' @export
sn_get_pixi_config_path <- function(environment = NULL) {
  environment <- .sn_normalize_pixi_environment(environment %||% "scvi")
  installed <- system.file("pixi", environment, "pixi.toml", package = "Shennong")
  if (nzchar(installed) && file.exists(installed)) {
    return(normalizePath(installed, winslash = "/", mustWork = TRUE))
  }
  source_path <- file.path(getwd(), "inst", "pixi", environment, "pixi.toml")
  if (file.exists(source_path)) {
    return(normalizePath(source_path, winslash = "/", mustWork = TRUE))
  }
  stop("Could not locate the bundled pixi config for environment: ", environment, call. = FALSE)
}

#' Prepare or call a Shennong pixi environment
#'
#' \code{sn_prepare_pixi_environment()} materializes a package-bundled
#' \code{inst/pixi/<family>/pixi.toml} template into the user-level
#' \code{~/.shennong/pixi/<family>/} workspace. \code{sn_call_pixi_environment()}
#' runs a command inside one of these environments.
#'
#' @details
#' The environment-specific aliases \code{sn_call_scvi()},
#' \code{sn_call_scanvi()}, \code{sn_call_mmochi()}, \code{sn_call_scarches()},
#' \code{sn_call_scpoli()}, \code{sn_call_infercnvpy()},
#' \code{sn_call_trajectory()}, \code{sn_call_cellphonedb()},
#' \code{sn_call_cell2location()}, \code{sn_call_tangram()},
#' \code{sn_call_squidpy()}, \code{sn_call_spatialdata()}, and
#' \code{sn_call_stlearn()} are deprecated compatibility wrappers that only
#' forward to \code{sn_call_pixi_environment()}. Call
#' \code{sn_call_pixi_environment("<environment>", command = ..., args = ...)}
#' directly; the aliases emit a deprecation warning and will be removed in a
#' future major release.
#'
#' @param environment Python environment name.
#' @param pixi_environment Pixi environment inside the manifest, for example
#'   \code{"cpu"}, \code{"gpu"}, or \code{"default"}. \code{"auto"} uses
#'   CUDA when available for GPU-aware configs and otherwise CPU/default.
#' @param runtime_dir Optional Shennong runtime directory.
#' @param project_dir Optional explicit pixi workspace directory.
#' @param manifest_path Optional explicit materialized manifest path.
#' @param overwrite Whether to overwrite an existing materialized manifest.
#' @param cuda_version CUDA runtime version used when rendering templates.
#' @param platforms Platforms that the bundled lock must cover. This does not
#'   narrow the platform inventory declared by the bundled manifest, because
#'   doing so would make its multi-platform lock stale. Defaults to the current
#'   platform.
#' @param mirror Mirror setting passed to \code{sn_configure_pixi_mirror()}.
#' @param install_pixi Ensure the standalone pixi binary is available.
#' @param install_environment Run \code{pixi install} for the selected
#'   environment after materializing the manifest.
#' @param pixi Optional pixi executable path.
#' @param pixi_version Pixi version used if installation is needed.
#' @param pixi_download_url Optional custom pixi binary download URL.
#' @param pixi_sha256 Optional expected SHA-256 digest for
#'   \code{pixi_download_url}. Custom URLs require this value; official pinned
#'   release downloads obtain their checksum sidecar automatically.
#' @param command Command to run inside the pixi environment.
#' @param args Character vector of command arguments.
#' @param quiet Logical; suppress status messages where possible.
#' @param ... Additional arguments passed from environment-specific helpers to
#'   \code{sn_call_pixi_environment()}.
#'
#' @return \code{sn_prepare_pixi_environment()} returns a named list of paths
#'   (including \code{manifest_path} and \code{lock_path}) and selected
#'   environment metadata. \code{sn_call_pixi_environment()} invisibly returns
#'   command output.
#'
#' @examples
#' sn_prepare_pixi_environment("scvi", runtime_dir = tempfile("shennong-home-"))
#' \dontrun{
#' sn_call_pixi_environment("scvi", command = "python", args = "--version")
#' }
#'
#' @export
sn_prepare_pixi_environment <- function(environment = NULL,
                                        pixi_environment = c("auto", "default", "cpu", "gpu"),
                                        runtime_dir = NULL,
                                        project_dir = NULL,
                                        manifest_path = NULL,
                                        overwrite = FALSE,
                                        cuda_version = NULL,
                                        platforms = NULL,
                                        mirror = c("default", "auto", "china", "tuna", "ustc", "bfsu"),
                                        install_pixi = FALSE,
                                        install_environment = FALSE,
                                        pixi = NULL,
                                        pixi_version = "0.69.0",
                                        pixi_download_url = NULL,
                                        pixi_sha256 = NULL,
                                        quiet = FALSE) {
  requested_environment <- tolower(as.character(environment %||% "scvi"))
  environment <- .sn_normalize_pixi_environment(requested_environment)
  pixi_environment <- match.arg(pixi_environment)
  mirror <- match.arg(mirror)
  selected <- .sn_select_pixi_environment(environment = environment, pixi_environment = pixi_environment, cuda_version = cuda_version)
  paths <- sn_get_pixi_paths(environment = environment, runtime_dir = runtime_dir)
  if (!is.null(project_dir) && nzchar(project_dir)) {
    paths$project_dir <- normalizePath(path.expand(project_dir), winslash = "/", mustWork = FALSE)
    paths$manifest_path <- normalizePath(file.path(paths$project_dir, "pixi.toml"), winslash = "/", mustWork = FALSE)
    paths$workspace_env_dir <- normalizePath(file.path(paths$project_dir, ".pixi", "envs"), winslash = "/", mustWork = FALSE)
  }
  if (!is.null(manifest_path) && nzchar(manifest_path)) {
    paths$manifest_path <- normalizePath(path.expand(manifest_path), winslash = "/", mustWork = FALSE)
  }

  dir.create(dirname(paths$manifest_path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(paths$manifest_path) || isTRUE(overwrite)) {
    rendered <- .sn_render_pixi_config(
      environment = environment,
      platforms = platforms,
      cuda_version = selected$cuda_version
    )
    writeLines(rendered, con = paths$manifest_path, useBytes = TRUE)
  }
  paths$manifest_path <- normalizePath(paths$manifest_path, winslash = "/", mustWork = TRUE)
  requested_platforms <- unique(as.character(platforms %||% .sn_current_pixi_platform()))
  paths$lock_path <- .sn_prepare_pixi_lock(
    environment = environment,
    manifest_path = paths$manifest_path,
    platforms = requested_platforms,
    cuda_version = selected$cuda_version,
    overwrite = overwrite
  )

  resolved_mirror <- .sn_resolve_pixi_mirror(mirror)
  if (!identical(resolved_mirror, "default")) {
    sn_configure_pixi_mirror(mirror = mirror, pixi_home = paths$pixi_home, append_original = TRUE)
  }

  pixi_info <- NULL
  if (isTRUE(install_pixi) || isTRUE(install_environment)) {
    pixi_info <- sn_ensure_pixi(
      pixi = pixi,
      install = TRUE,
      version = pixi_version,
      pixi_home = "~/.pixi",
      no_path_update = TRUE,
      download_url = pixi_download_url,
      sha256 = pixi_sha256,
      quiet = quiet
    )
  }
  if (isTRUE(install_environment)) {
    .sn_pixi_install_environment(
      pixi = pixi_info$path,
      manifest_path = paths$manifest_path,
      pixi_environment = selected$pixi_environment,
      pixi_home = paths$pixi_home
    )
  }

  invisible(c(
    paths,
    list(
      requested_environment = requested_environment,
      pixi_environment = selected$pixi_environment,
      accelerator = selected$accelerator,
      cuda_version = selected$cuda_version
    )
  ))
}

#' @rdname sn_prepare_pixi_environment
#' @export
sn_call_pixi_environment <- function(environment = NULL,
                                     command,
                                     args = character(),
                                     pixi_environment = c("auto", "default", "cpu", "gpu"),
                                     runtime_dir = NULL,
                                     project_dir = NULL,
                                     manifest_path = NULL,
                                     overwrite = FALSE,
                                     cuda_version = NULL,
                                     platforms = NULL,
                                     mirror = c("default", "auto", "china", "tuna", "ustc", "bfsu"),
                                     install_pixi = TRUE,
                                     pixi = NULL,
                                     pixi_version = "0.69.0",
                                     pixi_download_url = NULL,
                                     pixi_sha256 = NULL,
                                     quiet = FALSE) {
  if (missing(command) || is.null(command) || !nzchar(command)) {
    stop("`command` must be a non-empty command name.", call. = FALSE)
  }
  prepared <- sn_prepare_pixi_environment(
    environment = environment,
    pixi_environment = pixi_environment,
    runtime_dir = runtime_dir,
    project_dir = project_dir,
    manifest_path = manifest_path,
    overwrite = overwrite,
    cuda_version = cuda_version,
    platforms = platforms,
    mirror = mirror,
    install_pixi = install_pixi,
    install_environment = FALSE,
    pixi = pixi,
    pixi_version = pixi_version,
    pixi_download_url = pixi_download_url,
    pixi_sha256 = pixi_sha256,
    quiet = quiet
  )
  pixi_info <- sn_ensure_pixi(
    pixi = pixi,
    install = install_pixi,
    version = pixi_version,
    pixi_home = "~/.pixi",
    no_path_update = TRUE,
    download_url = pixi_download_url,
    sha256 = pixi_sha256,
    quiet = quiet
  )
  out <- .sn_pixi_run_command(
    pixi = pixi_info$path,
    manifest_path = prepared$manifest_path,
    pixi_environment = prepared$pixi_environment,
    command = command,
    args = args,
    pixi_home = prepared$pixi_home
  )
  invisible(out)
}

#' @rdname sn_prepare_pixi_environment
#' @export
sn_call_scvi <- function(command, args = character(), ...) {
  .Deprecated("sn_call_pixi_environment", package = "Shennong", old = "sn_call_scvi")
  sn_call_pixi_environment("scvi", command = command, args = args, ...)
}

#' @rdname sn_prepare_pixi_environment
#' @export
sn_call_scanvi <- function(command, args = character(), ...) {
  .Deprecated("sn_call_pixi_environment", package = "Shennong", old = "sn_call_scanvi")
  sn_call_pixi_environment("scanvi", command = command, args = args, ...)
}

#' @rdname sn_prepare_pixi_environment
#' @export
sn_call_mmochi <- function(command, args = character(), ...) {
  .Deprecated("sn_call_pixi_environment", package = "Shennong", old = "sn_call_mmochi")
  sn_call_pixi_environment("mmochi", command = command, args = args, ...)
}

#' @rdname sn_prepare_pixi_environment
#' @export
sn_call_scarches <- function(command, args = character(), ...) {
  .Deprecated("sn_call_pixi_environment", package = "Shennong", old = "sn_call_scarches")
  sn_call_pixi_environment("scarches", command = command, args = args, ...)
}

#' @rdname sn_prepare_pixi_environment
#' @export
sn_call_scpoli <- function(command, args = character(), ...) {
  .Deprecated("sn_call_pixi_environment", package = "Shennong", old = "sn_call_scpoli")
  sn_call_pixi_environment("scpoli", command = command, args = args, ...)
}

#' @rdname sn_prepare_pixi_environment
#' @export
sn_call_infercnvpy <- function(command, args = character(), ...) {
  .Deprecated("sn_call_pixi_environment", package = "Shennong", old = "sn_call_infercnvpy")
  sn_call_pixi_environment("infercnvpy", command = command, args = args, ...)
}

#' @rdname sn_prepare_pixi_environment
#' @export
sn_call_trajectory <- function(command, args = character(), ...) {
  .Deprecated("sn_call_pixi_environment", package = "Shennong", old = "sn_call_trajectory")
  sn_call_pixi_environment("trajectory", command = command, args = args, ...)
}

#' @rdname sn_prepare_pixi_environment
#' @export
sn_call_cellphonedb <- function(command, args = character(), ...) {
  .Deprecated("sn_call_pixi_environment", package = "Shennong", old = "sn_call_cellphonedb")
  sn_call_pixi_environment("cellphonedb", command = command, args = args, ...)
}

#' @rdname sn_prepare_pixi_environment
#' @export
sn_call_cell2location <- function(command, args = character(), ...) {
  .Deprecated("sn_call_pixi_environment", package = "Shennong", old = "sn_call_cell2location")
  sn_call_pixi_environment("cell2location", command = command, args = args, ...)
}

#' @rdname sn_prepare_pixi_environment
#' @export
sn_call_tangram <- function(command, args = character(), ...) {
  .Deprecated("sn_call_pixi_environment", package = "Shennong", old = "sn_call_tangram")
  sn_call_pixi_environment("tangram", command = command, args = args, ...)
}

#' @rdname sn_prepare_pixi_environment
#' @export
sn_call_squidpy <- function(command, args = character(), ...) {
  .Deprecated("sn_call_pixi_environment", package = "Shennong", old = "sn_call_squidpy")
  sn_call_pixi_environment("squidpy", command = command, args = args, ...)
}

#' @rdname sn_prepare_pixi_environment
#' @export
sn_call_spatialdata <- function(command, args = character(), ...) {
  .Deprecated("sn_call_pixi_environment", package = "Shennong", old = "sn_call_spatialdata")
  sn_call_pixi_environment("spatialdata", command = command, args = args, ...)
}

#' @rdname sn_prepare_pixi_environment
#' @export
sn_call_stlearn <- function(command, args = character(), ...) {
  .Deprecated("sn_call_pixi_environment", package = "Shennong", old = "sn_call_stlearn")
  sn_call_pixi_environment("stlearn", command = command, args = args, ...)
}

#' Configure pixi mirrors for Shennong runtime environments
#'
#' Writes a pixi \code{config.toml} under the selected \code{PIXI_HOME}. By
#' default Shennong uses a user-level pixi home such as
#' \code{~/.shennong/pixi/home}, so this does not mutate pixi's global
#' \code{~/.pixi/config.toml} unless \code{pixi_home} points there explicitly.
#'
#' @param mirror One of \code{"default"}, \code{"auto"}, \code{"china"},
#'   \code{"tuna"}, \code{"ustc"}, or \code{"bfsu"}.
#' @param pixi_home Pixi home directory. When \code{NULL}, it is derived from
#'   \code{runtime_dir}.
#' @param runtime_dir Optional Shennong runtime directory used to derive
#'   \code{pixi_home}.
#' @param config_path Optional explicit config path.
#' @param append_original Whether to keep the original conda-forge URL as a
#'   fallback after mirror URLs.
#'
#' @return Invisibly returns the config path, or \code{NA_character_} when no
#'   mirror config was written.
#'
#' @examples
#' tmp <- tempfile("pixi-home-")
#' sn_configure_pixi_mirror("tuna", pixi_home = tmp)
#'
#' @export
sn_configure_pixi_mirror <- function(mirror = c("default", "auto", "china", "tuna", "ustc", "bfsu"),
                                     pixi_home = NULL,
                                     runtime_dir = NULL,
                                     config_path = NULL,
                                     append_original = TRUE) {
  mirror <- match.arg(mirror)
  mirror <- .sn_resolve_pixi_mirror(mirror)
  if (identical(mirror, "default")) {
    return(invisible(NA_character_))
  }

  pixi_home <- pixi_home %||% file.path(.sn_shennong_runtime_dir(runtime_dir), "pixi", "home")
  pixi_home <- path.expand(pixi_home)
  config_path <- config_path %||% file.path(pixi_home, "config.toml")
  config_path <- path.expand(config_path)
  dir.create(dirname(config_path), recursive = TRUE, showWarnings = FALSE)

  lines <- .sn_pixi_mirror_config_lines(mirror = mirror, append_original = append_original)
  writeLines(lines, con = config_path, useBytes = TRUE)
  invisible(normalizePath(config_path, winslash = "/", mustWork = TRUE))
}

.sn_pixi_environment_names <- function() {
  installed_root <- system.file("pixi", package = "Shennong")
  roots <- c(
    if (nzchar(installed_root) && dir.exists(installed_root)) installed_root else character(0),
    file.path(getwd(), "inst", "pixi")
  )
  roots <- roots[dir.exists(roots)]
  names <- unique(unlist(lapply(roots, function(root) {
    dirs <- list.dirs(root, full.names = FALSE, recursive = FALSE)
    dirs[file.exists(file.path(root, dirs, "pixi.toml"))]
  }), use.names = FALSE))
  names[nzchar(names)]
}

.sn_normalize_pixi_environment <- function(environment) {
  environment <- tolower(as.character(environment %||% "scvi"))
  environment <- switch(
    environment,
    "scvi-tools" = "scvi",
    "scverse" = "scvi",
    "scanvi" = "scvi",
    "totalvi" = "scvi",
    "total-vi" = "scvi",
    "mmochi" = "mmochi",
    "mmochi-landmark" = "mmochi",
    "mmochi_landmark" = "mmochi",
    "scpoli" = "scarches",
    "sc_poli" = "scarches",
    "scarches" = "scarches",
    "scib_metrics" = "scib-metrics",
    "scibmetrics" = "scib-metrics",
    "scib-metrics" = "scib-metrics",
    "infercnv" = "infercnvpy",
    "infercnvpy" = "infercnvpy",
    "cellphone" = "cellphonedb",
    "cellphonedb" = "cellphonedb",
    "cell_phone_db" = "cellphonedb",
    "cell2location" = "cell2location",
    "cell_2_location" = "cell2location",
    "tarngram" = "tangram",
    "tangram" = "tangram",
    "tangram-sc" = "tangram",
    "squidpy" = "squidpy",
    "spatialdata" = "spatialdata",
    "spatial-data" = "spatialdata",
    "stlearn" = "stlearn",
    "stlearnr" = "stlearn",
    environment
  )
  known <- .sn_pixi_environment_names()
  if (!environment %in% known) {
    stop(
      "Unknown pixi environment '", environment, "'. Available environments: ",
      paste(known, collapse = ", "),
      call. = FALSE
    )
  }
  environment
}

.sn_pixi_gpu_aware_environment <- function(environment) {
  environment %in% c("scvi", "scarches", "scib-metrics", "cell2location", "tangram")
}

.sn_select_pixi_environment <- function(environment, pixi_environment = "auto", cuda_version = NULL) {
  pixi_environment <- match.arg(pixi_environment, c("auto", "default", "cpu", "gpu"))
  if (!.sn_pixi_gpu_aware_environment(environment)) {
    selected <- "default"
    accelerator <- "cpu"
    detected_cuda <- cuda_version %||% "12.6"
  } else {
    accelerator_info <- .sn_resolve_scvi_accelerator(if (identical(pixi_environment, "gpu")) "cuda" else pixi_environment)
    selected <- if (identical(pixi_environment, "auto")) accelerator_info$environment else pixi_environment
    accelerator <- accelerator_info$requested
    detected_cuda <- cuda_version %||% accelerator_info$cuda_version
  }
  list(
    pixi_environment = selected,
    accelerator = accelerator,
    cuda_version = .sn_default_scvi_cuda_version(detected_cuda)
  )
}

.sn_render_pixi_config <- function(environment, platforms = NULL, cuda_version = "12.6") {
  template <- readLines(sn_get_pixi_config_path(environment), warn = FALSE)
  cuda_version <- .sn_normalize_cuda_requirement(cuda_version)
  cuda_major <- sub("\\..*$", "", cuda_version)
  platforms <- platforms %||% .sn_current_pixi_platform()
  platform <- paste(unique(as.character(platforms)), collapse = '", "')
  replacements <- c(
    "{{ platform }}" = platform,
    "{{ cuda_version }}" = cuda_version,
    "{{ cuda_major }}" = cuda_major
  )
  for (pattern in names(replacements)) {
    template <- gsub(pattern, replacements[[pattern]], template, fixed = TRUE)
  }
  template
}

.sn_bundled_pixi_lock_path <- function(environment) {
  environment <- .sn_normalize_pixi_environment(environment)
  installed <- system.file("pixi", environment, "pixi.lock", package = "Shennong")
  if (nzchar(installed) && file.exists(installed)) {
    return(normalizePath(installed, winslash = "/", mustWork = TRUE))
  }
  source_path <- file.path(getwd(), "inst", "pixi", environment, "pixi.lock")
  if (file.exists(source_path)) {
    return(normalizePath(source_path, winslash = "/", mustWork = TRUE))
  }
  stop("Could not locate the bundled pixi lock for environment: ", environment, call. = FALSE)
}

.sn_pixi_lock_platforms <- function(lock_path) {
  lines <- readLines(lock_path, warn = FALSE)
  platform_start <- which(trimws(lines) == "platforms:")
  environment_start <- which(trimws(lines) == "environments:")
  if (length(platform_start) != 1L || length(environment_start) == 0L ||
      environment_start[[1L]] <= platform_start[[1L]]) {
    stop("Pixi lock is missing a parseable platform inventory: ", lock_path, call. = FALSE)
  }
  section <- lines[seq.int(platform_start[[1L]] + 1L, environment_start[[1L]] - 1L)]
  values <- sub("^[[:space:]]*-[[:space:]]*name:[[:space:]]*", "", section)
  values <- values[grepl("^[[:space:]]*-[[:space:]]*name:", section)]
  values <- gsub("['\"]", "", trimws(values))
  values[nzchar(values)]
}

.sn_prepare_pixi_lock <- function(environment,
                                  manifest_path,
                                  platforms,
                                  cuda_version = "12.6",
                                  overwrite = FALSE) {
  lock_path <- file.path(dirname(manifest_path), "pixi.lock")
  bundled_manifest <- .sn_render_pixi_config(
    environment = environment,
    # Bundled locks are solved against the package-tested CUDA requirement.
    # A caller changing that requirement owns a matching adjacent lock.
    cuda_version = "12.6"
  )
  current_manifest <- readLines(manifest_path, warn = FALSE)
  uses_bundled_manifest <- identical(current_manifest, bundled_manifest)
  if (uses_bundled_manifest) {
    source_lock <- .sn_bundled_pixi_lock_path(environment)
    same_path <- identical(
      normalizePath(source_lock, winslash = "/", mustWork = TRUE),
      normalizePath(lock_path, winslash = "/", mustWork = FALSE)
    )
    if (!same_path && (!file.exists(lock_path) || isTRUE(overwrite))) {
      if (!file.copy(source_lock, lock_path, overwrite = TRUE, copy.mode = TRUE)) {
        stop("Could not copy the bundled pixi lock to: ", lock_path, call. = FALSE)
      }
    }
  } else if (file.exists(lock_path)) {
    source_lock <- .sn_bundled_pixi_lock_path(environment)
    if (identical(
      digest::digest(lock_path, algo = "sha256", file = TRUE, serialize = FALSE),
      digest::digest(source_lock, algo = "sha256", file = TRUE, serialize = FALSE)
    )) {
      stop(
        "The custom pixi manifest differs from the package-tested manifest, ",
        "but its adjacent lock is still the bundled lock. Generate and review ",
        "a matching `pixi.lock` before running the backend.",
        call. = FALSE
      )
    }
  }
  if (!file.exists(lock_path)) {
    stop(
      "A custom pixi manifest requires a matching `pixi.lock` beside it; refusing an unpinned solve: ",
      manifest_path, ".",
      call. = FALSE
    )
  }
  available <- .sn_pixi_lock_platforms(lock_path)
  missing <- setdiff(platforms, available)
  if (length(missing) > 0L) {
    stop(
      "The pixi lock does not cover requested platform(s): ",
      paste(missing, collapse = ", "), ". Available: ",
      paste(available, collapse = ", "), ".",
      call. = FALSE
    )
  }
  normalizePath(lock_path, winslash = "/", mustWork = TRUE)
}

.sn_pixi_command_env <- function(pixi_home = NULL) {
  if (is.null(pixi_home) || !nzchar(pixi_home)) {
    return(character(0))
  }
  dir.create(pixi_home, recursive = TRUE, showWarnings = FALSE)
  paste0("PIXI_HOME=", normalizePath(pixi_home, winslash = "/", mustWork = TRUE))
}

.sn_pixi_install_environment <- function(pixi, manifest_path, pixi_environment = "default", pixi_home = NULL) {
  args <- c(
    "install",
    "--locked",
    "--manifest-path", manifest_path,
    if (!is.null(pixi_environment) && nzchar(pixi_environment)) c("--environment", pixi_environment) else character(0)
  )
  status <- .sn_system2_capture(
    command = pixi,
    args = args,
    env = .sn_pixi_command_env(pixi_home)
  )
  exit_code <- attr(status, "status") %||% 0L
  if (!identical(exit_code, 0L)) {
    stop("pixi environment installation failed.\n", paste(status, collapse = "\n"), call. = FALSE)
  }
  invisible(status)
}

.sn_pixi_run_command <- function(pixi, manifest_path, pixi_environment = "default", command, args = character(), pixi_home = NULL) {
  pixi_args <- c(
    "run",
    "--locked",
    "--manifest-path", manifest_path,
    if (!is.null(pixi_environment) && nzchar(pixi_environment)) c("--environment", pixi_environment) else character(0),
    command,
    args
  )
  status <- .sn_system2_capture(
    command = pixi,
    args = pixi_args,
    env = .sn_pixi_command_env(pixi_home)
  )
  exit_code <- attr(status, "status") %||% 0L
  if (!identical(exit_code, 0L)) {
    stop("pixi command failed.\n", paste(status, collapse = "\n"), call. = FALSE)
  }
  status
}

.sn_pixi_binary_name <- function() {
  if (.Platform$OS.type == "windows") "pixi.exe" else "pixi"
}

.sn_normalize_requested_pixi_version <- function(version) {
  if (length(version) != 1L || is.na(version)) {
    stop("`version` must be one explicit semantic Pixi release.", call. = FALSE)
  }
  version <- sub("^v", "", trimws(as.character(version)))
  if (!grepl("^[0-9]+\\.[0-9]+\\.[0-9]+$", version)) {
    stop(
      "`version` must be an explicit semantic Pixi release such as '0.69.0'; ",
      "mutable 'latest' downloads are not permitted.",
      call. = FALSE
    )
  }
  version
}

.sn_pixi_version_matches <- function(installed, requested) {
  if (is.null(installed) || length(installed) != 1L || is.na(installed)) {
    return(FALSE)
  }
  requested <- .sn_normalize_requested_pixi_version(requested)
  match <- regmatches(
    as.character(installed),
    regexpr("[0-9]+\\.[0-9]+\\.[0-9]+", as.character(installed))
  )
  length(match) == 1L && nzchar(match) && identical(match, requested)
}

.sn_resolve_pixi_path <- function(pixi = NULL) {
  if (!is.null(pixi)) {
    explicit <- path.expand(as.character(pixi)[[1L]])
    if (!nzchar(explicit) || !file.exists(explicit)) return(NULL)
    return(explicit)
  }
  candidates <- c(
    getOption("shennong.pixi", NULL),
    Sys.getenv("SHENNONG_PIXI", unset = ""),
    Sys.which("pixi")[["pixi"]],
    file.path("~", ".pixi", "bin", .sn_pixi_binary_name())
  )
  candidates <- path.expand(candidates[!is.na(candidates) & nzchar(candidates)])
  candidates <- candidates[file.exists(candidates)]
  if (length(candidates) == 0L) {
    return(NULL)
  }
  candidates[[1]]
}

.sn_pixi_release_asset <- function(version,
                                   sysname = Sys.info()[["sysname"]],
                                   machine = Sys.info()[["machine"]]) {
  version <- .sn_normalize_requested_pixi_version(version)
  machine <- tolower(as.character(machine %||% ""))
  architecture <- if (machine %in% c("x86_64", "x86-64", "amd64", "x64")) {
    "x86_64"
  } else if (machine %in% c("aarch64", "arm64")) {
    "aarch64"
  } else {
    stop("Unsupported Pixi installation architecture: ", machine, ".", call. = FALSE)
  }
  os <- tolower(as.character(sysname %||% ""))
  target <- if (identical(os, "linux")) {
    paste0(architecture, "-unknown-linux-musl")
  } else if (identical(os, "darwin")) {
    paste0(architecture, "-apple-darwin")
  } else if (grepl("windows|mingw", os)) {
    paste0(architecture, "-pc-windows-msvc")
  } else {
    stop("Unsupported Pixi installation operating system: ", os, ".", call. = FALSE)
  }
  extension <- if (grepl("windows", target, fixed = TRUE)) ".zip" else ".tar.gz"
  asset <- paste0("pixi-", target, extension)
  list(
    version = version,
    asset = asset,
    url = paste0(
      "https://github.com/prefix-dev/pixi/releases/download/v", version, "/", asset
    )
  )
}

.sn_download_file <- function(url, destfile, quiet = FALSE) {
  tryCatch(
    {
      utils::download.file(url = url, destfile = destfile, mode = "wb", quiet = quiet)
      TRUE
    },
    error = function(e) {
      stop("Download failed: ", conditionMessage(e), call. = FALSE)
    }
  )
}

.sn_verify_file_sha256 <- function(path, expected) {
  expected <- tolower(trimws(as.character(expected)[[1L]]))
  if (!grepl("^[0-9a-f]{64}$", expected)) {
    stop("`sha256` must be exactly 64 hexadecimal characters.", call. = FALSE)
  }
  actual <- tolower(digest::digest(path, algo = "sha256", file = TRUE, serialize = FALSE))
  if (!identical(actual, expected)) {
    stop(
      "Pixi download SHA-256 mismatch: expected ", expected,
      ", received ", actual, ".",
      call. = FALSE
    )
  }
  invisible(actual)
}

.sn_download_verified_pixi <- function(version, download_url = NULL,
                                       sha256 = NULL, quiet = FALSE) {
  release <- .sn_pixi_release_asset(version)
  custom_url <- !is.null(download_url) && nzchar(download_url)
  url <- if (custom_url) download_url else release$url
  asset <- basename(sub("[?#].*$", "", url))
  if (!nzchar(asset)) asset <- release$asset
  if (custom_url && is.null(sha256)) {
    stop("Custom `download_url` requires an explicit `sha256` digest.", call. = FALSE)
  }
  if (is.null(sha256)) {
    checksum_file <- tempfile("pixi-checksum-")
    on.exit(unlink(checksum_file, force = TRUE), add = TRUE)
    .sn_download_file(paste0(url, ".sha256"), checksum_file, quiet = quiet)
    checksum_text <- readLines(checksum_file, warn = FALSE)
    checksum_match <- regmatches(
      checksum_text,
      regexpr("[0-9A-Fa-f]{64}", checksum_text)
    )
    checksum_match <- checksum_match[nzchar(checksum_match)]
    if (length(checksum_match) == 0L) {
      stop("The Pixi release checksum sidecar did not contain a SHA-256 digest.", call. = FALSE)
    }
    sha256 <- checksum_match[[1L]]
  }
  suffix <- if (grepl("\\.tar\\.gz$", asset)) ".tar.gz" else paste0(".", tools::file_ext(asset))
  archive <- tempfile("pixi-release-", fileext = suffix)
  ok <- .sn_download_file(url = url, destfile = archive, quiet = quiet)
  if (!isTRUE(ok) || !file.exists(archive) || file.info(archive)$size <= 0) {
    stop("Failed to download the Pixi release asset from ", url, ".", call. = FALSE)
  }
  .sn_verify_file_sha256(archive, sha256)
  list(path = archive, asset = asset, sha256 = tolower(sha256), url = url)
}

.sn_install_verified_pixi_archive <- function(archive, asset, install_dir) {
  extract_dir <- tempfile("pixi-extract-")
  dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(extract_dir, recursive = TRUE, force = TRUE), add = TRUE)
  if (grepl("\\.zip$", asset, ignore.case = TRUE)) {
    utils::unzip(archive, exdir = extract_dir)
  } else if (grepl("\\.tar\\.gz$", asset, ignore.case = TRUE)) {
    utils::untar(archive, exdir = extract_dir)
  } else {
    file.copy(archive, file.path(extract_dir, .sn_pixi_binary_name()), overwrite = TRUE)
  }
  candidates <- list.files(extract_dir, recursive = TRUE, full.names = TRUE)
  candidates <- candidates[basename(candidates) == .sn_pixi_binary_name()]
  if (length(candidates) != 1L) {
    stop("Verified Pixi release did not contain exactly one Pixi executable.", call. = FALSE)
  }
  dir.create(install_dir, recursive = TRUE, showWarnings = FALSE)
  target <- file.path(install_dir, .sn_pixi_binary_name())
  if (!file.copy(candidates[[1L]], target, overwrite = TRUE, copy.mode = TRUE)) {
    stop("Could not install the verified Pixi executable at ", target, ".", call. = FALSE)
  }
  if (.Platform$OS.type != "windows") Sys.chmod(target, mode = "0755")
  normalizePath(target, winslash = "/", mustWork = TRUE)
}

.sn_detect_nvidia_gpu <- function() {
  nvidia_smi <- Sys.which("nvidia-smi")[["nvidia-smi"]]
  if (is.na(nvidia_smi) || !nzchar(nvidia_smi)) {
    return(list(available = FALSE, devices = character(0), cuda_version = NA_character_, raw = character(0)))
  }

  query <- tryCatch(
    suppressWarnings(system2(
      nvidia_smi,
      args = c("--query-gpu=name", "--format=csv,noheader"),
      stdout = TRUE,
      stderr = TRUE
    )),
    error = function(e) character(0)
  )
  status <- attr(query, "status") %||% 0L
  devices <- trimws(query)
  devices <- devices[nzchar(devices) & !grepl("failed|error", devices, ignore.case = TRUE)]
  raw <- tryCatch(
    suppressWarnings(system2(nvidia_smi, stdout = TRUE, stderr = TRUE)),
    error = function(e) character(0)
  )
  cuda_version <- .sn_parse_nvidia_cuda_version(raw)

  list(
    available = identical(status, 0L) && length(devices) > 0L,
    devices = devices,
    cuda_version = cuda_version,
    raw = raw
  )
}

.sn_parse_nvidia_cuda_version <- function(raw) {
  if (length(raw) == 0L) {
    return(NA_character_)
  }
  hit <- regmatches(raw, regexpr("CUDA Version:\\s*[0-9.]+", raw))
  hit <- hit[nzchar(hit)]
  if (length(hit) == 0L) {
    return(NA_character_)
  }
  sub("CUDA Version:\\s*", "", hit[[1]])
}

.sn_resolve_pixi_mirror <- function(mirror = c("default", "auto", "china", "tuna", "ustc", "bfsu")) {
  mirror <- match.arg(mirror)
  env_mirror <- Sys.getenv("SHENNONG_PIXI_MIRROR", unset = "")
  if (!identical(mirror, "auto")) {
    return(mirror)
  }
  if (nzchar(env_mirror)) {
    env_mirror <- tolower(env_mirror)
    if (env_mirror %in% c("default", "china", "tuna", "ustc", "bfsu")) {
      return(env_mirror)
    }
  }
  locale <- paste(Sys.getlocale(), Sys.getenv(c("LANG", "LC_ALL", "TZ"), unset = ""), collapse = " ")
  if (grepl("zh_CN|China|Asia/Shanghai|Asia/Chongqing|Asia/Urumqi|CN", locale, ignore.case = TRUE)) {
    return("china")
  }
  "default"
}

.sn_pixi_mirror_urls <- function(mirror) {
  switch(
    mirror,
    china = c(
      "https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge",
      "https://mirrors.ustc.edu.cn/anaconda/cloud/conda-forge",
      "https://mirrors.bfsu.edu.cn/anaconda/cloud/conda-forge"
    ),
    tuna = "https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge",
    ustc = "https://mirrors.ustc.edu.cn/anaconda/cloud/conda-forge",
    bfsu = "https://mirrors.bfsu.edu.cn/anaconda/cloud/conda-forge",
    character(0)
  )
}

.sn_pixi_mirror_pypi_urls <- function(mirror) {
  switch(
    mirror,
    china = list(
      simple = c(
        "https://pypi.tuna.tsinghua.edu.cn/simple",
        "https://mirrors.ustc.edu.cn/pypi/simple",
        "https://mirrors.bfsu.edu.cn/pypi/web/simple"
      ),
      packages = c(
        "https://pypi.tuna.tsinghua.edu.cn/packages",
        "https://mirrors.ustc.edu.cn/pypi/packages",
        "https://mirrors.bfsu.edu.cn/pypi/web/packages"
      )
    ),
    tuna = list(
      simple = "https://pypi.tuna.tsinghua.edu.cn/simple",
      packages = "https://pypi.tuna.tsinghua.edu.cn/packages"
    ),
    ustc = list(
      simple = "https://mirrors.ustc.edu.cn/pypi/simple",
      packages = "https://mirrors.ustc.edu.cn/pypi/packages"
    ),
    bfsu = list(
      simple = "https://mirrors.bfsu.edu.cn/pypi/web/simple",
      packages = "https://mirrors.bfsu.edu.cn/pypi/web/packages"
    ),
    list(simple = character(0), packages = character(0))
  )
}

.sn_pixi_mirror_config_lines <- function(mirror, append_original = TRUE) {
  conda_original <- "https://conda.anaconda.org/conda-forge"
  pypi_simple <- "https://pypi.org/simple"
  pypi_packages <- "https://files.pythonhosted.org/packages"
  conda_urls <- .sn_pixi_mirror_urls(mirror)
  pypi_urls <- .sn_pixi_mirror_pypi_urls(mirror)
  if (isTRUE(append_original)) {
    conda_urls <- unique(c(conda_urls, conda_original))
    pypi_urls$simple <- unique(c(pypi_urls$simple, pypi_simple))
    pypi_urls$packages <- unique(c(pypi_urls$packages, pypi_packages))
  }
  conda_quoted <- paste0('"', conda_urls, '"', collapse = ",\n  ")
  pypi_simple_quoted <- paste0('"', pypi_urls$simple, '"', collapse = ",\n  ")
  pypi_packages_quoted <- paste0('"', pypi_urls$packages, '"', collapse = ",\n  ")
  c(
    "# Generated by Shennong. Edit or remove this file to change pixi mirrors.",
    "[mirrors]",
    paste0('"', conda_original, '" = ['),
    paste0("  ", conda_quoted),
    "]",
    paste0('"', pypi_simple, '" = ['),
    paste0("  ", pypi_simple_quoted),
    "]",
    paste0('"', pypi_packages, '" = ['),
    paste0("  ", pypi_packages_quoted),
    "]",
    ""
  )
}

#' Deprecated alias of `sn_get_pixi_paths()`
#'
#' `sn_pixi_paths()` is a deprecated compatibility alias. Use [sn_get_pixi_paths()] directly;
#' the alias will be removed in a future release.
#'
#' @param ... Named arguments passed on to [sn_get_pixi_paths()].
#'
#' @return Result of `sn_get_pixi_paths(...)`.
#'
#' @export
sn_pixi_paths <- function(...) {
  .Deprecated("sn_get_pixi_paths", package = "Shennong")
  sn_get_pixi_paths(...)
}

#' Deprecated alias of `sn_get_pixi_config_path()`
#'
#' `sn_pixi_config_path()` is a deprecated compatibility alias. Use [sn_get_pixi_config_path()] directly;
#' the alias will be removed in a future release.
#'
#' @param ... Named arguments passed on to [sn_get_pixi_config_path()].
#'
#' @return Result of `sn_get_pixi_config_path(...)`.
#'
#' @export
sn_pixi_config_path <- function(...) {
  .Deprecated("sn_get_pixi_config_path", package = "Shennong")
  sn_get_pixi_config_path(...)
}
