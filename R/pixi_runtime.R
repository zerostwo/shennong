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
#' @param install Logical used by \code{sn_ensure_pixi()}; install pixi when it
#'   is not already available.
#' @param version Pixi version passed to the official installer script.
#'   Defaults to \code{"latest"}.
#' @param pixi_home Pixi home directory. The installer defaults to
#'   \code{"~/.pixi"}; Shennong runtime workflows set \code{PIXI_HOME} to a
#'   user-level path such as \code{"~/.shennong/pixi/home"}.
#' @param bin_dir Optional directory where the standalone pixi binary should be
#'   installed.
#' @param no_path_update Logical; when \code{TRUE}, ask the installer not to
#'   edit shell startup files.
#' @param download_url Optional custom pixi binary download URL. This is useful
#'   for institutional mirrors.
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
      suppressWarnings(system2(pixi_path, "--version", stdout = TRUE, stderr = TRUE)),
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
sn_install_pixi <- function(version = "latest",
                            pixi_home = "~/.pixi",
                            bin_dir = NULL,
                            no_path_update = TRUE,
                            download_url = NULL,
                            force = FALSE,
                            quiet = FALSE) {
  current <- sn_check_pixi(quiet = TRUE)
  if (isTRUE(current$installed) && !isTRUE(force)) {
    if (!quiet) {
      .sn_log_info("pixi is already installed at {current$path}.")
    }
    return(invisible(current))
  }

  installer <- .sn_download_pixi_installer(quiet = quiet)
  on.exit(unlink(installer), add = TRUE)

  env <- .sn_pixi_installer_env(
    version = version,
    pixi_home = pixi_home,
    bin_dir = bin_dir,
    no_path_update = no_path_update,
    download_url = download_url
  )
  status <- .sn_run_pixi_installer(installer = installer, env = env, quiet = quiet)
  if (!identical(status, 0L)) {
    stop("pixi installation failed with exit status ", status, ".", call. = FALSE)
  }

  installed_binary <- if (!is.null(bin_dir) && nzchar(bin_dir)) {
    file.path(path.expand(bin_dir), .sn_pixi_binary_name())
  } else {
    file.path(path.expand(pixi_home), "bin", .sn_pixi_binary_name())
  }
  refreshed <- sn_check_pixi(pixi = installed_binary, quiet = quiet)
  if (!isTRUE(refreshed$installed)) {
    refreshed <- sn_check_pixi(quiet = quiet)
  }
  if (!isTRUE(refreshed$installed)) {
    stop("pixi installer completed, but the pixi executable could not be found.", call. = FALSE)
  }

  invisible(refreshed)
}

#' @rdname sn_check_pixi
#' @export
sn_ensure_pixi <- function(pixi = NULL,
                           install = TRUE,
                           version = "latest",
                           pixi_home = "~/.pixi",
                           bin_dir = NULL,
                           no_path_update = TRUE,
                           download_url = NULL,
                           quiet = FALSE) {
  info <- sn_check_pixi(pixi = pixi, quiet = TRUE)
  if (isTRUE(info$installed)) {
    if (!quiet) {
      .sn_log_info("pixi detected at {info$path}.")
    }
    return(invisible(info))
  }

  if (!isTRUE(install)) {
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
#' @param environment Python environment name.
#' @param pixi_environment Pixi environment inside the manifest, for example
#'   \code{"cpu"}, \code{"gpu"}, or \code{"default"}. \code{"auto"} uses
#'   CUDA when available for GPU-aware configs and otherwise CPU/default.
#' @param runtime_dir Optional Shennong runtime directory.
#' @param project_dir Optional explicit pixi workspace directory.
#' @param manifest_path Optional explicit materialized manifest path.
#' @param overwrite Whether to overwrite an existing materialized manifest.
#' @param cuda_version CUDA runtime version used when rendering templates.
#' @param platforms Pixi platform vector. Defaults to the current platform.
#' @param mirror Mirror setting passed to \code{sn_configure_pixi_mirror()}.
#' @param install_pixi Ensure the standalone pixi binary is available.
#' @param install_environment Run \code{pixi install} for the selected
#'   environment after materializing the manifest.
#' @param pixi Optional pixi executable path.
#' @param pixi_version Pixi version used if installation is needed.
#' @param pixi_download_url Optional custom pixi binary download URL.
#' @param command Command to run inside the pixi environment.
#' @param args Character vector of command arguments.
#' @param quiet Logical; suppress status messages where possible.
#' @param ... Additional arguments passed from environment-specific helpers to
#'   \code{sn_call_pixi_environment()}.
#'
#' @return \code{sn_prepare_pixi_environment()} returns a named list of paths
#'   and selected environment metadata. \code{sn_call_pixi_environment()}
#'   invisibly returns command output.
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
                                        pixi_version = "latest",
                                        pixi_download_url = NULL,
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
                                     pixi_version = "latest",
                                     pixi_download_url = NULL,
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
    quiet = quiet
  )
  pixi_info <- sn_ensure_pixi(
    pixi = pixi,
    install = install_pixi,
    version = pixi_version,
    pixi_home = "~/.pixi",
    no_path_update = TRUE,
    download_url = pixi_download_url,
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
  sn_call_pixi_environment("scvi", command = command, args = args, ...)
}

#' @rdname sn_prepare_pixi_environment
#' @export
sn_call_scanvi <- function(command, args = character(), ...) {
  sn_call_pixi_environment("scanvi", command = command, args = args, ...)
}

#' @rdname sn_prepare_pixi_environment
#' @export
sn_call_mmochi <- function(command, args = character(), ...) {
  sn_call_pixi_environment("mmochi", command = command, args = args, ...)
}

#' @rdname sn_prepare_pixi_environment
#' @export
sn_call_scarches <- function(command, args = character(), ...) {
  sn_call_pixi_environment("scarches", command = command, args = args, ...)
}

#' @rdname sn_prepare_pixi_environment
#' @export
sn_call_scpoli <- function(command, args = character(), ...) {
  sn_call_pixi_environment("scpoli", command = command, args = args, ...)
}

#' @rdname sn_prepare_pixi_environment
#' @export
sn_call_infercnvpy <- function(command, args = character(), ...) {
  sn_call_pixi_environment("infercnvpy", command = command, args = args, ...)
}

#' @rdname sn_prepare_pixi_environment
#' @export
sn_call_trajectory <- function(command, args = character(), ...) {
  sn_call_pixi_environment("trajectory", command = command, args = args, ...)
}

#' @rdname sn_prepare_pixi_environment
#' @export
sn_call_cellphonedb <- function(command, args = character(), ...) {
  sn_call_pixi_environment("cellphonedb", command = command, args = args, ...)
}

#' @rdname sn_prepare_pixi_environment
#' @export
sn_call_cell2location <- function(command, args = character(), ...) {
  sn_call_pixi_environment("cell2location", command = command, args = args, ...)
}

#' @rdname sn_prepare_pixi_environment
#' @export
sn_call_tangram <- function(command, args = character(), ...) {
  sn_call_pixi_environment("tangram", command = command, args = args, ...)
}

#' @rdname sn_prepare_pixi_environment
#' @export
sn_call_squidpy <- function(command, args = character(), ...) {
  sn_call_pixi_environment("squidpy", command = command, args = args, ...)
}

#' @rdname sn_prepare_pixi_environment
#' @export
sn_call_spatialdata <- function(command, args = character(), ...) {
  sn_call_pixi_environment("spatialdata", command = command, args = args, ...)
}

#' @rdname sn_prepare_pixi_environment
#' @export
sn_call_stlearn <- function(command, args = character(), ...) {
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
    "--manifest-path", manifest_path,
    if (!is.null(pixi_environment) && nzchar(pixi_environment)) c("--environment", pixi_environment) else character(0)
  )
  status <- system2(
    command = pixi,
    args = args,
    env = .sn_pixi_command_env(pixi_home),
    stdout = TRUE,
    stderr = TRUE
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
    "--manifest-path", manifest_path,
    if (!is.null(pixi_environment) && nzchar(pixi_environment)) c("--environment", pixi_environment) else character(0),
    command,
    args
  )
  status <- system2(
    command = pixi,
    args = pixi_args,
    env = .sn_pixi_command_env(pixi_home),
    stdout = TRUE,
    stderr = TRUE
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

.sn_resolve_pixi_path <- function(pixi = NULL) {
  candidates <- c(
    pixi,
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

.sn_download_pixi_installer <- function(quiet = FALSE) {
  ext <- if (.Platform$OS.type == "windows") ".ps1" else ".sh"
  url <- if (.Platform$OS.type == "windows") "https://pixi.sh/install.ps1" else "https://pixi.sh/install.sh"
  installer <- tempfile("pixi-install-", fileext = ext)
  ok <- .sn_download_file(url = url, destfile = installer, quiet = quiet)
  if (!isTRUE(ok) || !file.exists(installer)) {
    stop("Failed to download the pixi installer from ", url, ".", call. = FALSE)
  }
  installer
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

.sn_pixi_installer_env <- function(version = "latest",
                                   pixi_home = "~/.pixi",
                                   bin_dir = NULL,
                                   no_path_update = TRUE,
                                   download_url = NULL) {
  env <- c(
    paste0("PIXI_VERSION=", version),
    paste0("PIXI_HOME=", path.expand(pixi_home))
  )
  if (!is.null(bin_dir) && nzchar(bin_dir)) {
    env <- c(env, paste0("PIXI_BIN_DIR=", path.expand(bin_dir)))
  }
  if (isTRUE(no_path_update)) {
    env <- c(env, "PIXI_NO_PATH_UPDATE=1")
  }
  if (!is.null(download_url) && nzchar(download_url)) {
    env <- c(env, paste0("PIXI_DOWNLOAD_URL=", download_url))
  }
  env
}

.sn_run_pixi_installer <- function(installer, env, quiet = FALSE) {
  output <- if (quiet) FALSE else ""
  if (.Platform$OS.type == "windows") {
    status <- system2(
      "powershell",
      args = c("-ExecutionPolicy", "Bypass", "-File", shQuote(installer)),
      env = env,
      stdout = output,
      stderr = output
    )
  } else {
    status <- system2(
      "sh",
      args = shQuote(installer),
      env = env,
      stdout = output,
      stderr = output
    )
  }
  status %||% 0L
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
