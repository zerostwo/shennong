#' Return installed Shennong Codex asset paths
#'
#' This helper exposes the Codex-related assets bundled with the installed
#' \pkg{Shennong} package. It can return the package-level Codex asset root, the
#' packaged project template, or the package and project skill directories.
#'
#' @param component One of \code{"codex_root"}, \code{"package_skills"},
#'   \code{"project_template"}, or \code{"project_template_skills"}.
#'
#' @return A character scalar path to the requested bundled asset directory.
#'
#' @examples
#' if (requireNamespace("Shennong", quietly = TRUE)) {
#'   sn_get_codex_skill_path("codex_root")
#'   sn_get_codex_skill_path("package_skills")
#' }
#' @export
sn_get_codex_skill_path <- function(
  component = c("codex_root", "package_skills", "project_template", "project_template_skills")
) {
  .sn_codex_component_path(component = component)
}

#' Install the bundled Shennong Codex skill for end-user agents
#'
#' This helper installs the packaged Shennong Codex skills into a user skill
#' directory, for example \code{~/.agents/skills}. Package usage skills are
#' distinct from initialized-project governance skills.
#'
#' @param path Destination directory that will contain the \code{shennong}
#'   skill folders. Defaults to \code{"~/.agents/skills"}.
#' @param type One of \code{"package_skills"}, \code{"project_skills"}, or
#'   \code{"both"}.
#' @param overwrite Whether to replace an existing installed copy of the target
#'   skill directories. Defaults to \code{TRUE}.
#'
#' @return Invisibly returns a named character vector of installed directories.
#'
#' @examples
#' target_dir <- file.path(tempdir(), "skills")
#' installed <- sn_install_codex_skill(
#'   path = target_dir,
#'   type = "package_skills",
#'   overwrite = TRUE
#' )
#' all(dir.exists(installed))
#' @export
sn_install_codex_skill <- function(
  path = "~/.agents/skills",
  type = c("package_skills", "project_skills", "both"),
  overwrite = TRUE
) {
  type <- match.arg(type)
  dest_root <- path.expand(path)
  if (!dir.exists(dest_root)) {
    dir.create(dest_root, recursive = TRUE, showWarnings = FALSE)
  }

  source_dirs <- switch(
    type,
    package_skills = .sn_skill_source_directories("package_skills"),
    project_skills = .sn_skill_source_directories("project_template_skills"),
    both = c(
      .sn_skill_source_directories("package_skills"),
      .sn_skill_source_directories("project_template_skills")
    )
  )

  installed <- character(0)
  for (source_dir in source_dirs) {
    dest_dir <- file.path(dest_root, basename(source_dir))
    if (dir.exists(dest_dir)) {
      if (!isTRUE(overwrite)) {
        stop(glue("Destination '{dest_dir}' already exists. Set `overwrite = TRUE` to replace it."))
      }
      unlink(dest_dir, recursive = TRUE, force = TRUE)
    }

    ok <- file.copy(from = source_dir, to = dest_root, recursive = TRUE)
    if (!isTRUE(ok)) {
      stop(glue("Failed to copy the bundled Codex skill '{basename(source_dir)}' to '{dest_root}'."))
    }

    installed <- c(installed, stats::setNames(dest_dir, basename(source_dir)))
  }

  invisible(installed)
}

#' Initialize a Shennong analysis project
#'
#' This helper bootstraps a clean analysis-project structure for end users who
#' want to analyze their own data with \pkg{Shennong}. When
#' \code{with_agent = TRUE}, it materializes the packaged project-governance
#' template; otherwise it creates the same project skeleton without the
#' governance layer.
#'
#' @param path Project directory to initialize. Defaults to the current working
#'   directory.
#' @param project_name Optional human-readable project name. If \code{NULL}, the
#'   final path component of \code{path} is used.
#' @param objective Short project objective written into the generated project
#'   README, memory, and config files.
#' @param with_agent Logical; if \code{TRUE}, materialize the full Codex-facing
#'   governance template. Defaults to \code{TRUE}.
#' @param overwrite Whether to replace existing managed files copied from the
#'   template. Defaults to \code{FALSE}. Managed files include the generated
#'   project \code{.gitignore} and \code{.Rproj} file.
#'
#' @return Invisibly returns a named list containing the initialized project
#'   directory plus the main created directory and file paths.
#'
#' @examples
#' project_dir <- file.path(tempdir(), "analysis-project")
#' created <- sn_initialize_project(
#'   path = project_dir,
#'   project_name = "PBMC pilot study",
#'   objective = "Build a reproducible Shennong-based PBMC analysis workflow.",
#'   with_agent = TRUE,
#'   overwrite = TRUE
#' )
#' file.exists(created$readme)
#' file.exists(created$gitignore)
#' file.exists(created$rproj)
#' @export
sn_initialize_project <- function(
  path = ".",
  project_name = NULL,
  objective = "Build a reproducible Shennong-based single-cell analysis workflow for this project.",
  with_agent = TRUE,
  overwrite = FALSE
) {
  if (isTRUE(with_agent)) {
    return(.sn_initialize_from_project_template(
      path = path,
      project_name = project_name,
      objective = objective,
      overwrite = overwrite,
      include_governance = TRUE
    ))
  }

  .sn_initialize_from_project_template(
    path = path,
    project_name = project_name,
    objective = objective,
    overwrite = overwrite,
    include_governance = FALSE
  )
}

#' List Shennong runtime and recommended R package dependencies
#'
#' This helper reads the package dependency declaration and returns a tidy table
#' covering required imports plus recommended optional packages from
#' \code{Suggests}. It also annotates the expected installation source, GitHub
#' remote when relevant, and whether each package is already installed.
#'
#' @param scope One of \code{"all"}, \code{"required"}, or
#'   \code{"recommended"}.
#'
#' @return A tibble with package names, requirement class, declared field,
#'   expected source, GitHub remote when relevant, and installed-version
#'   metadata.
#'
#' @examples
#' deps <- sn_list_dependencies()
#' head(deps)
#' subset(deps, !installed & requirement == "recommended")
#' @export
sn_list_dependencies <- function(scope = c("all", "required", "recommended")) {
  scope <- match.arg(scope)
  deps <- .sn_dependency_table()

  if (identical(scope, "required")) {
    deps <- deps[deps$requirement == "required", , drop = FALSE]
  } else if (identical(scope, "recommended")) {
    deps <- deps[deps$requirement == "recommended", , drop = FALSE]
  }

  tibble::as_tibble(deps)
}

#' Install missing Shennong dependencies in one step
#'
#' This helper installs missing required and/or recommended packages using the
#' declared source for each dependency. CRAN packages are installed with
#' \code{install.packages()}, Bioconductor packages with
#' \code{BiocManager::install()}, and GitHub packages with
#' \code{remotes::install_github()}.
#' Legacy \code{.qs} files remain supported when the archived \pkg{qs} package
#' is already available, but new installations should use \pkg{qs2}.
#'
#' @param scope One of \code{"all"}, \code{"required"}, or
#'   \code{"recommended"}.
#' @param packages Optional character vector restricting installation to a
#'   subset of packages returned by \code{sn_list_dependencies()}.
#' @param missing_only Logical; when \code{TRUE} (default), install only missing
#'   packages.
#' @param repos CRAN-like repositories used for CRAN installs and bootstrap
#'   installation of helper installers such as \pkg{BiocManager} and
#'   \pkg{remotes}.
#' @param ask Passed to \code{BiocManager::install()}. Defaults to interactive
#'   behavior.
#' @param upgrade Logical; when \code{TRUE}, allow updating already installed
#'   GitHub and Bioconductor packages during installation.
#' @param github_dependencies Dependency policy passed to
#'   \code{remotes::install_github()} for GitHub-hosted packages. The default
#'   installs required dependencies without pulling optional suggested packages.
#' @param ... Additional arguments forwarded to the underlying installer calls.
#'
#' @return Invisibly returns the refreshed dependency table from
#'   \code{\link{sn_list_dependencies}()} after installation.
#'
#' @examples
#' \dontrun{
#' sn_install_dependencies(scope = "required")
#' sn_install_dependencies(packages = c("Seurat", "clusterProfiler"))
#' }
#' @export
sn_install_dependencies <- function(scope = c("all", "required", "recommended"),
                                    packages = NULL,
                                    missing_only = TRUE,
                                    repos = getOption("repos"),
                                    ask = interactive(),
                                    upgrade = FALSE,
                                    github_dependencies = NA,
                                    ...) {
  scope <- match.arg(scope)
  deps <- .sn_dependency_table()

  if (identical(scope, "required")) {
    deps <- deps[deps$requirement == "required", , drop = FALSE]
  } else if (identical(scope, "recommended")) {
    deps <- deps[deps$requirement == "recommended", , drop = FALSE]
  }

  if (!is.null(packages)) {
    unknown <- setdiff(packages, deps$package)
    if (length(unknown) > 0) {
      stop(
        "Unknown package(s): ", paste(unknown, collapse = ", "),
        ". Use `sn_list_dependencies()` to see supported names.",
        call. = FALSE
      )
    }
    deps <- deps[deps$package %in% packages, , drop = FALSE]
  }

  if (isTRUE(missing_only)) {
    deps <- deps[!deps$installed, , drop = FALSE]
  }

  if (nrow(deps) == 0) {
    .sn_log_info("No dependency installation is required for the selected scope.")
    return(invisible(sn_list_dependencies(scope = "all")))
  }

  cran_pkgs <- deps$package[deps$source == "CRAN"]
  bioc_pkgs <- deps$package[deps$source == "Bioconductor"]
  github_remotes <- deps$remote[deps$source == "GitHub"]

  .sn_install_cran_packages(packages = cran_pkgs, repos = repos, ...)
  .sn_install_bioc_packages(packages = bioc_pkgs, ask = ask, update = upgrade, repos = repos, ...)
  .sn_install_github_packages(
    remotes = github_remotes,
    upgrade = upgrade,
    repos = repos,
    dependencies = github_dependencies,
    ...
  )

  failed <- .sn_find_missing_packages(deps$package)
  if (length(failed) > 0) {
    stop(
      "Failed to install package(s): ", paste(failed, collapse = ", "),
      ". Review the installer output above for the first compilation or dependency error.",
      call. = FALSE
    )
  }

  invisible(sn_list_dependencies(scope = "all"))
}

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
#' sn_pixi_paths("scvi", runtime_dir = tempfile("shennong-home-"))
#'
#' @export
sn_pixi_paths <- function(environment = NULL,
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
    source_config_path = sn_pixi_config_path(environment),
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
#' sn_pixi_config_path("scvi")
#'
#' @export
sn_pixi_config_path <- function(environment = NULL) {
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
  paths <- sn_pixi_paths(environment = environment, runtime_dir = runtime_dir)
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

#' Run a Python analysis command through a managed Shennong pixi environment
#'
#' These are analysis-oriented wrappers around \code{sn_call_pixi_environment()}.
#' They prepare the corresponding package-bundled environment and run the
#' requested command. When \code{object} is supplied, method wrappers use a
#' Seurat object-level contract: export the object, run the packaged pixi
#' runner script, and import method outputs back into the object when the
#' backend produces cell-level metadata or embeddings.
#'
#' @param object Seurat object. Shennong writes the object to a Python
#'   interchange directory, runs the corresponding pixi script, and imports
#'   supported results.
#' @param reference_object Optional reference Seurat object for tools that map
#'   a query/spatial object against a single-cell reference, such as Tangram.
#' @param reference_assay,reference_layer Assay and layer used when exporting
#'   \code{reference_object}.
#' @param reference_signatures Optional file path or data frame of reference
#'   cell-state signatures for cell2location.
#' @param group_by Metadata column used by CellPhoneDB cell groups.
#' @param batch_by,label_by Metadata columns used by scArches/scPoli-style
#'   object workflows.
#' @param spatial_cols Two metadata columns containing spatial coordinates for
#'   spatial tools.
#' @param cell_type_by Reference metadata column containing cell-type labels
#'   for Tangram projection.
#' @param cluster_by Metadata column used by Squidpy neighborhood enrichment.
#' @param method_control Optional named list of backend-specific settings passed
#'   to the Python runner config.
#' @param assay Assay used for object-level infercnvpy input.
#' @param layer Assay layer used for object-level infercnvpy input. Defaults to
#'   \code{"data"} when present and otherwise \code{"counts"}.
#' @param species Species used to match bundled gene positions when
#'   \code{gene_order} and \code{gtf_file} are not supplied.
#' @param reference_by Metadata column containing normal/tumor annotations.
#' @param reference_cat One or more values in \code{reference_by} denoting
#'   normal reference cells.
#' @param gene_order Optional data frame with gene positions. It must contain a
#'   gene identifier column such as \code{feature}, \code{gene},
#'   \code{gene_name}, or \code{gene_id}, plus chromosome/start/end columns.
#' @param gtf_file Optional GTF file used by infercnvpy to annotate genomic
#'   positions instead of Shennong's bundled GENCODE table.
#' @param gtf_gene_id GTF attribute used by infercnvpy for matching.
#' @param adata_gene_id Optional AnnData var column used for matching a GTF.
#' @param output_dir Optional run directory. Defaults to
#'   \code{~/.shennong/runs/infercnvpy_*}.
#' @param runtime_dir Optional Shennong runtime directory.
#' @param key_added infercnvpy key used for the CNV representation.
#' @param window_size,step,dynamic_threshold,exclude_chromosomes,chunksize,n_jobs,calculate_gene_values,lfc_clip
#'   Parameters forwarded to \code{infercnvpy.tl.infercnv()}.
#' @param run_pca,run_neighbors,run_leiden,run_umap,score Logical flags for
#'   downstream infercnvpy analysis steps.
#' @param leiden_resolution Resolution passed to infercnvpy Leiden clustering.
#' @param cnv_score_group_by Optional grouping column for infercnvpy CNV scores.
#' @param metadata_prefix Prefix added to imported infercnvpy metadata columns.
#' @param result_name Name used under \code{object@misc$infercnvpy}.
#' @param return_object Whether to return the updated object. If \code{FALSE},
#'   return a run manifest list.
#' @param ... Additional arguments passed to \code{sn_call_pixi_environment()}.
#'
#' @return A Seurat object or a run manifest.
#'
#' @examples
#' \dontrun{
#' object <- sn_run_infercnvpy(
#'   object = object,
#'   reference_by = "cell_type",
#'   reference_cat = c("T cell", "Myeloid")
#' )
#' spatial <- sn_run_tangram(
#'   object = spatial,
#'   reference_object = reference,
#'   cell_type_by = "cell_type",
#'   spatial_cols = c("x", "y")
#' )
#' }
#'
#' @export
sn_run_scarches <- function(object,
                            assay = NULL,
                            layer = NULL,
                            batch_by = NULL,
                            label_by = NULL,
                            output_dir = NULL,
                            runtime_dir = NULL,
                            metadata_prefix = "scarches_",
                            result_name = "scarches",
                            return_object = TRUE,
                            method_control = list(),
                            ...) {
  .sn_run_python_object_method(
    object = object,
    environment = "scarches",
    script_name = "scarches_run.py",
    method = "scarches",
    assay = assay,
    layer = layer,
    output_dir = output_dir,
    runtime_dir = runtime_dir,
    metadata_prefix = metadata_prefix,
    result_name = result_name,
    return_object = return_object,
    config = c(list(batch_key = batch_by, labels_key = label_by), method_control),
    ...
  )
}

#' @rdname sn_run_scarches
#' @export
sn_run_scpoli <- function(object,
                          assay = NULL,
                          layer = NULL,
                          batch_by = NULL,
                          label_by = NULL,
                          output_dir = NULL,
                          runtime_dir = NULL,
                          metadata_prefix = "scpoli_",
                          result_name = "scpoli",
                          return_object = TRUE,
                          method_control = list(),
                          ...) {
  .sn_run_python_object_method(
    object = object,
    environment = "scarches",
    script_name = "scarches_run.py",
    method = "scpoli",
    assay = assay,
    layer = layer,
    output_dir = output_dir,
    runtime_dir = runtime_dir,
    metadata_prefix = metadata_prefix,
    result_name = result_name,
    return_object = return_object,
    config = c(list(batch_key = batch_by, labels_key = label_by), method_control),
    ...
  )
}

#' @rdname sn_run_scarches
#' @export
sn_run_infercnvpy <- function(object,
                              assay = NULL,
                              layer = NULL,
                              species = NULL,
                              reference_by = NULL,
                              reference_cat = NULL,
                              gene_order = NULL,
                              gtf_file = NULL,
                              gtf_gene_id = c("gene_name", "gene_id"),
                              adata_gene_id = NULL,
                              output_dir = NULL,
                              runtime_dir = NULL,
                              key_added = "cnv",
                              window_size = 100,
                              step = 10,
                              dynamic_threshold = 1.5,
                              exclude_chromosomes = c("chrX", "chrY"),
                              chunksize = 5000,
                              n_jobs = NULL,
                              calculate_gene_values = FALSE,
                              lfc_clip = 3,
                              run_pca = TRUE,
                              run_neighbors = TRUE,
                              run_leiden = TRUE,
                              run_umap = FALSE,
                              score = TRUE,
                              leiden_resolution = 1,
                              cnv_score_group_by = NULL,
                              metadata_prefix = "infercnvpy_",
                              result_name = "infercnvpy",
                              return_object = TRUE,
                              ...) {
  gtf_gene_id <- match.arg(gtf_gene_id)

  .sn_run_infercnvpy_object(
    object = object,
    assay = assay,
    layer = layer,
    species = species,
    reference_key = reference_by,
    reference_cat = reference_cat,
    gene_order = gene_order,
    gtf_file = gtf_file,
    gtf_gene_id = gtf_gene_id,
    adata_gene_id = adata_gene_id,
    output_dir = output_dir,
    runtime_dir = runtime_dir,
    key_added = key_added,
    window_size = window_size,
    step = step,
    dynamic_threshold = dynamic_threshold,
    exclude_chromosomes = exclude_chromosomes,
    chunksize = chunksize,
    n_jobs = n_jobs,
    calculate_gene_values = calculate_gene_values,
    lfc_clip = lfc_clip,
    run_pca = run_pca,
    run_neighbors = run_neighbors,
    run_leiden = run_leiden,
    run_umap = run_umap,
    score = score,
    leiden_resolution = leiden_resolution,
    cnv_score_group_by = cnv_score_group_by,
    metadata_prefix = metadata_prefix,
    result_name = result_name,
    return_object = return_object,
    ...
  )
}

#' @rdname sn_run_scarches
#' @export
sn_run_cellphonedb <- function(object,
                               assay = NULL,
                               layer = "counts",
                               group_by = NULL,
                               output_dir = NULL,
                               runtime_dir = NULL,
                               result_name = "cellphonedb",
                               return_object = TRUE,
                               method_control = list(),
                               ...) {
  if (is.null(group_by) || !nzchar(group_by)) {
    stop("`group_by` is required for `sn_run_cellphonedb()`.", call. = FALSE)
  }
  .sn_run_python_object_method(
    object = object,
    environment = "cellphonedb",
    script_name = "cellphonedb_run.py",
    method = "cellphonedb",
    assay = assay,
    layer = layer,
    output_dir = output_dir,
    runtime_dir = runtime_dir,
    metadata_prefix = "cellphonedb_",
    result_name = result_name,
    return_object = return_object,
    config = c(list(groupby = group_by), method_control),
    ...
  )
}

#' @rdname sn_run_scarches
#' @export
sn_run_cell2location <- function(object,
                                 assay = NULL,
                                 layer = "counts",
                                 reference_signatures = NULL,
                                 spatial_cols = NULL,
                                 output_dir = NULL,
                                 runtime_dir = NULL,
                                 metadata_prefix = "cell2location_",
                                 result_name = "cell2location",
                                 return_object = TRUE,
                                 method_control = list(),
                                 ...) {
  .sn_run_python_object_method(
    object = object,
    environment = "cell2location",
    script_name = "cell2location_run.py",
    method = "cell2location",
    assay = assay,
    layer = layer,
    spatial_cols = spatial_cols,
    output_dir = output_dir,
    runtime_dir = runtime_dir,
    metadata_prefix = metadata_prefix,
    result_name = result_name,
    return_object = return_object,
    config = c(list(reference_signatures = reference_signatures), method_control),
    ...
  )
}

#' @rdname sn_run_scarches
#' @export
sn_run_tangram <- function(object,
                           reference_object = NULL,
                           assay = NULL,
                           layer = NULL,
                           reference_assay = NULL,
                           reference_layer = NULL,
                           spatial_cols = NULL,
                           cell_type_by = NULL,
                           output_dir = NULL,
                           runtime_dir = NULL,
                           metadata_prefix = "tangram_",
                           result_name = "tangram",
                           return_object = TRUE,
                           method_control = list(),
                           ...) {
  if (is.null(reference_object)) {
    stop("`reference_object` is required for `sn_run_tangram()`.", call. = FALSE)
  }
  .sn_run_python_object_method(
    object = object,
    reference_object = reference_object,
    environment = "tangram",
    script_name = "tangram_run.py",
    method = "tangram",
    assay = assay,
    layer = layer,
    reference_assay = reference_assay,
    reference_layer = reference_layer,
    spatial_cols = spatial_cols,
    output_dir = output_dir,
    runtime_dir = runtime_dir,
    metadata_prefix = metadata_prefix,
    result_name = result_name,
    return_object = return_object,
    config = c(list(cell_type_key = cell_type_by), method_control),
    ...
  )
}

#' @rdname sn_run_scarches
#' @export
sn_run_squidpy <- function(object,
                           assay = NULL,
                           layer = NULL,
                           spatial_cols = NULL,
                           cluster_by = NULL,
                           output_dir = NULL,
                           runtime_dir = NULL,
                           metadata_prefix = "squidpy_",
                           result_name = "squidpy",
                           return_object = TRUE,
                           method_control = list(),
                           ...) {
  .sn_run_python_object_method(
    object = object,
    environment = "squidpy",
    script_name = "squidpy_run.py",
    method = "squidpy",
    assay = assay,
    layer = layer,
    spatial_cols = spatial_cols,
    output_dir = output_dir,
    runtime_dir = runtime_dir,
    metadata_prefix = metadata_prefix,
    result_name = result_name,
    return_object = return_object,
    config = c(list(cluster_key = cluster_by), method_control),
    ...
  )
}

#' @rdname sn_run_scarches
#' @export
sn_run_spatialdata <- function(object,
                               assay = NULL,
                               layer = NULL,
                               spatial_cols = NULL,
                               output_dir = NULL,
                               runtime_dir = NULL,
                               metadata_prefix = "spatialdata_",
                               result_name = "spatialdata",
                               return_object = TRUE,
                               method_control = list(),
                               ...) {
  .sn_run_python_object_method(
    object = object,
    environment = "spatialdata",
    script_name = "spatialdata_run.py",
    method = "spatialdata",
    assay = assay,
    layer = layer,
    spatial_cols = spatial_cols,
    output_dir = output_dir,
    runtime_dir = runtime_dir,
    metadata_prefix = metadata_prefix,
    result_name = result_name,
    return_object = return_object,
    config = method_control,
    ...
  )
}

#' @rdname sn_run_scarches
#' @export
sn_run_stlearn <- function(object,
                           assay = NULL,
                           layer = NULL,
                           spatial_cols = NULL,
                           output_dir = NULL,
                           runtime_dir = NULL,
                           metadata_prefix = "stlearn_",
                           result_name = "stlearn",
                           return_object = TRUE,
                           method_control = list(),
                           ...) {
  .sn_run_python_object_method(
    object = object,
    environment = "stlearn",
    script_name = "stlearn_run.py",
    method = "stlearn",
    assay = assay,
    layer = layer,
    spatial_cols = spatial_cols,
    output_dir = output_dir,
    runtime_dir = runtime_dir,
    metadata_prefix = metadata_prefix,
    result_name = result_name,
    return_object = return_object,
    config = method_control,
    ...
  )
}

.sn_run_python_object_method <- function(object,
                                         environment,
                                         script_name,
                                         method,
                                         assay = NULL,
                                         layer = NULL,
                                         reference_object = NULL,
                                         reference_assay = NULL,
                                         reference_layer = NULL,
                                         spatial_cols = NULL,
                                         output_dir = NULL,
                                         runtime_dir = NULL,
                                         metadata_prefix = paste0(method, "_"),
                                         result_name = method,
                                         return_object = TRUE,
                                         config = list(),
                                         ...) {
  check_installed(pkg = "Seurat", reason = glue::glue("to run {method} on a Seurat object."))
  assay <- assay %||% SeuratObject::DefaultAssay(object = object)
  layer <- layer %||% .sn_select_python_object_layer(object = object, assay = assay)
  output_dir <- output_dir %||% .sn_default_python_run_dir(method, runtime_dir = runtime_dir)
  input_dir <- file.path(output_dir, "input")
  result_dir <- file.path(output_dir, "output")
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  spatial_cols <- .sn_resolve_spatial_cols(object = object, spatial_cols = spatial_cols, required = environment %in% c("cell2location", "tangram", "squidpy", "spatialdata", "stlearn"))
  query <- .sn_write_python_object_input(
    object = object,
    input_dir = file.path(input_dir, "query"),
    assay = assay,
    layer = layer,
    spatial_cols = spatial_cols
  )

  reference <- NULL
  if (!is.null(reference_object)) {
    reference_assay <- reference_assay %||% SeuratObject::DefaultAssay(object = reference_object)
    reference_layer <- reference_layer %||% .sn_select_python_object_layer(object = reference_object, assay = reference_assay)
    reference <- .sn_write_python_object_input(
      object = reference_object,
      input_dir = file.path(input_dir, "reference"),
      assay = reference_assay,
      layer = reference_layer,
      spatial_cols = NULL
    )
  }

  config <- .sn_prepare_python_object_config(config = config, input_dir = input_dir)
  config <- c(
    list(
      method = method,
      assay = assay,
      layer = layer,
      reference_assay = reference_assay,
      reference_layer = reference_layer
    ),
    config
  )
  config_path <- .sn_write_json_file(config, file.path(output_dir, paste0(method, "_config.json")))
  script <- .sn_pixi_script_path(environment = environment, script_name = script_name)
  .sn_execute_python_object_pixi(
    environment = environment,
    script = script,
    input_dir = input_dir,
    output_dir = result_dir,
    config_path = config_path,
    ...
  )
  .sn_import_python_object_results(
    object = object,
    method = method,
    result_name = result_name,
    output_dir = result_dir,
    run_dir = output_dir,
    assay = assay,
    query = query,
    reference = reference,
    config = config,
    metadata_prefix = metadata_prefix,
    return_object = return_object
  )
}

.sn_select_python_object_layer <- function(object, assay) {
  layers <- SeuratObject::Layers(object[[assay]])
  if ("data" %in% layers || any(grepl("^data\\.", layers))) {
    return("data")
  }
  "counts"
}

.sn_pixi_script_path <- function(environment, script_name) {
  installed <- system.file("pixi", environment, "scripts", script_name, package = "Shennong")
  if (nzchar(installed) && file.exists(installed)) {
    return(normalizePath(installed, winslash = "/", mustWork = TRUE))
  }
  source_path <- file.path(getwd(), "inst", "pixi", environment, "scripts", script_name)
  if (file.exists(source_path)) {
    return(normalizePath(source_path, winslash = "/", mustWork = TRUE))
  }
  stop("Could not locate Shennong pixi runner script: ", file.path(environment, "scripts", script_name), call. = FALSE)
}

.sn_write_python_object_input <- function(object,
                                          input_dir,
                                          assay,
                                          layer,
                                          spatial_cols = NULL) {
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  expr <- .sn_get_seurat_layer_data(object = object, assay = assay, layer = layer)
  expr <- .sn_as_sparse_matrix(expr)
  expr <- expr[, colnames(object), drop = FALSE]

  matrix_path <- file.path(input_dir, "matrix.mtx")
  obs_path <- file.path(input_dir, "obs.csv")
  var_path <- file.path(input_dir, "var.csv")
  Matrix::writeMM(obj = expr, file = matrix_path)

  obs <- object[[]][colnames(expr), , drop = FALSE]
  obs <- data.frame(cell_id = rownames(obs), obs, check.names = FALSE)
  utils::write.csv(obs, file = obs_path, row.names = FALSE, quote = TRUE)
  var <- data.frame(feature_id = rownames(expr), stringsAsFactors = FALSE)
  utils::write.csv(var, file = var_path, row.names = FALSE, quote = TRUE)

  spatial_path <- NULL
  if (!is.null(spatial_cols)) {
    spatial <- obs[, spatial_cols, drop = FALSE]
    rownames(spatial) <- obs$cell_id
    spatial_path <- file.path(input_dir, "spatial.csv")
    utils::write.csv(spatial, file = spatial_path, row.names = TRUE, quote = TRUE)
  }

  list(
    input_dir = normalizePath(input_dir, winslash = "/", mustWork = TRUE),
    matrix_path = normalizePath(matrix_path, winslash = "/", mustWork = TRUE),
    obs_path = normalizePath(obs_path, winslash = "/", mustWork = TRUE),
    var_path = normalizePath(var_path, winslash = "/", mustWork = TRUE),
    spatial_path = if (!is.null(spatial_path)) normalizePath(spatial_path, winslash = "/", mustWork = TRUE) else NULL,
    assay = assay,
    layer = layer,
    features = rownames(expr),
    n_features = nrow(expr),
    n_cells = ncol(expr)
  )
}

.sn_resolve_spatial_cols <- function(object, spatial_cols = NULL, required = FALSE) {
  metadata_cols <- colnames(object[[]])
  if (!is.null(spatial_cols)) {
    if (length(spatial_cols) != 2L || !all(spatial_cols %in% metadata_cols)) {
      stop("`spatial_cols` must contain two metadata columns present in `object`.", call. = FALSE)
    }
    return(spatial_cols)
  }
  candidates <- list(
    c("x", "y"),
    c("spatial_x", "spatial_y"),
    c("imagecol", "imagerow"),
    c("pxl_col_in_fullres", "pxl_row_in_fullres"),
    c("array_col", "array_row")
  )
  for (candidate in candidates) {
    if (all(candidate %in% metadata_cols)) {
      return(candidate)
    }
  }
  if (isTRUE(required)) {
    stop("Spatial object workflow requires `spatial_cols` or recognized coordinate metadata columns.", call. = FALSE)
  }
  NULL
}

.sn_prepare_python_object_config <- function(config, input_dir) {
  if (!is.null(config$reference_signatures)) {
    signatures <- config$reference_signatures
    if (is.data.frame(signatures) || is.matrix(signatures)) {
      signatures_path <- file.path(input_dir, "reference_signatures.csv")
      utils::write.csv(as.data.frame(signatures), file = signatures_path, row.names = TRUE, quote = TRUE)
      config$reference_signatures <- normalizePath(signatures_path, winslash = "/", mustWork = TRUE)
    } else if (is.character(signatures) && length(signatures) == 1L && file.exists(path.expand(signatures))) {
      config$reference_signatures <- normalizePath(path.expand(signatures), winslash = "/", mustWork = TRUE)
    }
  }
  config
}

.sn_execute_python_object_pixi <- function(environment,
                                          script,
                                          input_dir,
                                          output_dir,
                                          config_path,
                                          ...) {
  sn_call_pixi_environment(
    environment = environment,
    command = "python",
    args = c(
      shQuote(script),
      "--input-dir", shQuote(input_dir),
      "--output-dir", shQuote(output_dir),
      "--config", shQuote(config_path)
    ),
    ...
  )
}

.sn_import_python_object_results <- function(object,
                                            method,
                                            result_name,
                                            output_dir,
                                            run_dir,
                                            assay,
                                            query,
                                            reference = NULL,
                                            config = list(),
                                            metadata_prefix = paste0(method, "_"),
                                            return_object = TRUE) {
  obs_path <- file.path(output_dir, "obs.csv")
  imported_metadata <- character(0)
  if (file.exists(obs_path)) {
    metadata <- utils::read.csv(obs_path, row.names = 1, check.names = FALSE)
    shared_cells <- intersect(colnames(object), rownames(metadata))
    if (length(shared_cells) > 0L && ncol(metadata) > 0L) {
      metadata <- metadata[shared_cells, , drop = FALSE]
      colnames(metadata) <- .sn_prefix_metadata_columns(colnames(metadata), metadata_prefix)
      imported_metadata <- colnames(metadata)
      object <- Seurat::AddMetaData(object = object, metadata = metadata)
    }
  }

  reductions <- character(0)
  embedding_files <- list.files(output_dir, pattern = "\\.csv$", full.names = TRUE)
  embedding_files <- embedding_files[grepl("(latent|pca|umap|embedding)", basename(embedding_files), ignore.case = TRUE)]
  for (embedding_file in embedding_files) {
    embedding <- tryCatch(.sn_read_embedding_csv(embedding_file, cells = colnames(object)), error = function(e) NULL)
    if (is.null(embedding)) {
      next
    }
    reduction_name <- paste0(metadata_prefix, tools::file_path_sans_ext(basename(embedding_file)))
    reduction_key <- paste0(gsub("[^A-Za-z0-9]", "", toupper(reduction_name)), "_")
    colnames(embedding) <- paste0(reduction_key, seq_len(ncol(embedding)))
    object[[reduction_name]] <- Seurat::CreateDimReducObject(
      embeddings = embedding,
      key = reduction_key,
      assay = assay
    )
    reductions <- c(reductions, reduction_name)
  }

  manifest_path <- file.path(output_dir, "manifest.json")
  backend_manifest <- if (file.exists(manifest_path)) {
    jsonlite::read_json(manifest_path, simplifyVector = TRUE)
  } else {
    list()
  }
  run_manifest <- c(
    list(
      method = method,
      result_name = result_name,
      assay = assay,
      source_layer = query$layer,
      run_dir = normalizePath(run_dir, winslash = "/", mustWork = TRUE),
      output_dir = normalizePath(output_dir, winslash = "/", mustWork = TRUE),
      input_features = query$features,
      n_features = query$n_features,
      n_cells = query$n_cells,
      reference = reference,
      imported_metadata = imported_metadata,
      imported_reductions = reductions,
      config = config
    ),
    backend_manifest
  )
  object@misc[[method]] <- object@misc[[method]] %||% list()
  object@misc[[method]][[result_name]] <- run_manifest
  object <- .sn_log_seurat_command(object = object, name = paste0("sn_run_", method))
  if (isTRUE(return_object)) {
    return(object)
  }
  run_manifest
}

.sn_run_infercnvpy_object <- function(object,
                                      assay = NULL,
                                      layer = NULL,
                                      species = NULL,
                                      reference_key = NULL,
                                      reference_cat = NULL,
                                      gene_order = NULL,
                                      gtf_file = NULL,
                                      gtf_gene_id = "gene_name",
                                      adata_gene_id = NULL,
                                      output_dir = NULL,
                                      runtime_dir = NULL,
                                      key_added = "cnv",
                                      window_size = 100,
                                      step = 10,
                                      dynamic_threshold = 1.5,
                                      exclude_chromosomes = c("chrX", "chrY"),
                                      chunksize = 5000,
                                      n_jobs = NULL,
                                      calculate_gene_values = FALSE,
                                      lfc_clip = 3,
                                      run_pca = TRUE,
                                      run_neighbors = TRUE,
                                      run_leiden = TRUE,
                                      run_umap = FALSE,
                                      score = TRUE,
                                      leiden_resolution = 1,
                                      cnv_score_group_by = NULL,
                                      metadata_prefix = "infercnvpy_",
                                      result_name = "infercnvpy",
                                      return_object = TRUE,
                                      ...) {
  check_installed(pkg = "Seurat", reason = "to run infercnvpy on a Seurat object.")
  assay <- assay %||% SeuratObject::DefaultAssay(object = object)
  layer <- layer %||% .sn_select_infercnvpy_layer(object = object, assay = assay)
  if (!is.null(reference_key) && (!nzchar(reference_key) || !reference_key %in% colnames(object[[]]))) {
    stop("`reference_key` must name a metadata column in `object`.", call. = FALSE)
  }

  output_dir <- output_dir %||% .sn_default_python_run_dir("infercnvpy", runtime_dir = runtime_dir)
  input_dir <- file.path(output_dir, "input")
  result_dir <- file.path(output_dir, "output")
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  input <- .sn_write_infercnvpy_input(
    object = object,
    input_dir = input_dir,
    assay = assay,
    layer = layer,
    species = species,
    gene_order = gene_order,
    gtf_file = gtf_file
  )
  script <- .sn_infercnvpy_script_path()
  config <- list(
    method = "infercnvpy",
    assay = assay,
    layer = NULL,
    source_layer = layer,
    reference_key = reference_key,
    reference_cat = reference_cat,
    gtf_file = if (!is.null(gtf_file)) normalizePath(path.expand(gtf_file), winslash = "/", mustWork = TRUE) else NULL,
    gtf_gene_id = gtf_gene_id,
    adata_gene_id = adata_gene_id,
    key_added = key_added,
    window_size = window_size,
    step = step,
    dynamic_threshold = dynamic_threshold,
    exclude_chromosomes = exclude_chromosomes,
    chunksize = chunksize,
    n_jobs = n_jobs,
    calculate_gene_values = isTRUE(calculate_gene_values),
    lfc_clip = lfc_clip,
    run_pca = isTRUE(run_pca),
    run_neighbors = isTRUE(run_neighbors),
    run_leiden = isTRUE(run_leiden),
    run_umap = isTRUE(run_umap),
    score = isTRUE(score),
    leiden_resolution = leiden_resolution,
    cnv_score_groupby = cnv_score_group_by,
    write_h5ad = TRUE
  )
  config_path <- .sn_write_json_file(config, file.path(output_dir, "infercnvpy_config.json"))

  .sn_execute_infercnvpy_pixi(
    script = script,
    input_dir = input$input_dir,
    output_dir = result_dir,
    config_path = config_path,
    ...
  )

  manifest <- .sn_import_infercnvpy_results(
    object = object,
    output_dir = result_dir,
    run_dir = output_dir,
    assay = assay,
    input = input,
    config = config,
    metadata_prefix = metadata_prefix,
    result_name = result_name,
    return_object = return_object
  )
  manifest
}

.sn_default_python_run_dir <- function(method, runtime_dir = NULL) {
  runtime_dir <- .sn_shennong_runtime_dir(runtime_dir)
  run_id <- paste0(
    method,
    "_",
    format(Sys.time(), "%Y%m%d_%H%M%S"),
    "_",
    Sys.getpid()
  )
  file.path(runtime_dir, "runs", run_id)
}

.sn_select_infercnvpy_layer <- function(object, assay) {
  layers <- SeuratObject::Layers(object[[assay]])
  if ("data" %in% layers || any(grepl("^data\\.", layers))) {
    return("data")
  }
  .sn_log_warn(
    "No normalized `data` layer was found for assay '{assay}'; using `counts`. ",
    "infercnvpy expects normalized log-transformed expression."
  )
  "counts"
}

.sn_infercnvpy_script_path <- function(script = NULL) {
  if (!is.null(script) && nzchar(script)) {
    script <- path.expand(script)
    if (!file.exists(script)) {
      stop("infercnvpy runner script does not exist: ", script, call. = FALSE)
    }
    return(normalizePath(script, winslash = "/", mustWork = TRUE))
  }

  installed <- system.file("pixi", "infercnvpy", "scripts", "infercnvpy_run.py", package = "Shennong")
  if (nzchar(installed) && file.exists(installed)) {
    return(normalizePath(installed, winslash = "/", mustWork = TRUE))
  }
  source_path <- file.path(getwd(), "inst", "pixi", "infercnvpy", "scripts", "infercnvpy_run.py")
  if (file.exists(source_path)) {
    return(normalizePath(source_path, winslash = "/", mustWork = TRUE))
  }
  stop("Could not locate Shennong's infercnvpy Python runner.", call. = FALSE)
}

.sn_write_infercnvpy_input <- function(object,
                                       input_dir,
                                       assay,
                                       layer,
                                       species = NULL,
                                       gene_order = NULL,
                                       gtf_file = NULL) {
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  expr <- .sn_get_seurat_layer_data(object = object, assay = assay, layer = layer)
  expr <- .sn_as_sparse_matrix(expr)
  expr <- expr[, colnames(object), drop = FALSE]

  gene_positions <- .sn_prepare_infercnvpy_gene_positions(
    features = rownames(expr),
    object = object,
    species = species,
    gene_order = gene_order,
    gtf_file = gtf_file
  )
  keep_features <- gene_positions$feature_id
  expr <- expr[keep_features, , drop = FALSE]
  gene_positions <- gene_positions[match(rownames(expr), gene_positions$feature_id), , drop = FALSE]

  matrix_path <- file.path(input_dir, "matrix.mtx")
  obs_path <- file.path(input_dir, "obs.csv")
  var_path <- file.path(input_dir, "var.csv")
  Matrix::writeMM(obj = expr, file = matrix_path)

  obs <- object[[]][colnames(expr), , drop = FALSE]
  obs <- data.frame(cell_id = rownames(obs), obs, check.names = FALSE)
  utils::write.csv(obs, file = obs_path, row.names = FALSE, quote = TRUE)

  utils::write.csv(gene_positions, file = var_path, row.names = FALSE, quote = TRUE)
  list(
    input_dir = normalizePath(input_dir, winslash = "/", mustWork = TRUE),
    matrix_path = normalizePath(matrix_path, winslash = "/", mustWork = TRUE),
    obs_path = normalizePath(obs_path, winslash = "/", mustWork = TRUE),
    var_path = normalizePath(var_path, winslash = "/", mustWork = TRUE),
    assay = assay,
    layer = layer,
    features = rownames(expr),
    n_features = nrow(expr),
    n_cells = ncol(expr)
  )
}

.sn_prepare_infercnvpy_gene_positions <- function(features,
                                                  object,
                                                  species = NULL,
                                                  gene_order = NULL,
                                                  gtf_file = NULL) {
  features <- as.character(features)
  if (!is.null(gene_order)) {
    positions <- .sn_match_user_gene_order(features = features, gene_order = gene_order)
  } else if (!is.null(gtf_file)) {
    if (!file.exists(path.expand(gtf_file))) {
      stop("`gtf_file` does not exist: ", gtf_file, call. = FALSE)
    }
    positions <- data.frame(feature_id = features, stringsAsFactors = FALSE)
  } else {
    species <- sn_get_species(object = object, species = species)
    annotations <- .sn_get_gene_annotation_table(species = species)
    positions <- .sn_match_annotation_gene_order(features = features, annotations = annotations)
  }

  if (!all(c("feature_id") %in% colnames(positions))) {
    stop("Internal error: infercnvpy gene positions must include `feature_id`.", call. = FALSE)
  }
  if (is.null(gtf_file) && !all(c("chromosome", "start", "end") %in% colnames(positions))) {
    stop("Gene positions must include `chromosome`, `start`, and `end`.", call. = FALSE)
  }
  positions <- positions[positions$feature_id %in% features, , drop = FALSE]
  positions <- positions[!duplicated(positions$feature_id), , drop = FALSE]
  positions <- positions[match(intersect(features, positions$feature_id), positions$feature_id), , drop = FALSE]
  if (nrow(positions) == 0L) {
    stop("No input genes could be matched to genomic positions for infercnvpy.", call. = FALSE)
  }
  missing <- setdiff(features, positions$feature_id)
  if (length(missing) > 0L) {
    .sn_log_warn("Dropping {length(missing)} gene(s) without genomic positions before infercnvpy.")
  }
  positions
}

.sn_match_annotation_gene_order <- function(features, annotations) {
  feature_base <- sub("\\..*$", "", features)
  by_name <- annotations[!is.na(annotations$gene_name) & nzchar(annotations$gene_name), , drop = FALSE]
  by_name <- by_name[!duplicated(by_name$gene_name), , drop = FALSE]
  rownames(by_name) <- by_name$gene_name
  by_id <- annotations[!is.na(annotations$gene_id) & nzchar(annotations$gene_id), , drop = FALSE]
  by_id <- by_id[!duplicated(by_id$gene_id), , drop = FALSE]
  rownames(by_id) <- by_id$gene_id
  by_base <- annotations[!is.na(annotations$gene_id_base) & nzchar(annotations$gene_id_base), , drop = FALSE]
  by_base <- by_base[!duplicated(by_base$gene_id_base), , drop = FALSE]
  rownames(by_base) <- by_base$gene_id_base

  matched <- vector("list", length(features))
  for (i in seq_along(features)) {
    feature <- features[[i]]
    row <- NULL
    matched_by <- NA_character_
    if (feature %in% rownames(by_name)) {
      row <- by_name[feature, , drop = FALSE]
      matched_by <- "gene_name"
    } else if (feature %in% rownames(by_id)) {
      row <- by_id[feature, , drop = FALSE]
      matched_by <- "gene_id"
    } else if (feature_base[[i]] %in% rownames(by_base)) {
      row <- by_base[feature_base[[i]], , drop = FALSE]
      matched_by <- "gene_id_base"
    }
    if (!is.null(row)) {
      matched[[i]] <- data.frame(
        feature_id = feature,
        chromosome = row$seqname,
        start = row$start,
        end = row$end,
        gene_id = row$gene_id,
        gene_id_base = row$gene_id_base,
        gene_name = row$gene_name,
        matched_by = matched_by,
        stringsAsFactors = FALSE
      )
    }
  }
  matched <- matched[!vapply(matched, is.null, logical(1))]
  if (length(matched) == 0L) {
    return(data.frame(feature_id = character(0), stringsAsFactors = FALSE))
  }
  do.call(rbind, matched)
}

.sn_match_user_gene_order <- function(features, gene_order) {
  gene_order <- as.data.frame(gene_order, stringsAsFactors = FALSE)
  id_col <- intersect(c("feature_id", "feature", "gene", "gene_name", "gene_id"), colnames(gene_order))
  if (length(id_col) == 0L) {
    stop(
      "`gene_order` must contain one of: feature_id, feature, gene, gene_name, gene_id.",
      call. = FALSE
    )
  }
  id_col <- id_col[[1]]
  chr_col <- intersect(c("chromosome", "seqname", "chr"), colnames(gene_order))
  if (length(chr_col) == 0L || !"start" %in% colnames(gene_order) || !"end" %in% colnames(gene_order)) {
    stop("`gene_order` must contain chromosome/seqname, start, and end columns.", call. = FALSE)
  }
  chr_col <- chr_col[[1]]
  gene_order <- gene_order[!duplicated(gene_order[[id_col]]), , drop = FALSE]
  rownames(gene_order) <- as.character(gene_order[[id_col]])
  keep <- features[features %in% rownames(gene_order)]
  out <- gene_order[keep, , drop = FALSE]
  data.frame(
    feature_id = keep,
    chromosome = out[[chr_col]],
    start = out[["start"]],
    end = out[["end"]],
    stringsAsFactors = FALSE
  )
}

.sn_execute_infercnvpy_pixi <- function(script,
                                        input_dir,
                                        output_dir,
                                        config_path,
                                        ...) {
  sn_call_infercnvpy(
    command = "python",
    args = c(
      shQuote(script),
      "--input-dir", shQuote(input_dir),
      "--output-dir", shQuote(output_dir),
      "--config", shQuote(config_path)
    ),
    ...
  )
}

.sn_import_infercnvpy_results <- function(object,
                                          output_dir,
                                          run_dir,
                                          assay,
                                          input,
                                          config,
                                          metadata_prefix = "infercnvpy_",
                                          result_name = "infercnvpy",
                                          return_object = TRUE) {
  obs_path <- file.path(output_dir, "obs.csv")
  manifest_path <- file.path(output_dir, "manifest.json")
  if (!file.exists(obs_path)) {
    stop("infercnvpy output is missing `obs.csv`: ", obs_path, call. = FALSE)
  }

  metadata <- utils::read.csv(obs_path, row.names = 1, check.names = FALSE)
  shared_cells <- intersect(colnames(object), rownames(metadata))
  if (length(shared_cells) == 0L) {
    stop("infercnvpy metadata output does not contain cells from the input object.", call. = FALSE)
  }
  metadata <- metadata[shared_cells, , drop = FALSE]
  colnames(metadata) <- .sn_prefix_metadata_columns(colnames(metadata), metadata_prefix)
  object <- Seurat::AddMetaData(object = object, metadata = metadata)

  pca_path <- file.path(output_dir, "cnv_pca.csv")
  umap_path <- file.path(output_dir, "cnv_umap.csv")
  reductions <- character(0)
  if (file.exists(pca_path)) {
    pca <- .sn_read_embedding_csv(pca_path, cells = colnames(object))
    reduction_name <- paste0(metadata_prefix, "cnv_pca")
    colnames(pca) <- paste0("CNVPCA_", seq_len(ncol(pca)))
    object[[reduction_name]] <- Seurat::CreateDimReducObject(
      embeddings = pca,
      key = "CNVPCA_",
      assay = assay
    )
    reductions <- c(reductions, reduction_name)
  }
  if (file.exists(umap_path)) {
    umap <- .sn_read_embedding_csv(umap_path, cells = colnames(object))
    reduction_name <- paste0(metadata_prefix, "cnv_umap")
    colnames(umap) <- paste0("CNVUMAP_", seq_len(ncol(umap)))
    object[[reduction_name]] <- Seurat::CreateDimReducObject(
      embeddings = umap,
      key = "CNVUMAP_",
      assay = assay
    )
    reductions <- c(reductions, reduction_name)
  }

  backend_manifest <- if (file.exists(manifest_path)) {
    jsonlite::read_json(manifest_path, simplifyVector = TRUE)
  } else {
    list()
  }
  run_manifest <- c(
    list(
      method = "infercnvpy",
      result_name = result_name,
      assay = assay,
      source_layer = input$layer,
      run_dir = normalizePath(run_dir, winslash = "/", mustWork = TRUE),
      output_dir = normalizePath(output_dir, winslash = "/", mustWork = TRUE),
      input_features = input$features,
      n_features = input$n_features,
      n_cells = input$n_cells,
      imported_metadata = colnames(metadata),
      imported_reductions = reductions,
      config = config
    ),
    backend_manifest
  )

  object@misc$infercnvpy <- object@misc$infercnvpy %||% list()
  object@misc$infercnvpy[[result_name]] <- run_manifest
  object <- .sn_log_seurat_command(object = object, name = "sn_run_infercnvpy")
  if (isTRUE(return_object)) {
    return(object)
  }
  run_manifest
}

.sn_prefix_metadata_columns <- function(columns, prefix) {
  if (is.null(prefix) || !nzchar(prefix)) {
    return(columns)
  }
  ifelse(startsWith(columns, prefix), columns, paste0(prefix, columns))
}

.sn_read_embedding_csv <- function(path, cells) {
  embedding <- utils::read.csv(path, row.names = 1, check.names = FALSE)
  missing_cells <- setdiff(cells, rownames(embedding))
  if (length(missing_cells) > 0L) {
    stop("Embedding output is missing cells from the input object: ", path, call. = FALSE)
  }
  embedding <- as.matrix(embedding[cells, , drop = FALSE])
  storage.mode(embedding) <- "numeric"
  embedding
}

#' Detect local accelerator support for pixi-managed Python methods
#'
#' This helper performs lightweight command-line checks for GPUs. CUDA/NVIDIA is
#' currently the only accelerator profile Shennong uses to select a pixi GPU
#' environment automatically; other detected accelerators are reported for the
#' user's information and fall back to CPU unless explicitly handled later.
#'
#' @param quiet Logical; suppress status messages.
#'
#' @return A named list with \code{has_gpu}, \code{backend}, \code{devices},
#'   and detected CUDA version when available.
#'
#' @examples
#' accel <- sn_detect_accelerator(quiet = TRUE)
#' accel$backend
#'
#' @export
sn_detect_accelerator <- function(quiet = FALSE) {
  nvidia <- .sn_detect_nvidia_gpu()
  if (isTRUE(nvidia$available)) {
    result <- list(
      has_gpu = TRUE,
      backend = "cuda",
      devices = nvidia$devices,
      cuda_version = nvidia$cuda_version,
      raw = nvidia$raw
    )
  } else if (.sn_command_available("rocminfo") || .sn_command_available("rocm-smi")) {
    result <- list(
      has_gpu = TRUE,
      backend = "rocm",
      devices = character(0),
      cuda_version = NA_character_,
      raw = character(0)
    )
  } else {
    result <- list(
      has_gpu = FALSE,
      backend = "cpu",
      devices = character(0),
      cuda_version = NA_character_,
      raw = character(0)
    )
  }

  if (!quiet) {
    .sn_log_info(
      "Accelerator detection: backend = {result$backend}; ",
      "has_gpu = {result$has_gpu}; cuda = {result$cuda_version %||% NA_character_}."
    )
  }

  invisible(result)
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

#' Check whether Shennong is up to date
#'
#' This helper compares the installed package version against the latest
#' available CRAN or GitHub development version. When `channel = "auto"`, it
#' prefers CRAN if a release exists there and otherwise falls back to GitHub.
#'
#' @param channel One of \code{"auto"}, \code{"cran"}, or \code{"github"}.
#' @param package Package name to check. Defaults to \code{"Shennong"}.
#' @param source GitHub repository in \code{"owner/repo"} format when checking
#'   the development channel.
#' @param ref GitHub ref to inspect. Defaults to \code{"main"}.
#' @param repos CRAN-like repositories used for version lookup.
#' @param quiet Logical; if \code{TRUE}, suppress the summary message.
#'
#' @return A named list with the installed version, remote version, selected
#'   channel, whether the package is up to date, and the recommended install
#'   command.
#'
#' @examples
#' \dontrun{
#' sn_check_version()
#' sn_check_version(channel = "github")
#' }
#'
#' @export
sn_check_version <- function(
  channel = c("auto", "cran", "github"),
  package = "Shennong",
  source = "zerostwo/shennong",
  ref = "main",
  repos = getOption("repos"),
  quiet = FALSE
) {
  channel <- match.arg(channel)
  installed_version <- .sn_get_installed_version(package = package)
  cran_version <- .sn_get_cran_version(package = package, repos = repos)
  github_version <- .sn_get_github_version(repo = source, ref = ref)
  resolved_channel <- .sn_resolve_release_channel(
    channel = channel,
    cran_version = cran_version,
    github_version = github_version
  )
  remote_version <- switch(
    resolved_channel,
    cran = cran_version,
    github = github_version
  )
  status <- .sn_compare_version_status(
    installed_version = installed_version,
    remote_version = remote_version
  )
  install_command <- switch(
    resolved_channel,
    cran = sprintf('install.packages("%s")', package),
    github = sprintf('remotes::install_github("%s", ref = "%s")', source, ref)
  )

  result <- list(
    package = package,
    channel = resolved_channel,
    installed_version = installed_version,
    remote_version = remote_version,
    remote_available = !is.null(remote_version),
    up_to_date = status$up_to_date,
    status = status$status,
    install_command = install_command
  )

  if (!quiet) {
    installed_label <- if (is.null(installed_version)) "not installed" else as.character(installed_version)
    remote_label <- if (is.null(remote_version)) "unavailable" else as.character(remote_version)
    .sn_log_info(
      "Shennong ({resolved_channel}): installed = {installed_label}; ",
      "latest = {remote_label}; status = {status$status}."
    )
    if (!isTRUE(status$up_to_date)) {
      .sn_log_info("Install or update with: {install_command}.")
    }
  }

  invisible(result)
}

#' Install Shennong from CRAN, GitHub, or a local source
#'
#' This helper installs the stable CRAN release when available, or the GitHub
#' development version when requested, or installs from a local source tree or
#' tarball. When `channel = "auto"`, it prefers CRAN and falls back to GitHub
#' if no CRAN release is available.
#'
#' @param channel One of \code{"auto"}, \code{"cran"}, \code{"github"}, or
#'   \code{"local"}.
#' @param package Package name. Defaults to \code{"Shennong"}.
#' @param source Installation source. For \code{channel = "github"}, supply an
#'   \code{"owner/repo"} string; when omitted, Shennong uses
#'   \code{"zerostwo/shennong"}. For \code{channel = "local"}, supply a local
#'   package directory or source tarball path.
#' @param ref GitHub ref used for \code{channel = "github"}. Defaults to
#'   \code{"main"}.
#' @param repos CRAN-like repositories used by \code{install.packages()}.
#' @param ... Additional arguments passed to \code{utils::install.packages()} or
#'   \code{remotes::install_github()} / \code{remotes::install_local()}. For
#'   GitHub installs, Shennong defaults to \code{dependencies = FALSE} and
#'   \code{upgrade = "never"} unless you override them explicitly.
#'
#' @return Invisibly returns the chosen installation channel.
#'
#' @examples
#' \dontrun{
#' sn_install_shennong(channel = "github")
#' sn_install_shennong(channel = "github", source = "zerostwo/shennong", ref = "main")
#' sn_install_shennong(channel = "local", source = "~/personal/packages/shennong")
#' }
#'
#' @export
sn_install_shennong <- function(
  channel = c("auto", "cran", "github", "local"),
  package = "Shennong",
  source = NULL,
  ref = "main",
  repos = getOption("repos"),
  ...
) {
  channel <- match.arg(channel)

  if (identical(channel, "local")) {
    if (is.null(source) || !nzchar(source)) {
      stop("`source` must be supplied when `channel = \"local\"`.", call. = FALSE)
    }
    check_installed("remotes", reason = "to install Shennong from a local path.")
    .sn_install_local_release(
      path = source,
      args = list(...)
    )
    return(invisible(channel))
  }

  source <- source %||% "zerostwo/shennong"
  cran_version <- .sn_get_cran_version(package = package, repos = repos)
  github_version <- .sn_get_github_version(repo = source, ref = ref)
  resolved_channel <- .sn_resolve_release_channel(
    channel = channel,
    cran_version = cran_version,
    github_version = github_version
  )

  if (resolved_channel == "cran") {
    utils::install.packages(package, repos = repos, ...)
    return(invisible(resolved_channel))
  }

  check_installed("remotes", reason = "to install the GitHub development version.")
  github_args <- list(...)
  if (is.null(github_args$dependencies)) {
    github_args$dependencies <- FALSE
  }
  if (is.null(github_args$upgrade)) {
    github_args$upgrade <- "never"
  }

  .sn_install_github_release(
    repo = source,
    ref = ref,
    args = github_args
  )
  invisible(resolved_channel)
}

.sn_install_github_release <- function(repo, ref, args = list()) {
  do.call(
    remotes::install_github,
    c(list(repo = repo, ref = ref), args)
  )
}

.sn_install_local_release <- function(path, args = list()) {
  do.call(
    remotes::install_local,
    c(list(path = path), args)
  )
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
  environment %in% c("scvi", "scarches", "cell2location", "tangram")
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
  template <- readLines(sn_pixi_config_path(environment), warn = FALSE)
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

.sn_command_available <- function(command) {
  path <- Sys.which(command)[[command]]
  !is.na(path) && nzchar(path)
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

.sn_description_root <- function() {
  namespace_path <- .sn_namespace_path()
  if (nzchar(namespace_path) && file.exists(file.path(namespace_path, "DESCRIPTION"))) {
    return(namespace_path)
  }

  system_path <- system.file(package = "Shennong")
  if (nzchar(system_path) && file.exists(file.path(system_path, "DESCRIPTION"))) {
    return(system_path)
  }

  getwd()
}

.sn_read_package_description <- function() {
  desc_path <- file.path(.sn_description_root(), "DESCRIPTION")
  if (!file.exists(desc_path)) {
    stop("Could not locate the package DESCRIPTION file.", call. = FALSE)
  }

  as.list(read.dcf(desc_path, all = TRUE)[1, , drop = FALSE])
}

.sn_split_description_packages <- function(field_value) {
  if (is.null(field_value) || !nzchar(field_value)) {
    return(character(0))
  }

  entries <- trimws(unlist(strsplit(field_value, ",", fixed = TRUE), use.names = FALSE))
  entries <- entries[nzchar(entries)]
  entries <- sub("\\s*\\(.*\\)$", "", entries)
  entries <- trimws(entries)
  entries[entries != "R"]
}

.sn_dependency_source_overrides <- function() {
  data.frame(
    package = c(
      "anndataR", "BayesPrism", "BPCells", "COSG", "CellChat", "GapClust",
      "ROGUE", "SignatuR", "catplot", "harmony", "liana", "lisi",
      "nichenetr", "tidytemplate"
    ),
    source = rep("GitHub", 14),
    remote = c(
      "scverse/anndataR",
      "Danko-Lab/BayesPrism/BayesPrism",
      "bnprks/BPCells/r",
      "genecell/COSGR",
      "jinworks/CellChat",
      "fabotao/GapClust",
      "PaulingLiu/ROGUE",
      "carmonalab/SignatuR",
      "catplot/catplot",
      "immunogenomics/harmony@harmony2",
      "saezlab/liana",
      "immunogenomics/lisi",
      "saeyslab/nichenetr",
      "tidyverse/tidytemplate"
    ),
    stringsAsFactors = FALSE
  )
}

.sn_bioconductor_packages <- function() {
  c(
    "BiocParallel", "celda", "clusterProfiler", "Coralysis", "decoupleR", "DESeq2", "dorothea", "edgeR",
    "glmGamPoi", "miloR", "org.Hs.eg.db", "org.Mm.eg.db", "progeny", "rhdf5",
    "rtracklayer", "S4Vectors", "scDblFinder", "scDesign3", "scran", "SingleCellExperiment", "Nebulosa",
    "SummarizedExperiment", "limma"
  )
}

.sn_dependency_table <- function() {
  desc <- .sn_read_package_description()
  required <- .sn_split_description_packages(desc$Imports)
  recommended <- .sn_split_description_packages(desc$Suggests)
  pkg_names <- c(required, recommended)
  declared_in <- c(rep("Imports", length(required)), rep("Suggests", length(recommended)))
  requirement <- c(rep("required", length(required)), rep("recommended", length(recommended)))

  deps <- data.frame(
    package = pkg_names,
    requirement = requirement,
    declared_in = declared_in,
    stringsAsFactors = FALSE
  )
  deps <- deps[!duplicated(deps$package), , drop = FALSE]

  overrides <- .sn_dependency_source_overrides()
  deps$source <- ifelse(deps$package %in% .sn_bioconductor_packages(), "Bioconductor", "CRAN")
  deps$remote <- NA_character_

  matched_override <- match(deps$package, overrides$package)
  has_override <- !is.na(matched_override)
  deps$source[has_override] <- overrides$source[matched_override[has_override]]
  deps$remote[has_override] <- overrides$remote[matched_override[has_override]]

  deps$installed <- vapply(
    deps$package,
    function(package) suppressWarnings(rlang::is_installed(package)),
    logical(1)
  )
  deps$version <- vapply(
    seq_len(nrow(deps)),
    function(i) {
      if (!deps$installed[[i]]) {
        return(NA_character_)
      }
      as.character(utils::packageVersion(deps$package[[i]]))
    },
    character(1)
  )

  deps[order(deps$requirement, deps$package), , drop = FALSE]
}

.sn_install_cran_packages <- function(packages, repos = getOption("repos"), ...) {
  packages <- unique(stats::na.omit(packages))
  if (length(packages) == 0) {
    return(invisible(NULL))
  }

  utils::install.packages(packages, repos = repos, ...)
  invisible(packages)
}

.sn_install_bioc_packages <- function(packages,
                                      ask = interactive(),
                                      update = FALSE,
                                      repos = getOption("repos"),
                                      ...) {
  packages <- unique(stats::na.omit(packages))
  if (length(packages) == 0) {
    return(invisible(NULL))
  }

  if (!rlang::is_installed("BiocManager")) {
    utils::install.packages("BiocManager", repos = repos)
  }

  BiocManager::install(packages, ask = ask, update = update, ...)
  invisible(packages)
}

.sn_install_github_packages <- function(remotes,
                                        upgrade = FALSE,
                                        repos = getOption("repos"),
                                        dependencies = NA,
                                        ...) {
  remotes <- unique(stats::na.omit(remotes))
  if (length(remotes) == 0) {
    return(invisible(NULL))
  }

  if (!rlang::is_installed("remotes")) {
    utils::install.packages("remotes", repos = repos)
  }

  upgrade_policy <- if (isTRUE(upgrade)) "always" else "never"
  for (remote in remotes) {
    remotes::install_github(remote, upgrade = upgrade_policy, dependencies = dependencies, ...)
  }

  invisible(remotes)
}

.sn_find_missing_packages <- function(packages) {
  packages <- unique(stats::na.omit(packages))
  packages[!vapply(packages, rlang::is_installed, logical(1))]
}

.sn_get_installed_version <- function(package = "Shennong") {
  if (!rlang::is_installed(package)) {
    return(NULL)
  }

  utils::packageVersion(package)
}

.sn_get_cran_version <- function(package = "Shennong", repos = getOption("repos")) {
  available <- tryCatch(
    utils::available.packages(repos = repos),
    error = function(e) NULL
  )

  if (is.null(available) || !package %in% rownames(available)) {
    return(NULL)
  }

  package_version(available[package, "Version"])
}

.sn_get_github_version <- function(repo = "zerostwo/shennong", ref = "main") {
  if (length(strsplit(repo, "/", fixed = TRUE)[[1]]) != 2) {
    stop("`source` must use the form 'owner/repo'.", call. = FALSE)
  }

  url <- sprintf(
    "https://raw.githubusercontent.com/%s/%s/DESCRIPTION",
    repo,
    ref
  )

  desc_text <- tryCatch(
    rawToChar(curl::curl_fetch_memory(
      url,
      handle = curl::new_handle(timeout = 20)
    )$content),
    error = function(e) NULL
  )

  if (is.null(desc_text)) {
    return(NULL)
  }

  con <- textConnection(desc_text)
  on.exit(close(con), add = TRUE)
  parsed <- tryCatch(
    read.dcf(con),
    error = function(e) NULL
  )

  if (is.null(parsed) || !"Version" %in% colnames(parsed)) {
    return(NULL)
  }

  package_version(parsed[1, "Version"])
}

.sn_resolve_release_channel <- function(channel = c("auto", "cran", "github"),
                                        cran_version = NULL,
                                        github_version = NULL) {
  channel <- match.arg(channel)

  if (channel == "auto") {
    if (!is.null(cran_version)) {
      return("cran")
    }
    if (!is.null(github_version)) {
      return("github")
    }
    stop("Could not determine a remote version from CRAN or GitHub.", call. = FALSE)
  }

  if (channel == "cran" && is.null(cran_version)) {
    stop("Shennong is not currently available on CRAN.", call. = FALSE)
  }

  if (channel == "github" && is.null(github_version)) {
    stop("Could not retrieve the GitHub development version.", call. = FALSE)
  }

  channel
}

.sn_namespace_path <- function() {
  tryCatch(getNamespaceInfo(asNamespace("Shennong"), "path"), error = function(e) "")
}

.sn_template_path <- function(relative_path) {
  installed_path <- system.file(file.path("templates", relative_path), package = "Shennong")
  if (nzchar(installed_path) && file.exists(installed_path)) {
    return(installed_path)
  }

  ns_path <- .sn_namespace_path()
  candidates <- c(
    file.path(ns_path, "inst", "templates", relative_path),
    file.path(ns_path, "templates", relative_path)
  )
  existing <- candidates[file.exists(candidates)]
  if (length(existing) > 0) {
    return(existing[[1]])
  }

  stop(glue("Template file '{relative_path}' could not be found in the installed package."), call. = FALSE)
}

.sn_codex_component_path <- function(
  component = c("codex_root", "package_skills", "project_template", "project_template_skills")
) {
  component <- match.arg(component)
  relative_path <- switch(
    component,
    codex_root = "codex",
    package_skills = file.path("codex", "package-skills"),
    project_template = file.path("codex", "project-template"),
    project_template_skills = file.path("codex", "project-template", "skills")
  )

  installed_path <- system.file(relative_path, package = "Shennong")
  if (nzchar(installed_path) && dir.exists(installed_path)) {
    return(installed_path)
  }

  ns_path <- .sn_namespace_path()
  candidates <- c(
    file.path(ns_path, "inst", relative_path),
    file.path(ns_path, relative_path)
  )
  existing <- candidates[dir.exists(candidates)]
  if (length(existing) > 0) {
    return(existing[[1]])
  }

  stop(glue("Bundled Codex asset '{component}' could not be found in the installed package."), call. = FALSE)
}

.sn_render_template <- function(relative_path, context = list()) {
  template_path <- .sn_template_path(relative_path)
  lines <- readLines(template_path, warn = FALSE)
  if (length(lines) == 0) {
    return(character(0))
  }

  rendered <- lines
  for (name in names(context)) {
    rendered <- gsub(
      pattern = paste0("\\{\\{", name, "\\}\\}"),
      replacement = context[[name]],
      x = rendered,
      fixed = FALSE
    )
  }

  rendered
}

.sn_render_text_file <- function(file_path, context = list()) {
  lines <- readLines(file_path, warn = FALSE)
  if (length(lines) == 0) {
    return(character(0))
  }

  rendered <- lines
  for (name in names(context)) {
    rendered <- gsub(
      pattern = paste0("\\{\\{", name, "\\}\\}"),
      replacement = context[[name]],
      x = rendered,
      fixed = FALSE
    )
  }

  rendered
}

.sn_skill_source_directories <- function(component = c("package_skills", "project_template_skills")) {
  component <- match.arg(component)
  source_root <- .sn_codex_component_path(component = component)
  source_dirs <- list.dirs(source_root, recursive = FALSE, full.names = TRUE)
  source_dirs[dir.exists(source_dirs)]
}

.sn_project_template_text_files <- function() {
  c(
    "AGENTS.md",
    "README.md",
    "directories.txt",
    "gitignore.template",
    "project.Rproj.template",
    file.path("memory", "Decisions.md"),
    file.path("memory", "Plan.md"),
    file.path("memory", "Prompt.md"),
    file.path("memory", "Status.md"),
    file.path("docs", "standards", "BioinformaticsAnalysisConventions.md"),
    file.path("config", "default.yaml")
  )
}

.sn_project_template_directory_manifest <- function() {
  "directories.txt"
}

.sn_project_template_directories <- function() {
  manifest_path <- file.path(
    .sn_codex_component_path("project_template"),
    .sn_project_template_directory_manifest()
  )

  if (!file.exists(manifest_path)) {
    return(character(0))
  }

  dirs <- readLines(manifest_path, warn = FALSE)
  dirs <- trimws(dirs)
  unique(dirs[nzchar(dirs)])
}

.sn_create_directories <- function(root, relative_paths) {
  if (length(relative_paths) == 0) {
    return(invisible(character(0)))
  }

  created_paths <- file.path(root, relative_paths)
  for (dir_path in created_paths) {
    dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
  }

  invisible(created_paths)
}

.sn_project_template_skip_patterns <- function(include_governance = TRUE) {
  if (isTRUE(include_governance)) {
    return(character(0))
  }

  c(
    "^AGENTS\\.md$",
    "^memory/",
    "^docs/standards/",
    "^skills/"
  )
}

.sn_should_skip_template_path <- function(relative_path, include_governance = TRUE) {
  patterns <- .sn_project_template_skip_patterns(include_governance = include_governance)
  if (length(patterns) == 0) {
    return(FALSE)
  }

  any(vapply(patterns, function(pattern) grepl(pattern, relative_path), logical(1)))
}

.sn_project_rproj_filename <- function(project_dir) {
  paste0(basename(project_dir), ".Rproj")
}

.sn_project_template_destination_path <- function(relative_path, project_dir) {
  if (identical(relative_path, "gitignore.template")) {
    return(".gitignore")
  }

  if (identical(relative_path, "project.Rproj.template")) {
    return(.sn_project_rproj_filename(project_dir))
  }

  relative_path
}

.sn_initialize_from_project_template <- function(
  path = ".",
  project_name = NULL,
  objective = "",
  overwrite = FALSE,
  include_governance = TRUE
) {
  project_dir <- path.expand(path)
  if (!dir.exists(project_dir)) {
    dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
  }
  project_dir <- normalizePath(project_dir, winslash = "/", mustWork = TRUE)

  if (is.null(project_name) || !nzchar(project_name)) {
    project_name <- basename(project_dir)
  }

  template_dir <- .sn_codex_component_path("project_template")
  rel_dirs <- .sn_project_template_directories()
  rel_dirs <- rel_dirs[
    !vapply(rel_dirs, .sn_should_skip_template_path, logical(1), include_governance = include_governance)
  ]
  rel_files <- list.files(
    template_dir,
    recursive = TRUE,
    all.files = TRUE,
    no.. = TRUE,
    include.dirs = FALSE
  )
  rel_files <- setdiff(rel_files, .sn_project_template_directory_manifest())
  rel_files <- rel_files[!vapply(rel_files, .sn_should_skip_template_path, logical(1), include_governance = include_governance)]

  .sn_create_directories(project_dir, rel_dirs)

  context <- list(
    project_name = project_name,
    objective = objective,
    date = as.character(Sys.Date()),
    rproj_file = .sn_project_rproj_filename(project_dir)
  )
  text_files <- .sn_project_template_text_files()

  for (relative_path in rel_files) {
    source_path <- file.path(template_dir, relative_path)
    dest_relative_path <- .sn_project_template_destination_path(relative_path, project_dir)
    dest_path <- file.path(project_dir, dest_relative_path)
    dest_parent <- dirname(dest_path)
    if (!dir.exists(dest_parent)) {
      dir.create(dest_parent, recursive = TRUE, showWarnings = FALSE)
    }

    if (file.exists(dest_path) && !isTRUE(overwrite)) {
      next
    }

    if (relative_path %in% text_files) {
      writeLines(.sn_render_text_file(source_path, context = context), con = dest_path, useBytes = TRUE)
    } else {
      file.copy(source_path, dest_path, overwrite = isTRUE(overwrite))
    }
  }

  invisible(.sn_initialized_project_paths(project_dir = project_dir, include_governance = include_governance))
}

.sn_initialized_project_paths <- function(project_dir, include_governance = TRUE) {
  out <- list(
    project_dir = project_dir,
    readme = file.path(project_dir, "README.md"),
    gitignore = file.path(project_dir, ".gitignore"),
    rproj = file.path(project_dir, .sn_project_rproj_filename(project_dir)),
    config = file.path(project_dir, "config"),
    config_default = file.path(project_dir, "config", "default.yaml"),
    data = file.path(project_dir, "data"),
    data_raw = file.path(project_dir, "data", "raw"),
    data_processed = file.path(project_dir, "data", "processed"),
    data_metadata = file.path(project_dir, "data", "metadata"),
    scripts = file.path(project_dir, "scripts"),
    notebooks = file.path(project_dir, "notebooks"),
    runs = file.path(project_dir, "runs"),
    results = file.path(project_dir, "results"),
    results_figures = file.path(project_dir, "results", "figures"),
    results_tables = file.path(project_dir, "results", "tables"),
    results_reports = file.path(project_dir, "results", "reports")
  )

  if (isTRUE(include_governance)) {
    out <- c(out, list(
      agents_md = file.path(project_dir, "AGENTS.md"),
      agents = file.path(project_dir, "AGENTS.md"),
      memory = file.path(project_dir, "memory"),
      memory_decisions = file.path(project_dir, "memory", "Decisions.md"),
      decisions = file.path(project_dir, "memory", "Decisions.md"),
      memory_plan = file.path(project_dir, "memory", "Plan.md"),
      plan = file.path(project_dir, "memory", "Plan.md"),
      memory_prompt = file.path(project_dir, "memory", "Prompt.md"),
      prompt = file.path(project_dir, "memory", "Prompt.md"),
      memory_status = file.path(project_dir, "memory", "Status.md"),
      status = file.path(project_dir, "memory", "Status.md"),
      standards = file.path(project_dir, "docs", "standards"),
      conventions = file.path(project_dir, "docs", "standards", "BioinformaticsAnalysisConventions.md"),
      skills = file.path(project_dir, "skills")
    ))
  }

  out
}

.sn_compare_version_status <- function(installed_version = NULL, remote_version = NULL) {
  if (is.null(remote_version)) {
    return(list(status = "remote unavailable", up_to_date = FALSE))
  }

  if (is.null(installed_version)) {
    return(list(status = "not installed", up_to_date = FALSE))
  }

  if (installed_version < remote_version) {
    return(list(status = "update available", up_to_date = FALSE))
  }

  if (installed_version > remote_version) {
    return(list(status = "ahead of remote", up_to_date = TRUE))
  }

  list(status = "up to date", up_to_date = TRUE)
}
