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
#'   template. Defaults to \code{FALSE}.
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
#' file.exists(created$agents_md)
#' @export
sn_initialize_project <- function(
  path = ".",
  project_name = NULL,
  objective = "Build a reproducible Shennong-based single-cell analysis workflow for this project.",
  with_agent = TRUE,
  overwrite = FALSE
) {
  if (isTRUE(with_agent)) {
    return(sn_initialize_codex_project(
      path = path,
      project_name = project_name,
      objective = objective,
      overwrite = overwrite
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

#' Initialize Codex-style project guidance for a Shennong analysis
#'
#' This helper scaffolds a governed analysis project from the packaged
#' \code{inst/codex/project-template/} assets.
#'
#' @inheritParams sn_initialize_project
#'
#' @return Invisibly returns the same structure as
#'   \code{\link{sn_initialize_project}()}.
#'
#' @examples
#' project_dir <- file.path(tempdir(), "codex-analysis-project")
#' created <- sn_initialize_codex_project(
#'   path = project_dir,
#'   project_name = "PBMC pilot study",
#'   objective = "Build a reproducible Shennong-based PBMC analysis workflow.",
#'   overwrite = TRUE
#' )
#' file.exists(created$agents)
#' file.exists(created$prompt)
#' @export
sn_initialize_codex_project <- function(
  path = ".",
  project_name = NULL,
  objective = "Build a reproducible Shennong-based single-cell analysis workflow for this project.",
  overwrite = FALSE
) {
  .sn_initialize_from_project_template(
    path = path,
    project_name = project_name,
    objective = objective,
    overwrite = overwrite,
    include_governance = TRUE
  )
}

#' Check whether Shennong is up to date
#'
#' This helper compares the installed package version against the latest
#' available CRAN or GitHub development version. When `channel = "auto"`, it
#' prefers CRAN if a release exists there and otherwise falls back to GitHub.
#'
#' @param channel One of \code{"auto"}, \code{"cran"}, or \code{"github"}.
#' @param package Package name to check. Defaults to \code{"Shennong"}.
#' @param github_repo GitHub repository in \code{"owner/repo"} format.
#' @param github_ref GitHub ref to inspect. Defaults to \code{"main"}.
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
  github_repo = "zerostwo/shennong",
  github_ref = "main",
  repos = getOption("repos"),
  quiet = FALSE
) {
  channel <- match.arg(channel)
  installed_version <- .sn_get_installed_version(package = package)
  cran_version <- .sn_get_cran_version(package = package, repos = repos)
  github_version <- .sn_get_github_version(repo = github_repo, ref = github_ref)
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
    github = sprintf('remotes::install_github("%s", ref = "%s")', github_repo, github_ref)
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

#' Install Shennong from CRAN or GitHub
#'
#' This helper installs the stable CRAN release when available, or the GitHub
#' development version when requested, or installs from a local source tree or
#' tarball. When `channel = "auto"`, it prefers CRAN and falls back to GitHub
#' if no CRAN release is available.
#'
#' @param channel One of \code{"auto"}, \code{"cran"}, \code{"github"}, or
#'   \code{"local"}.
#' @param package Package name. Defaults to \code{"Shennong"}.
#' @param github_repo GitHub repository in \code{"owner/repo"} format.
#' @param github_ref GitHub ref to install from. Defaults to \code{"main"}.
#' @param local_path Local package directory or source tarball used when
#'   \code{channel = "local"}.
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
#' }
#'
#' @export
sn_install_shennong <- function(
  channel = c("auto", "cran", "github", "local"),
  package = "Shennong",
  github_repo = "zerostwo/shennong",
  github_ref = "main",
  local_path = NULL,
  repos = getOption("repos"),
  ...
) {
  channel <- match.arg(channel)
  if (identical(channel, "local")) {
    if (is.null(local_path) || !nzchar(local_path)) {
      stop("`local_path` must be supplied when `channel = \"local\"`.", call. = FALSE)
    }
    check_installed("remotes", reason = "to install Shennong from a local path.")
    .sn_install_local_release(
      path = local_path,
      args = list(...)
    )
    return(invisible(channel))
  }

  cran_version <- .sn_get_cran_version(package = package, repos = repos)
  github_version <- .sn_get_github_version(repo = github_repo, ref = github_ref)
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
    repo = github_repo,
    ref = github_ref,
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
    stop("`github_repo` must use the form 'owner/repo'.", call. = FALSE)
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
    date = as.character(Sys.Date())
  )
  text_files <- .sn_project_template_text_files()

  for (relative_path in rel_files) {
    source_path <- file.path(template_dir, relative_path)
    dest_path <- file.path(project_dir, relative_path)
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
