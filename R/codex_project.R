# Codex assets and project initialization.
#
# Extracted from package_tools.R: bundled Codex asset path resolution, skill installation, and project initialization from the packaged template.
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
