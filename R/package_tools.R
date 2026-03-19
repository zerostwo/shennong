#' Return the installed Shennong Codex skill path
#'
#' This helper returns the skill directory bundled with the installed
#' \pkg{Shennong} package.
#'
#' @return A character scalar path to the bundled skill directory.
#'
#' @examples
#' if (requireNamespace("Shennong", quietly = TRUE)) {
#'   path <- sn_get_codex_skill_path()
#'   dir.exists(path)
#' }
#' @export
sn_get_codex_skill_path <- function() {
  path <- system.file("codex/skills/shennong", package = "Shennong")
  if (!nzchar(path)) {
    stop("The bundled Shennong Codex skill could not be found in the installed package.")
  }
  path
}

#' Install the bundled Shennong Codex skill for end-user agents
#'
#' This helper copies the packaged Shennong skill into a user skill directory,
#' for example \code{~/.agents/skills}. The installed skill is intended for
#' analysis environments where the package is already available.
#'
#' @param path Destination directory that will contain the \code{shennong}
#'   skill folder. Defaults to \code{"~/.agents/skills"}.
#' @param overwrite Whether to replace an existing installed copy of the skill.
#'   Defaults to \code{TRUE}.
#'
#' @return Invisibly returns the installed skill directory path.
#'
#' @examples
#' target_dir <- file.path(tempdir(), "skills")
#' installed <- sn_install_codex_skill(path = target_dir, overwrite = TRUE)
#' dir.exists(installed)
#' @export
sn_install_codex_skill <- function(path = "~/.agents/skills", overwrite = TRUE) {
  source_dir <- sn_get_codex_skill_path()
  dest_root <- path.expand(path)
  if (!dir.exists(dest_root)) {
    dir.create(dest_root, recursive = TRUE, showWarnings = FALSE)
  }

  dest_dir <- file.path(dest_root, "shennong")
  if (dir.exists(dest_dir)) {
    if (!isTRUE(overwrite)) {
      stop(glue("Destination '{dest_dir}' already exists. Set `overwrite = TRUE` to replace it."))
    }
    unlink(dest_dir, recursive = TRUE, force = TRUE)
  }

  ok <- file.copy(from = source_dir, to = dest_root, recursive = TRUE)
  if (!isTRUE(ok)) {
    stop(glue("Failed to copy the bundled Shennong Codex skill to '{dest_root}'."))
  }

  invisible(dest_dir)
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
    message(
      glue(
        "Shennong ({resolved_channel}): installed = {installed_label}; ",
        "latest = {remote_label}; status = {status$status}"
      )
    )
    if (!isTRUE(status$up_to_date)) {
      message("Install/update with: ", install_command)
    }
  }

  invisible(result)
}

#' Install Shennong from CRAN or GitHub
#'
#' This helper installs the stable CRAN release when available, or the GitHub
#' development version when requested. When `channel = "auto"`, it prefers CRAN
#' and falls back to GitHub if no CRAN release is available.
#'
#' @param channel One of \code{"auto"}, \code{"cran"}, or \code{"github"}.
#' @param package Package name. Defaults to \code{"Shennong"}.
#' @param github_repo GitHub repository in \code{"owner/repo"} format.
#' @param github_ref GitHub ref to install from. Defaults to \code{"main"}.
#' @param repos CRAN-like repositories used by \code{install.packages()}.
#' @param ... Additional arguments passed to \code{utils::install.packages()} or
#'   \code{remotes::install_github()}.
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
  channel = c("auto", "cran", "github"),
  package = "Shennong",
  github_repo = "zerostwo/shennong",
  github_ref = "main",
  repos = getOption("repos"),
  ...
) {
  channel <- match.arg(channel)
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
  remotes::install_github(repo = github_repo, ref = github_ref, ...)
  invisible(resolved_channel)
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
