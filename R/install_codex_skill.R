#' Return the installed Shennong Codex skill path
#'
#' This helper returns the skill directory bundled with the installed
#' \pkg{Shennong} package.
#'
#' @return A character scalar path to the bundled skill directory.
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
