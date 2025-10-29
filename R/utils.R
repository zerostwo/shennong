#' Check if files exist
#'
#' This function takes a vector of file paths as input and checks if each file exists.
#' If a file does not exist, the function returns a message indicating which file(s) do(es) not exist.
#'
#' @param x A vector of file paths.
#' @param stop A logical value indicating whether to stop the function if a file does not exist.
#' Default is TRUE.
#'
#' @return If stop is TRUE, the function stops and returns a message indicating which file(s) do(es) not exist.
#' If stop is FALSE, the function returns a vector of file paths that do not exist.
#'
#' @examples
#' sn_check_file(c("file1.txt", "file2.txt", "file3.txt"))
#' sn_check_file(c("file1.txt", "file2.txt", "file3.txt"), stop = FALSE)
#'
#' @export
#' @keywords file, existence
#' @seealso \code{\link{file.exists}}, \code{\link{stop}}
sn_check_file <- function(x, stop = TRUE) {
  not_exist <- x[!file.exists(x)]
  if (length(x = not_exist) > 0) {
    if (stop) {
      stop(
        paste0(
          "Please check the file path, the file listed below does not exist:\n"
        ),
        paste0(not_exist, collapse = "\n")
      )
    } else {
      return(not_exist)
    }
  }
}

#' @export
sn_set_path <- function(path) {
  path <- glue(path)
  if (!file.exists(path)) {
    dir.create(path, recursive = TRUE)
  }
  return(path)
}

check_installed_github <- function(pkg, repo, reason = NULL) {
  check_installed(pkg = "remotes", reason = "to install packages from GitHub")
  check_installed(
    pkg,
    action = function(pkg, ...) {
        remotes::install_github(repo = repo)
    },
    reason = reason
  )
}
