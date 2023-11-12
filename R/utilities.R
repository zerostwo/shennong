#' Check if packages are installed and load them
#'
#' This function checks if a list of packages are installed and loads them if
#' they are.
#'
#' @param ... A list of package names to check and load.
#' @param error A logical value indicating whether to stop the function with an
#'   error message if any of the packages are not installed. Defaults to TRUE.
#'
#' @return A logical vector indicating whether each package is installed or not.
#'
#' @examples
#' # Check if packages are installed and load them
#' sn_check_package("dplyr", "ggplot2")
#' @export
sn_check_package <- function(..., error = TRUE) {
  pkgs <- unlist(list(...))
  package_installed <-
    suppressWarnings(sapply(pkgs, requireNamespace))
  if (error && any(!package_installed)) {
    missing_pkgs <- pkgs[!package_installed]
    stop(
      paste(
        "Cannot find the following packages:",
        paste(missing_pkgs, collapse = ", ")
      ),
      ". Please install."
    )
  }
  invisible(package_installed)
}
