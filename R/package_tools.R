# Package tools.
#
# Dependency listing/installation, version checks, self-installation channels, DESCRIPTION parsing, and shared template/command utilities.
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
#' tarball. When `channel = "auto"`, an explicit local `source` is installed
#' directly. Otherwise it prefers CRAN, falls back to GitHub, and, if neither
#' remote can be reached, uses the current working directory when it is a
#' source tree for `package`.
#'
#' @param channel One of \code{"auto"}, \code{"cran"}, \code{"github"}, or
#'   \code{"local"}.
#' @param package Package name. Defaults to \code{"Shennong"}.
#' @param source Installation source. For \code{channel = "github"}, supply an
#'   \code{"owner/repo"} string; when omitted, Shennong uses
#'   \code{"zerostwo/shennong"}. For \code{channel = "local"}, supply a local
#'   package directory or source tarball path. An existing local package source
#'   also selects the local channel automatically when `channel = "auto"`.
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
#' sn_install_shennong(source = ".")
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

  local_source <- .sn_find_local_package_source(
    source = source,
    package = package,
    discover = FALSE
  )
  if (identical(channel, "auto") && !is.null(local_source)) {
    check_installed("remotes", reason = "to install Shennong from a local path.")
    .sn_install_local_release(path = local_source, args = list(...))
    return(invisible("local"))
  }

  source <- source %||% "zerostwo/shennong"
  cran_version <- .sn_get_cran_version(package = package, repos = repos)
  github_version <- .sn_get_github_version(repo = source, ref = ref)

  if (identical(channel, "auto") &&
      is.null(cran_version) && is.null(github_version)) {
    local_source <- .sn_find_local_package_source(
      package = package,
      discover = TRUE
    )
    if (!is.null(local_source)) {
      warning(
        sprintf(
          "CRAN and GitHub versions are unavailable; falling back to the local %s source at '%s'.",
          package,
          local_source
        ),
        call. = FALSE
      )
      check_installed("remotes", reason = "to install Shennong from a local path.")
      .sn_install_local_release(path = local_source, args = list(...))
      return(invisible("local"))
    }
  }

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

.sn_find_local_package_source <- function(source = NULL,
                                          package = "Shennong",
                                          discover = FALSE) {
  candidate <- source
  if (is.null(candidate) && isTRUE(discover)) {
    candidate <- getwd()
  }
  if (is.null(candidate) || length(candidate) != 1L || !nzchar(candidate)) {
    return(NULL)
  }

  candidate <- path.expand(candidate)
  if (!dir.exists(candidate)) {
    return(NULL)
  }

  description <- file.path(candidate, "DESCRIPTION")
  if (!file.exists(description)) {
    return(NULL)
  }

  package_name <- tryCatch(
    read.dcf(description, fields = "Package")[[1L]],
    error = function(e) NA_character_
  )
  if (is.na(package_name) || !identical(package_name, package)) {
    return(NULL)
  }

  normalizePath(candidate, winslash = "/", mustWork = TRUE)
}

.sn_command_available <- function(command) {
  path <- Sys.which(command)[[command]]
  !is.na(path) && nzchar(path)
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
      "ROGUE", "SignatuR", "catplot", "copykat", "harmony", "liana", "lisi",
      "nichenetr", "multinichenetr", "scMetabolism"
    ),
    source = rep("GitHub", 16),
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
      "navinlabcode/copykat",
      "immunogenomics/harmony@harmony2",
      "saezlab/liana",
      "immunogenomics/lisi",
      "saeyslab/nichenetr",
      "saeyslab/multinichenetr",
      "wu-yc/scMetabolism"
    ),
    stringsAsFactors = FALSE
  )
}

.sn_bioconductor_packages <- function() {
  c(
    "apeglm", "Banksy", "BiocParallel", "clusterProfiler", "Coralysis", "decoupleR", "decontX", "DESeq2", "dorothea", "edgeR", "GENIE3", "GSEABase",
    "glmGamPoi", "miloR", "nnSVG", "org.Hs.eg.db", "org.Mm.eg.db", "progeny", "rhdf5",
    "rtracklayer", "S4Vectors", "scDblFinder", "scDesign3", "scran", "SingleCellExperiment", "SpatialExperiment", "Nebulosa",
    "SummarizedExperiment", "limma", "variancePartition"
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
