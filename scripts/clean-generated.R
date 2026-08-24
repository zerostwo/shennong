#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", all_args, value = TRUE)
script_file <- if (length(script_arg)) {
  sub("^--file=", "", script_arg[[1L]])
} else {
  file.path(getwd(), "scripts", "clean-generated.R")
}
repo_root <- normalizePath(
  file.path(dirname(normalizePath(script_file, mustWork = TRUE)), ".."),
  winslash = "/",
  mustWork = TRUE
)

if ("--help" %in% args) {
  cat(
    "Usage: Rscript scripts/clean-generated.R [--apply] [--include-site]\n",
    "\n",
    "Without --apply, print the exact ignored build/cache targets that would\n",
    "be removed. Scientific fixtures, dev outputs, user history, project\n",
    "settings, CodeGraph state, Git state, and benchmark inputs are never\n",
    "selected. --include-site also selects the reproducible pkgdown site.\n",
    sep = ""
  )
  quit(status = 0L)
}

unknown <- setdiff(args, c("--apply", "--include-site"))
if (length(unknown)) {
  stop("Unknown option(s): ", paste(unknown, collapse = ", "), call. = FALSE)
}
apply <- "--apply" %in% args
include_site <- "--include-site" %in% args

relative_targets <- c(
  "Shennong.Rcheck",
  "README.html",
  "tests/testthat/Rplots.pdf",
  "tests/testthat/omnipathr-log",
  "omnipathr-log",
  "vignettes/omnipathr-log"
)
tarballs <- list.files(
  repo_root,
  pattern = "^Shennong_[0-9][A-Za-z0-9_.-]*[.]tar[.]gz$",
  full.names = FALSE
)
python_caches <- list.dirs(
  file.path(repo_root, "inst", "pixi"),
  recursive = TRUE,
  full.names = TRUE
)
python_caches <- python_caches[basename(python_caches) == "__pycache__"]
python_caches <- substring(python_caches, nchar(repo_root) + 2L)
relative_targets <- unique(c(
  relative_targets,
  tarballs,
  python_caches,
  if (include_site) "site" else character()
))

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || is.na(x)) y else x
}

.inside_repo <- function(path) {
  normalized <- normalizePath(path, winslash = "/", mustWork = FALSE)
  startsWith(normalized, paste0(repo_root, "/")) && !identical(normalized, repo_root)
}

.path_size <- function(path) {
  if (!file.exists(path) && !dir.exists(path)) return(0)
  if (!dir.exists(path)) return(unname(file.info(path)$size %||% 0))
  files <- list.files(path, recursive = TRUE, all.files = TRUE, full.names = TRUE)
  files <- files[file.exists(files) & !dir.exists(files)]
  if (!length(files)) return(0)
  sum(file.info(files)$size, na.rm = TRUE)
}

rows <- lapply(relative_targets, function(relative) {
  path <- file.path(repo_root, relative)
  if (!.inside_repo(path)) stop("Refusing target outside repository: ", path)
  exists <- file.exists(path) || dir.exists(path)
  link <- if (exists) Sys.readlink(path) else ""
  if (!is.na(link) && nzchar(link)) {
    stop("Refusing symbolic-link target: ", relative)
  }
  tracked <- system2(
    "git",
    c("-c", paste0("safe.directory=", repo_root), "ls-files", "--", relative),
    stdout = TRUE,
    stderr = FALSE
  )
  if (length(tracked)) {
    stop("Refusing target containing tracked files: ", relative, call. = FALSE)
  }
  data.frame(
    path = relative,
    exists = exists,
    bytes = .path_size(path),
    stringsAsFactors = FALSE
  )
})
report <- do.call(rbind, rows)
report <- report[report$exists, , drop = FALSE]
report$mib <- round(report$bytes / 1024^2, 3)
print(report[, c("path", "mib")], row.names = FALSE)
cat(sprintf("Total selected: %.2f MiB\n", sum(report$bytes) / 1024^2))

if (!apply) {
  cat("Dry run only. Re-run with --apply to remove exactly these targets.\n")
  quit(status = 0L)
}

for (relative in report$path) {
  path <- file.path(repo_root, relative)
  status <- unlink(path, recursive = TRUE, force = FALSE)
  if (!identical(status, 0L) || file.exists(path) || dir.exists(path)) {
    stop("Could not remove generated target: ", relative, call. = FALSE)
  }
  cat("Removed ", relative, "\n", sep = "")
}
