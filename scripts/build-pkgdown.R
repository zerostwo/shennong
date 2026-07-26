#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", all_args, value = TRUE)
script_file <- if (length(script_arg)) {
  sub("^--file=", "", script_arg[[1]])
} else {
  file.path(getwd(), "scripts", "build-pkgdown.R")
}
repo_root <- normalizePath(
  dirname(dirname(script_file)),
  winslash = "/",
  mustWork = TRUE
)

.usage <- function() {
  cat(
    "Usage: Rscript scripts/build-pkgdown.R [--full] [--real [--extended]] [--data-root PATH]\n",
    "\n",
    "Build modes:\n",
    "  default           Incremental build; skip unchanged pages and examples.\n",
    "  --full            Clean-quality build; rebuild all pages and run examples.\n",
    "  --real            Validate all four local real-data bundles, then run\n",
    "                    vignette chunks with SHENNONG_REAL_PROFILE=core.\n",
    "  --extended        With --real, use SHENNONG_REAL_PROFILE=all so articles\n",
    "                    may run their extended external-backend coverage.\n",
    "  --data-root PATH  Real-data root (also accepts --data-root=PATH). Defaults\n",
    "                    to SHENNONG_REAL_DATA_DIR or data-local/pkgdown-real.\n",
    "  --help            Show this message.\n",
    "\n",
    "Real-data modes require existing local artifacts and never download data.\n",
    "Combine --full with --real for a clean real-data rebuild.\n",
    sep = ""
  )
}

if ("--help" %in% args) {
  .usage()
  quit(status = 0L)
}

.parse_args <- function(args) {
  parsed <- list(full = FALSE, real = FALSE, extended = FALSE, data_root = NULL)
  index <- 1L
  while (index <= length(args)) {
    current <- args[[index]]
    if (identical(current, "--full")) {
      parsed$full <- TRUE
    } else if (identical(current, "--real")) {
      parsed$real <- TRUE
    } else if (identical(current, "--extended")) {
      parsed$extended <- TRUE
    } else if (identical(current, "--data-root")) {
      if (index == length(args) || startsWith(args[[index + 1L]], "--")) {
        stop("`--data-root` requires a path.", call. = FALSE)
      }
      index <- index + 1L
      parsed$data_root <- args[[index]]
    } else if (startsWith(current, "--data-root=")) {
      parsed$data_root <- sub("^--data-root=", "", current)
      if (!nzchar(parsed$data_root)) {
        stop("`--data-root` requires a path.", call. = FALSE)
      }
    } else {
      stop("Unknown option `", current, "`; use `--help` for usage.", call. = FALSE)
    }
    index <- index + 1L
  }
  if (parsed$extended && !parsed$real) {
    stop("`--extended` requires `--real`.", call. = FALSE)
  }
  if (!is.null(parsed$data_root) && !parsed$real) {
    stop("`--data-root` requires `--real`.", call. = FALSE)
  }
  parsed
}

.absolute_path <- function(path) {
  path <- path.expand(path)
  if (!grepl("^/", path)) path <- file.path(repo_root, path)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

config <- .parse_args(args)
full <- config$full

if (config$real) {
  default_data_root <- Sys.getenv(
    "SHENNONG_REAL_DATA_DIR",
    unset = file.path(repo_root, "data-local", "pkgdown-real")
  )
  selected_data_root <- if (is.null(config$data_root)) {
    default_data_root
  } else {
    config$data_root
  }
  data_root <- .absolute_path(selected_data_root)
  profile <- if (config$extended) "all" else "core"
  validator <- file.path(repo_root, "scripts", "real-data", "validate-real-data.R")
  if (!file.exists(validator)) {
    stop("Real-data validator was not found at `", validator, "`.", call. = FALSE)
  }

  message("Validating all four local real-data bundles before the pkgdown build.")
  validation_status <- system2(
    file.path(R.home("bin"), "Rscript"),
    c(
      shQuote(validator),
      "--bundle", "all",
      "--root", shQuote(data_root)
    )
  )
  if (!identical(as.integer(validation_status), 0L)) {
    stop(
      "Real-data validation failed; pkgdown was not built and no data were downloaded.",
      call. = FALSE
    )
  }

  Sys.setenv(
    SHENNONG_REAL_DATA_DIR = data_root,
    SHENNONG_RUN_VIGNETTES = "true",
    SHENNONG_REAL_PROFILE = profile
  )
}

if (!requireNamespace("pkgdown", quietly = TRUE)) {
  stop("Install `pkgdown` before building the site.", call. = FALSE)
}

start <- Sys.time()
build_kind <- if (full) "full" else "incremental"
if (config$real) {
  message(
    "Building the ", build_kind, " pkgdown site with real-data profile '",
    Sys.getenv("SHENNONG_REAL_PROFILE"), "'."
  )
} else {
  message("Building the ", build_kind, " pkgdown site.")
}
pkgdown::build_site(
  lazy = !full,
  examples = full,
  new_process = full,
  install = full
)
elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
message(sprintf("pkgdown build completed in %.1fs.", elapsed))

if (config$real) {
  manifest_helper <- file.path(
    repo_root,
    "scripts",
    "real-data",
    "pkgdown-build-manifest.R"
  )
  if (!file.exists(manifest_helper)) {
    stop("Pkgdown build-manifest helper was not found at `", manifest_helper, "`.", call. = FALSE)
  }
  sys.source(manifest_helper, envir = environment())

  site_dir <- normalizePath(
    file.path(repo_root, "site", "dev"),
    winslash = "/",
    mustWork = TRUE
  )
  removed_logs <- .pkgdown_remove_generated_logs(site_dir)
  if (length(removed_logs)) {
    message(
      "Removed ", length(removed_logs),
      " generated log file(s) from the static site before manifesting."
    )
  }
  runtime_report <- file.path(
    data_root,
    "results",
    "runtime-coverage-core-summary.json"
  )
  if (!file.exists(runtime_report)) {
    stop(
      "The core runtime coverage report is missing; run ",
      "scripts/real-data/run-runtime-coverage.R before the real pkgdown build.",
      call. = FALSE
    )
  }
  manifest_path <- file.path(
    data_root,
    "results",
    .pkgdown_manifest_filename
  )
  .pkgdown_write_build_manifest(
    site_dir = site_dir,
    repo_root = repo_root,
    runtime_report = runtime_report,
    manifest_path = manifest_path
  )
  message("Wrote same-build pkgdown manifest: ", manifest_path)
}
