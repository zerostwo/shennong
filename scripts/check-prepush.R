#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

.parse_flag <- function(name) {
  any(args == paste0("--", name))
}

.parse_value <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  matched <- args[startsWith(args, prefix)]
  if (length(matched) == 0) {
    return(default)
  }
  sub(prefix, "", matched[[1]], fixed = TRUE)
}

.usage <- function() {
  cat(
    paste(
      "Usage: Rscript scripts/check-prepush.R [options]",
      "",
      "Options:",
      "  --filter=<regex>       Run a targeted test pass before the full suite.",
      "  --skip-document        Skip devtools::document().",
      "  --skip-tests           Skip all testthat runs.",
      "  --skip-build           Skip R CMD build.",
      "  --skip-check           Skip R CMD check --no-manual.",
      "  --help                 Show this message.",
      sep = "\n"
    ),
    "\n"
  )
}

if (.parse_flag("help")) {
  .usage()
  quit(status = 0)
}

filter <- .parse_value("filter")
skip_document <- .parse_flag("skip-document")
skip_tests <- .parse_flag("skip-tests")
skip_build <- .parse_flag("skip-build")
skip_check <- .parse_flag("skip-check")

run_step <- function(label, expr) {
  message("==> ", label)
  force(expr)
  invisible(TRUE)
}

run_cmd <- function(args, wd = getwd()) {
  r_bin <- file.path(R.home("bin"), "R")
  oldwd <- setwd(wd)
  on.exit(setwd(oldwd), add = TRUE)
  status <- system2(r_bin, args = args, stdout = "", stderr = "", wait = TRUE, env = character())
  if (!identical(status, 0L)) {
    stop("Command failed: ", paste(c(r_bin, args), collapse = " "), call. = FALSE)
  }
  invisible(TRUE)
}

desc <- read.dcf("DESCRIPTION", fields = c("Package", "Version"))
pkg_name <- desc[1, "Package"]
pkg_version <- desc[1, "Version"]
tarball <- file.path(getwd(), sprintf("%s_%s.tar.gz", pkg_name, pkg_version))

if (!skip_document) {
  run_step("Regenerating documentation", {
    if (!requireNamespace("devtools", quietly = TRUE)) {
      stop("`devtools` must be installed to run the pre-push check.", call. = FALSE)
    }
    devtools::document(quiet = TRUE)
  })
}

if (!skip_tests) {
  run_step("Running targeted tests", {
    if (!is.null(filter) && nzchar(filter)) {
      testthat::test_local(filter = filter, stop_on_failure = TRUE)
    } else {
      message("Skipping targeted tests because no `--filter=` value was supplied.")
    }
  })

  run_step("Running full test suite", {
    testthat::test_local(stop_on_failure = TRUE)
  })
}

if (!skip_build) {
  run_step("Building source tarball", {
    if (file.exists(tarball)) {
      unlink(tarball)
    }
    run_cmd(c("CMD", "build", "."))
    if (!file.exists(tarball)) {
      stop("Expected tarball was not created: ", tarball, call. = FALSE)
    }
  })
}

if (!skip_check) {
  run_step("Running R CMD check --no-manual", {
    if (!file.exists(tarball)) {
      stop("Tarball not found. Run the build step before check.", call. = FALSE)
    }
    run_cmd(c("CMD", "check", "--no-manual", basename(tarball)))
  })
}

message("Pre-push check completed successfully.")
