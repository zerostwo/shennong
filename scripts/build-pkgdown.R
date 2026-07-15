#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
full <- "--full" %in% args

if ("--help" %in% args) {
  message(
    "Usage: Rscript scripts/build-pkgdown.R [--full]\n",
    "  default  Incremental build; skip unchanged pages and example execution.\n",
    "  --full   Clean-quality build; rebuild all pages and run examples."
  )
  quit(status = 0L)
}

if (!requireNamespace("pkgdown", quietly = TRUE)) {
  stop("Install `pkgdown` before building the site.", call. = FALSE)
}

start <- Sys.time()
message(if (full) "Building the full pkgdown site." else "Building the incremental pkgdown site.")
pkgdown::build_site(
  lazy = !full,
  examples = full,
  new_process = full,
  install = full
)
elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))
message(sprintf("pkgdown build completed in %.1fs.", elapsed))
