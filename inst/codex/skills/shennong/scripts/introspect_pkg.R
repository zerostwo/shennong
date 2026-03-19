#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
pkg <- if (length(args) >= 1) args[[1]] else "Shennong"

if (!requireNamespace(pkg, quietly = TRUE)) {
  stop(sprintf("Package '%s' is not installed.", pkg))
}

ns <- asNamespace(pkg)
exports <- getNamespaceExports(ns)
exports <- exports[grepl("^sn_", exports)]
exports <- sort(exports)

cat(sprintf("Package: %s\n", pkg))
cat(sprintf("Version: %s\n\n", as.character(utils::packageVersion(pkg))))
cat("Exported user-facing functions:\n")
cat(paste0("- ", exports, collapse = "\n"))
cat("\n\n")

family_counts <- sort(table(sub("^((sn_[^_]+)_?).*$", "\\1", exports)), decreasing = TRUE)
cat("Function families:\n")
for (family in names(family_counts)) {
  cat(sprintf("- %s: %d\n", family, family_counts[[family]]))
}
cat("\n")

cat("Representative help topics:\n")
for (topic in utils::head(exports, 10)) {
  help_path <- utils::help(topic, package = pkg)
  cat(sprintf("- %s: %s\n", topic, as.character(help_path)))
}
