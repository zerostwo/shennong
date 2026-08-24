#!/usr/bin/env Rscript

all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", all_args, value = TRUE)
script_file <- if (length(script_arg)) {
  sub("^--file=", "", script_arg[[1L]])
} else {
  file.path(getwd(), "scripts", "conformance", "build-public-api-matrix.R")
}
repo_root <- normalizePath(
  file.path(dirname(normalizePath(script_file, mustWork = TRUE)), "..", ".."),
  winslash = "/",
  mustWork = TRUE
)
trailing_args <- commandArgs(trailingOnly = TRUE)
output_args <- grep("^--output=", trailing_args, value = TRUE)
if (length(output_args) > 1L) {
  stop("Supply at most one `--output=PATH` argument.", call. = FALSE)
}
output <- if (length(output_args) == 1L) {
  path.expand(sub("^--output=", "", output_args[[1L]]))
} else {
  file.path(repo_root, "inst", "conformance", "public-api-parameters.json")
}

required <- c("digest", "jsonlite", "pkgload")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing matrix package(s): ", paste(missing, collapse = ", "), ".", call. = FALSE)
}
pkgload::load_all(repo_root, quiet = TRUE)

.default_text <- function(value) {
  if (identical(value, quote(expr = ))) return("<required>")
  paste(deparse(value, width.cutoff = 500L), collapse = "")
}

.selector_values <- function(value) {
  if (!is.call(value) || !identical(value[[1L]], as.name("c"))) return(NULL)
  values <- as.list(value)[-1L]
  if (!length(values) || !all(vapply(values, is.character, logical(1)))) return(NULL)
  unname(unlist(values, use.names = FALSE))
}

.parameter_class <- function(name, value) {
  if (identical(name, "...")) return("dots")
  if (identical(value, quote(expr = ))) return("required")
  if (!is.null(.selector_values(value))) return("selector")
  if (is.null(value) || identical(value, quote(NULL))) return("nullable")
  if (is.logical(value) && length(value) == 1L) return("boolean")
  if (is.numeric(value)) return("numeric")
  if (is.character(value)) return("character")
  "expression"
}

coverage <- utils::read.csv(
  file.path(repo_root, "scripts", "real-data", "coverage.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
exclusions <- utils::read.csv(
  file.path(repo_root, "scripts", "real-data", "coverage-exclusions.csv"),
  stringsAsFactors = FALSE,
  check.names = FALSE
)
classified <- c(coverage[["function"]], exclusions[["function"]])

exports <- sort(
  grep("^sn_", getNamespaceExports("Shennong"), value = TRUE),
  method = "radix"
)
subjects <- lapply(exports, function(name) {
  fun <- get(name, envir = asNamespace("Shennong"), inherits = FALSE)
  defaults <- formals(fun)
  parameters <- Map(function(parameter, value) {
    selector <- .selector_values(value)
    list(
      name = parameter,
      class = .parameter_class(parameter, value),
      default = .default_text(value),
      selector_values = selector,
      minimum_cases = switch(
        .parameter_class(parameter, value),
        selector = length(selector),
        boolean = 2L,
        numeric = 4L,
        nullable = 3L,
        required = 2L,
        dots = 1L,
        1L
      )
    )
  }, names(defaults), as.list(defaults))
  list(
    function_name = name,
    runtime_classified = name %in% classified,
    formal_digest = digest::digest(
      vapply(defaults, .default_text, character(1)),
      algo = "sha256",
      serialize = TRUE
    ),
    parameters = unname(parameters)
  )
})

payload <- list(
  schema_version = "shennong.public-api-parameters/v1",
  package_version = unname(read.dcf(
    file.path(repo_root, "DESCRIPTION"),
    fields = "Version"
  )[[1L]]),
  function_count = length(subjects),
  generation_rule = paste(
    "Every exported sn_* formal is inventoried. Selector values require one",
    "case each; other minimum_cases are planning requirements, not evidence",
    "that those cases have already passed."
  ),
  functions = subjects
)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(
  payload,
  output,
  auto_unbox = TRUE,
  pretty = TRUE,
  null = "null",
  na = "null",
  digits = NA
)
cat("Wrote ", normalizePath(output, winslash = "/", mustWork = TRUE), "\n", sep = "")
