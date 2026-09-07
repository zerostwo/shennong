#!/usr/bin/env Rscript

# Provision the same locked interpreters used by the candidate and oracle.
# Contract versions are verified below against the package-shipped lock; this
# script must never rewrite the manifest after the lock has been materialized.
conformance_runtime_dir <- Sys.getenv("SHENNONG_RUNTIME_DIR", unset = "")
if (!nzchar(conformance_runtime_dir)) {
  stop("Set SHENNONG_RUNTIME_DIR to a dedicated CI directory before provisioning.", call. = FALSE)
}
pkgload::load_all(".", quiet = TRUE)
options(shennong.runtime_dir = conformance_runtime_dir)
contracts <- lapply(list.files("inst/conformance/contracts", pattern = "[.]json$",
                              full.names = TRUE), jsonlite::read_json)
contracts <- Filter(function(x) identical(tolower(x$upstream$language), "python"), contracts)
pixi <- sn_ensure_pixi(install = TRUE, version = "0.69.0")
for (contract in contracts) {
  environment <- sub("^pixi environment `([^`]+)`.*$", "\\1", contract$upstream$runtime)
  if (identical(environment, contract$upstream$runtime)) stop("Missing pixi environment in contract")
  paths <- sn_prepare_pixi_environment(
    environment,
    overwrite = TRUE,
    install_environment = TRUE,
    pixi = pixi$path
  )
  # Fail during provisioning if metadata disagrees, before running the tests.
  for (spec in unlist(contract$versions$upstream, use.names = FALSE)) {
    package <- sub(" .*", "", spec)
    version <- sub("^[^ ]+ ", "", spec)
    python <- file.path(paths$workspace_env_dir, "default", "bin", "python")
    actual <- system2(python, c("-c", shQuote(
      "import importlib.metadata, sys; print(importlib.metadata.version(sys.argv[1]))"),
      shQuote(package)), stdout = TRUE)
    stopifnot(is.null(attr(actual, "status")), identical(actual, version))
    message(environment, ": ", package, " ", actual)
  }
}
