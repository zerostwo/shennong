#!/usr/bin/env Rscript

# Provision the same managed interpreters used by the candidate and oracle.
# Pin the recorded Python distributions only in this CI runtime, leaving the
# package's user-facing manifests unchanged.
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
  paths <- sn_prepare_pixi_environment(environment, overwrite = TRUE)
  manifest <- readLines(paths$manifest_path)
  for (spec in unlist(contract$versions$upstream, use.names = FALSE)) {
    package <- sub(" .*", "", spec)
    version <- sub("^[^ ]+ ", "", spec)
    entry <- which(startsWith(manifest, paste0(package, " = ")))
    pin <- paste0(package, ' = "==', version, '"')
    if (length(entry)) {
      manifest[entry] <- pin
    } else {
      section <- match("[pypi-dependencies]", manifest)
      if (is.na(section)) stop("No Python dependency section for ", package)
      manifest <- append(manifest, pin, after = section)
    }
  }
  writeLines(manifest, paths$manifest_path)
  sn_prepare_pixi_environment(environment, install_environment = TRUE, pixi = pixi$path)
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
