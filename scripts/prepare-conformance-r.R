#!/usr/bin/env Rscript
# Install the exact direct R dependencies declared by the executable contracts.
# Never rewrite expected versions to match whatever the resolver happened to install.
contracts <- lapply(list.files("inst/conformance/contracts", "[.]json$", full.names=TRUE), jsonlite::read_json)
contracts <- Filter(function(x) identical(tolower(x$upstream$language), "r"), contracts)
specs <- unique(unlist(lapply(contracts, function(x) x$versions$upstream)))
repos <- BiocManager::repositories()
for (spec in specs) {
  package <- sub(" .*", "", spec)
  version <- sub("^[^ ]+ ", "", spec)
  installed <- tryCatch(as.character(utils::packageVersion(package)), error=function(e) "")
  if (!identical(installed, version)) {
    remotes::install_version(package, version=version, repos=repos, upgrade="never")
  }
  actual <- as.character(utils::packageVersion(package))
  if (!identical(actual, version)) stop(package, ": expected ", version, ", got ", actual)
  message(package, " ", actual)
}
