.conformance_project_root <- function() {
  normalizePath(
    testthat::test_path("..", ".."),
    winslash = "/",
    mustWork = TRUE
  )
}

.conformance_contract_dir <- function() {
  installed <- system.file("conformance", "contracts", package = "Shennong")
  source <- file.path(
    .conformance_project_root(),
    "inst", "conformance", "contracts"
  )
  candidates <- unique(c(if (nzchar(installed)) installed else character(), source))
  candidates <- candidates[dir.exists(candidates)]
  if (length(candidates) == 0L) {
    stop("Backend conformance contracts were not found.", call. = FALSE)
  }
  candidates[[1]]
}

.conformance_read_contracts <- function() {
  paths <- sort(list.files(
    .conformance_contract_dir(),
    pattern = "[.]json$",
    full.names = TRUE
  ))
  if (length(paths) == 0L) {
    stop("No backend conformance contracts were found.", call. = FALSE)
  }
  contracts <- lapply(paths, jsonlite::read_json, simplifyVector = FALSE)
  names(contracts) <- vapply(contracts, `[[`, character(1), "id")
  contracts
}

.conformance_contract <- function(id) {
  contracts <- .conformance_read_contracts()
  if (!id %in% names(contracts)) {
    stop("Unknown backend conformance contract: ", id, call. = FALSE)
  }
  contracts[[id]]
}

.conformance_read_method_inventory <- function(filename) {
  path <- file.path(
    .conformance_project_root(),
    "tests", "conformance", filename
  )
  lines <- trimws(readLines(path, warn = FALSE))
  sort(lines[nzchar(lines) & !startsWith(lines, "#")])
}

.conformance_read_legacy_methods <- function() {
  .conformance_read_method_inventory("legacy-methods.txt")
}

.conformance_read_pending_methods <- function() {
  .conformance_read_method_inventory("legacy-pending-methods.txt")
}

.conformance_strict <- function() {
  tolower(Sys.getenv("SHENNONG_CONFORMANCE_STRICT", unset = "false")) %in%
    c("1", "true", "yes", "on")
}

.conformance_require_package <- function(package) {
  if (requireNamespace(package, quietly = TRUE)) {
    return(invisible(TRUE))
  }
  message <- paste0("Required conformance dependency `", package, "` is unavailable.")
  if (.conformance_strict()) {
    testthat::fail(message)
  }
  testthat::skip(message)
}

.conformance_fingerprint <- function(object) {
  digest::digest(object, algo = "sha256", serialize = TRUE)
}

.conformance_expect_unchanged <- function(object, before, label = "input") {
  testthat::expect_identical(
    .conformance_fingerprint(object),
    before,
    info = paste(label, "was modified by the candidate call")
  )
}

.conformance_without_acceleration <- function(code) {
  expression <- substitute(code)
  old_options <- options(shennong.acceleration = FALSE)
  old_disabled <- Sys.getenv("SHENNONG_ACCELERATION_DISABLED", unset = NA_character_)
  on.exit({
    options(old_options)
    if (is.na(old_disabled)) {
      Sys.unsetenv("SHENNONG_ACCELERATION_DISABLED")
    } else {
      Sys.setenv(SHENNONG_ACCELERATION_DISABLED = old_disabled)
    }
  }, add = TRUE)
  Sys.setenv(SHENNONG_ACCELERATION_DISABLED = "true")

  if (base::isNamespaceLoaded("ShennongOpt")) {
    with_disabled <- Shennong:::.sn_acceleration_get_export(
      "sn_with_acceleration_disabled"
    )
    return(with_disabled(eval.parent(expression)))
  }
  eval.parent(expression)
}

.conformance_expect_equal <- function(actual, expected, contract, info = NULL) {
  tolerance <- max(
    as.numeric(contract$equivalence$absolute_tolerance %||% 0),
    as.numeric(contract$equivalence$relative_tolerance %||% 0)
  )
  testthat::expect_equal(
    actual,
    expected,
    tolerance = tolerance,
    info = info
  )
}

.conformance_scrublet_env_installed <- function() {
  tryCatch(
    {
      Shennong:::.sn_pixi_script_path(environment = "scrublet", script_name = "scrublet_run.py")
      paths <- Shennong::sn_pixi_paths(environment = "scrublet")
      file.exists(file.path(paths$workspace_env_dir, "default", "bin", "python"))
    },
    error = function(e) FALSE
  )
}
