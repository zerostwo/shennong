make_autozyme_mock_state <- function(active = character(0),
                                     versions = NULL,
                                     fail_patch = NULL,
                                     fail_deactivate_patch = NULL) {
  manifest <- Shennong:::.sn_autozyme_patch_manifest
  state <- new.env(parent = emptyenv())
  state$status <- stats::setNames(
    rep("inactive", length(manifest)),
    names(manifest)
  )
  state$status[active] <- "active"
  state$activated <- character(0)
  state$deactivated <- character(0)
  state$fail_patch <- fail_patch
  state$fail_deactivate_patch <- fail_deactivate_patch

  default_versions <- vapply(manifest, function(spec) spec$versions[[1]], character(1))
  names(default_versions) <- vapply(manifest, `[[`, character(1), "upstream")
  if (!is.null(versions)) {
    default_versions[names(versions)] <- versions
  }
  state$versions <- default_versions

  state$call <- local({
    .state <- state
    function(fun, ...) {
      args <- list(...)
      if (identical(fun, "list_patches")) {
        return(names(.state$status))
      }
      if (identical(fun, "activate")) {
        patch <- args[[1]]
        .state$activated <- c(.state$activated, patch)
        if (identical(patch, .state$fail_patch)) {
          return(FALSE)
        }
        .state$status[[patch]] <- "active"
        return(TRUE)
      }
      if (identical(fun, "deactivate")) {
        patch <- args[[1]]
        .state$deactivated <- c(.state$deactivated, patch)
        if (identical(patch, .state$fail_deactivate_patch)) {
          stop("mock deactivation failure")
        }
        .state$status[[patch]] <- "inactive"
        return(invisible(NULL))
      }
      stop("Unexpected mocked AutoZyme call: ", fun)
    }
  })

  state
}

mock_autozyme_bindings <- function(state,
                                   autozyme_version = "0.3.1",
                                   remote_sha = Shennong:::.sn_autozyme_expected_sha,
                                   .local_envir = parent.frame()) {
  testthat::local_mocked_bindings(
    .sn_autozyme_is_installed = function(package = "autozyme") TRUE,
    .sn_autozyme_installed_version = function(package) {
      value <- state$versions[package]
      if (length(value) == 0L || is.na(value[[1]])) NA_character_ else unname(value[[1]])
    },
    .sn_autozyme_description = function() {
      list(
        version = autozyme_version,
        remote_sha = remote_sha
      )
    },
    .sn_autozyme_status_vector = function() state$status,
    .sn_autozyme_call = state$call,
    .package = "Shennong",
    .env = .local_envir
  )
}

test_that("AutoZyme APIs use conservative defaults and validate inputs", {
  expect_identical(
    eval(formals(Shennong:::sn_check_autozyme)$patches),
    c("cellchat", "nichenetr")
  )
  expect_identical(
    eval(formals(Shennong:::sn_enable_autozyme)$patches),
    c("cellchat", "nichenetr")
  )
  expect_true(all(vapply(
    c(
      "sn_check_autozyme", "sn_enable_autozyme",
      "sn_disable_autozyme", "sn_with_autozyme"
    ),
    function(name) is.function(get(name, envir = asNamespace("Shennong"))),
    logical(1)
  )))
  expect_false(any(c("threads", "install", "python") %in%
    names(formals(Shennong:::sn_enable_autozyme))))

  expect_error(
    Shennong:::sn_check_autozyme(character(0)),
    "non-empty character vector"
  )
  expect_error(
    Shennong:::sn_check_autozyme("not-a-patch"),
    "Unsupported AutoZyme patch"
  )
  expect_error(Shennong:::sn_check_autozyme(strict = 1), "must be TRUE or FALSE")
  expect_error(
    Shennong:::sn_check_autozyme(allow_approximate = NA),
    "must be TRUE or FALSE"
  )
})

test_that("AutoZyme absence is reported without installation or evaluation", {
  testthat::local_mocked_bindings(
    .sn_autozyme_is_installed = function(package = "autozyme") FALSE,
    .sn_autozyme_installed_version = function(package) NA_character_,
    .package = "Shennong"
  )

  status <- Shennong:::sn_check_autozyme()
  expect_identical(status$patch, c("cellchat", "nichenetr"))
  expect_true(all(!status$autozyme_installed))
  expect_true(all(!status$eligible))
  expect_true(all(status$reason == "autozyme is not installed"))

  expect_error(Shennong:::sn_enable_autozyme(), "optional package `autozyme` is not installed")
  expect_false(Shennong:::sn_disable_autozyme())

  evaluated <- FALSE
  expect_error(
    Shennong:::sn_with_autozyme({
      evaluated <- TRUE
    }),
    "optional package `autozyme` is not installed"
  )
  expect_false(evaluated)

  provenance <- Shennong:::.sn_autozyme_provenance()
  expect_identical(provenance, list())
})

test_that("strict gating requires the pinned AutoZyme build", {
  state <- make_autozyme_mock_state()
  mock_autozyme_bindings(state, autozyme_version = "0.3.2")
  version_drift <- Shennong:::sn_check_autozyme("cellchat")
  expect_false(version_drift$eligible)
  expect_match(version_drift$reason, "does not match pinned version")

  state <- make_autozyme_mock_state()
  mock_autozyme_bindings(state, remote_sha = NA_character_)
  source_drift <- Shennong:::sn_check_autozyme("cellchat")
  expect_false(source_drift$eligible)
  expect_match(source_drift$reason, "source does not match")

  relaxed <- Shennong:::sn_check_autozyme("cellchat", strict = FALSE)
  expect_true(relaxed$eligible)
})

test_that("strict gating blocks upstream version drift", {
  state <- make_autozyme_mock_state(
    versions = c(CellChat = "2.2.0")
  )
  mock_autozyme_bindings(state)

  strict <- Shennong:::sn_check_autozyme("cellchat")
  expect_false(strict$version_match)
  expect_false(strict$eligible)
  expect_match(strict$reason, "not exactly validated")
  expect_error(
    Shennong:::sn_enable_autozyme("cellchat"),
    "not exactly validated"
  )
  expect_length(state$activated, 0L)

  relaxed <- Shennong:::sn_check_autozyme("cellchat", strict = FALSE)
  expect_true(relaxed$eligible)
  enabled <- Shennong:::sn_enable_autozyme("cellchat", strict = FALSE)
  expect_identical(state$status[["cellchat"]], "active")
  expect_true(enabled$active)
  expect_true(Shennong:::sn_disable_autozyme("cellchat"))
  expect_identical(state$status[["cellchat"]], "inactive")
})

test_that("approximate patches require a second explicit opt-in", {
  state <- make_autozyme_mock_state()
  mock_autozyme_bindings(state)

  status <- Shennong:::sn_check_autozyme("slingshot")
  expect_true(status$approximate)
  expect_false(status$eligible)
  expect_match(status$reason, "allow_approximate = TRUE")
  expect_error(
    Shennong:::sn_enable_autozyme("slingshot"),
    "allow_approximate = TRUE"
  )

  expect_warning(
    Shennong:::sn_enable_autozyme("slingshot", allow_approximate = TRUE),
    "approximate AutoZyme patch"
  )
  expect_identical(state$status[["slingshot"]], "active")

  provenance <- Shennong:::.sn_autozyme_provenance()
  expect_identical(provenance$active_patches, "slingshot")
  expect_true(provenance$patches$slingshot$version_match)
  expect_true(provenance$patches$slingshot$approximate)
})

test_that("WGCNA remains opt-in and reports numeric-tolerance equivalence", {
  state <- make_autozyme_mock_state()
  mock_autozyme_bindings(state)

  status <- Shennong:::sn_check_autozyme("wgcna")

  expect_true(status$eligible)
  expect_false(status$default)
  expect_false(status$approximate)
  expect_identical(status$equivalence, "numeric_tolerance_scoped")
})

test_that("temporary activation restores only patches owned by its scope", {
  state <- make_autozyme_mock_state(active = "cellchat")
  mock_autozyme_bindings(state)

  value <- Shennong:::sn_with_autozyme(
    {
      expect_identical(state$status[["cellchat"]], "active")
      expect_identical(state$status[["nichenetr"]], "active")
      42L
    },
    patches = c("cellchat", "nichenetr")
  )
  expect_identical(value, 42L)
  expect_identical(state$status[["cellchat"]], "active")
  expect_identical(state$status[["nichenetr"]], "inactive")
  expect_identical(state$deactivated, "nichenetr")

  expect_error(
    Shennong:::sn_with_autozyme(
      stop("scope failure"),
      patches = "nichenetr"
    ),
    "scope failure"
  )
  expect_identical(state$status[["cellchat"]], "active")
  expect_identical(state$status[["nichenetr"]], "inactive")
  expect_identical(tail(state$deactivated, 1L), "nichenetr")
})

test_that("failed multi-patch activation rolls back newly activated patches", {
  state <- make_autozyme_mock_state(fail_patch = "nichenetr")
  mock_autozyme_bindings(state)

  expect_error(
    Shennong:::sn_enable_autozyme(),
    "did not activate patch `nichenetr`"
  )
  expect_identical(state$status[["cellchat"]], "inactive")
  expect_identical(state$status[["nichenetr"]], "inactive")
  expect_identical(state$deactivated, "cellchat")
})

test_that("failed activation reports patches that remain active after rollback", {
  state <- make_autozyme_mock_state(
    fail_patch = "nichenetr",
    fail_deactivate_patch = "cellchat"
  )
  mock_autozyme_bindings(state)

  expect_error(
    Shennong:::sn_enable_autozyme(),
    paste0(
      "did not activate patch `nichenetr`.*",
      "deactivation error.*cellchat.*",
      "residual active patch\\(es\\): cellchat"
    )
  )
  expect_identical(state$status[["cellchat"]], "active")
  expect_identical(state$status[["nichenetr"]], "inactive")
})

test_that("reading or upgrading results does not add ambient acceleration provenance", {
  original <- Shennong:::.sn_new_analysis_result(
    analysis_type = "test_result",
    name = "before_acceleration",
    method = "test",
    tables = list(primary = tibble::tibble(value = 1))
  )
  expect_null(original$provenance$acceleration)

  testthat::local_mocked_bindings(
    .sn_autozyme_provenance = function() list(active_patches = "cellchat"),
    .package = "Shennong"
  )
  upgraded <- Shennong:::.sn_upgrade_analysis_result(
    original,
    analysis_type = "test_result",
    name = "before_acceleration"
  )
  expect_null(upgraded$provenance$acceleration)

  created_while_active <- Shennong:::.sn_new_analysis_result(
    analysis_type = "test_result",
    name = "during_acceleration",
    method = "test",
    tables = list(primary = tibble::tibble(value = 1))
  )
  expect_identical(
    created_while_active$provenance$acceleration$active_patches,
    "cellchat"
  )
})

test_that("DESCRIPTION pins AutoZyme without making it a required dependency", {
  desc <- utils::packageDescription("Shennong")
  expect_match(desc[["Suggests"]], "(^|,)[[:space:]]*autozyme([[:space:]]|,|$)")
  expect_match(
    desc[["Remotes"]],
    paste0(
      "ElliotXie/autozyme/autozyme_r@",
      Shennong:::.sn_autozyme_expected_sha
    ),
    fixed = TRUE
  )
})
