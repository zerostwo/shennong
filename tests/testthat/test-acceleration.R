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

test_that("AutoZyme APIs expose automatic defaults and validate inputs", {
  expect_identical(
    eval(formals(Shennong:::sn_check_autozyme)$patches),
    Shennong:::.sn_autozyme_default_patches
  )
  expect_identical(
    eval(formals(Shennong:::sn_enable_autozyme)$patches),
    c("cellchat", "nichenetr")
  )
  expect_identical(
    eval(formals(Shennong:::sn_disable_autozyme)$patches),
    Shennong:::.sn_autozyme_default_patches
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
  expect_identical(status$patch, Shennong:::.sn_autozyme_default_patches)
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

test_that("automatic patch set covers every integrated non-approximate backend", {
  expect_identical(
    Shennong:::.sn_autozyme_default_patches,
    c(
      "cellchat", "clusterprofiler", "fgsea", "nichenetr", "seurat",
      "tradeseq", "wgcna"
    )
  )

  state <- make_autozyme_mock_state()
  mock_autozyme_bindings(state)

  status <- Shennong:::sn_check_autozyme("wgcna")

  expect_true(status$eligible)
  expect_true(status$default)
  expect_false(status$approximate)
  expect_identical(status$equivalence, "numeric_tolerance_scoped")
})

test_that("automatic AutoZyme activation is strict, scoped, and provenance-aware", {
  state <- make_autozyme_mock_state()
  mock_autozyme_bindings(state)

  result <- Shennong:::.sn_with_default_autozyme(
    Shennong:::.sn_new_analysis_result(
      analysis_type = "test_result",
      name = "default_acceleration",
      method = "test",
      tables = list(primary = tibble::tibble(value = 1))
    ),
    patches = c("cellchat", "wgcna")
  )
  expect_identical(state$activated, c("cellchat", "wgcna"))
  expect_identical(state$deactivated, c("wgcna", "cellchat"))
  expect_identical(state$status[["cellchat"]], "inactive")
  expect_identical(state$status[["wgcna"]], "inactive")
  expect_identical(
    result$provenance$acceleration$active_patches,
    c("cellchat", "wgcna")
  )
})

test_that("automatic AutoZyme activation honors session opt-outs", {
  state <- make_autozyme_mock_state()
  mock_autozyme_bindings(state)
  withr::local_options(list(shennong.autozyme = FALSE))

  expect_identical(
    Shennong:::.sn_with_default_autozyme(41L, "cellchat"),
    41L
  )
  expect_length(state$activated, 0L)

  withr::local_options(list(shennong.autozyme = TRUE))
  withr::local_envvar(c(AUTOZYME_DISABLED = "true"))
  expect_identical(
    Shennong:::.sn_with_default_autozyme(42L, "cellchat"),
    42L
  )
  expect_length(state$activated, 0L)

  withr::local_envvar(c(AUTOZYME_DISABLED = NA, AUTOZYME_DISABLE = "yes"))
  expect_identical(
    Shennong:::.sn_with_default_autozyme(43L, "cellchat"),
    43L
  )
  expect_length(state$activated, 0L)

  # The automatic opt-out does not change explicit helper behavior.
  withr::local_options(list(shennong.autozyme = FALSE))
  withr::local_envvar(c(AUTOZYME_DISABLED = NA, AUTOZYME_DISABLE = NA))
  Shennong:::sn_enable_autozyme("cellchat")
  expect_identical(state$status[["cellchat"]], "active")
})

test_that("automatic activation restores AutoZyme's future option side effect", {
  state <- make_autozyme_mock_state()
  original_call <- state$call
  state$call <- function(fun, ...) {
    if (identical(fun, "list_patches")) {
      options(future.globals.maxSize = 16 * 1024^3)
    }
    original_call(fun, ...)
  }
  mock_autozyme_bindings(state)
  withr::local_options(list(future.globals.maxSize = 123))

  value <- Shennong:::.sn_with_default_autozyme(
    {
      expect_identical(getOption("future.globals.maxSize"), 123)
      options(future.globals.maxSize = 456)
      42L
    },
    patches = "cellchat"
  )
  expect_identical(value, 42L)
  expect_identical(getOption("future.globals.maxSize"), 456)

  previous <- options()
  previous_present <- "future.globals.maxSize" %in% names(previous)
  previous_value <- previous[["future.globals.maxSize"]]
  on.exit({
    value <- if (previous_present) previous_value else NULL
    options(future.globals.maxSize = value)
  }, add = TRUE)
  options(future.globals.maxSize = NULL)
  state <- make_autozyme_mock_state()
  original_call <- state$call
  state$call <- function(fun, ...) {
    if (identical(fun, "list_patches")) {
      options(future.globals.maxSize = 16 * 1024^3)
    }
    original_call(fun, ...)
  }
  mock_autozyme_bindings(state)
  expect_identical(
    Shennong:::.sn_with_default_autozyme(
      {
        expect_null(getOption("future.globals.maxSize"))
        43L
      },
      patches = "cellchat"
    ),
    43L
  )
  expect_null(getOption("future.globals.maxSize"))
})

test_that("automatic activation safely falls back after a complete rollback", {
  state <- make_autozyme_mock_state(fail_patch = "nichenetr")
  original_call <- state$call
  state$call <- function(fun, ...) {
    if (identical(fun, "list_patches")) {
      options(future.globals.maxSize = 16 * 1024^3)
    }
    original_call(fun, ...)
  }
  mock_autozyme_bindings(state)
  withr::local_options(list(future.globals.maxSize = 123))
  executions <- 0L

  expect_warning(
    value <- Shennong:::.sn_with_default_autozyme(
      {
        expect_identical(getOption("future.globals.maxSize"), 123)
        executions <- executions + 1L
        42L
      },
      patches = c("cellchat", "nichenetr")
    ),
    "continuing without it"
  )
  expect_identical(value, 42L)
  expect_identical(executions, 1L)
  expect_identical(state$status[["cellchat"]], "inactive")
  expect_identical(state$status[["nichenetr"]], "inactive")
  expect_identical(state$deactivated, "cellchat")
  expect_identical(getOption("future.globals.maxSize"), 123)
})

test_that("automatic activation stops when rollback leaves residual state", {
  state <- make_autozyme_mock_state(
    fail_patch = "nichenetr",
    fail_deactivate_patch = "cellchat"
  )
  mock_autozyme_bindings(state)
  executions <- 0L

  expect_error(
    Shennong:::.sn_with_default_autozyme(
      {
        executions <- executions + 1L
        42L
      },
      patches = c("cellchat", "nichenetr")
    ),
    "could not be safely restored.*residual active patch.*cellchat"
  )
  expect_identical(executions, 0L)
  expect_identical(state$status[["cellchat"]], "active")
  expect_identical(state$status[["nichenetr"]], "inactive")
})

test_that("automatic scope restores state after analytical and cleanup failures", {
  state <- make_autozyme_mock_state(active = "cellchat")
  mock_autozyme_bindings(state)

  expect_error(
    Shennong:::.sn_with_default_autozyme(
      stop("analysis failure"),
      patches = c("cellchat", "nichenetr")
    ),
    "analysis failure"
  )
  expect_identical(state$status[["cellchat"]], "active")
  expect_identical(state$status[["nichenetr"]], "inactive")
  expect_identical(state$deactivated, "nichenetr")

  state <- make_autozyme_mock_state(fail_deactivate_patch = "cellchat")
  mock_autozyme_bindings(state)
  executions <- 0L
  expect_error(
    Shennong:::.sn_with_default_autozyme(
      {
        executions <- executions + 1L
        42L
      },
      patches = "cellchat"
    ),
    "scope could not restore.*residual active patch.*cellchat"
  )
  expect_identical(executions, 1L)
  expect_identical(state$status[["cellchat"]], "active")
})

test_that("automatic scope restores state after a non-local return", {
  state <- make_autozyme_mock_state()
  mock_autozyme_bindings(state)

  runner <- function() {
    Shennong:::.sn_with_default_autozyme(
      return(42L),
      patches = "cellchat"
    )
    0L
  }

  expect_identical(runner(), 42L)
  expect_identical(state$status[["cellchat"]], "inactive")
  expect_identical(state$deactivated, "cellchat")
})

test_that("automatic provenance includes only patches relevant to the workflow", {
  state <- make_autozyme_mock_state(active = c("cellchat", "nichenetr"))
  mock_autozyme_bindings(state)

  result <- Shennong:::.sn_with_default_autozyme(
    Shennong:::.sn_new_analysis_result(
      analysis_type = "test_result",
      name = "filtered_acceleration",
      method = "test",
      tables = list(primary = tibble::tibble(value = 1))
    ),
    patches = "cellchat"
  )
  expect_identical(
    result$provenance$acceleration$active_patches,
    "cellchat"
  )

  excluded <- Shennong:::.sn_with_default_autozyme(
    Shennong:::.sn_new_analysis_result(
      analysis_type = "test_result",
      name = "excluded_acceleration",
      method = "test",
      tables = list(primary = tibble::tibble(value = 1))
    ),
    patches = character(0)
  )
  expect_null(excluded$provenance$acceleration)
})

test_that("workflow provenance retains used patches after scoped rollback", {
  state <- make_autozyme_mock_state(active = "slingshot")
  mock_autozyme_bindings(state)

  result <- Shennong:::.sn_with_autozyme_provenance_context(
    {
      expect_identical(
        Shennong:::.sn_with_default_autozyme(42L, patches = "tradeseq"),
        42L
      )
      Shennong:::.sn_new_analysis_result(
        analysis_type = "test_result",
        name = "trajectory_acceleration",
        method = "test",
        tables = list(primary = tibble::tibble(value = 1))
      )
    },
    patches = c("slingshot", "tradeseq")
  )

  expect_identical(state$status[["slingshot"]], "active")
  expect_identical(state$status[["tradeseq"]], "inactive")
  expect_identical(
    result$provenance$acceleration$active_patches,
    c("slingshot", "tradeseq")
  )
})

test_that("automatic AutoZyme silently skips missing and version-drifted patches", {
  testthat::local_mocked_bindings(
    .sn_autozyme_is_installed = function(package = "autozyme") FALSE,
    .package = "Shennong"
  )
  expect_identical(
    Shennong:::.sn_with_default_autozyme(41L, "cellchat"),
    41L
  )

  state <- make_autozyme_mock_state(versions = c(CellChat = "2.2.0"))
  mock_autozyme_bindings(state)
  expect_identical(
    Shennong:::.sn_with_default_autozyme(42L, "cellchat"),
    42L
  )
  expect_length(state$activated, 0L)
})

test_that("unsafe backend calls can temporarily suspend active patches", {
  suspended <- FALSE
  testthat::local_mocked_bindings(
    .sn_autozyme_is_installed = function(package = "autozyme") TRUE,
    .sn_autozyme_namespace_loaded = function() TRUE,
    .sn_autozyme_with_disabled = function() {
      function(expr) {
        suspended <<- TRUE
        force(expr)
      }
    },
    .package = "Shennong"
  )

  expect_identical(Shennong:::.sn_with_autozyme_disabled(42L), 42L)
  expect_true(suspended)

  expect_error(
    Shennong:::.sn_with_autozyme_disabled(stop("analysis failure")),
    "analysis failure"
  )
})

test_that("Seurat automatic acceleration bypasses BPCells-backed objects", {
  state <- make_autozyme_mock_state()
  mock_autozyme_bindings(state)
  suspended <- FALSE
  fake_object <- structure(list(), class = "Seurat")
  testthat::local_mocked_bindings(
    .sn_seurat_uses_bpcells = function(object, assay = NULL) TRUE,
    .sn_with_autozyme_disabled = function(expr) {
      suspended <<- TRUE
      force(expr)
    },
    .package = "Shennong"
  )

  expect_identical(
    Shennong:::.sn_with_default_seurat_autozyme(42L, fake_object),
    42L
  )
  expect_true(suspended)
  expect_length(state$activated, 0L)

  suspended <- FALSE
  testthat::local_mocked_bindings(
    .sn_seurat_uses_bpcells = function(object, assay = NULL) FALSE,
    .package = "Shennong"
  )
  expect_identical(
    Shennong:::.sn_with_default_seurat_autozyme(43L, fake_object),
    43L
  )
  expect_false(suspended)
  expect_identical(state$activated, "seurat")
  expect_identical(state$deactivated, "seurat")
  expect_identical(state$status[["seurat"]], "inactive")
})

test_that("BPCells suppression prevents ambient Seurat provenance claims", {
  state <- make_autozyme_mock_state(active = "seurat")
  mock_autozyme_bindings(state)
  fake_object <- structure(list(), class = "Seurat")
  testthat::local_mocked_bindings(
    .sn_seurat_uses_bpcells = function(object, assay = NULL) TRUE,
    .sn_with_autozyme_disabled = function(expr) force(expr),
    .package = "Shennong"
  )

  result <- Shennong:::.sn_with_autozyme_provenance_context(
    {
      expect_identical(
        Shennong:::.sn_with_default_seurat_autozyme(42L, fake_object),
        42L
      )
      Shennong:::.sn_new_analysis_result(
        analysis_type = "test_result",
        name = "bpcells_without_seurat_acceleration",
        method = "test",
        tables = list(primary = tibble::tibble(value = 1))
      )
    },
    patches = "seurat"
  )

  expect_identical(state$status[["seurat"]], "active")
  expect_null(result$provenance$acceleration)
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
