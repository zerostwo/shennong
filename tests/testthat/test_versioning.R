library(testthat)

test_that("release channel resolution prefers CRAN for auto when available", {
  expect_equal(
    Shennong:::.sn_resolve_release_channel(
      channel = "auto",
      cran_version = package_version("1.0.0"),
      github_version = package_version("1.1.0")
    ),
    "cran"
  )

  expect_equal(
    Shennong:::.sn_resolve_release_channel(
      channel = "auto",
      cran_version = NULL,
      github_version = package_version("1.1.0")
    ),
    "github"
  )
})

test_that("release channel resolution errors when requested channel is unavailable", {
  expect_error(
    Shennong:::.sn_resolve_release_channel(
      channel = "cran",
      cran_version = NULL,
      github_version = package_version("1.1.0")
    ),
    "not currently available on CRAN"
  )

  expect_error(
    Shennong:::.sn_resolve_release_channel(
      channel = "github",
      cran_version = package_version("1.0.0"),
      github_version = NULL
    ),
    "Could not retrieve the GitHub development version"
  )
})

test_that("version comparison reports update states correctly", {
  expect_equal(
    Shennong:::.sn_compare_version_status(
      installed_version = package_version("1.0.0"),
      remote_version = package_version("1.0.1")
    )$status,
    "update available"
  )

  expect_equal(
    Shennong:::.sn_compare_version_status(
      installed_version = package_version("1.0.1"),
      remote_version = package_version("1.0.1")
    )$status,
    "up to date"
  )

  expect_equal(
    Shennong:::.sn_compare_version_status(
      installed_version = package_version("1.0.2"),
      remote_version = package_version("1.0.1")
    )$status,
    "ahead of remote"
  )
})

test_that("sn_install_shennong uses conservative defaults for GitHub installs", {
  local_mocked_bindings(
    .sn_get_cran_version = function(...) NULL,
    .sn_get_github_version = function(...) package_version("1.0.0"),
    check_installed = function(...) invisible(TRUE),
    .env = asNamespace("Shennong")
  )

  captured <- NULL
  local_mocked_bindings(
    .sn_install_github_release = function(repo, ref, args = list()) {
      captured <<- c(list(repo = repo, ref = ref), args)
      invisible(TRUE)
    },
    .env = asNamespace("Shennong")
  )

  expect_invisible(
    Shennong::sn_install_shennong(channel = "github")
  )

  expect_equal(captured$repo, "zerostwo/shennong")
  expect_equal(captured$ref, "main")
  expect_false(isTRUE(captured$dependencies))
  expect_equal(captured$upgrade, "never")
})

test_that("sn_install_shennong respects explicit GitHub install overrides", {
  local_mocked_bindings(
    .sn_get_cran_version = function(...) NULL,
    .sn_get_github_version = function(...) package_version("1.0.0"),
    check_installed = function(...) invisible(TRUE),
    .env = asNamespace("Shennong")
  )

  captured <- NULL
  local_mocked_bindings(
    .sn_install_github_release = function(repo, ref, args = list()) {
      captured <<- c(list(repo = repo, ref = ref), args)
      invisible(TRUE)
    },
    .env = asNamespace("Shennong")
  )

  expect_invisible(
    Shennong::sn_install_shennong(
      channel = "github",
      source = "acme/shennong",
      ref = "dev",
      dependencies = TRUE,
      upgrade = "always"
    )
  )

  expect_equal(captured$repo, "acme/shennong")
  expect_equal(captured$ref, "dev")
  expect_true(isTRUE(captured$dependencies))
  expect_equal(captured$upgrade, "always")
})

test_that("sn_install_shennong supports local installs", {
  local_mocked_bindings(
    check_installed = function(...) invisible(TRUE),
    .env = asNamespace("Shennong")
  )

  captured <- NULL
  local_mocked_bindings(
    .sn_install_local_release = function(path, args = list()) {
      captured <<- list(path = path, args = args)
      invisible(TRUE)
    },
    .env = asNamespace("Shennong")
  )

  expect_invisible(
    Shennong::sn_install_shennong(
      channel = "local",
      source = "/tmp/Shennong"
    )
  )

  expect_equal(captured$path, "/tmp/Shennong")
})

test_that("sn_install_shennong rejects conflicting source arguments", {
  expect_error(
    Shennong::sn_install_shennong(
      channel = "github",
      source = "acme/shennong",
      github_repo = "zerostwo/shennong"
    ),
    "conflicting values"
  )

  expect_error(
    Shennong::sn_install_shennong(
      channel = "local",
      source = "/tmp/new",
      local_path = "/tmp/old"
    ),
    "conflicting values"
  )

  expect_error(
    Shennong::sn_install_shennong(
      channel = "github",
      ref = "dev",
      github_ref = "main-release"
    ),
    "conflicting values"
  )
})
