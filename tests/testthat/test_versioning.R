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
