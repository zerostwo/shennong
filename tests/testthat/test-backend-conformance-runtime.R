test_that("optional Python probes honor the configured runtime without a personal directory", {
  root <- tempfile("conformance-runtime-")
  withr::local_options(list(shennong.runtime_dir = root))
  withr::defer(unlink(root, recursive = TRUE))
  withr::local_envvar(SHENNONG_CONFORMANCE_STRICT = "true")
  expect_false(.conformance_pixi_env_installed("popv"))
  expect_error(.conformance_require_pixi_environment("popv"), "configured Shennong runtime")
  python <- .conformance_pixi_python("popv")
  expect_true(startsWith(python, normalizePath(root, winslash = "/")))
  dir.create(dirname(python), recursive = TRUE)
  file.create(python)
  expect_true(.conformance_pixi_env_installed("popv"))
  expect_true(.conformance_require_pixi_environment("popv"))
})

test_that("conformance version checks distinguish R and Python dependencies", {
  expect_identical(
    .conformance_upstream_version(list(upstream = list(language = "r")), "stats"),
    as.character(utils::packageVersion("stats"))
  )
  root <- tempfile("conformance-version-runtime-")
  withr::local_options(list(shennong.runtime_dir = root))
  withr::defer(unlink(root, recursive = TRUE))
  withr::local_envvar(SHENNONG_CONFORMANCE_STRICT = "true")
  expect_error(.conformance_upstream_version(.conformance_contract("annotation::popv"), "popv"),
               "popv pixi environment is not installed")
  expect_error(.conformance_upstream_version(list(upstream = list(language = "other")), "stats"),
               "Unsupported conformance language")
})
