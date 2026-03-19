library(testthat)

test_that("codex_skill_path points to the bundled skill directory", {
  path <- sn_get_codex_skill_path()

  expect_true(dir.exists(path))
  expect_true(file.exists(file.path(path, "SKILL.md")))
})

test_that("install_codex_skill installs the bundled skill into a target directory", {
  target_root <- tempfile("skill-root-")
  dir.create(target_root)

  installed_path <- sn_install_codex_skill(path = target_root, overwrite = TRUE)

  expect_true(dir.exists(installed_path))
  expect_true(file.exists(file.path(installed_path, "SKILL.md")))
  expect_true(file.exists(file.path(installed_path, "references", "package_overview.md")))
})

test_that("sn_install_codex_skill refuses to overwrite when disabled", {
  target_root <- tempfile("skill-root-")
  dir.create(target_root)
  dir.create(file.path(target_root, "shennong"))

  expect_error(
    sn_install_codex_skill(path = target_root, overwrite = FALSE),
    "overwrite = TRUE"
  )
})
