python_autozyme_environments <- c(
  "bbknn", "cell2location", "cellphonedb", "scarches", "stlearn", "trajectory"
)

test_that("managed Python acceleration pins the same AutoZyme fork revision", {
  expected <- paste0(
    'autozyme = { git = "https://github.com/zerostwo/autozyme.git", ',
    'rev = "8fc2e9c3a7f70302f97589aaa9b0395dcf86f9bc", ',
    'subdirectory = "autozyme_py" }'
  )
  for (environment in python_autozyme_environments) {
    lines <- readLines(sn_pixi_config_path(environment), warn = FALSE)
    expect_identical(grep("^autozyme[[:space:]]*=", lines, value = TRUE), expected)
  }

  scanpy_environments <- c("bbknn", "cell2location", "scarches", "stlearn", "trajectory")
  for (environment in scanpy_environments) {
    lines <- readLines(sn_pixi_config_path(environment), warn = FALSE)
    expect_true('scanpy = "==1.11.5"' %in% lines)
  }
  cell2location <- readLines(sn_pixi_config_path("cell2location"), warn = FALSE)
  expect_true('cell2location = "==0.1.5"' %in% cell2location)
  expect_true('scvi-tools = "==1.3.3"' %in% cell2location)
  expect_true('pyro-ppl = "==1.9.1"' %in% cell2location)
  expect_true('cellphonedb = "==5.0.1"' %in%
    readLines(sn_pixi_config_path("cellphonedb"), warn = FALSE))
  trajectory <- readLines(sn_pixi_config_path("trajectory"), warn = FALSE)
  expect_true(any(grepl("f63c0e70596ced2f1bee8cf07e8ab66037cf86b2", trajectory, fixed = TRUE)))
})

test_that("only Python backends with a matching patch request AutoZyme", {
  unaffected <- c(
    "infercnvpy", "mmochi", "scib-metrics", "scvi", "spatialdata", "squidpy",
    "tangram"
  )
  for (environment in unaffected) {
    lines <- readLines(sn_pixi_config_path(environment), warn = FALSE)
    expect_false(any(grepl("^autozyme[[:space:]]*=", lines)))
  }

  scripts <- c(
    bbknn = "bbknn_integration.py",
    cell2location = "cell2location_run.py",
    cellphonedb = "cellphonedb_run.py",
    scarches = "scarches_run.py",
    stlearn = "stlearn_run.py",
    trajectory = "trajectory_run.py"
  )
  for (environment in names(scripts)) {
    script <- system.file(
      "pixi", environment, "scripts", scripts[[environment]],
      package = "Shennong"
    )
    lines <- readLines(script, warn = FALSE)
    expect_true(any(grepl("shennong_autozyme", lines, fixed = TRUE)))
    expect_true(any(grepl('"autozyme"', lines, fixed = TRUE)))
  }
})

test_that("existing Shennong-managed manifests gain the pinned dependency once", {
  manifest <- tempfile(fileext = ".toml")
  writeLines(
    c(
      "[workspace]",
      'name = "shennong-trajectory"',
      "",
      "[pypi-dependencies]",
      'regvelo = ">=0.1.0"',
      "",
      "[environments]",
      "default = []"
    ),
    manifest
  )
  expect_true(Shennong:::.sn_refresh_managed_pixi_autozyme(manifest, "trajectory"))
  refreshed <- readLines(manifest, warn = FALSE)
  expect_length(grep("^autozyme[[:space:]]*=", refreshed), 1L)
  expect_false(Shennong:::.sn_refresh_managed_pixi_autozyme(manifest, "trajectory"))
  expect_length(grep("^autozyme[[:space:]]*=", readLines(manifest, warn = FALSE)), 1L)

  custom <- tempfile(fileext = ".toml")
  writeLines(c("[workspace]", 'name = "custom"', "[pypi-dependencies]"), custom)
  expect_false(Shennong:::.sn_refresh_managed_pixi_autozyme(custom, "trajectory"))
})

test_that("bundled Python AutoZyme helpers and runners compile", {
  python <- Sys.which("python3")
  skip_if(!nzchar(python), "python3 is not available")
  files <- c(
    system.file("pixi", "_shared", "shennong_autozyme.py", package = "Shennong"),
    system.file("pixi", "_shared", "sitecustomize.py", package = "Shennong"),
    vapply(names(c(
      bbknn = "bbknn_integration.py",
      cell2location = "cell2location_run.py",
      cellphonedb = "cellphonedb_run.py",
      scarches = "scarches_run.py",
      stlearn = "stlearn_run.py",
      trajectory = "trajectory_run.py"
    )), function(environment) {
      script <- c(
        bbknn = "bbknn_integration.py",
        cell2location = "cell2location_run.py",
        cellphonedb = "cellphonedb_run.py",
        scarches = "scarches_run.py",
        stlearn = "stlearn_run.py",
        trajectory = "trajectory_run.py"
      )[[environment]]
      system.file("pixi", environment, "scripts", script, package = "Shennong")
    }, character(1))
  )
  status <- system2(python, c("-m", "py_compile", files), stdout = TRUE, stderr = TRUE)
  expect_identical(attr(status, "status") %||% 0L, 0L, info = paste(status, collapse = "\n"))
})

test_that("Python activation is fail-closed on unvalidated versions and unsafe cell2location modes", {
  helper <- readLines(
    system.file("pixi", "_shared", "shennong_autozyme.py", package = "Shennong"),
    warn = FALSE
  )
  expect_true(any(grepl("VALIDATED_UPSTREAM_VERSIONS", helper, fixed = TRUE)))
  expect_true(any(grepl("compatibility_errors", helper, fixed = TRUE)))
  expect_true(any(grepl('if report["disabled"]:', helper, fixed = TRUE)))

  runner <- readLines(
    system.file(
      "pixi", "cell2location", "scripts", "cell2location_run.py",
      package = "Shennong"
    ),
    warn = FALSE
  )
  expect_true(any(grepl("safe_autozyme_scope", runner, fixed = TRUE)))
  expect_true(any(grepl("torch.cuda.is_available", runner, fixed = TRUE)))
  expect_true(any(grepl('config.get("batch_size") is None', runner, fixed = TRUE)))
})
