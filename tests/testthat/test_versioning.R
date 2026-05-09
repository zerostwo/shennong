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

test_that("sn_list_dependencies reports declared package metadata", {
  deps <- Shennong::sn_list_dependencies()

  expect_s3_class(deps, "tbl_df")
  expect_true(all(c(
    "package", "requirement", "declared_in", "source",
    "remote", "installed", "version"
  ) %in% colnames(deps)))
  expect_true("cluster" %in% deps$package)
  expect_true("Seurat" %in% deps$package)
  expect_equal(
    deps$requirement[deps$package == "cluster"][[1]],
    "required"
  )
  expect_equal(
    deps$declared_in[deps$package == "Seurat"][[1]],
    "Suggests"
  )
  expect_equal(
    deps$source[deps$package == "harmony"][[1]],
    "GitHub"
  )
  expect_equal(
    deps$source[deps$package == "anndataR"][[1]],
    "GitHub"
  )
  expect_equal(
    deps$remote[deps$package == "anndataR"][[1]],
    "scverse/anndataR"
  )
  expect_equal(
    deps$source[deps$package == "tidytemplate"][[1]],
    "GitHub"
  )
  expect_equal(
    deps$remote[deps$package == "tidytemplate"][[1]],
    "tidyverse/tidytemplate"
  )
  expect_equal(
    deps$source[deps$package == "Nebulosa"][[1]],
    "Bioconductor"
  )
  expect_false("qs" %in% deps$package)
  expect_true("qs2" %in% deps$package)
})

test_that("sn_install_dependencies dispatches installs by declared source", {
  fake_deps <- data.frame(
    package = c("cli", "scran", "harmony"),
    requirement = c("required", "recommended", "recommended"),
    declared_in = c("Imports", "Suggests", "Suggests"),
    source = c("CRAN", "Bioconductor", "GitHub"),
    remote = c(NA_character_, NA_character_, "immunogenomics/harmony@harmony2"),
    installed = c(FALSE, FALSE, FALSE),
    version = c(NA_character_, NA_character_, NA_character_),
    stringsAsFactors = FALSE
  )

  captured <- list(cran = NULL, bioc = NULL, github = NULL)
  local_mocked_bindings(
    .sn_dependency_table = function() fake_deps,
    .sn_install_cran_packages = function(packages, repos = getOption("repos"), ...) {
      captured$cran <<- packages
      invisible(packages)
    },
    .sn_install_bioc_packages = function(packages, ask = interactive(), update = FALSE, repos = getOption("repos"), ...) {
      captured$bioc <<- packages
      invisible(packages)
    },
    .sn_install_github_packages = function(remotes,
                                           upgrade = FALSE,
                                           repos = getOption("repos"),
                                           dependencies = NA,
                                           ...) {
      captured$github <<- remotes
      captured$github_dependencies <<- dependencies
      invisible(remotes)
    },
    .sn_find_missing_packages = function(packages) character(0),
    .env = asNamespace("Shennong")
  )

  expect_invisible(
    Shennong::sn_install_dependencies()
  )

  expect_equal(captured$cran, "cli")
  expect_equal(captured$bioc, "scran")
  expect_equal(captured$github, "immunogenomics/harmony@harmony2")
  expect_equal(captured$github_dependencies, NA)
})

test_that("sn_install_dependencies errors when installer leaves packages missing", {
  fake_deps <- data.frame(
    package = "cli",
    requirement = "required",
    declared_in = "Imports",
    source = "CRAN",
    remote = NA_character_,
    installed = FALSE,
    version = NA_character_,
    stringsAsFactors = FALSE
  )

  local_mocked_bindings(
    .sn_dependency_table = function() fake_deps,
    .sn_install_cran_packages = function(packages, repos = getOption("repos"), ...) {
      invisible(packages)
    },
    .sn_install_bioc_packages = function(packages, ask = interactive(), update = FALSE, repos = getOption("repos"), ...) {
      invisible(packages)
    },
    .sn_install_github_packages = function(remotes, upgrade = FALSE, repos = getOption("repos"), dependencies = NA, ...) {
      invisible(remotes)
    },
    .sn_find_missing_packages = function(packages) "cli",
    .env = asNamespace("Shennong")
  )

  expect_error(
    Shennong::sn_install_dependencies(),
    "Failed to install package"
  )
})

test_that("sn_install_dependencies validates requested package names", {
  local_mocked_bindings(
    .sn_dependency_table = function() {
      data.frame(
        package = "cli",
        requirement = "required",
        declared_in = "Imports",
        source = "CRAN",
        remote = NA_character_,
        installed = FALSE,
        version = NA_character_,
        stringsAsFactors = FALSE
      )
    },
    .env = asNamespace("Shennong")
  )

  expect_error(
    Shennong::sn_install_dependencies(packages = "not_real_pkg"),
    "Unknown package"
  )
})

test_that("pixi helpers detect executables, expose runtime paths, and write mirror config", {
  old_home <- Sys.getenv("HOME", unset = NA_character_)
  old_path <- Sys.getenv("PATH", unset = NA_character_)
  old_env <- Sys.getenv("SHENNONG_PIXI", unset = NA_character_)
  old_option <- getOption("shennong.pixi", NULL)
  fake_home <- tempfile("pixi-home-empty-")
  fake_path <- tempfile("pixi-path-empty-")
  dir.create(fake_home, recursive = TRUE)
  dir.create(fake_path, recursive = TRUE)
  on.exit({
    if (is.na(old_home)) Sys.unsetenv("HOME") else Sys.setenv(HOME = old_home)
    if (is.na(old_path)) Sys.unsetenv("PATH") else Sys.setenv(PATH = old_path)
    if (is.na(old_env)) Sys.unsetenv("SHENNONG_PIXI") else Sys.setenv(SHENNONG_PIXI = old_env)
    options(shennong.pixi = old_option)
  }, add = TRUE)

  Sys.setenv(HOME = fake_home, PATH = fake_path, SHENNONG_PIXI = "")
  options(shennong.pixi = NULL)
  missing_info <- sn_check_pixi(quiet = TRUE)
  expect_false(missing_info$installed)
  expect_true(is.na(missing_info$path))

  fake_dir <- tempfile("pixi-bin-")
  dir.create(fake_dir, recursive = TRUE)
  fake_pixi <- file.path(fake_dir, if (.Platform$OS.type == "windows") "pixi.exe" else "pixi")
  writeLines(c("#!/bin/sh", "echo 'pixi 0.99.0'"), fake_pixi)
  Sys.chmod(fake_pixi, mode = "755")

  info <- sn_check_pixi(pixi = fake_pixi, quiet = TRUE)
  expect_true(info$installed)
  expect_equal(info$version, "0.99.0")

  ensured <- sn_ensure_pixi(pixi = fake_pixi, install = FALSE, quiet = TRUE)
  expect_equal(ensured$path, fake_pixi)

  runtime_dir <- tempfile("shennong-home-")
  paths <- sn_pixi_paths("scanvi", runtime_dir = runtime_dir)
  expect_equal(paths$environment, "scanvi")
  expect_equal(paths$family, "scvi")
  expect_true(grepl("/\\.shennong|shennong-home-", paths$runtime_dir))
  expect_true(grepl("/pixi/scvi/pixi\\.toml$", paths$manifest_path))
  expect_true(grepl("/pixi/home$", paths$pixi_home))
  expect_true(file.exists(paths$source_config_path))
  expect_true("scvi" %in% sn_list_pixi_environments())
  expect_true("scarches" %in% sn_list_pixi_environments())
  expect_true("infercnvpy" %in% sn_list_pixi_environments())
  expect_true("cell2location" %in% sn_list_pixi_environments())
  expect_true("tangram" %in% sn_list_pixi_environments())
  expect_true("squidpy" %in% sn_list_pixi_environments())
  expect_true("spatialdata" %in% sn_list_pixi_environments())
  expect_true("stlearn" %in% sn_list_pixi_environments())
  expect_true("mmochi" %in% sn_list_pixi_environments())
  expect_false("scanvi" %in% sn_list_pixi_environments())
  expect_false("scpoli" %in% sn_list_pixi_environments())
  expect_false("spatial" %in% sn_list_pixi_environments())
  expect_true(file.exists(sn_pixi_config_path("cellphonedb")))
  expect_equal(dirname(sn_pixi_config_path("scanvi")), dirname(sn_pixi_config_path("scvi")))
  expect_equal(dirname(sn_pixi_config_path("scpoli")), dirname(sn_pixi_config_path("scarches")))
  expect_equal(dirname(sn_pixi_config_path("tarngram")), dirname(sn_pixi_config_path("tangram")))
  expect_equal(basename(dirname(sn_pixi_config_path("mmochi-landmark"))), "mmochi")

  prepared <- sn_prepare_pixi_environment(
    "scpoli",
    runtime_dir = runtime_dir,
    pixi_environment = "cpu",
    install_environment = FALSE,
    overwrite = TRUE,
    platforms = "linux-64"
  )
  expect_true(file.exists(prepared$manifest_path))
  expect_equal(prepared$family, "scarches")
  prepared_manifest <- readLines(prepared$manifest_path)
  expect_true(any(grepl('name = "shennong-scarches"', prepared_manifest, fixed = TRUE)))
  expect_true(any(grepl('platforms = ["linux-64"]', prepared_manifest, fixed = TRUE)))

  pixi_home <- tempfile("pixi-home-")
  config_path <- sn_configure_pixi_mirror("tuna", pixi_home = pixi_home)
  config <- readLines(config_path)
  expect_true(any(grepl("\\[mirrors\\]", config)))
  expect_true(any(grepl("mirrors.tuna.tsinghua.edu.cn", config, fixed = TRUE)))
  expect_true(any(grepl("pypi.tuna.tsinghua.edu.cn", config, fixed = TRUE)))
  expect_true(any(grepl("files.pythonhosted.org/packages", config, fixed = TRUE)))
})

test_that("scVI pixi manifest includes CPU and GPU environments", {
  manifest <- Shennong:::.sn_scvi_pixi_manifest_lines(cuda_version = "12.6", platforms = "linux-64")

  expect_true(any(grepl('platforms = \\["linux-64"\\]', manifest)))
  expect_true(any(grepl("\\[feature.cpu.dependencies\\]", manifest)))
  expect_true(any(grepl("pytorch-cpu", manifest, fixed = TRUE)))
  expect_true(any(grepl("\\[feature.gpu.system-requirements\\]", manifest)))
  expect_true(any(grepl('cuda = "12"', manifest, fixed = TRUE)))
  expect_true(any(grepl("pytorch-gpu", manifest, fixed = TRUE)))
  expect_true(any(grepl('cuda-version = "12.6.*"', manifest, fixed = TRUE)))
  expect_true(any(grepl("gpu = \\[\"gpu\"\\]", manifest)))
  expect_equal(Shennong:::.sn_normalize_cuda_requirement("12.6"), "12.6")
  expect_equal(Shennong:::.sn_normalize_cuda_requirement("12"), "12.0")
  expect_equal(Shennong:::.sn_default_scvi_cuda_version("13.0"), "12.6")
  expect_equal(Shennong:::.sn_default_scvi_cuda_version("11.8"), "11.8")
})
