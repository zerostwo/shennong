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

test_that("version checks use current source/ref arguments", {
  expect_true(all(c("source", "ref") %in% names(formals(Shennong::sn_check_version))))
  expect_false(any(c("github_repo", "github_ref") %in% names(formals(Shennong::sn_check_version))))
})

test_that("sn_install_shennong uses conservative defaults for GitHub installs", {
  local_mocked_bindings(
    .sn_get_cran_version = function(...) NULL,
    .sn_get_github_version = function(...) package_version("1.0.0"),
    check_installed = function(...) invisible(TRUE),
    .package = "Shennong"
  )

  captured <- NULL
  local_mocked_bindings(
    .sn_install_github_release = function(repo, ref, args = list()) {
      captured <<- c(list(repo = repo, ref = ref), args)
      invisible(TRUE)
    },
    .package = "Shennong"
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
    .package = "Shennong"
  )

  captured <- NULL
  local_mocked_bindings(
    .sn_install_github_release = function(repo, ref, args = list()) {
      captured <<- c(list(repo = repo, ref = ref), args)
      invisible(TRUE)
    },
    .package = "Shennong"
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
    .package = "Shennong"
  )

  captured <- NULL
  local_mocked_bindings(
    .sn_install_local_release = function(path, args = list()) {
      captured <<- list(path = path, args = args)
      invisible(TRUE)
    },
    .package = "Shennong"
  )

  expect_invisible(
    Shennong::sn_install_shennong(
      channel = "local",
      source = "/tmp/Shennong"
    )
  )

  expect_equal(captured$path, "/tmp/Shennong")
})

test_that("sn_install_shennong auto installs an explicit local source without remote probes", {
  source_dir <- withr::local_tempdir()
  writeLines(c("Package: Shennong", "Version: 0.3.0.9000"), file.path(source_dir, "DESCRIPTION"))

  captured <- NULL
  local_mocked_bindings(
    .sn_get_cran_version = function(...) stop("CRAN should not be queried"),
    .sn_get_github_version = function(...) stop("GitHub should not be queried"),
    check_installed = function(...) invisible(TRUE),
    .sn_install_local_release = function(path, args = list()) {
      captured <<- list(path = path, args = args)
      invisible(TRUE)
    },
    .package = "Shennong"
  )

  result <- Shennong::sn_install_shennong(source = source_dir)

  expect_identical(result, "local")
  expect_equal(captured$path, normalizePath(source_dir))
})

test_that("sn_install_shennong auto falls back to the current source tree", {
  source_dir <- withr::local_tempdir()
  writeLines(c("Package: Shennong", "Version: 0.3.0.9000"), file.path(source_dir, "DESCRIPTION"))
  withr::local_dir(source_dir)

  captured <- NULL
  local_mocked_bindings(
    .sn_get_cran_version = function(...) NULL,
    .sn_get_github_version = function(...) NULL,
    check_installed = function(...) invisible(TRUE),
    .sn_install_local_release = function(path, args = list()) {
      captured <<- list(path = path, args = args)
      invisible(TRUE)
    },
    .package = "Shennong"
  )

  expect_warning(
    result <- Shennong::sn_install_shennong(),
    "falling back to the local Shennong source"
  )

  expect_identical(result, "local")
  expect_equal(captured$path, normalizePath(source_dir))
})

test_that("sn_install_shennong auto still errors without any available source", {
  withr::local_dir(withr::local_tempdir())
  local_mocked_bindings(
    .sn_get_cran_version = function(...) NULL,
    .sn_get_github_version = function(...) NULL,
    .package = "Shennong"
  )

  expect_error(
    Shennong::sn_install_shennong(),
    "Could not determine a remote version"
  )
})

test_that("sn_install_shennong requires source for local installs", {
  expect_error(
    Shennong::sn_install_shennong(
      channel = "local"
    ),
    "`source` must be supplied"
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
  expected_github <- c(
    ShennongOpt = "zerostwo/shennong-opt",
    monocle3 = "cole-trapnell-lab/monocle3",
    RareQ = "xiaolab-xjtu/RareQ",
    Scissor = "sunduanchen/Scissor"
  )
  for (package in names(expected_github)) {
    expect_identical(deps$source[deps$package == package][[1]], "GitHub")
    expect_identical(deps$remote[deps$package == package][[1]], expected_github[[package]])
  }
  expect_false("tidytemplate" %in% deps$package)
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
    .package = "Shennong"
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
    .package = "Shennong"
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
    .package = "Shennong"
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

  ensured <- sn_ensure_pixi(
    pixi = fake_pixi, install = FALSE, version = "0.99.0", quiet = TRUE
  )
  expect_equal(ensured$path, fake_pixi)
  expect_error(
    sn_ensure_pixi(
      pixi = fake_pixi, install = FALSE, version = "0.69.0", quiet = TRUE
    ),
    "does not match required version"
  )

  runtime_dir <- tempfile("shennong-home-")
  paths <- sn_get_pixi_paths("scanvi", runtime_dir = runtime_dir)
  expect_equal(paths$environment, "scanvi")
  expect_equal(paths$family, "scvi")
  expect_true(grepl("/\\.shennong|shennong-home-", paths$runtime_dir))
  expect_true(grepl("/pixi/scvi/pixi\\.toml$", paths$manifest_path))
  expect_true(grepl("/pixi/home$", paths$pixi_home))
  expect_true(file.exists(paths$source_config_path))
  expect_true("scvi" %in% sn_list_pixi_environments())
  expect_true("scarches" %in% sn_list_pixi_environments())
  expect_true("bbknn" %in% sn_list_pixi_environments())
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
  expect_true(file.exists(sn_get_pixi_config_path("cellphonedb")))
  expect_equal(dirname(sn_get_pixi_config_path("scanvi")), dirname(sn_get_pixi_config_path("scvi")))
  expect_equal(dirname(sn_get_pixi_config_path("scpoli")), dirname(sn_get_pixi_config_path("scarches")))
  expect_true(file.exists(sn_get_pixi_config_path("bbknn")))
  expect_equal(dirname(sn_get_pixi_config_path("tarngram")), dirname(sn_get_pixi_config_path("tangram")))
  expect_equal(basename(dirname(sn_get_pixi_config_path("mmochi-landmark"))), "mmochi")

  prepared <- sn_prepare_pixi_environment(
    "scpoli",
    runtime_dir = runtime_dir,
    pixi_environment = "cpu",
    install_environment = FALSE,
    overwrite = TRUE,
    platforms = "linux-64"
  )
  expect_true(file.exists(prepared$manifest_path))
  expect_true(file.exists(prepared$lock_path))
  expect_equal(prepared$family, "scarches")
  prepared_manifest <- readLines(prepared$manifest_path)
  expect_true(any(grepl('name = "shennong-scarches"', prepared_manifest, fixed = TRUE)))
  expect_true(any(grepl(
    'platforms = ["linux-64", "osx-64", "osx-arm64", "win-64"]',
    prepared_manifest,
    fixed = TRUE
  )))
  expect_true("linux-64" %in% Shennong:::.sn_pixi_lock_platforms(prepared$lock_path))

  python_backends <- c(
    system.file("pixi", "scvi", "scripts", "scvi_integration.py", package = "Shennong"),
    system.file("pixi", "scarches", "scripts", "scpoli_integration.py", package = "Shennong")
  )
  expect_true(all(nzchar(python_backends) & file.exists(python_backends)))
  backend_source <- paste(unlist(lapply(python_backends, readLines, warn = FALSE)), collapse = "\n")
  expect_false(grepl("adata.X.toarray()", backend_source, fixed = TRUE))
  expect_false(grepl("protein_counts.toarray()", backend_source, fixed = TRUE))
  expect_true(grepl("adata.X[start:stop, :]", backend_source, fixed = TRUE))

  pixi_home <- tempfile("pixi-home-")
  config_path <- sn_configure_pixi_mirror("tuna", pixi_home = pixi_home)
  config <- readLines(config_path)
  expect_true(any(grepl("\\[mirrors\\]", config)))
  expect_true(any(grepl("mirrors.tuna.tsinghua.edu.cn", config, fixed = TRUE)))
  expect_true(any(grepl("pypi.tuna.tsinghua.edu.cn", config, fixed = TRUE)))
  expect_true(any(grepl("files.pythonhosted.org/packages", config, fixed = TRUE)))
})

test_that("Pixi bootstrap uses explicit versioned and verified release assets", {
  linux <- Shennong:::.sn_pixi_release_asset(
    "0.69.0", sysname = "Linux", machine = "x86_64"
  )
  mac <- Shennong:::.sn_pixi_release_asset(
    "v0.69.0", sysname = "Darwin", machine = "arm64"
  )
  windows <- Shennong:::.sn_pixi_release_asset(
    "0.69.0", sysname = "Windows", machine = "AMD64"
  )
  windows_x86_hyphen <- Shennong:::.sn_pixi_release_asset(
    "0.69.0", sysname = "Windows", machine = "x86-64"
  )
  expect_match(linux$url, "/v0.69.0/pixi-x86_64-unknown-linux-musl\\.tar\\.gz$")
  expect_match(mac$asset, "pixi-aarch64-apple-darwin\\.tar\\.gz$")
  expect_match(windows$asset, "pixi-x86_64-pc-windows-msvc\\.zip$")
  expect_identical(windows_x86_hyphen$asset, windows$asset)
  expect_error(Shennong:::.sn_pixi_release_asset("latest"), "explicit semantic")
  expect_true(Shennong:::.sn_pixi_version_matches("pixi 0.69.0", "0.69.0"))
  expect_false(Shennong:::.sn_pixi_version_matches("0.68.0", "0.69.0"))

  payload <- tempfile("pixi-digest-")
  writeBin(charToRaw("verified payload"), payload)
  expected <- digest::digest(payload, algo = "sha256", file = TRUE, serialize = FALSE)
  expect_invisible(Shennong:::.sn_verify_file_sha256(payload, expected))
  expect_error(
    Shennong:::.sn_verify_file_sha256(payload, paste(rep("0", 64), collapse = "")),
    "mismatch"
  )
})

test_that("Pixi environment helpers forward custom download SHA-256", {
  runtime_dir <- tempfile("pixi-sha-runtime-")
  withr::defer(unlink(runtime_dir, recursive = TRUE, force = TRUE))
  custom_url <- "https://mirror.example/pixi-x86_64.tar.gz"
  custom_sha <- paste(rep("a", 64L), collapse = "")
  ensure_calls <- list()

  output <- testthat::with_mocked_bindings(
    sn_call_pixi_environment(
      environment = "scvi",
      command = "python",
      args = "--version",
      runtime_dir = runtime_dir,
      overwrite = TRUE,
      platforms = "linux-64",
      install_pixi = TRUE,
      pixi_download_url = custom_url,
      pixi_sha256 = custom_sha,
      quiet = TRUE
    ),
    sn_ensure_pixi = function(download_url = NULL, sha256 = NULL, ...) {
      ensure_calls[[length(ensure_calls) + 1L]] <<- list(
        download_url = download_url,
        sha256 = sha256
      )
      list(path = "/mock/pixi")
    },
    .sn_pixi_run_command = function(...) "pixi 0.69.0",
    .package = "Shennong"
  )

  expect_identical(output, "pixi 0.69.0")
  expect_length(ensure_calls, 2L)
  expect_true(all(vapply(
    ensure_calls,
    function(call) identical(call$download_url, custom_url),
    logical(1)
  )))
  expect_true(all(vapply(
    ensure_calls,
    function(call) identical(call$sha256, custom_sha),
    logical(1)
  )))
  expect_true("pixi_sha256" %in% names(formals(sn_prepare_pixi_environment)))
  expect_true("pixi_sha256" %in% names(formals(sn_call_pixi_environment)))
})

test_that("scVI pixi manifest includes CPU and GPU environments", {
  manifest <- Shennong:::.sn_scvi_pixi_manifest_lines(cuda_version = "12.6", platforms = "linux-64")

  expect_true(any(grepl(
    'platforms = ["linux-64", "osx-64", "osx-arm64", "win-64"]',
    manifest,
    fixed = TRUE
  )))
  expect_true(any(grepl("\\[feature.cpu.dependencies\\]", manifest)))
  expect_true(any(grepl("pytorch-cpu", manifest, fixed = TRUE)))
  expect_true(any(grepl("\\[feature.gpu.system-requirements\\]", manifest)))
  expect_true(any(grepl('cuda = "12"', manifest, fixed = TRUE)))
  expect_true(any(grepl("pytorch-gpu", manifest, fixed = TRUE)))
  expect_true(any(grepl('cuda-version = "==12.6"', manifest, fixed = TRUE)))
  expect_true(any(grepl("gpu = \\[\"gpu\"\\]", manifest)))
  expect_equal(Shennong:::.sn_normalize_cuda_requirement("12.6"), "12.6")
  expect_equal(Shennong:::.sn_normalize_cuda_requirement("12"), "12.0")
  expect_equal(Shennong:::.sn_default_scvi_cuda_version("13.0"), "12.6")
  expect_equal(Shennong:::.sn_default_scvi_cuda_version("11.8"), "11.8")
})

test_that("integration Pixi projects always materialize a matching lock", {
  project <- tempfile("scvi-locked-project-")
  manifest_path <- Shennong:::.sn_prepare_scvi_pixi_project(
    project_dir = project,
    environment = "scvi",
    overwrite = TRUE,
    platforms = "linux-64"
  )
  lock_path <- file.path(dirname(manifest_path), "pixi.lock")
  expect_true(file.exists(lock_path))
  expect_true("linux-64" %in% Shennong:::.sn_pixi_lock_platforms(lock_path))

  custom <- tempfile("scvi-custom-project-")
  expect_error(
    Shennong:::.sn_prepare_scvi_pixi_project(
      project_dir = custom,
      environment = "scvi",
      manifest_lines = c(
        "[workspace]",
        'name = "custom"',
        'channels = ["conda-forge"]',
        'platforms = ["linux-64"]'
      ),
      overwrite = TRUE,
      platforms = "linux-64"
    ),
    "custom pixi manifest requires a matching"
  )
})

test_that("deprecated environment-specific pixi call aliases warn and forward", {
  env <- tempfile("shennong-pixi-alias-")
  on.exit(unlink(env, recursive = TRUE, force = TRUE), add = TRUE)

  expect_warning(
    resolved <- tryCatch(
      sn_call_stlearn("definitely-not-a-real-command", args = character(), runtime_dir = env, quiet = TRUE),
      error = function(e) e
    ),
    "deprecated",
    ignore.case = TRUE
  )

  # The alias must forward to sn_call_pixi_environment(): the failure mode is
  # the generic unknown-command error from the managed environment, not a
  # missing-function error.
  expect_s3_class(resolved, "error")
  expect_false(grepl("could not find function|no applicable method", conditionMessage(resolved)))
})
