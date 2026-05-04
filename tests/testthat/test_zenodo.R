test_that("sn_upload_zenodo dry-run builds a versioned manifest", {
  skip_if_not_installed("zen4R")

  path <- tempfile(fileext = ".csv")
  write.csv(data.frame(x = 1:2), path, row.names = FALSE)
  manifest_dir <- tempfile("zenodo-manifest-")

  result <- sn_upload_zenodo(
    files = path,
    title = "Reusable test dataset",
    creators = "Duan, Songqi",
    version = "0.1.0",
    sandbox = TRUE,
    dry_run = TRUE,
    manifest_dir = manifest_dir
  )

  expect_true(result$dry_run)
  expect_false(result$published)
  expect_equal(result$version, "0.1.0")
  expect_true(file.exists(result$manifest_path))
  expect_equal(result$files$file, basename(path))
  expect_true(all(c("md5", "sha256", "size") %in% colnames(result$files)))

  manifest <- jsonlite::read_json(result$manifest_path, simplifyVector = TRUE)
  expect_equal(manifest$title, "Reusable test dataset")
  expect_equal(manifest$dataset_version, "0.1.0")
  expect_equal(manifest$package, "Shennong")
  expect_equal(manifest$files$file, basename(path))
})

test_that("sn_upload_zenodo validates file and versioning inputs", {
  skip_if_not_installed("zen4R")

  expect_error(
    sn_upload_zenodo(files = tempfile(), title = "missing", dry_run = TRUE),
    "not found"
  )
  path <- tempfile(fileext = ".txt")
  writeLines("payload", path)
  expect_error(
    sn_upload_zenodo(files = path, title = "new version", new_version = TRUE, dry_run = TRUE),
    "`record_id` is required"
  )
})

test_that("sn_upload_zenodo uploads through zen4R manager without publishing by default", {
  skip_if_not_installed("zen4R")

  path <- tempfile(fileext = ".txt")
  writeLines("payload", path)
  uploaded <- character()
  deposited_record <- NULL

  fake_manager <- new.env(parent = emptyenv())
  fake_manager$depositRecord <- function(record, reserveDOI = TRUE, publish = FALSE) {
    record$id <- "12345"
    record$pids <- list(doi = list(identifier = "10.5281/zenodo.12345"))
    record$links <- list(record_html = "https://zenodo.org/records/12345")
    deposited_record <<- record
    record
  }
  fake_manager$uploadFile <- function(path, record = NULL) {
    uploaded <<- c(uploaded, basename(path))
    list(key = basename(path), status = "completed")
  }
  fake_manager$publishRecord <- function(recordId) {
    stop("publishRecord should not be called")
  }

  result <- testthat::with_mocked_bindings(
    sn_upload_zenodo(
      files = path,
      title = "Draft upload",
      creators = data.frame(name = "Duan, Songqi", orcid = "0000-0002-0822-5883"),
      version = "2026.05.04",
      token = "fake-token",
      dry_run = FALSE
    ),
    .sn_zenodo_manager = function(...) fake_manager,
    .package = "Shennong"
  )

  expect_equal(result$record_id, "12345")
  expect_equal(result$doi, "10.5281/zenodo.12345")
  expect_false(result$published)
  expect_true(basename(path) %in% uploaded)
  expect_true("shennong_zenodo_manifest.json" %in% uploaded)
  expect_equal(deposited_record$metadata$title, "Draft upload")
  expect_equal(deposited_record$metadata$version, "2026.05.04")
})

test_that("sn_download_zenodo downloads named files without requiring an environment token", {
  old_token <- Sys.getenv("ZENODO_TOKEN", unset = NA_character_)
  Sys.setenv(ZENODO_TOKEN = "env-token-should-not-be-used")
  on.exit({
    if (is.na(old_token)) {
      Sys.unsetenv("ZENODO_TOKEN")
    } else {
      Sys.setenv(ZENODO_TOKEN = old_token)
    }
  }, add = TRUE)

  save_dir <- tempfile("zenodo-cache-")
  captured <- NULL

  paths <- testthat::with_mocked_bindings(
    sn_download_zenodo(
      record_id = "14884845",
      files = c("pbmc1k_filtered_feature_bc_matrix.h5", "pbmc3k_filtered_feature_bc_matrix.h5"),
      save_dir = save_dir,
      quiet = TRUE
    ),
    .sn_download_zenodo_file = function(url, destfile, token = NULL, quiet = FALSE) {
      captured <<- rbind(
        captured,
        data.frame(
          url = url,
          destfile = destfile,
          token = if (is.null(token)) NA_character_ else token,
          quiet = quiet,
          stringsAsFactors = FALSE
        )
      )
      dir.create(dirname(destfile), recursive = TRUE, showWarnings = FALSE)
      writeLines("payload", destfile)
      normalizePath(destfile, winslash = "/", mustWork = TRUE)
    },
    .package = "Shennong"
  )

  expect_named(paths, c("pbmc1k_filtered_feature_bc_matrix.h5", "pbmc3k_filtered_feature_bc_matrix.h5"))
  expect_true(all(file.exists(paths)))
  expect_true(all(is.na(captured$token)))
  expect_true(all(captured$quiet))
  expect_match(captured$url[[1]], "zenodo\\.org/records/14884845/files/pbmc1k_filtered_feature_bc_matrix\\.h5\\?download=1")
})

test_that("sn_download_zenodo can query record metadata and pass an explicit token", {
  save_dir <- tempfile("zenodo-cache-")
  captured <- NULL

  paths <- testthat::with_mocked_bindings(
    sn_download_zenodo(
      record_id = "12345",
      files = NULL,
      save_dir = save_dir,
      token = "restricted-token",
      sandbox = TRUE
    ),
    .sn_list_zenodo_record_files = function(record_id, token = NULL, sandbox = FALSE) {
      expect_equal(record_id, "12345")
      expect_equal(token, "restricted-token")
      expect_true(sandbox)
      c("dataset-a.qs", "dataset-b.csv")
    },
    .sn_download_zenodo_file = function(url, destfile, token = NULL, quiet = FALSE) {
      captured <<- c(captured, token)
      dir.create(dirname(destfile), recursive = TRUE, showWarnings = FALSE)
      writeLines("payload", destfile)
      normalizePath(destfile, winslash = "/", mustWork = TRUE)
    },
    .package = "Shennong"
  )

  expect_named(paths, c("dataset-a.qs", "dataset-b.csv"))
  expect_equal(captured, c("restricted-token", "restricted-token"))
  expect_true(all(file.exists(paths)))
})

test_that("sn_download_zenodo reuses cached files unless overwrite is requested", {
  save_dir <- tempfile("zenodo-cache-")
  dir.create(save_dir)
  cached <- file.path(save_dir, "dataset.qs")
  writeLines("cached", cached)

  reused <- testthat::with_mocked_bindings(
    sn_download_zenodo(record_id = "12345", files = "dataset.qs", save_dir = save_dir),
    .sn_download_zenodo_file = function(...) {
      stop("cached files should not be downloaded")
    },
    .package = "Shennong"
  )

  expect_equal(unname(reused), normalizePath(cached, winslash = "/", mustWork = TRUE))

  overwritten <- FALSE
  testthat::with_mocked_bindings(
    sn_download_zenodo(record_id = "12345", files = "dataset.qs", save_dir = save_dir, overwrite = TRUE),
    .sn_download_zenodo_file = function(url, destfile, token = NULL, quiet = FALSE) {
      overwritten <<- TRUE
      writeLines("new", destfile)
      normalizePath(destfile, winslash = "/", mustWork = TRUE)
    },
    .package = "Shennong"
  )

  expect_true(overwritten)
  expect_equal(readLines(cached), "new")
})
