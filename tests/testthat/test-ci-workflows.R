test_that("routine package checks avoid duplicate full tests", {
  path <- test_path("..", "..", ".github", "workflows", "R-CMD-check.yaml")
  skip_if_not(file.exists(path), "Repository-only workflow is excluded from source packages.")
  workflow <- readLines(path, warn = FALSE)
  text <- paste(workflow, collapse = "\n")

  expect_match(text, "workflow_dispatch:", fixed = TRUE)
  expect_match(text, "--no-tests", fixed = TRUE)
  expect_match(text, "cancel-in-progress: true", fixed = TRUE)
})

test_that("pkgdown CI checks builds without deployment permissions", {
  workflow_path <- test_path("..", "..", ".github", "workflows", "pkgdown.yaml")
  helper_path <- test_path("..", "..", "scripts", "build-pkgdown.R")
  skip_if_not(file.exists(workflow_path) && file.exists(helper_path), "Repository-only CI files are excluded from source packages.")
  workflow <- readLines(workflow_path, warn = FALSE)
  helper <- readLines(helper_path, warn = FALSE)
  workflow_text <- paste(workflow, collapse = "\n")
  helper_text <- paste(helper, collapse = "\n")

  expect_match(workflow_text, "contents: read", fixed = TRUE)
  expect_false(grepl("contents: write", workflow_text, fixed = TRUE))
  expect_false(grepl("deploy_to_branch", workflow_text, fixed = TRUE))
  expect_false(grepl("git push", workflow_text, fixed = TRUE))
  expect_false(grepl("gh-pages", workflow_text, fixed = TRUE))
  expect_match(workflow_text, "Rscript scripts/build-pkgdown.R --full", fixed = TRUE)
  expect_match(workflow_text, "Rscript scripts/build-pkgdown.R", fixed = TRUE)
  expect_match(workflow_text, "PKGDOWN_FULL", fixed = TRUE)
  expect_match(helper_text, "lazy = !full", fixed = TRUE)
  expect_match(helper_text, "examples = full", fixed = TRUE)
  expect_match(helper_text, "--full", fixed = TRUE)
  expect_match(helper_text, "--real", fixed = TRUE)
  expect_match(helper_text, "--extended", fixed = TRUE)
  expect_match(helper_text, "SHENNONG_REAL_DATA_DIR", fixed = TRUE)
  expect_match(helper_text, "SHENNONG_RUN_VIGNETTES", fixed = TRUE)
  expect_match(helper_text, "SHENNONG_REAL_PROFILE", fixed = TRUE)
  expect_match(helper_text, '"--bundle", "all"', fixed = TRUE)
  expect_match(helper_text, "pkgdown-build-manifest.R", fixed = TRUE)
  expect_match(helper_text, ".pkgdown_remove_generated_logs", fixed = TRUE)
  expect_match(helper_text, ".pkgdown_write_build_manifest", fixed = TRUE)
})

test_that("local pkgdown real-data CLI documents modes and fails before building", {
  helper_path <- test_path("..", "..", "scripts", "build-pkgdown.R")
  skip_if_not(file.exists(helper_path), "Repository-only helper is excluded from source packages.")
  helper_path <- normalizePath(helper_path, winslash = "/", mustWork = TRUE)

  run_helper <- function(arguments) {
    output <- suppressWarnings(system2(
      file.path(R.home("bin"), "Rscript"),
      c(shQuote(helper_path), vapply(arguments, shQuote, character(1))),
      stdout = TRUE,
      stderr = TRUE
    ))
    status <- attr(output, "status")
    list(
      status = if (is.null(status)) 0L else as.integer(status),
      output = paste(output, collapse = "\n")
    )
  }

  help <- run_helper("--help")
  expect_identical(help$status, 0L)
  expect_match(help$output, "--real", fixed = TRUE)
  expect_match(help$output, "--extended", fixed = TRUE)
  expect_match(help$output, "--data-root PATH", fixed = TRUE)
  expect_match(help$output, "never download data", fixed = TRUE)

  invalid_mode <- run_helper("--extended")
  expect_gt(invalid_mode$status, 0L)
  expect_match(invalid_mode$output, "`--extended` requires `--real`", fixed = TRUE)

  missing_root <- tempfile("shennong-missing-real-data-")
  unlink(missing_root, recursive = TRUE, force = TRUE)
  missing <- run_helper(c("--real", "--data-root", missing_root))
  expect_gt(missing$status, 0L)
  expect_match(
    missing$output,
    "Validating all four local real-data bundles before the pkgdown build.",
    fixed = TRUE
  )
  expect_match(
    missing$output,
    "Real-data validation failed; pkgdown was not built and no data were downloaded.",
    fixed = TRUE
  )
  expect_false(dir.exists(missing_root))
})

run_pkgdown_publish_cli <- function(script, arguments) {
  output <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c(shQuote(script), vapply(arguments, shQuote, character(1))),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(output, "status")
  list(
    status = if (is.null(status)) 0L else as.integer(status),
    output = paste(output, collapse = "\n")
  )
}

test_that("real pkgdown publisher is audit-only unless explicitly published", {
  root <- test_path("..", "..")
  publisher <- file.path(root, "scripts", "publish-pkgdown-real.R")
  skip_if_not(file.exists(publisher), "Repository-only publisher is excluded from source packages.")
  skip_if_not_installed("jsonlite")
  publisher <- normalizePath(publisher, winslash = "/", mustWork = TRUE)

  help <- run_pkgdown_publish_cli(publisher, "--help")
  expect_identical(help$status, 0L)
  expect_match(help$output, "dry run", fixed = TRUE)
  expect_match(help$output, "--publish", fixed = TRUE)
  expect_match(help$output, "--data-root PATH", fixed = TRUE)
  expect_match(help$output, "--site-dir PATH", fixed = TRUE)
  expect_match(help$output, "replace only dev/", fixed = TRUE)

  source <- paste(readLines(publisher, warn = FALSE), collapse = "\n")
  expect_false(grepl("pkgdown::build", source, fixed = TRUE))
  expect_false(grepl("rmarkdown::render", source, fixed = TRUE))
  expect_false(grepl("download.file", source, fixed = TRUE))
  expect_match(source, "runtime-coverage-core-summary.json", fixed = TRUE)
  expect_match(source, "pkgdown-build-manifest.R", fixed = TRUE)
  expect_match(source, ".pkgdown_assert_static_allowlist", fixed = TRUE)
  expect_match(source, ".pkgdown_verify_build_manifest", fixed = TRUE)
  expect_match(source, "worktree", fixed = TRUE)
  expect_match(source, 'c("push", "origin", "HEAD:gh-pages")', fixed = TRUE)
})

test_that("real pkgdown publisher audits an existing static site before Git", {
  source_root <- normalizePath(test_path("..", ".."), winslash = "/", mustWork = TRUE)
  publisher_source <- file.path(source_root, "scripts", "publish-pkgdown-real.R")
  manifest_helper_source <- file.path(
    source_root,
    "scripts",
    "real-data",
    "pkgdown-build-manifest.R"
  )
  skip_if_not(file.exists(publisher_source), "Repository-only publisher is excluded from source packages.")
  skip_if_not(file.exists(manifest_helper_source), "Repository-only manifest helper is excluded from source packages.")
  skip_if_not_installed("jsonlite")

  sandbox <- tempfile("shennong-pkgdown-publish-audit-")
  root <- file.path(sandbox, "repo")
  dir.create(file.path(root, "scripts", "real-data"), recursive = TRUE)
  dir.create(file.path(root, "R"), recursive = TRUE)
  dir.create(file.path(root, "man"), recursive = TRUE)
  dir.create(file.path(root, "vignettes"), recursive = TRUE)
  on.exit(unlink(sandbox, recursive = TRUE, force = TRUE), add = TRUE)

  publisher <- file.path(root, "scripts", "publish-pkgdown-real.R")
  manifest_helper <- file.path(root, "scripts", "real-data", "pkgdown-build-manifest.R")
  expect_true(file.copy(publisher_source, publisher))
  expect_true(file.copy(manifest_helper_source, manifest_helper))
  version <- unname(read.dcf(
    file.path(source_root, "DESCRIPTION"),
    fields = "Version"
  )[[1]])
  writeLines(
    c("Package: Shennong", paste0("Version: ", version)),
    file.path(root, "DESCRIPTION")
  )
  writeLines("export(sn_plot_example)", file.path(root, "NAMESPACE"))
  writeLines("sn_plot_example <- function() NULL", file.path(root, "R", "example.R"))
  writeLines("# Example", file.path(root, "vignettes", "example.Rmd"))
  writeLines(
    c(
      "\\name{sn_plot_example}",
      "\\alias{sn_plot_example}",
      "\\title{Example plot}",
      "\\usage{sn_plot_example()}",
      "\\description{A test-only plot helper.}"
    ),
    file.path(root, "man", "sn_plot_example.Rd")
  )
  run_git <- function(arguments) {
    output <- suppressWarnings(system2(
      "git",
      c("-C", shQuote(root), vapply(arguments, shQuote, character(1))),
      stdout = TRUE,
      stderr = TRUE
    ))
    status <- attr(output, "status")
    if (is.null(status)) 0L else as.integer(status)
  }
  expect_identical(run_git("init"), 0L)
  expect_identical(run_git(c("config", "user.name", "Shennong Test")), 0L)
  expect_identical(run_git(c("config", "user.email", "shennong-test@example.invalid")), 0L)
  expect_identical(run_git(c("add", "-A")), 0L)
  expect_identical(
    run_git(c("commit", "--no-gpg-sign", "--no-verify", "-m", "test fixture")),
    0L
  )

  helper <- new.env(parent = globalenv())
  sys.source(manifest_helper, envir = helper)

  data_root <- file.path(sandbox, "real-data")
  results_root <- file.path(data_root, "results")
  runtime_article_root <- file.path(results_root, "runtime-articles", "core")
  site_dir <- file.path(sandbox, "site-dev")
  article_assets <- file.path(
    site_dir,
    "articles",
    "example_files",
    "figure-html"
  )
  dir.create(runtime_article_root, recursive = TRUE)
  dir.create(article_assets, recursive = TRUE)

  runtime_output <- file.path(runtime_article_root, "example.html")
  writeLines("<html><body>runtime output</body></html>", runtime_output)
  writeLines(
    '<html><body><a href="https://github.com/zerostwo/shennong">source</a></body></html>',
    file.path(site_dir, "index.html")
  )
  article_html <- file.path(site_dir, "articles", "example.html")
  writeLines(
    paste0(
      "<html><body><pre><code>#&gt; real output</code></pre>",
      '<img src="example_files/figure-html/plot-1.png"></body></html>'
    ),
    article_html
  )
  writeBin(charToRaw("static figure"), file.path(article_assets, "plot-1.png"))

  expected_reference <- helper$.pkgdown_expected_reference_pages(root)
  reference_paths <- file.path(site_dir, expected_reference)
  invisible(lapply(reference_paths, function(path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    writeLines("<html><body>current reference</body></html>", path)
  }))

  report <- list(
    schema_version = "1.0.0",
    profile = "core",
    data_root = normalizePath(data_root, winslash = "/", mustWork = TRUE),
    package = list(name = "Shennong", version = version),
    validation = list(status = "success", bundles = "all"),
    network_policy = list(downloads_allowed = FALSE, attempts = list()),
    totals = list(core_missing = 0L, articles_failed = 0L),
    core_missing = list(),
    coverage = list(list(
      "function" = "sn_plot_example",
      category = "visualization",
      level = "core",
      article = "vignettes/example.Rmd",
      runtime_status = "observed",
      article_status = "success"
    )),
    articles = list(list(
      article = "vignettes/example.Rmd",
      output = normalizePath(runtime_output, winslash = "/", mustWork = TRUE),
      status = "success",
      observed_function_count = 1L
    ))
  )
  jsonlite::write_json(
    report,
    file.path(results_root, "runtime-coverage-core-summary.json"),
    auto_unbox = TRUE,
    pretty = TRUE
  )
  runtime_report <- file.path(results_root, "runtime-coverage-core-summary.json")
  manifest_path <- file.path(results_root, helper$.pkgdown_manifest_filename)
  write_manifest <- function() {
    helper$.pkgdown_write_build_manifest(
      site_dir = site_dir,
      repo_root = root,
      runtime_report = runtime_report,
      manifest_path = manifest_path
    )
  }
  write_manifest()

  dry_run <- run_pkgdown_publish_cli(
    publisher,
    c("--data-root", data_root, "--site-dir", site_dir)
  )
  expect_identical(dry_run$status, 0L, info = dry_run$output)
  expect_match(dry_run$output, "Audit passed", fixed = TRUE)
  expect_match(dry_run$output, "Dry run only", fixed = TRUE)
  expect_match(dry_run$output, "no worktree, commit, or push", fixed = TRUE)

  writeLines("<html><body><p>fixture missing</p></body></html>", article_html)
  write_manifest()
  missing_fixture <- run_pkgdown_publish_cli(
    publisher,
    c("--data-root", data_root, "--site-dir", site_dir)
  )
  expect_gt(missing_fixture$status, 0L)
  expect_match(missing_fixture$output, "reports a missing fixture", fixed = TRUE)

  writeLines(
    paste0(
      '<html><body><div class="sourceCode"><pre><code>',
      "#&gt; text that appears only inside source code",
      "</code></pre></div></body></html>"
    ),
    article_html
  )
  unlink(file.path(article_assets, "plot-1.png"))
  write_manifest()
  code_only <- run_pkgdown_publish_cli(
    publisher,
    c("--data-root", data_root, "--site-dir", site_dir)
  )
  expect_gt(code_only$status, 0L)
  expect_match(code_only$output, "has no rendered table, console output, or figure", fixed = TRUE)

  writeLines(
    "<html><body><pre><code>#&gt; real output</code></pre></body></html>",
    article_html
  )
  writeBin(charToRaw("static figure"), file.path(article_assets, "plot-1.png"))
  writeBin(charToRaw("not allowed"), file.path(site_dir, "raw-fixture.qs2"))
  forbidden_data <- run_pkgdown_publish_cli(
    publisher,
    c("--data-root", data_root, "--site-dir", site_dir)
  )
  expect_gt(forbidden_data$status, 0L)
  expect_match(
    forbidden_data$output,
    "forbidden log/raw/data/table/database/archive files",
    fixed = TRUE
  )
  unlink(file.path(site_dir, "raw-fixture.qs2"))

  writeLines("OmniPath runtime details", file.path(site_dir, "omnipath.log"))
  forbidden_log <- run_pkgdown_publish_cli(
    publisher,
    c("--data-root", data_root, "--site-dir", site_dir)
  )
  expect_gt(forbidden_log$status, 0L)
  expect_match(forbidden_log$output, "forbidden log/raw/data", fixed = TRUE)
  unlink(file.path(site_dir, "omnipath.log"))

  writeLines("suspicious", file.path(site_dir, "fixture"))
  extensionless <- run_pkgdown_publish_cli(
    publisher,
    c("--data-root", data_root, "--site-dir", site_dir)
  )
  expect_gt(extensionless$status, 0L)
  expect_match(extensionless$output, "suspicious extensionless files", fixed = TRUE)
  unlink(file.path(site_dir, "fixture"))

  writeBin(charToRaw("binary"), file.path(site_dir, "payload.exe"))
  unknown_static <- run_pkgdown_publish_cli(
    publisher,
    c("--data-root", data_root, "--site-dir", site_dir)
  )
  expect_gt(unknown_static$status, 0L)
  expect_match(unknown_static$output, "outside the conservative static-site allowlist", fixed = TRUE)
  unlink(file.path(site_dir, "payload.exe"))

  stale_reference <- file.path(site_dir, "reference", "sn_removed_api.html")
  writeLines("<html><body>stale</body></html>", stale_reference)
  extra_reference <- run_pkgdown_publish_cli(
    publisher,
    c("--data-root", data_root, "--site-dir", site_dir)
  )
  expect_gt(extra_reference$status, 0L)
  expect_match(extra_reference$output, "stale/extra", fixed = TRUE)
  expect_match(extra_reference$output, "sn_removed_api.html", fixed = TRUE)
  unlink(stale_reference)

  removed_reference <- reference_paths[[1]]
  unlink(removed_reference)
  missing_reference <- run_pkgdown_publish_cli(
    publisher,
    c("--data-root", data_root, "--site-dir", site_dir)
  )
  expect_gt(missing_reference$status, 0L)
  expect_match(missing_reference$output, "missing:", fixed = TRUE)
  writeLines("<html><body>current reference</body></html>", removed_reference)

  write_manifest()
  writeLines(
    "<html><body><pre><code>#&gt; changed after build</code></pre></body></html>",
    article_html
  )
  changed_site <- run_pkgdown_publish_cli(
    publisher,
    c("--data-root", data_root, "--site-dir", site_dir)
  )
  expect_gt(changed_site$status, 0L)
  expect_match(changed_site$output, "changed after the manifest was created", fixed = TRUE)

  writeLines(
    "<html><body><pre><code>#&gt; real output</code></pre></body></html>",
    article_html
  )
  write_manifest()
  cat(" ", file = runtime_report, append = TRUE)
  changed_runtime <- run_pkgdown_publish_cli(
    publisher,
    c("--data-root", data_root, "--site-dir", site_dir)
  )
  expect_gt(changed_runtime$status, 0L)
  expect_match(changed_runtime$output, "runtime-report SHA-256 does not match", fixed = TRUE)
})
