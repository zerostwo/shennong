#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", all_args, value = TRUE)
script_file <- if (length(script_arg)) {
  sub("^--file=", "", script_arg[[1]])
} else {
  file.path(getwd(), "scripts", "publish-pkgdown-real.R")
}
repo_root <- normalizePath(
  dirname(dirname(script_file)),
  winslash = "/",
  mustWork = TRUE
)
manifest_helper <- file.path(
  repo_root,
  "scripts",
  "real-data",
  "pkgdown-build-manifest.R"
)
if (!file.exists(manifest_helper)) {
  stop("Pkgdown build-manifest helper was not found: ", manifest_helper, call. = FALSE)
}
sys.source(manifest_helper, envir = environment())

.usage <- function() {
  cat(
    "Usage: Rscript scripts/publish-pkgdown-real.R [options]\n",
    "\n",
    "Audit an already-rendered real-data pkgdown development site. The script\n",
    "never knits articles, downloads data, or installs packages. By default it\n",
    "is a dry run and does not create a worktree, commit, or push.\n",
    "\n",
    "Options:\n",
    "  --data-root PATH  Local real-data root (default: SHENNONG_REAL_DATA_DIR\n",
    "                    or data-local/pkgdown-real). The required report is\n",
    "                    results/runtime-coverage-core-summary.json plus the\n",
    "                    same-build pkgdown manifest.\n",
    "  --site-dir PATH   Existing rendered development site (default: site/dev).\n",
    "  --publish         After every audit passes, replace only dev/ on the\n",
    "                    existing origin/gh-pages branch, commit, and push.\n",
    "  --help            Show this message.\n",
    "\n",
    "The stable gh-pages root is preserved. Only conservative static-site file\n",
    "types are accepted, and every file must match the same-build manifest.\n",
    sep = ""
  )
}

.parse_args <- function(values) {
  parsed <- list(publish = FALSE, data_root = NULL, site_dir = NULL)
  index <- 1L
  while (index <= length(values)) {
    value <- values[[index]]
    if (identical(value, "--help")) {
      .usage()
      quit(status = 0L)
    } else if (identical(value, "--publish")) {
      parsed$publish <- TRUE
    } else if (startsWith(value, "--data-root=")) {
      parsed$data_root <- sub("^--data-root=", "", value)
      if (!nzchar(parsed$data_root)) {
        stop("`--data-root` requires a non-empty path.", call. = FALSE)
      }
    } else if (identical(value, "--data-root")) {
      if (index == length(values) || startsWith(values[[index + 1L]], "--")) {
        stop("`--data-root` requires a path.", call. = FALSE)
      }
      index <- index + 1L
      parsed$data_root <- values[[index]]
    } else if (startsWith(value, "--site-dir=")) {
      parsed$site_dir <- sub("^--site-dir=", "", value)
      if (!nzchar(parsed$site_dir)) {
        stop("`--site-dir` requires a non-empty path.", call. = FALSE)
      }
    } else if (identical(value, "--site-dir")) {
      if (index == length(values) || startsWith(values[[index + 1L]], "--")) {
        stop("`--site-dir` requires a path.", call. = FALSE)
      }
      index <- index + 1L
      parsed$site_dir <- values[[index]]
    } else {
      stop("Unknown option `", value, "`. Use `--help` for usage.", call. = FALSE)
    }
    index <- index + 1L
  }
  parsed
}

.absolute_path <- function(path) {
  path <- path.expand(path)
  if (!grepl("^/", path)) path <- file.path(repo_root, path)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

.canonical_existing <- function(path, type = c("directory", "file")) {
  type <- match.arg(type)
  exists <- if (identical(type, "directory")) dir.exists(path) else file.exists(path)
  if (!exists) stop("Required ", type, " does not exist: ", path, call. = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.inside <- function(path, parent) {
  path <- sub("/+$", "", normalizePath(path, winslash = "/", mustWork = FALSE))
  parent <- sub("/+$", "", normalizePath(parent, winslash = "/", mustWork = FALSE))
  identical(path, parent) || startsWith(paste0(path, "/"), paste0(parent, "/"))
}

.scalar_character <- function(value, label) {
  if (length(value) != 1L || !is.character(value) || is.na(value) || !nzchar(value)) {
    stop("Runtime coverage JSON has invalid `", label, "`.", call. = FALSE)
  }
  value
}

.scalar_number <- function(value, label) {
  if (length(value) != 1L || !is.numeric(value) || is.na(value) || !is.finite(value)) {
    stop("Runtime coverage JSON has invalid `", label, "`.", call. = FALSE)
  }
  unname(value)
}

.latest_analysis_source_mtime <- function() {
  candidates <- c(
    list.files(file.path(repo_root, "R"), recursive = TRUE, full.names = TRUE),
    list.files(
      file.path(repo_root, "vignettes"),
      pattern = "[.]Rmd$",
      full.names = TRUE
    ),
    file.path(repo_root, "scripts", "real-data", "coverage.csv"),
    file.path(repo_root, "DESCRIPTION"),
    file.path(repo_root, "NAMESPACE")
  )
  candidates <- candidates[file.exists(candidates) & !dir.exists(candidates)]
  if (!length(candidates)) {
    stop("Unable to locate analysis sources for freshness validation.", call. = FALSE)
  }
  max(file.info(candidates)$mtime, na.rm = TRUE)
}

.read_runtime_report <- function(data_root) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Install `jsonlite` to audit the runtime coverage report.", call. = FALSE)
  }
  report_path <- file.path(
    data_root,
    "results",
    "runtime-coverage-core-summary.json"
  )
  report_path <- .canonical_existing(report_path, "file")
  if (file.info(report_path)$mtime < .latest_analysis_source_mtime()) {
    stop(
      "Runtime coverage JSON is older than the current analysis or vignette sources; ",
      "rerun scripts/real-data/run-runtime-coverage.R.",
      call. = FALSE
    )
  }
  report <- jsonlite::read_json(report_path, simplifyVector = FALSE)

  if (!identical(report$schema_version, "1.0.0")) {
    stop("Runtime coverage schema_version must be 1.0.0.", call. = FALSE)
  }
  if (!identical(report$profile, "core")) {
    stop("Runtime coverage report must use the core profile.", call. = FALSE)
  }
  if (!identical(report$validation$status, "success")) {
    stop("Runtime coverage fixture validation did not succeed.", call. = FALSE)
  }
  if (!identical(report$network_policy$downloads_allowed, FALSE)) {
    stop("Runtime coverage must explicitly disable downloads.", call. = FALSE)
  }
  attempts <- report$network_policy$attempts
  if (is.null(attempts)) {
    stop("Runtime coverage JSON is missing `network_policy.attempts`.", call. = FALSE)
  }
  if (length(attempts) != 0L) {
    stop(
      "Runtime coverage recorded ", length(attempts), " download attempt(s).",
      call. = FALSE
    )
  }

  core_missing <- .scalar_number(
    report$totals$core_missing,
    "totals.core_missing"
  )
  articles_failed <- .scalar_number(
    report$totals$articles_failed,
    "totals.articles_failed"
  )
  if (core_missing != 0) {
    stop("Runtime coverage has ", core_missing, " missing core function(s).", call. = FALSE)
  }
  if (articles_failed != 0) {
    stop("Runtime coverage has ", articles_failed, " failed article(s).", call. = FALSE)
  }
  if (!is.null(report$core_missing) && length(report$core_missing) != 0L) {
    stop("Runtime coverage `core_missing` must be empty.", call. = FALSE)
  }

  reported_root <- .canonical_existing(
    .scalar_character(report$data_root, "data_root"),
    "directory"
  )
  if (!identical(reported_root, data_root)) {
    stop(
      "Runtime coverage was produced for a different data root: ",
      reported_root,
      call. = FALSE
    )
  }
  expected_version <- unname(read.dcf(
    file.path(repo_root, "DESCRIPTION"),
    fields = "Version"
  )[[1]])
  reported_version <- .scalar_character(
    report$package$version,
    "package.version"
  )
  if (!identical(reported_version, expected_version)) {
    stop(
      "Runtime coverage package version ", reported_version,
      " does not match DESCRIPTION version ", expected_version, ".",
      call. = FALSE
    )
  }

  coverage <- report$coverage
  if (!is.list(coverage) || !length(coverage)) {
    stop("Runtime coverage report has no coverage declarations.", call. = FALSE)
  }
  core_rows <- Filter(function(row) identical(row$level, "core"), coverage)
  if (!length(core_rows)) {
    stop("Runtime coverage report has no core declarations.", call. = FALSE)
  }
  unobserved <- vapply(
    core_rows,
    function(row) !identical(row$runtime_status, "observed"),
    logical(1)
  )
  if (any(unobserved)) {
    names <- vapply(
      core_rows[unobserved],
      function(row) row[["function"]] %||% "<unknown>",
      character(1)
    )
    stop(
      "Core coverage rows are not observed: ",
      paste(names, collapse = ", "),
      call. = FALSE
    )
  }
  mapped_rows <- Filter(
    function(row) {
      is.character(row$article) && length(row$article) == 1L && nzchar(row$article)
    },
    core_rows
  )
  if (!length(mapped_rows)) {
    stop("Runtime coverage report has no mapped core articles.", call. = FALSE)
  }
  mapped_articles <- sort(unique(vapply(
    mapped_rows,
    function(row) row$article,
    character(1)
  )))

  article_runs <- report$articles
  if (!is.list(article_runs) || !length(article_runs)) {
    stop("Runtime coverage report has no article runs.", call. = FALSE)
  }
  article_names <- vapply(
    article_runs,
    function(run) run$article %||% "",
    character(1)
  )
  if (any(!nzchar(article_names)) || anyDuplicated(article_names)) {
    stop("Runtime coverage article runs must have unique article paths.", call. = FALSE)
  }
  names(article_runs) <- article_names
  missing_runs <- setdiff(mapped_articles, article_names)
  if (length(missing_runs)) {
    stop(
      "Mapped articles are absent from runtime coverage: ",
      paste(missing_runs, collapse = ", "),
      call. = FALSE
    )
  }
  results_root <- .canonical_existing(file.path(data_root, "results"), "directory")
  for (article in mapped_articles) {
    run <- article_runs[[article]]
    if (!identical(run$status, "success")) {
      stop("Mapped article did not render successfully: ", article, call. = FALSE)
    }
    observed_count <- .scalar_number(
      run$observed_function_count,
      paste0("articles[", article, "].observed_function_count")
    )
    if (observed_count <= 0) {
      stop("Mapped article observed no Shennong functions: ", article, call. = FALSE)
    }
    output <- .canonical_existing(
      .scalar_character(run$output, paste0("articles[", article, "].output")),
      "file"
    )
    if (!.inside(output, results_root)) {
      stop("Runtime article output is outside the local results root: ", output, call. = FALSE)
    }
  }

  list(
    path = report_path,
    report = report,
    core_rows = core_rows,
    mapped_rows = mapped_rows,
    mapped_articles = mapped_articles
  )
}

.site_files <- function(site_dir, include_directories = FALSE) {
  entries <- list.files(
    site_dir,
    recursive = TRUE,
    all.files = TRUE,
    full.names = TRUE,
    include.dirs = include_directories,
    no.. = TRUE
  )
  entries[nzchar(entries)]
}

.read_text <- function(path) {
  paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
}

.without_source_blocks <- function(html) {
  gsub(
    "(?s)<div class=\"sourceCode\"[^>]*>.*?</div>",
    " ",
    html,
    perl = TRUE
  )
}

.audit_site <- function(site_dir, data_root, runtime) {
  index_path <- file.path(site_dir, "index.html")
  if (!file.exists(index_path)) {
    stop("Rendered site is missing index.html: ", site_dir, call. = FALSE)
  }
  if (.inside(site_dir, data_root) || .inside(data_root, site_dir)) {
    stop("Rendered site and real-data root must not contain one another.", call. = FALSE)
  }

  all_entries <- .site_files(site_dir, include_directories = TRUE)
  symlinks <- all_entries[nzchar(Sys.readlink(all_entries))]
  if (length(symlinks)) {
    stop(
      "Rendered site contains symbolic links: ",
      paste(symlinks, collapse = ", "),
      call. = FALSE
    )
  }
  files <- all_entries[!dir.exists(all_entries)]
  if (!length(files)) stop("Rendered site contains no static files.", call. = FALSE)
  .pkgdown_assert_static_allowlist(site_dir)
  .pkgdown_assert_reference_pages(site_dir, repo_root)
  manifest_path <- file.path(
    data_root,
    "results",
    .pkgdown_manifest_filename
  )
  manifest <- .pkgdown_verify_build_manifest(
    site_dir = site_dir,
    repo_root = repo_root,
    runtime_report = runtime$path,
    manifest_path = manifest_path
  )

  text_files <- files[grepl(
    "[.](html?|md|json|ya?ml|xml|txt|css|js|map|svg|webmanifest)$",
    files,
    ignore.case = TRUE,
    perl = TRUE
  )]
  local_roots <- unique(c(
    data_root,
    repo_root,
    site_dir,
    normalizePath(path.expand("~"), winslash = "/", mustWork = TRUE)
  ))
  generic_local_path <- paste0(
    "(?i)(file:///|(?<![:A-Za-z0-9.])/(home|Users)/[^/[:space:]<]+/|",
    "(?<![A-Za-z0-9])[A-Z]:[\\\\/]Users[\\\\/][^\\\\/[:space:]<]+[\\\\/])"
  )
  for (path in text_files) {
    text <- .read_text(path)
    leaked_root <- local_roots[vapply(
      local_roots,
      function(root) grepl(root, text, fixed = TRUE),
      logical(1)
    )]
    if (length(leaked_root) || grepl(generic_local_path, text, perl = TRUE)) {
      stop("Rendered site exposes a local absolute path in: ", path, call. = FALSE)
    }
    if (grepl("[.]html?$", path, ignore.case = TRUE)) {
      visible_text <- .without_source_blocks(text)
      visible_text <- gsub("<[^>]+>", " ", visible_text, perl = TRUE)
      if (grepl(
        "fixture[[:space:]_-]+(?:or[[:space:]]+[^.]{0,60}[[:space:]]+)?missing",
        visible_text,
        ignore.case = TRUE,
        perl = TRUE
      )) {
        stop("Rendered site reports a missing fixture in: ", path, call. = FALSE)
      }
    }
  }

  article_results <- vector("list", length(runtime$mapped_articles))
  names(article_results) <- runtime$mapped_articles
  for (article in runtime$mapped_articles) {
    slug <- tools::file_path_sans_ext(basename(article))
    html_path <- file.path(site_dir, "articles", paste0(slug, ".html"))
    if (!file.exists(html_path) || file.info(html_path)$size <= 0) {
      stop("Mapped article is absent or empty in the static site: ", article, call. = FALSE)
    }
    html <- .read_text(html_path)
    asset_dir <- file.path(site_dir, "articles", paste0(slug, "_files"))
    plot_assets <- if (dir.exists(asset_dir)) {
      list.files(
        asset_dir,
        recursive = TRUE,
        full.names = TRUE,
        pattern = "[.](png|svg|pdf|jpe?g|webp)$",
        ignore.case = TRUE
      )
    } else {
      character()
    }
    rendered_html <- .without_source_blocks(html)
    has_output <- length(plot_assets) > 0L ||
      grepl("<table[ >]", rendered_html, ignore.case = TRUE, perl = TRUE) ||
      grepl("#&gt;|class=\"r-output\"", rendered_html, perl = TRUE)
    if (!has_output) {
      stop("Mapped article has no rendered table, console output, or figure: ", article, call. = FALSE)
    }
    article_rows <- Filter(
      function(row) identical(row$article, article),
      runtime$mapped_rows
    )
    requires_plot <- any(vapply(
      article_rows,
      function(row) identical(row$category, "visualization"),
      logical(1)
    ))
    if (requires_plot && !length(plot_assets)) {
      stop(
        "Article mapped to visualization functions has no figure asset: ",
        article,
        call. = FALSE
      )
    }
    article_results[[article]] <- list(
      html = html_path,
      plot_assets = length(plot_assets),
      requires_plot = requires_plot
    )
  }

  list(files = length(files), articles = article_results, manifest = manifest)
}

.command <- function(command, arguments, fail = TRUE) {
  output <- suppressWarnings(system2(
    command,
    vapply(arguments, shQuote, character(1)),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(output, "status")
  status <- if (is.null(status)) 0L else as.integer(status)
  if (fail && status != 0L) {
    stop(
      "Command failed (", status, "): ", command, " ",
      paste(arguments, collapse = " "),
      if (length(output)) paste0("\n", paste(output, collapse = "\n")) else "",
      call. = FALSE
    )
  }
  list(status = status, output = output)
}

.git <- function(path, arguments, fail = TRUE) {
  .command("git", c("-C", path, arguments), fail = fail)
}

.copy_site <- function(source, destination) {
  if (file.exists(destination) || dir.exists(destination)) {
    unlink(destination, recursive = TRUE, force = TRUE)
  }
  if (!dir.create(destination, recursive = TRUE, showWarnings = FALSE)) {
    stop("Unable to create gh-pages dev directory: ", destination, call. = FALSE)
  }
  entries <- list.files(
    source,
    all.files = TRUE,
    full.names = TRUE,
    no.. = TRUE
  )
  if (!length(entries)) stop("Rendered site directory is empty.", call. = FALSE)
  copied <- file.copy(
    entries,
    destination,
    recursive = TRUE,
    copy.mode = TRUE,
    copy.date = TRUE
  )
  if (!all(copied)) {
    stop(
      "Failed to copy rendered site entries: ",
      paste(entries[!copied], collapse = ", "),
      call. = FALSE
    )
  }
  invisible(destination)
}

.publish <- function(site_dir, data_root, runtime) {
  status <- .git(
    repo_root,
    c("status", "--porcelain", "--untracked-files=normal")
  )$output
  if (length(status)) {
    stop(
      "Refusing to publish from a dirty source worktree. Commit or stash source changes first.",
      call. = FALSE
    )
  }

  .git(
    repo_root,
    c("fetch", "origin", "gh-pages:refs/remotes/origin/gh-pages")
  )
  worktree <- tempfile("shennong-pkgdown-gh-pages-")
  worktree_added <- FALSE
  cleanup <- function() {
    if (worktree_added) {
      .git(
        repo_root,
        c("worktree", "remove", "--force", worktree),
        fail = FALSE
      )
    }
    if (dir.exists(worktree)) unlink(worktree, recursive = TRUE, force = TRUE)
    .git(repo_root, c("worktree", "prune"), fail = FALSE)
    invisible(NULL)
  }
  on.exit(cleanup(), add = TRUE)

  .git(
    repo_root,
    c("worktree", "add", "--detach", worktree, "origin/gh-pages")
  )
  worktree_added <- TRUE
  target_dev <- file.path(worktree, "dev")
  .copy_site(site_dir, target_dev)
  .audit_site(target_dev, data_root, runtime)

  .git(worktree, c("add", "-A", "--", "dev"))
  staged <- .git(
    worktree,
    c("diff", "--cached", "--name-only", "--", "dev")
  )$output
  staged_outside_dev <- staged[!startsWith(staged, "dev/")]
  if (length(staged_outside_dev)) {
    stop("Staged files escaped gh-pages dev/: ", paste(staged_outside_dev, collapse = ", "), call. = FALSE)
  }
  if (!length(staged)) {
    message("Audited site is identical to origin/gh-pages dev/; nothing to publish.")
    return(invisible(NULL))
  }

  source_sha <- .git(repo_root, c("rev-parse", "--short=12", "HEAD"))$output[[1]]
  version <- unname(read.dcf(
    file.path(repo_root, "DESCRIPTION"),
    fields = "Version"
  )[[1]])
  commit_message <- paste0(
    "docs(pkgdown): publish audited real-data dev site for Shennong@",
    version,
    " ",
    source_sha
  )
  .git(worktree, c("commit", "-m", commit_message))
  committed <- .git(
    worktree,
    c("diff-tree", "--no-commit-id", "--name-only", "-r", "HEAD")
  )$output
  committed_outside_dev <- committed[!startsWith(committed, "dev/")]
  if (length(committed_outside_dev)) {
    stop(
      "Publication commit escaped gh-pages dev/: ",
      paste(committed_outside_dev, collapse = ", "),
      call. = FALSE
    )
  }
  .git(worktree, c("push", "origin", "HEAD:gh-pages"))
  message("Published audited static site to origin/gh-pages dev/.")
  invisible(NULL)
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

.main <- function() {
  options <- .parse_args(args)
  default_data_root <- Sys.getenv(
    "SHENNONG_REAL_DATA_DIR",
    unset = file.path(repo_root, "data-local", "pkgdown-real")
  )
  data_root <- .canonical_existing(
    .absolute_path(options$data_root %||% default_data_root),
    "directory"
  )
  site_dir <- .canonical_existing(
    .absolute_path(options$site_dir %||% file.path(repo_root, "site", "dev")),
    "directory"
  )

  runtime <- .read_runtime_report(data_root)
  audit <- .audit_site(site_dir, data_root, runtime)
  figure_count <- sum(vapply(
    audit$articles,
    function(article) article$plot_assets,
    integer(1)
  ))
  message(
    "Audit passed: ", audit$files, " static files; ",
    length(audit$articles), " mapped core articles; ",
    figure_count, " article figure assets; zero missing core functions, ",
    "failed articles, or download attempts."
  )

  if (!isTRUE(options$publish)) {
    message(
      "Dry run only: no worktree, commit, or push was created. ",
      "Re-run with --publish to replace only gh-pages dev/."
    )
    return(invisible(NULL))
  }
  .publish(site_dir, data_root, runtime)
  invisible(NULL)
}

status <- tryCatch(
  {
    .main()
    0L
  },
  error = function(error) {
    message("pkgdown publication audit failed: ", conditionMessage(error))
    1L
  }
)
quit(status = status)
