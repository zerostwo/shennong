#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", all_args, value = TRUE)
script_file <- if (length(script_arg)) {
  sub("^--file=", "", script_arg[[1]])
} else {
  file.path(getwd(), "scripts", "real-data", "run-runtime-coverage.R")
}
repo_root <- normalizePath(
  file.path(dirname(script_file), "..", ".."),
  winslash = "/",
  mustWork = TRUE
)

.usage <- function() {
  cat(
    "Usage: Rscript scripts/real-data/run-runtime-coverage.R [options]\n",
    "\n",
    "Render the mapped real-data articles in this R process while tracing the\n",
    "current Shennong checkout namespace. The command validates all four local fixtures\n",
    "first and never downloads data or installs optional backends.\n",
    "\n",
    "Options:\n",
    "  --data-root PATH   Local fixture root (default: SHENNONG_REAL_DATA_DIR or\n",
    "                     data-local/pkgdown-real)\n",
    "  --extended         Render core and extended article branches (default: core)\n",
    "  --help             Show this message\n",
    sep = ""
  )
}

.parse_args <- function(values) {
  if ("--help" %in% values) {
    .usage()
    quit(status = 0L)
  }
  parsed <- list(data_root = NULL, extended = FALSE)
  index <- 1L
  while (index <= length(values)) {
    value <- values[[index]]
    if (identical(value, "--extended")) {
      parsed$extended <- TRUE
    } else if (startsWith(value, "--data-root=")) {
      parsed$data_root <- sub("^--data-root=", "", value)
      if (!nzchar(parsed$data_root)) {
        stop("`--data-root` requires a non-empty path.", call. = FALSE)
      }
    } else if (identical(value, "--data-root")) {
      if (index == length(values)) {
        stop("`--data-root` requires a path.", call. = FALSE)
      }
      index <- index + 1L
      parsed$data_root <- values[[index]]
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

.command_status <- function(output) {
  status <- attr(output, "status")
  if (is.null(status)) 0L else as.integer(status)
}

.run_validator <- function(root) {
  validator <- file.path(repo_root, "scripts", "real-data", "validate-real-data.R")
  output <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c(
      shQuote(validator),
      "--bundle", "all",
      "--root", shQuote(root)
    ),
    stdout = TRUE,
    stderr = TRUE
  ))
  if (length(output)) cat(paste(output, collapse = "\n"), "\n", sep = "")
  if (.command_status(output) != 0L) {
    stop(
      "Real-data validation failed; runtime coverage was not rendered and no data were downloaded.",
      call. = FALSE
    )
  }
  invisible(output)
}

.iso_time <- function(time = Sys.time()) {
  format(time, "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC")
}

.restore_environment <- function(previous) {
  missing <- is.na(previous)
  if (any(missing)) Sys.unsetenv(names(previous)[missing])
  present <- !missing
  if (any(present)) do.call(Sys.setenv, as.list(previous[present]))
}

.main <- function() {
  options <- .parse_args(args)
  default_root <- Sys.getenv(
    "SHENNONG_REAL_DATA_DIR",
    unset = file.path(repo_root, "data-local", "pkgdown-real")
  )
  root <- .absolute_path(options$data_root %||% default_root)

  message("Validating all four local real-data bundles before runtime tracing.")
  .run_validator(root)

  required_packages <- c("jsonlite", "pkgload", "rmarkdown")
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
  ]
  if (length(missing_packages)) {
    stop(
      "Install runtime-coverage dependencies: ",
      paste(missing_packages, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  pkgload::load_all(
    path = repo_root,
    reset = TRUE,
    recompile = FALSE,
    export_all = FALSE,
    helpers = FALSE,
    quiet = TRUE
  )
  if ("package:Shennong" %in% search()) {
    detach("package:Shennong", unload = FALSE, character.only = TRUE)
  }

  coverage_path <- file.path(repo_root, "scripts", "real-data", "coverage.csv")
  coverage <- utils::read.csv(
    coverage_path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = character()
  )
  required_columns <- c(
    "function", "category", "bundle", "level", "backend", "notes", "article"
  )
  if (!all(required_columns %in% names(coverage))) {
    stop(
      "coverage.csv must contain: ",
      paste(required_columns, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  coverage[["article"]] <- trimws(coverage[["article"]])
  profile <- if (isTRUE(options$extended)) "all" else "core"
  profile_label <- if (isTRUE(options$extended)) "extended" else "core"
  selected <- coverage$level == "core" | isTRUE(options$extended)

  environment_names <- c(
    "SHENNONG_REAL_DATA_DIR", "SHENNONG_RUN_VIGNETTES",
    "SHENNONG_REAL_PROFILE", "SHENNONG_RUN_EXTERNAL_BACKENDS"
  )
  previous_environment <- Sys.getenv(environment_names, unset = NA_character_)
  names(previous_environment) <- environment_names
  on.exit(.restore_environment(previous_environment), add = TRUE)
  Sys.setenv(
    SHENNONG_REAL_DATA_DIR = root,
    SHENNONG_RUN_VIGNETTES = "true",
    SHENNONG_REAL_PROFILE = profile,
    SHENNONG_RUN_EXTERNAL_BACKENDS = "false"
  )

  output_root <- file.path(root, "results")
  article_output_root <- file.path(
    output_root,
    "runtime-articles",
    profile_label
  )
  dir.create(article_output_root, recursive = TRUE, showWarnings = FALSE)

  namespace <- asNamespace("Shennong")
  namespace_exports <- getNamespaceExports("Shennong")
  namespace_path <- getNamespaceInfo(namespace, "path")
  trace_function <- get("trace", envir = asNamespace("methods"))
  untrace_function <- get("untrace", envir = asNamespace("methods"))
  trace_state <- new.env(parent = emptyenv())
  trace_state$current_article <- ""
  trace_state$calls <- new.env(parent = emptyenv())
  trace_state$download_attempts <- list()
  traced <- list()

  .call_key <- function(article, function_name) {
    paste(article, function_name, sep = "\034")
  }
  .call_count <- function(article, function_name) {
    get0(
      .call_key(article, function_name),
      envir = trace_state$calls,
      inherits = FALSE,
      ifnotfound = 0L
    )
  }
  .record_call <- function(function_name) {
    article <- trace_state$current_article
    if (!nzchar(article)) return(invisible(NULL))
    key <- .call_key(article, function_name)
    assign(
      key,
      get0(key, envir = trace_state$calls, inherits = FALSE, ifnotfound = 0L) + 1L,
      envir = trace_state$calls
    )
    invisible(NULL)
  }
  .record_download <- function(api) {
    attempt <- list(
      article = trace_state$current_article,
      api = api,
      observed_at = .iso_time()
    )
    trace_state$download_attempts[[length(trace_state$download_attempts) + 1L]] <- attempt
    stop(
      "Network/download API `", api,
      "` is disabled during real-data runtime coverage.",
      call. = FALSE
    )
  }

  call_helper_name <- paste0(
    ".shennong_runtime_coverage_hit_",
    Sys.getpid()
  )
  download_helper_name <- paste0(
    ".shennong_runtime_download_guard_",
    Sys.getpid()
  )
  assign(call_helper_name, .record_call, envir = globalenv())
  assign(download_helper_name, .record_download, envir = globalenv())

  .install_trace <- function(package, function_name, tracer_expression, kind) {
    package_namespace <- tryCatch(
      asNamespace(package),
      error = function(error) NULL
    )
    if (is.null(package_namespace) ||
        !exists(function_name, envir = package_namespace, inherits = FALSE)) {
      return(list(status = "unavailable", error = "function is not installed"))
    }
    target <- get(function_name, envir = package_namespace, inherits = FALSE)
    if (!is.function(target)) {
      return(list(status = "unavailable", error = "namespace binding is not a function"))
    }
    if (inherits(target, "functionWithTrace")) {
      return(list(status = "pretraced", error = "existing trace was preserved"))
    }
    error <- tryCatch(
      {
        trace_function(
          what = function_name,
          tracer = tracer_expression,
          where = package_namespace,
          print = FALSE
        )
        NULL
      },
      error = identity
    )
    if (!is.null(error)) {
      return(list(status = "failed", error = conditionMessage(error)))
    }
    traced[[length(traced) + 1L]] <<- list(
      package = package,
      function_name = function_name,
      namespace = package_namespace,
      kind = kind
    )
    list(status = "traced", error = "")
  }

  .remove_traces <- function() {
    if (length(traced)) {
      for (entry in rev(traced)) {
        try(
          untrace_function(
            what = entry$function_name,
            where = entry$namespace
          ),
          silent = TRUE
        )
      }
    }
    if ("package:Shennong" %in% search()) {
      detach("package:Shennong", unload = FALSE, character.only = TRUE)
    }
    for (helper_name in c(call_helper_name, download_helper_name)) {
      if (exists(helper_name, envir = globalenv(), inherits = FALSE)) {
        rm(list = helper_name, envir = globalenv())
      }
    }
    invisible(NULL)
  }
  on.exit(.remove_traces(), add = TRUE)

  trace_status <- stats::setNames(
    rep("not_selected", nrow(coverage)),
    coverage[["function"]]
  )
  trace_error <- stats::setNames(
    rep("", nrow(coverage)),
    coverage[["function"]]
  )
  for (function_name in coverage[["function"]][selected]) {
    if (!function_name %in% namespace_exports) {
      trace_status[[function_name]] <- "not_exported_by_installed_package"
      trace_error[[function_name]] <- "function is absent from the installed namespace exports"
      next
    }
    tracer <- substitute(
      get(HELPER, envir = globalenv(), inherits = FALSE)(FUNCTION_NAME),
      list(HELPER = call_helper_name, FUNCTION_NAME = function_name)
    )
    installed <- .install_trace("Shennong", function_name, tracer, "coverage")
    trace_status[[function_name]] <- installed$status
    trace_error[[function_name]] <- installed$error
  }

  download_guards <- list(
    c("utils", "download.file"),
    c("curl", "curl_download"),
    c("curl", "curl_fetch_memory"),
    c("curl", "curl_fetch_disk"),
    c("httr", "GET"),
    c("httr", "POST"),
    c("httr2", "req_perform")
  )
  guarded_apis <- character()
  for (specification in download_guards) {
    package <- specification[[1]]
    function_name <- specification[[2]]
    api <- paste0(package, "::", function_name)
    tracer <- substitute(
      get(HELPER, envir = globalenv(), inherits = FALSE)(API_NAME),
      list(HELPER = download_helper_name, API_NAME = api)
    )
    installed <- .install_trace(package, function_name, tracer, "download_guard")
    if (identical(installed$status, "traced")) guarded_apis <- c(guarded_apis, api)
  }

  mapped_articles <- unique(
    coverage[["article"]][selected & nzchar(coverage[["article"]])]
  )
  mapped_articles <- sort(mapped_articles)
  article_runs <- vector("list", length(mapped_articles))
  names(article_runs) <- mapped_articles

  for (article in mapped_articles) {
    article_path <- file.path(repo_root, article)
    output_path <- file.path(
      article_output_root,
      paste0(tools::file_path_sans_ext(basename(article)), ".html")
    )
    started <- Sys.time()
    warnings <- character()
    error_message <- ""
    rendered_path <- ""
    trace_state$current_article <- article
    message("Rendering with runtime traces: ", article)
    render_error <- tryCatch(
      withCallingHandlers(
        {
          rendered_path <- rmarkdown::render(
            input = article_path,
            output_file = basename(output_path),
            output_dir = dirname(output_path),
            envir = new.env(parent = globalenv()),
            clean = TRUE,
            quiet = TRUE
          )
          NULL
        },
        warning = function(warning) {
          warnings <<- c(warnings, conditionMessage(warning))
          invokeRestart("muffleWarning")
        }
      ),
      error = identity
    )
    trace_state$current_article <- ""
    finished <- Sys.time()
    if (!is.null(render_error)) error_message <- conditionMessage(render_error)
    article_runs[[article]] <- list(
      article = article,
      output = if (nzchar(rendered_path)) {
        normalizePath(rendered_path, winslash = "/", mustWork = FALSE)
      } else {
        normalizePath(output_path, winslash = "/", mustWork = FALSE)
      },
      status = if (is.null(render_error)) "success" else "failed",
      started_at = .iso_time(started),
      finished_at = .iso_time(finished),
      elapsed_seconds = unname(as.numeric(difftime(finished, started, units = "secs"))),
      warnings = unique(warnings),
      error = error_message
    )
  }

  article_table <- if (length(article_runs)) {
    data.frame(
      article = vapply(article_runs, `[[`, character(1), "article"),
      article_status = vapply(article_runs, `[[`, character(1), "status"),
      article_started_at = vapply(article_runs, `[[`, character(1), "started_at"),
      article_finished_at = vapply(article_runs, `[[`, character(1), "finished_at"),
      article_elapsed_seconds = vapply(article_runs, `[[`, numeric(1), "elapsed_seconds"),
      article_error = vapply(article_runs, `[[`, character(1), "error"),
      article_output = vapply(article_runs, `[[`, character(1), "output"),
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(
      article = character(), article_status = character(),
      article_started_at = character(), article_finished_at = character(),
      article_elapsed_seconds = numeric(), article_error = character(),
      article_output = character(), stringsAsFactors = FALSE
    )
  }

  function_calls <- integer(nrow(coverage))
  observed_articles <- character(nrow(coverage))
  for (index in seq_len(nrow(coverage))) {
    function_name <- coverage[["function"]][[index]]
    per_article <- vapply(
      mapped_articles,
      .call_count,
      function_name = function_name,
      FUN.VALUE = integer(1)
    )
    function_calls[[index]] <- sum(per_article)
    observed_articles[[index]] <- paste(mapped_articles[per_article > 0L], collapse = ";")
  }

  runtime <- coverage
  runtime$profile <- profile
  runtime$selected <- selected
  runtime$trace_status <- unname(trace_status[coverage[["function"]]])
  runtime$trace_error <- unname(trace_error[coverage[["function"]]])
  runtime$observed <- function_calls > 0L
  runtime$call_count <- function_calls
  runtime$observed_articles <- observed_articles
  runtime$article_status <- ifelse(
    !selected,
    "not_run",
    ifelse(nzchar(runtime[["article"]]), "not_rendered", "not_mapped")
  )
  runtime$article_started_at <- ""
  runtime$article_finished_at <- ""
  runtime$article_elapsed_seconds <- NA_real_
  runtime$article_error <- ""
  runtime$article_output <- ""
  if (nrow(article_table)) {
    match_index <- match(runtime[["article"]], article_table$article)
    has_article_run <- !is.na(match_index) & selected
    for (column in setdiff(names(article_table), "article")) {
      runtime[[column]][has_article_run] <- article_table[[column]][match_index[has_article_run]]
    }
  }

  runtime$runtime_status <- ifelse(
    !selected,
    "not_run",
    ifelse(
      runtime$observed,
      "observed",
      ifelse(runtime$level == "core", "missing", "not_run")
    )
  )
  runtime$status_reason <- ""
  runtime$status_reason[!selected] <- "extended row excluded by the core profile"
  runtime$status_reason[
    selected & !runtime$observed & !nzchar(runtime[["article"]])
  ] <-
    "no article is mapped to this coverage row"
  runtime$status_reason[
    selected & !runtime$observed & nzchar(runtime[["article"]]) &
      runtime$article_status == "success" & runtime$level == "extended"
  ] <- "optional backend branch was not observed"
  runtime$status_reason[
    selected & !runtime$observed & nzchar(runtime[["article"]]) &
      runtime$article_status == "success" & runtime$level == "core"
  ] <- "core function was not observed while its mapped article rendered"
  runtime$status_reason[
    selected & !runtime$observed & runtime$article_status == "failed"
  ] <- "mapped article failed before this function was observed"
  runtime$status_reason[
    selected & !runtime$observed & runtime$trace_status != "traced"
  ] <- paste0("trace unavailable: ", runtime$trace_status[
    selected & !runtime$observed & runtime$trace_status != "traced"
  ])

  for (article in names(article_runs)) {
    per_function_calls <- vapply(
      runtime[["function"]],
      function(function_name) .call_count(article, function_name),
      integer(1)
    )
    observed <- runtime[["function"]][per_function_calls > 0L]
    article_runs[[article]]$observed_functions <- observed
    article_runs[[article]]$observed_function_count <- length(observed)
    article_runs[[article]]$observed_calls <- as.list(stats::setNames(
      per_function_calls[per_function_calls > 0L],
      observed
    ))
  }

  csv_path <- file.path(output_root, "runtime-coverage.csv")
  json_path <- file.path(output_root, "runtime-coverage-summary.json")
  profile_csv_path <- file.path(
    output_root,
    paste0("runtime-coverage-", profile_label, ".csv")
  )
  profile_json_path <- file.path(
    output_root,
    paste0("runtime-coverage-", profile_label, "-summary.json")
  )
  utils::write.csv(runtime, csv_path, row.names = FALSE, na = "")
  utils::write.csv(runtime, profile_csv_path, row.names = FALSE, na = "")

  core_missing <- runtime[["function"]][
    runtime$level == "core" & runtime$runtime_status == "missing"
  ]
  failed_articles <- vapply(
    article_runs,
    function(run) identical(run$status, "failed"),
    logical(1)
  )
  summary <- list(
    schema_version = "1.0.0",
    generated_at = .iso_time(),
    profile = profile,
    repo_root = repo_root,
    data_root = root,
    package = list(
      name = "Shennong",
      version = as.character(utils::packageVersion("Shennong")),
      loader = "pkgload::load_all",
      namespace_path = normalizePath(namespace_path, winslash = "/", mustWork = TRUE)
    ),
    validation = list(status = "success", bundles = "all"),
    network_policy = list(
      downloads_allowed = FALSE,
      guarded_apis = guarded_apis,
      attempts = trace_state$download_attempts
    ),
    totals = list(
      declarations = nrow(runtime),
      selected = sum(runtime$selected),
      observed = sum(runtime$runtime_status == "observed"),
      core_missing = length(core_missing),
      extended_not_run = sum(
        runtime$level == "extended" & runtime$runtime_status == "not_run"
      ),
      articles = length(article_runs),
      articles_succeeded = sum(!failed_articles),
      articles_failed = sum(failed_articles)
    ),
    core_missing = core_missing,
    coverage = runtime,
    articles = unname(article_runs),
    outputs = list(
      csv = normalizePath(csv_path, winslash = "/", mustWork = FALSE),
      json = normalizePath(json_path, winslash = "/", mustWork = FALSE),
      profile_csv = normalizePath(
        profile_csv_path, winslash = "/", mustWork = FALSE
      ),
      profile_json = normalizePath(
        profile_json_path, winslash = "/", mustWork = FALSE
      ),
      articles = normalizePath(
        article_output_root, winslash = "/", mustWork = FALSE
      )
    )
  )
  jsonlite::write_json(
    summary,
    json_path,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null",
    na = "null"
  )
  jsonlite::write_json(
    summary,
    profile_json_path,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null",
    na = "null"
  )

  message("Runtime coverage CSV: ", csv_path)
  message("Runtime coverage JSON: ", json_path)
  if (length(core_missing)) {
    message("Unobserved core functions (", length(core_missing), "): ", paste(core_missing, collapse = ", "))
  }
  if (any(failed_articles)) {
    message(
      "Failed articles (", sum(failed_articles), "): ",
      paste(names(article_runs)[failed_articles], collapse = ", ")
    )
  }
  if (length(trace_state$download_attempts)) {
    message("Blocked download attempts: ", length(trace_state$download_attempts))
  }

  if (length(core_missing) || any(failed_articles) || length(trace_state$download_attempts)) {
    return(1L)
  }
  0L
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

status <- tryCatch(
  .main(),
  error = function(error) {
    message("Runtime coverage failed: ", conditionMessage(error))
    1L
  }
)
quit(status = status)
