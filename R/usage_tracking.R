# Local workflow usage and timing --------------------------------------------

.sn_usage_schema_version <- 2L
.sn_usage_instrumentation_version <- 2L
.sn_usage_state <- new.env(parent = emptyenv())
.sn_usage_state$enabled <- FALSE
.sn_usage_state$config <- NULL
.sn_usage_state$session_id <- NULL
.sn_usage_state$pid <- NULL
.sn_usage_state$bindings <- list()
.sn_usage_state$stack <- list()
.sn_usage_state$warned <- character()

.sn_usage_now <- function() {
  format(Sys.time(), "%Y-%m-%dT%H:%M:%OS6Z", tz = "UTC")
}

.sn_usage_validate_flag <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop("`", arg, "` must be TRUE or FALSE.", call. = FALSE)
  }
  x
}

.sn_usage_require_database <- function() {
  missing <- c("DBI", "RSQLite")[!vapply(
    c("DBI", "RSQLite"),
    requireNamespace,
    logical(1),
    quietly = TRUE
  )]
  if (length(missing) > 0L) {
    stop(
      "Usage tracking requires optional package(s): ",
      paste(missing, collapse = ", "),
      ". Install them before enabling tracking.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.sn_usage_normalize_path <- function(path, must_exist = FALSE) {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !nzchar(trimws(path))) {
    stop("`path` must be one non-empty local file path.", call. = FALSE)
  }
  path <- path.expand(path)
  if (dir.exists(path)) {
    stop("`path` must name a SQLite file, not a directory.", call. = FALSE)
  }
  if (isTRUE(must_exist) && !file.exists(path)) {
    stop("Usage database does not exist: ", path, call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

.sn_usage_validate_table_prefix <- function(x) {
  if (!is.character(x) || length(x) != 1L || is.na(x) ||
      !grepl("^[A-Za-z][A-Za-z0-9_]{0,62}$", x)) {
    stop(
      "`table_prefix` must start with a letter and contain at most 63 letters, numbers, or underscores.",
      call. = FALSE
    )
  }
  x
}

.sn_usage_validate_store <- function(store) {
  if (!inherits(store, "sn_usage_store") || !is.list(store)) {
    stop("`store` must be created by `sn_create_usage_store()`.", call. = FALSE)
  }
  store
}

.sn_usage_allowed_categories <- c(
  "api_name", "timing", "safe_parameters", "package_versions",
  "error_class", "acceleration"
)

.sn_usage_consent_id <- function(
    scope,
    policy_id,
    policy_version,
    purposes,
    data_categories,
    expires_at,
    study_id) {
  expiry <- if (is.null(expires_at)) {
    NULL
  } else {
    format(
      as.POSIXct(expires_at, tz = "UTC"),
      "%Y-%m-%dT%H:%M:%OS6Z",
      tz = "UTC"
    )
  }
  normalized_study <- if (is.null(study_id) || length(study_id) == 0L ||
      is.na(study_id[[1L]])) {
    NA_character_
  } else {
    as.character(study_id[[1L]])
  }
  substr(digest::digest(
    list(
      scope = as.character(scope),
      policy_id = as.character(policy_id),
      policy_version = as.character(policy_version),
      purposes = sort(unique(trimws(as.character(purposes))), method = "radix"),
      data_categories = sort(
        unique(trimws(as.character(data_categories))),
        method = "radix"
      ),
      expires_at = expiry,
      study_id = normalized_study
    ),
    algo = "sha256",
    serialize = TRUE
  ), 1L, 32L)
}

.sn_usage_validate_consent <- function(consent, remote = FALSE) {
  if (!inherits(consent, "sn_usage_consent") || !is.list(consent)) {
    if (remote) {
      stop(
        "Remote usage research requires consent from `sn_confirm_usage_consent()`.",
        call. = FALSE
      )
    }
    return(NULL)
  }
  required <- c(
    "consent_id", "scope", "policy_id", "policy_version", "purposes",
    "data_categories", "granted_at", "expires_at", "study_id"
  )
  if (!all(required %in% names(consent)) ||
      !is.character(consent$consent_id) || length(consent$consent_id) != 1L ||
      is.na(consent$consent_id) ||
      !is.character(consent$data_categories) ||
      any(!consent$data_categories %in% .sn_usage_allowed_categories) ||
      (!is.null(consent$expires_at) &&
       (!inherits(consent$expires_at, "POSIXt") ||
        length(consent$expires_at) != 1L || is.na(consent$expires_at)))) {
    stop("The usage-research consent receipt is invalid or has been modified.", call. = FALSE)
  }
  expected_id <- .sn_usage_consent_id(
    scope = consent$scope,
    policy_id = consent$policy_id,
    policy_version = consent$policy_version,
    purposes = consent$purposes,
    data_categories = consent$data_categories,
    expires_at = consent$expires_at,
    study_id = consent$study_id
  )
  if (!identical(consent$consent_id, expected_id)) {
    stop("The usage-research consent receipt is invalid or has been modified.", call. = FALSE)
  }
  if (!is.null(consent$expires_at) && Sys.time() >= consent$expires_at) {
    stop("Usage-research consent has expired.", call. = FALSE)
  }
  if (remote && !identical(consent$scope, "remote_research")) {
    stop("Remote stores require `scope = \"remote_research\"` consent.", call. = FALSE)
  }
  consent
}

.sn_usage_resolve_store <- function(path = NULL, store = NULL) {
  if (is.null(store)) {
    path <- .sn_usage_normalize_path(path)
    return(structure(
      list(
        backend = "sqlite",
        path = path,
        connect = NULL,
        table_prefix = "shennong_usage",
        allow_remote_ddl = FALSE
      ),
      class = "sn_usage_store"
    ))
  }
  store <- .sn_usage_validate_store(store)
  if (!is.null(path)) {
    stop("Supply either `path` or `store`, not both.", call. = FALSE)
  }
  store
}

.sn_usage_consent_json <- function(consent) {
  if (is.null(consent)) return("{}")
  as.character(.sn_usage_json(list(
    consent_id = consent$consent_id,
    scope = consent$scope,
    policy_id = consent$policy_id,
    policy_version = consent$policy_version,
    purposes = consent$purposes,
    data_categories = consent$data_categories,
    granted_at = consent$granted_at,
    expires_at = if (is.null(consent$expires_at)) NULL else format(
      consent$expires_at,
      "%Y-%m-%dT%H:%M:%OS6Z",
      tz = "UTC"
    ),
    study_id = consent$study_id
  )))
}

.sn_usage_connect <- function(path) {
  connection <- DBI::dbConnect(RSQLite::SQLite(), dbname = path)
  DBI::dbExecute(connection, "PRAGMA foreign_keys = ON")
  DBI::dbExecute(connection, "PRAGMA busy_timeout = 5000")
  connection
}

.sn_usage_db_retry <- function(code, strict, operation) {
  if (!is.function(code)) {
    stop("Internal usage database retry code must be a function.", call. = FALSE)
  }
  last_error <- NULL
  for (attempt in seq_len(5L)) {
    value <- tryCatch(
      list(ok = TRUE, value = code()),
      error = function(error) {
        last_error <<- error
        list(ok = FALSE, value = NULL)
      }
    )
    if (isTRUE(value$ok)) {
      return(value$value)
    }
    if (attempt < 5L) {
      Sys.sleep(0.01 * (2 ^ (attempt - 1L)))
    }
  }

  message <- paste0(
    "Could not ", operation, " in the Shennong usage database: ",
    conditionMessage(last_error)
  )
  if (isTRUE(strict)) {
    stop(message, call. = FALSE)
  }
  warning_key <- paste(operation, conditionMessage(last_error), sep = "\r")
  if (!warning_key %in% .sn_usage_state$warned) {
    .sn_usage_state$warned <- c(.sn_usage_state$warned, warning_key)
    warning(message, call. = FALSE)
  }
  NULL
}

.sn_usage_initialize_database <- function(path) {
  parent <- dirname(path)
  if (!dir.exists(parent)) {
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(parent)) {
    stop("Could not create usage database directory: ", parent, call. = FALSE)
  }

  connection <- .sn_usage_connect(path)
  on.exit(DBI::dbDisconnect(connection), add = TRUE)
  DBI::dbExecute(connection, "PRAGMA journal_mode = WAL")
  DBI::dbExecute(connection, paste(
    "CREATE TABLE IF NOT EXISTS sessions (",
    "session_id TEXT PRIMARY KEY,",
    "started_at TEXT NOT NULL,",
    "finished_at TEXT,",
    "mode TEXT NOT NULL CHECK (mode IN",
    "('development','production','test','benchmark')),",
    "shennong_version TEXT NOT NULL,",
    "r_version TEXT NOT NULL,",
    "platform TEXT,",
    "packages_json TEXT NOT NULL DEFAULT '{}',",
    "consent_json TEXT NOT NULL DEFAULT '{}',",
    "destination_kind TEXT NOT NULL DEFAULT 'sqlite',",
    "study_id TEXT,",
    "config_json TEXT NOT NULL DEFAULT '{}')"
  ))
  if (!"packages_json" %in% DBI::dbListFields(connection, "sessions")) {
    DBI::dbExecute(
      connection,
      "ALTER TABLE sessions ADD COLUMN packages_json TEXT NOT NULL DEFAULT '{}'"
    )
  }
  session_fields <- DBI::dbListFields(connection, "sessions")
  if (!"consent_json" %in% session_fields) {
    DBI::dbExecute(
      connection,
      "ALTER TABLE sessions ADD COLUMN consent_json TEXT NOT NULL DEFAULT '{}'"
    )
  }
  if (!"destination_kind" %in% session_fields) {
    DBI::dbExecute(
      connection,
      "ALTER TABLE sessions ADD COLUMN destination_kind TEXT NOT NULL DEFAULT 'sqlite'"
    )
  }
  if (!"study_id" %in% session_fields) {
    DBI::dbExecute(connection, "ALTER TABLE sessions ADD COLUMN study_id TEXT")
  }
  DBI::dbExecute(connection, paste(
    "CREATE TABLE IF NOT EXISTS workflow_runs (",
    "run_id TEXT PRIMARY KEY,",
    "session_id TEXT NOT NULL REFERENCES sessions(session_id),",
    "parent_run_id TEXT REFERENCES workflow_runs(run_id),",
    "root_run_id TEXT NOT NULL,",
    "depth INTEGER NOT NULL DEFAULT 0 CHECK (depth >= 0),",
    "pid INTEGER NOT NULL,",
    "mode TEXT NOT NULL CHECK (mode IN",
    "('development','production','test','benchmark')),",
    "workflow TEXT NOT NULL,",
    "category TEXT NOT NULL,",
    "call_scope TEXT NOT NULL DEFAULT 'entrypoint',",
    "method TEXT,",
    "backend TEXT,",
    "status TEXT NOT NULL CHECK (status IN",
    "('running','ok','error','interrupted','abandoned')),",
    "started_at TEXT NOT NULL,",
    "finished_at TEXT,",
    "elapsed_ms REAL,",
    "cpu_user_ms REAL,",
    "cpu_system_ms REAL,",
    "params_json TEXT NOT NULL DEFAULT '{}',",
    "params_sha256 TEXT NOT NULL,",
    "invocation_number INTEGER NOT NULL,",
    "parameter_set_invocation INTEGER NOT NULL,",
    "warning_count INTEGER NOT NULL DEFAULT 0,",
    "acceleration_json TEXT NOT NULL DEFAULT '{}',",
    "error_class TEXT,",
    "error_message_redacted TEXT,",
    "remote_synced_at TEXT,",
    "remote_sync_error TEXT,",
    "instrumentation_version INTEGER NOT NULL DEFAULT 1)"
  ))
  run_fields <- DBI::dbListFields(connection, "workflow_runs")
  if (!"call_scope" %in% run_fields) {
    DBI::dbExecute(
      connection,
      "ALTER TABLE workflow_runs ADD COLUMN call_scope TEXT NOT NULL DEFAULT 'entrypoint'"
    )
  }
  if (!"remote_synced_at" %in% run_fields) {
    DBI::dbExecute(connection, "ALTER TABLE workflow_runs ADD COLUMN remote_synced_at TEXT")
  }
  if (!"remote_sync_error" %in% run_fields) {
    DBI::dbExecute(connection, "ALTER TABLE workflow_runs ADD COLUMN remote_sync_error TEXT")
  }
  DBI::dbExecute(
    connection,
    paste(
      "CREATE INDEX IF NOT EXISTS workflow_runs_workflow_started",
      "ON workflow_runs(workflow, started_at)"
    )
  )
  DBI::dbExecute(
    connection,
    paste(
      "CREATE INDEX IF NOT EXISTS workflow_runs_session",
      "ON workflow_runs(session_id)"
    )
  )
  DBI::dbExecute(
    connection,
    paste(
      "CREATE INDEX IF NOT EXISTS workflow_runs_parent",
      "ON workflow_runs(parent_run_id)"
    )
  )
  DBI::dbExecute(
    connection,
    paste(
      "CREATE INDEX IF NOT EXISTS workflow_runs_params",
      "ON workflow_runs(workflow, params_sha256)"
    )
  )
  DBI::dbExecute(
    connection,
    paste0("PRAGMA user_version = ", .sn_usage_schema_version)
  )
  invisible(path)
}

.sn_usage_sqlite_id <- function(connection) {
  as.character(DBI::dbGetQuery(
    connection,
    "SELECT lower(hex(randomblob(16))) AS id"
  )$id[[1L]])
}

.sn_usage_package_version <- function() {
  tryCatch(
    as.character(utils::packageVersion("Shennong")),
    error = function(error) "development"
  )
}

.sn_usage_package_versions <- function() {
  packages <- c(
    "Shennong", "ShennongOpt", "Matrix", "Seurat", "SeuratObject",
    "clusterProfiler", "enrichit", "scDblFinder", "lisi", "UCell",
    "CellChat", "nichenetr", "SoupX", "decontX", "tradeSeq", "WGCNA",
    "harmony", "edgeR", "DESeq2", "limma", "miloR", "SingleR",
    "slingshot", "BiocParallel"
  )
  installed <- packages[vapply(packages, function(package) {
    length(tryCatch(find.package(package, quiet = TRUE), error = function(error) character())) == 1L
  }, logical(1))]
  versions <- stats::setNames(
    lapply(installed, function(package) {
      list(version = as.character(utils::packageVersion(package)))
    }),
    installed
  )
  if ("ShennongOpt" %in% names(versions)) {
    description <- tryCatch(utils::packageDescription("ShennongOpt"), error = function(error) NULL)
    versions$ShennongOpt$remote_sha <- description$RemoteSha %||% NULL
  }
  versions
}

.sn_usage_json <- function(x) {
  jsonlite::toJSON(
    x,
    auto_unbox = TRUE,
    null = "null",
    na = "null",
    digits = NA,
    POSIXt = "ISO8601"
  )
}

.sn_usage_create_session <- function(path, mode, config, consent = NULL,
                                     destination_kind = "sqlite") {
  connection <- .sn_usage_connect(path)
  on.exit(DBI::dbDisconnect(connection), add = TRUE)
  session_id <- .sn_usage_sqlite_id(connection)
  DBI::dbExecute(
    connection,
    paste(
      "INSERT INTO sessions",
      "(session_id, started_at, mode, shennong_version, r_version,",
      "platform, packages_json, consent_json, destination_kind, study_id,",
      "config_json) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
    ),
    params = list(
      session_id,
      .sn_usage_now(),
      mode,
      .sn_usage_package_version(),
      R.version.string,
      R.version$platform,
      .sn_usage_json(.sn_usage_package_versions()),
      .sn_usage_consent_json(consent),
      destination_kind,
      consent$study_id %||% NA_character_,
      .sn_usage_json(config)
    )
  )
  session_id
}

.sn_usage_finish_session <- function(path, session_id, strict) {
  if (is.null(path) || is.null(session_id) || !file.exists(path)) {
    return(invisible(FALSE))
  }
  .sn_usage_db_retry(function() {
    connection <- .sn_usage_connect(path)
    on.exit(DBI::dbDisconnect(connection), add = TRUE)
    DBI::dbExecute(
      connection,
      "UPDATE sessions SET finished_at = ? WHERE session_id = ?",
      params = list(.sn_usage_now(), session_id)
    )
  }, strict = strict, operation = "finish the usage session")
  invisible(TRUE)
}

.sn_usage_control_functions <- c(
  "sn_check_usage_tracking",
  "sn_confirm_usage_consent",
  "sn_create_usage_store",
  "sn_disable_usage_tracking",
  "sn_enable_usage_tracking",
  "sn_flush_usage_tracking",
  "sn_list_usage_runs",
  "sn_summarize_usage",
  "sn_time_call",
  "sn_with_usage_tracking"
)

.sn_usage_registry <- function(functions = NULL) {
  exports <- tryCatch(
    getNamespaceExports("Shennong"),
    error = function(error) character()
  )
  included <- sort(setdiff(exports, .sn_usage_control_functions))
  if (!is.null(functions)) {
    if (!is.character(functions) || anyNA(functions) ||
        any(!nzchar(trimws(functions)))) {
      stop("`functions` must be NULL or a non-empty character vector.", call. = FALSE)
    }
    functions <- unique(trimws(functions))
    unknown <- setdiff(functions, included)
    if (length(unknown) > 0L) {
      stop(
        "Usage tracking cannot instrument function(s): ",
        paste(unknown, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    included <- sort(functions)
  }
  category <- vapply(included, function(workflow) {
    if (grepl("^[.]?(import|export)[.]rio_", workflow)) return("io")
    if (identical(workflow, "sn_enrich")) return("enrichment")
    if (grepl("^sn_plot_", workflow)) return("visualization")
    if (grepl("^sn_(get|list|check|method)_", workflow)) return("query")
    if (grepl("^sn_(read|write|export|convert)_", workflow)) return("io")
    if (grepl("^sn_(store|add|update|delete|set|upgrade)_", workflow)) {
      return("state_management")
    }
    if (grepl("^sn_(install|ensure|prepare_pixi|call|pixi|configure)_", workflow)) {
      return("runtime")
    }
    if (grepl("^sn_(normalize|filter|initialize_seurat|standardize|remove_)", workflow)) {
      return("preprocessing")
    }
    if (grepl("^sn_(calculate|assess|compare|sweep)_", workflow)) return("metrics")
    if (grepl("^sn_(interpret|annotate|annotation|review|transfer|map)_", workflow)) {
      return("interpretation")
    }
    if (grepl("^sn_(score|discover|test)_", workflow)) return("programs")
    if (grepl("^sn_(run|find|prioritize|deconvolve|integrate|simulate)_", workflow)) {
      return("analysis")
    }
    "utility"
  }, character(1))
  data.frame(
    workflow = included,
    category = unname(category),
    parameter_policy = ifelse(included == "sn_run_llm", "llm", "minimal"),
    stringsAsFactors = FALSE
  )
}

.sn_usage_sensitive_parameter <- function(name) {
  grepl(
    paste0(
      "(^|_)(token|password|secret|credential|cookie|authorization|api_?key|",
      "prompts?|messages?|responses?|queries|texts?|descriptions?|notes?|contexts?|",
      "paths?|files?|dirs?|urls?|uris?|endpoint|host|user(name)?|email|",
      "x|objects?|data|metadata|matri(x|ces)|counts?|patients?|subjects?|",
      "cells?|samples?|genes?|features?|signatures?|references?|names?|",
      "contrasts?|ident|idents|subsets?)(_|$)"
    ),
    tolower(name),
    perl = TRUE
  )
}

.sn_usage_selector_parameters <- c(
  "method", "workflow", "integration_method", "normalization_method",
  "cluster_algorithm", "database", "modality", "backend"
)

.sn_usage_safe_selector_value <- function(value, use_first = FALSE) {
  if (is.factor(value)) value <- as.character(value)
  if (!is.atomic(value) || length(value) == 0L || length(value) > 5L ||
      anyNA(value)) {
    return(list(ok = FALSE, value = NULL))
  }
  if (is.character(value)) {
    safe <- nchar(value, type = "bytes") <= 64L &
      grepl("^[A-Za-z0-9_.:+-]*$", value)
    if (!all(safe)) return(list(ok = FALSE, value = NULL))
  }
  if (isTRUE(use_first) && length(value) > 1L) value <- value[[1L]]
  list(ok = TRUE, value = unname(value))
}

.sn_usage_resolved_selectors <- function(call, evaluation_frame, selector_names) {
  selector_names <- intersect(
    .sn_usage_selector_parameters,
    as.character(selector_names %||% character())
  )
  if (!is.environment(evaluation_frame) || length(selector_names) == 0L) {
    return(list())
  }
  arguments <- .sn_usage_flatten_call(call)
  dots <- as.list(call)[["..."]]
  dots <- if (is.pairlist(dots)) as.list(dots) else list()
  dot_names <- names(dots) %||% rep("", length(dots))
  resolved <- list()
  for (name in selector_names) {
    explicit <- name %in% names(arguments)
    if (explicit && !is.symbol(arguments[[name]])) next
    dot_position <- match(name, dot_names, nomatch = 0L)
    value <- if (exists(name, envir = evaluation_frame, inherits = FALSE)) {
      tryCatch(
        get(name, envir = evaluation_frame, inherits = FALSE),
        error = function(error) structure(list(), class = "sn_usage_unresolved")
      )
    } else if (dot_position > 0L) {
      tryCatch(
        eval(
          as.call(list(as.name("...elt"), as.integer(dot_position))),
          envir = evaluation_frame
        ),
        error = function(error) structure(list(), class = "sn_usage_unresolved")
      )
    } else {
      structure(list(), class = "sn_usage_unresolved")
    }
    if (inherits(value, "sn_usage_unresolved")) next
    safe <- .sn_usage_safe_selector_value(value, use_first = !explicit)
    if (isTRUE(safe$ok)) resolved[[name]] <- safe$value
  }
  resolved
}

.sn_usage_expression_summary <- function(expression, name = "") {
  if (.sn_usage_sensitive_parameter(name)) {
    return(list(type = "redacted"))
  }
  if (is.null(expression)) return(NULL)
  if (is.atomic(expression)) {
    if (length(expression) > 100L) {
      return(list(type = typeof(expression), length = length(expression)))
    }
    if (is.character(expression)) {
      safe <- nchar(expression, type = "bytes") <= 64L &
        grepl("^[A-Za-z0-9_.:+-]*$", expression)
      if (!all(safe)) {
        return(list(type = "character", length = length(expression)))
      }
    }
    return(unname(expression))
  }
  if (is.symbol(expression)) {
    value <- as.character(expression)
    if (value %in% c("TRUE", "FALSE", "NA", "NULL", "Inf", "NaN")) {
      if (identical(value, "TRUE")) return(TRUE)
      if (identical(value, "FALSE")) return(FALSE)
      if (identical(value, "NA")) return(NA)
      if (identical(value, "NULL")) return(NULL)
      if (identical(value, "Inf")) return(Inf)
      return(NaN)
    }
    return(list(type = "expression"))
  }
  if (is.call(expression)) {
    head <- as.character(expression[[1L]])
    if (identical(head, "c")) {
      values <- lapply(as.list(expression)[-1L], .sn_usage_expression_summary)
      if (all(vapply(values, function(x) is.atomic(x) && length(x) <= 1L, logical(1)))) {
        value <- unlist(values, recursive = FALSE, use.names = FALSE)
        if (length(value) <= 100L) return(value)
      }
      return(list(type = "vector_expression", length = length(values)))
    }
    if (head %in% c("-", "+") && length(expression) == 2L &&
        is.numeric(expression[[2L]])) {
      return(if (head == "-") -expression[[2L]] else expression[[2L]])
    }
    if (identical(head, "~")) return(list(type = "formula"))
    if (identical(head, "list")) {
      return(list(type = "list", names = names(as.list(expression)[-1L])))
    }
    return(list(type = "expression", call = head[[1L]]))
  }
  list(type = class(expression)[[1L]], length = length(expression))
}

.sn_usage_flatten_call <- function(call) {
  arguments <- as.list(call)[-1L]
  dots <- arguments[["..."]]
  arguments[["..."]] <- NULL
  if (is.pairlist(dots)) {
    dots <- as.list(dots)
    named <- nzchar(names(dots) %||% rep("", length(dots)))
    arguments <- c(arguments, dots[named])
    if (any(!named)) {
      arguments[["unnamed_dots"]] <- list(
        type = "arguments",
        length = sum(!named)
      )
    }
  }
  arguments
}

.sn_usage_call_parameters <- function(
    workflow,
    call,
    evaluation_frame = NULL,
    selector_names = NULL) {
  arguments <- .sn_usage_flatten_call(call)
  result <- list()
  if (length(arguments) > 0L) {
    names(arguments) <- make.unique(
      names(arguments) %||% rep("argument", length(arguments))
    )
    result <- Map(
      function(expression, name) .sn_usage_expression_summary(expression, name),
      arguments,
      names(arguments)
    )
    names(result) <- names(arguments)
  }
  resolved <- .sn_usage_resolved_selectors(
    call = call,
    evaluation_frame = evaluation_frame,
    selector_names = selector_names
  )
  if (length(resolved) > 0L) result[names(resolved)] <- resolved
  if (identical(workflow, "sn_run_llm")) {
    keep <- intersect(names(result), c("provider", "model", "temperature"))
    result <- result[keep]
  }
  if (length(result) == 0L) return(list())
  result[order(names(result))]
}

.sn_usage_scalar_parameter <- function(parameters, names) {
  for (name in names) {
    value <- parameters[[name]]
    if (is.atomic(value) && length(value) > 0L && length(value) <= 5L) {
      return(paste(value, collapse = ","))
    }
  }
  NA_character_
}

.sn_usage_method <- function(workflow, parameters) {
  selectors <- intersect(
    c(
      "method", "workflow", "integration_method", "normalization_method",
      "cluster_algorithm", "database", "modality"
    ),
    names(parameters)
  )
  values <- vapply(selectors, function(name) {
    value <- parameters[[name]]
    if (is.atomic(value) && length(value) > 0L && length(value) <= 5L) {
      paste0(name, "=", paste(value, collapse = ","))
    } else {
      NA_character_
    }
  }, character(1))
  values <- values[!is.na(values)]
  if (length(values) == 0L) NA_character_ else paste(values, collapse = ";")
}

.sn_usage_current_token <- function() {
  stack <- .sn_usage_state$stack
  if (length(stack) == 0L) return(NULL)
  stack[[length(stack)]]
}

.sn_usage_push <- function(token) {
  .sn_usage_state$stack[[length(.sn_usage_state$stack) + 1L]] <- token
  invisible(token)
}

.sn_usage_pop <- function(token) {
  stack <- .sn_usage_state$stack
  if (length(stack) == 0L) return(invisible(FALSE))
  positions <- which(vapply(stack, identical, logical(1), y = token))
  if (length(positions) == 0L) return(invisible(FALSE))
  .sn_usage_state$stack <- stack[-positions[[length(positions)]]]
  invisible(TRUE)
}

.sn_usage_insert_run <- function(token) {
  config <- .sn_usage_state$config
  connection <- .sn_usage_connect(config$path)
  on.exit(DBI::dbDisconnect(connection), add = TRUE)
  committed <- FALSE
  DBI::dbExecute(connection, "BEGIN IMMEDIATE")
  on.exit({
    if (!committed) try(DBI::dbExecute(connection, "ROLLBACK"), silent = TRUE)
  }, add = TRUE)

  run_id <- .sn_usage_sqlite_id(connection)
  ordinal <- DBI::dbGetQuery(
    connection,
    paste(
      "SELECT COALESCE(MAX(invocation_number), 0) + 1 AS value",
      "FROM workflow_runs WHERE workflow = ?"
    ),
    params = list(token$workflow)
  )$value[[1L]]
  parameter_ordinal <- DBI::dbGetQuery(
    connection,
    paste(
      "SELECT COALESCE(MAX(parameter_set_invocation), 0) + 1 AS value",
      "FROM workflow_runs WHERE workflow = ? AND params_sha256 = ?"
    ),
    params = list(token$workflow, token$params_sha256)
  )$value[[1L]]
  parent <- token$parent
  parent_id <- if (is.environment(parent) && isTRUE(parent$recorded)) {
    parent$run_id
  } else {
    NA_character_
  }
  root_id <- if (is.environment(parent) && isTRUE(parent$recorded)) {
    parent$root_run_id
  } else {
    run_id
  }
  DBI::dbExecute(
    connection,
    paste(
      "INSERT INTO workflow_runs",
      "(run_id, session_id, parent_run_id, root_run_id, depth, pid, mode,",
      "workflow, category, call_scope, method, backend, status, started_at, params_json,",
      "params_sha256, invocation_number, parameter_set_invocation,",
      "instrumentation_version)",
      "VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 'running', ?, ?, ?, ?, ?, ?)"
    ),
    params = list(
      run_id,
      .sn_usage_state$session_id,
      parent_id,
      root_id,
      token$depth,
      as.integer(Sys.getpid()),
      config$mode,
      token$workflow,
      token$category,
      token$call_scope,
      token$method,
      token$backend,
      token$started_at,
      token$params_json,
      token$params_sha256,
      as.integer(ordinal),
      as.integer(parameter_ordinal),
      .sn_usage_instrumentation_version
    )
  )
  DBI::dbExecute(connection, "COMMIT")
  committed <- TRUE
  list(
    run_id = run_id,
    root_run_id = root_id,
    invocation_number = as.integer(ordinal),
    parameter_set_invocation = as.integer(parameter_ordinal)
  )
}

.sn_usage_begin <- function(
    workflow,
    call,
    category = NULL,
    evaluation_frame = NULL,
    selector_names = NULL) {
  token <- new.env(parent = emptyenv())
  token$workflow <- workflow
  token$started_at <- .sn_usage_now()
  token$started_proc <- proc.time()
  token$warning_count <- 0L
  token$acceleration_patches <- character()
  token$recorded <- FALSE
  token$tracking_active <- FALSE
  token$finished <- FALSE
  token$run_id <- NULL
  token$root_run_id <- NULL
  token$invocation_number <- NA_integer_
  token$parameter_set_invocation <- NA_integer_
  token$parent <- .sn_usage_current_token()
  token$depth <- if (is.environment(token$parent)) token$parent$depth + 1L else 0L
  token$call_scope <- if (token$depth == 0L) "entrypoint" else "component"

  if (!isTRUE(.sn_usage_state$enabled) ||
      !identical(as.integer(Sys.getpid()), as.integer(.sn_usage_state$pid))) {
    return(token)
  }
  config <- .sn_usage_state$config
  if (token$depth > 0L && !isTRUE(config$nested)) {
    return(token)
  }
  if (is.environment(token$parent) &&
      isTRUE(token$parent$tracking_active) &&
      !isTRUE(token$parent$recorded)) {
    return(token)
  }
  token$tracking_active <- TRUE

  registry <- .sn_usage_registry()
  row <- registry[registry$workflow == workflow, , drop = FALSE]
  token$category <- category %||%
    if (nrow(row) == 1L) row$category[[1L]] else "ad_hoc"
  parameters <- .sn_usage_call_parameters(
    workflow,
    call,
    evaluation_frame = evaluation_frame,
    selector_names = selector_names
  )
  token$params_json <- as.character(.sn_usage_json(parameters))
  token$params_sha256 <- digest::digest(token$params_json, algo = "sha256", serialize = FALSE)
  token$method <- .sn_usage_method(workflow, parameters)
  token$backend <- .sn_usage_scalar_parameter(parameters, c("backend"))

  inserted <- .sn_usage_db_retry(
    function() .sn_usage_insert_run(token),
    strict = config$strict,
    operation = paste0("start usage record for `", workflow, "`")
  )
  if (!is.null(inserted)) {
    token$recorded <- TRUE
    token$run_id <- inserted$run_id
    token$root_run_id <- inserted$root_run_id
    token$invocation_number <- inserted$invocation_number
    token$parameter_set_invocation <- inserted$parameter_set_invocation
  }
  # Exclude the initial SQLite insert from the function runtime. The matching
  # completion update is also performed only after the elapsed snapshot.
  token$started_proc <- proc.time()
  .sn_usage_push(token)
  token
}

.sn_usage_count_warning <- function(token) {
  current <- .sn_usage_current_token()
  if (is.environment(token) && is.environment(current) && identical(token, current)) {
    token$warning_count <- token$warning_count + 1L
  }
  invisible(NULL)
}

.sn_usage_record_acceleration <- function(patches) {
  token <- .sn_usage_current_token()
  if (is.environment(token)) {
    token$acceleration_patches <- union(
      token$acceleration_patches,
      unique(as.character(patches))
    )
  }
  invisible(NULL)
}

.sn_usage_redact_error <- function(condition) {
  if (is.null(condition)) return(NA_character_)
  message <- conditionMessage(condition)
  message <- gsub("https?://[^[:space:]]+", "<url>", message, perl = TRUE)
  message <- gsub(
    "[[:alnum:]._%+-]+@[[:alnum:].-]+[.][[:alpha:]]{2,}",
    "<email>",
    message,
    perl = TRUE
  )
  message <- gsub("(^|[[:space:]'\"])/[^[:space:]'\"]+", "\\1<path>", message, perl = TRUE)
  message <- gsub(
    "(?i)(token|password|secret|api[_-]?key)[[:space:]]*[:=][[:space:]]*[^,;[:space:]]+",
    "\\1=<redacted>",
    message,
    perl = TRUE
  )
  substr(message, 1L, 512L)
}

.sn_usage_finish_run <- function(token, status, condition, timing) {
  config <- .sn_usage_state$config
  acceleration <- list(
    patches_activated = sort(unique(token$acceleration_patches)),
    evidence = if (length(token$acceleration_patches) > 0L) {
      "scope_activation_only"
    } else {
      "none"
    },
    fast_path_hit = NULL
  )
  error_class <- if (is.null(condition)) {
    NA_character_
  } else {
    paste(class(condition), collapse = ",")
  }
  connection <- .sn_usage_connect(config$path)
  on.exit(DBI::dbDisconnect(connection), add = TRUE)
  DBI::dbExecute(
    connection,
    paste(
      "UPDATE workflow_runs SET status = ?, finished_at = ?, elapsed_ms = ?,",
      "cpu_user_ms = ?, cpu_system_ms = ?, warning_count = ?,",
      "acceleration_json = ?, error_class = ?, error_message_redacted = ?",
      "WHERE run_id = ?"
    ),
    params = list(
      status,
      .sn_usage_now(),
      unname(1000 * timing[["elapsed"]]),
      unname(1000 * timing[["user.self"]]),
      unname(1000 * timing[["sys.self"]]),
      as.integer(token$warning_count),
      as.character(.sn_usage_json(acceleration)),
      error_class,
      .sn_usage_redact_error(condition),
      token$run_id
    )
  )
}

.sn_usage_finish <- function(token, status = "ok", condition = NULL) {
  if (!is.environment(token) || isTRUE(token$finished)) return(invisible(NULL))
  token$finished <- TRUE
  timing <- proc.time() - token$started_proc
  .sn_usage_pop(token)
  config <- .sn_usage_state$config
  if (isTRUE(token$recorded) && !is.null(config)) {
    finish_strict <- isTRUE(config$strict) && is.null(condition)
    .sn_usage_db_retry(
      function() .sn_usage_finish_run(token, status, condition, timing),
      strict = finish_strict,
      operation = paste0("finish usage record for `", token$workflow, "`")
    )
  }
  if (!is.null(config) && isTRUE(token$tracking_active) && isTRUE(config$display)) {
    number <- if (isTRUE(token$recorded)) {
      paste0(" #", token$invocation_number)
    } else {
      ""
    }
    message(sprintf(
      "[Shennong timing] %s%s (%s): %.3f s [%s]",
      token$workflow,
      number,
      config$mode,
      unname(timing[["elapsed"]]),
      status
    ))
  }
  invisible(NULL)
}

.sn_usage_make_wrapper <- function(workflow, original) {
  wrapper <- original
  original_body <- body(original)
  selector_names <- intersect(
    names(formals(original)),
    .sn_usage_selector_parameters
  )
  if ("..." %in% names(formals(original))) {
    selector_names <- union(selector_names, .sn_usage_selector_parameters)
  }
  body(wrapper) <- substitute(
    {
      .sn_usage_token <- .sn_usage_begin(
        WORKFLOW,
        match.call(expand.dots = FALSE),
        evaluation_frame = environment(),
        selector_names = SELECTOR_NAMES
      )
      .sn_usage_status <- "ok"
      .sn_usage_condition <- NULL
      tryCatch(
        withCallingHandlers(
          ORIGINAL_BODY,
          warning = function(condition) {
            .sn_usage_count_warning(.sn_usage_token)
          }
        ),
        interrupt = function(condition) {
          .sn_usage_status <<- "interrupted"
          .sn_usage_condition <<- condition
          stop(condition)
        },
        error = function(condition) {
          .sn_usage_status <<- "error"
          .sn_usage_condition <<- condition
          stop(condition)
        },
        finally = .sn_usage_finish(
          .sn_usage_token,
          status = .sn_usage_status,
          condition = .sn_usage_condition
        )
      )
    },
    list(
      WORKFLOW = workflow,
      ORIGINAL_BODY = original_body,
      SELECTOR_NAMES = selector_names
    )
  )
  environment(wrapper) <- environment(original)
  wrapper
}

.sn_usage_replace_binding <- function(environment, name, value) {
  if (!exists(name, envir = environment, inherits = FALSE)) return(FALSE)
  active <- bindingIsActive(name, environment)
  if (isTRUE(active)) return(FALSE)
  locked <- unname(rlang::env_binding_are_locked(environment, name))
  if (locked) rlang::env_binding_unlock(environment, name)
  on.exit(if (locked) rlang::env_binding_lock(environment, name), add = TRUE)
  assign(name, value, envir = environment)
  TRUE
}

.sn_usage_instrument <- function(functions = NULL) {
  namespace <- asNamespace("Shennong")
  attached <- if ("package:Shennong" %in% search()) {
    as.environment("package:Shennong")
  } else {
    NULL
  }
  registry <- .sn_usage_registry(functions = functions)
  bindings <- list()
  complete <- FALSE
  on.exit({
    if (!complete && length(bindings) > 0L) {
      .sn_usage_state$bindings <- bindings
      .sn_usage_restore_instrumentation()
    }
  }, add = TRUE)
  for (workflow in registry$workflow) {
    original <- get(workflow, envir = namespace, inherits = FALSE)
    if (!is.function(original)) next
    wrapper <- .sn_usage_make_wrapper(workflow, original)
    namespace_changed <- .sn_usage_replace_binding(namespace, workflow, wrapper)
    if (!isTRUE(namespace_changed)) next
    attached_original <- NULL
    attached_changed <- FALSE
    if (is.environment(attached) &&
        exists(workflow, envir = attached, inherits = FALSE)) {
      attached_original <- get(workflow, envir = attached, inherits = FALSE)
      attached_changed <- .sn_usage_replace_binding(attached, workflow, wrapper)
    }
    bindings[[workflow]] <- list(
      original = original,
      wrapper = wrapper,
      attached_original = attached_original,
      attached_changed = attached_changed
    )
    .sn_usage_state$bindings <- bindings
  }
  .sn_usage_state$bindings <- bindings
  complete <- TRUE
  invisible(names(bindings))
}

.sn_usage_restore_instrumentation <- function() {
  bindings <- .sn_usage_state$bindings
  if (length(bindings) == 0L) return(invisible(character()))
  namespace <- asNamespace("Shennong")
  attached <- if ("package:Shennong" %in% search()) {
    as.environment("package:Shennong")
  } else {
    NULL
  }
  restored <- character()
  for (workflow in rev(names(bindings))) {
    binding <- bindings[[workflow]]
    current <- get(workflow, envir = namespace, inherits = FALSE)
    if (identical(current, binding$wrapper)) {
      .sn_usage_replace_binding(namespace, workflow, binding$original)
      restored <- c(restored, workflow)
    }
    if (is.environment(attached) && isTRUE(binding$attached_changed) &&
        exists(workflow, envir = attached, inherits = FALSE)) {
      current_attached <- get(workflow, envir = attached, inherits = FALSE)
      if (identical(current_attached, binding$wrapper)) {
        .sn_usage_replace_binding(
          attached,
          workflow,
          binding$attached_original
        )
      }
    }
  }
  .sn_usage_state$bindings <- list()
  invisible(restored)
}

.sn_usage_resolve_query_store <- function(path = NULL, store = NULL) {
  if (!is.null(store)) {
    if (!is.null(path)) stop("Supply either `path` or `store`, not both.", call. = FALSE)
    return(.sn_usage_validate_store(store))
  }
  if (is.null(path)) {
    if (!isTRUE(.sn_usage_state$enabled) || is.null(.sn_usage_state$config$store)) {
      stop(
        "Supply `path` or `store`, or enable usage tracking in this R process first.",
        call. = FALSE
      )
    }
    return(.sn_usage_state$config$store)
  }
  structure(
    list(
      backend = "sqlite",
      path = .sn_usage_normalize_path(path, must_exist = TRUE),
      connect = NULL,
      table_prefix = "shennong_usage",
      allow_remote_ddl = FALSE
    ),
    class = "sn_usage_store"
  )
}

.sn_usage_remote_connection <- function(store) {
  connection <- store$connect()
  if (!inherits(connection, "DBIConnection") || !DBI::dbIsValid(connection)) {
    try(DBI::dbDisconnect(connection), silent = TRUE)
    stop("The usage store connection factory did not return a valid DBI connection.", call. = FALSE)
  }
  connection
}

.sn_usage_remote_table_names <- function(store) {
  c(
    sessions = paste0(store$table_prefix, "_sessions"),
    runs = paste0(store$table_prefix, "_runs")
  )
}

.sn_usage_remote_session_prototype <- function() {
  data.frame(
    session_id = character(), started_at = character(), mode = character(),
    shennong_version = character(), r_version = character(),
    packages_json = character(), consent_id = character(),
    consent_scope = character(), policy_id = character(),
    policy_version = character(), study_id = character(),
    instrumentation_version = integer(), stringsAsFactors = FALSE
  )
}

.sn_usage_remote_run_prototype <- function() {
  data.frame(
    run_id = character(), session_id = character(), parent_run_id = character(),
    root_run_id = character(), depth = integer(), mode = character(),
    workflow = character(), category = character(), call_scope = character(),
    method = character(), backend = character(), status = character(),
    started_at = character(), finished_at = character(), elapsed_ms = numeric(),
    cpu_user_ms = numeric(), cpu_system_ms = numeric(), params_json = character(),
    params_sha256 = character(), warning_count = integer(),
    acceleration_json = character(), error_class = character(),
    instrumentation_version = integer(), study_id = character(),
    stringsAsFactors = FALSE
  )
}

.sn_usage_remote_existing_ids <- function(connection, table, column, ids) {
  ids <- unique(as.character(ids))
  ids <- ids[!is.na(ids) & nzchar(ids)]
  if (length(ids) == 0L) return(character())
  table_sql <- as.character(DBI::dbQuoteIdentifier(connection, table))
  column_sql <- as.character(DBI::dbQuoteIdentifier(connection, column))
  chunks <- split(ids, ceiling(seq_along(ids) / 500L))
  unique(unlist(lapply(chunks, function(current) {
    literals <- paste(
      as.character(DBI::dbQuoteString(connection, current)),
      collapse = ","
    )
    query <- paste0(
      "SELECT ", column_sql, " FROM ", table_sql,
      " WHERE ", column_sql, " IN (", literals, ")"
    )
    as.character(DBI::dbGetQuery(connection, query)[[column]])
  }), use.names = FALSE))
}

.sn_usage_initialize_remote_tables <- function(connection, store) {
  tables <- .sn_usage_remote_table_names(store)
  missing <- tables[!vapply(
    tables,
    function(table) DBI::dbExistsTable(connection, table),
    logical(1)
  )]
  if (length(missing) > 0L && !isTRUE(store$allow_remote_ddl)) {
    stop(
      "Remote usage tables are absent. Recreate the store with `allow_remote_ddl = TRUE` for an administrator-controlled initialization.",
      call. = FALSE
    )
  }
  if (tables[["sessions"]] %in% missing) {
    DBI::dbCreateTable(
      connection,
      tables[["sessions"]],
      .sn_usage_remote_session_prototype()
    )
  }
  if (tables[["runs"]] %in% missing) {
    DBI::dbCreateTable(
      connection,
      tables[["runs"]],
      .sn_usage_remote_run_prototype()
    )
  }
  if (isTRUE(store$allow_remote_ddl)) {
    unique_columns <- c(sessions = "session_id", runs = "run_id")
    for (kind in names(unique_columns)) {
      table <- tables[[kind]]
      column <- unique_columns[[kind]]
      index <- paste0(
        "sn_usage_",
        substr(
          digest::digest(paste(table, column, sep = "\r"), serialize = FALSE),
          1L,
          20L
        ),
        "_uidx"
      )
      DBI::dbExecute(
        connection,
        paste(
          "CREATE UNIQUE INDEX IF NOT EXISTS",
          as.character(DBI::dbQuoteIdentifier(connection, index)),
          "ON", as.character(DBI::dbQuoteIdentifier(connection, table)),
          paste0("(", as.character(DBI::dbQuoteIdentifier(connection, column)), ")")
        )
      )
    }
  }
  invisible(tables)
}

#' Create a local or managed-remote usage store
#'
#' A remote DBI store always retains a local SQLite outbox at `path`. Scientific
#' calls write only to the outbox; [sn_flush_usage_tracking()] performs the
#' explicit, consent-gated remote delivery. The `connect` closure and its
#' credentials are never serialized into the usage database.
#'
#' @param backend `"sqlite"` for local-only storage or `"dbi"` for a managed
#'   remote DBI destination backed by a local SQLite outbox.
#' @param path Local SQLite database/outbox path.
#' @param connect For `backend = "dbi"`, a zero-argument function returning a
#'   new valid `DBIConnection`. Do not return a connection opened before a fork.
#' @param table_prefix Remote table prefix.
#' @param allow_remote_ddl Permit an explicitly consented flush to create the
#'   two remote tables. Keep `FALSE` when migrations are administrator-owned.
#'
#' @return An `sn_usage_store` object. It contains a connection factory in
#'   memory and must not be serialized or committed.
#' @export
sn_create_usage_store <- function(
    backend = c("sqlite", "dbi"),
    path,
    connect = NULL,
    table_prefix = "shennong_usage",
    allow_remote_ddl = FALSE) {
  backend <- match.arg(backend)
  path <- .sn_usage_normalize_path(path)
  allow_remote_ddl <- .sn_usage_validate_flag(
    allow_remote_ddl,
    "allow_remote_ddl"
  )
  table_prefix <- .sn_usage_validate_table_prefix(table_prefix)
  if (identical(backend, "dbi") && !is.function(connect)) {
    stop("A remote DBI usage store requires a zero-argument `connect` function.", call. = FALSE)
  }
  if (identical(backend, "sqlite") && !is.null(connect)) {
    stop("`connect` is only used by `backend = \"dbi\"`.", call. = FALSE)
  }
  structure(
    list(
      backend = backend,
      path = path,
      connect = connect,
      table_prefix = table_prefix,
      allow_remote_ddl = allow_remote_ddl
    ),
    class = "sn_usage_store"
  )
}

#' Record explicit consent for usage research
#'
#' This constructor does not contact a server. It creates a bounded consent
#' receipt that is stored with the local session and required before any remote
#' DBI connection. The receipt binds its field categories and is checked for
#' modification before delivery. Shennong never creates a participant or
#' installation ID.
#'
#' @param scope `"local"` or `"remote_research"`.
#' @param policy_id,policy_version Research policy identifiers. Both are
#'   required for remote research.
#' @param purposes Non-empty research purposes.
#' @param data_categories Allowed remote fields. Supported values are
#'   `"api_name"`, `"timing"`, `"safe_parameters"`, `"package_versions"`,
#'   `"error_class"`, and `"acceleration"`.
#' @param expires_at Optional future `POSIXct` expiry.
#' @param study_id Optional short non-identifying study label. This is not a
#'   participant identifier.
#'
#' @return An immutable-style `sn_usage_consent` list.
#' @export
sn_confirm_usage_consent <- function(
    scope = c("local", "remote_research"),
    policy_id = NULL,
    policy_version = NULL,
    purposes = "performance_research",
    data_categories = c("api_name", "timing", "safe_parameters"),
    expires_at = NULL,
    study_id = NULL) {
  scope <- match.arg(scope)
  for (value in c(purposes, data_categories)) {
    if (!is.character(value) || length(value) == 0L || anyNA(value) ||
        any(!nzchar(trimws(value)))) {
      stop("Consent purposes and data categories must be non-empty character vectors.", call. = FALSE)
    }
  }
  data_categories <- unique(trimws(data_categories))
  unknown <- setdiff(data_categories, .sn_usage_allowed_categories)
  if (length(unknown) > 0L) {
    stop("Unsupported consent data category: ", paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  if (identical(scope, "remote_research")) {
    for (field in c("policy_id", "policy_version")) {
      value <- get(field)
      if (!is.character(value) || length(value) != 1L || is.na(value) ||
          !grepl("^[A-Za-z0-9_.-]{1,80}$", value)) {
        stop("Remote consent requires safe non-empty `", field, "`.", call. = FALSE)
      }
    }
  }
  if (!is.null(study_id) &&
      (!is.character(study_id) || length(study_id) != 1L || is.na(study_id) ||
       !grepl("^[A-Za-z0-9_.-]{1,80}$", study_id))) {
    stop("`study_id` must be NULL or a short non-identifying label.", call. = FALSE)
  }
  if (!is.null(expires_at)) {
    if (!inherits(expires_at, "POSIXt") || length(expires_at) != 1L ||
        is.na(expires_at) || expires_at <= Sys.time()) {
      stop("`expires_at` must be one future POSIXct value.", call. = FALSE)
    }
    expires_at <- as.POSIXct(expires_at, tz = "UTC")
  }
  granted_at <- .sn_usage_now()
  normalized_policy_id <- policy_id %||% "local"
  normalized_policy_version <- policy_version %||% "1"
  normalized_study_id <- study_id %||% NA_character_
  consent_id <- .sn_usage_consent_id(
    scope = scope,
    policy_id = normalized_policy_id,
    policy_version = normalized_policy_version,
    purposes = purposes,
    data_categories = data_categories,
    expires_at = expires_at,
    study_id = normalized_study_id
  )
  structure(
    list(
      consent_id = consent_id,
      scope = scope,
      policy_id = normalized_policy_id,
      policy_version = normalized_policy_version,
      purposes = unique(trimws(purposes)),
      data_categories = data_categories,
      granted_at = granted_at,
      expires_at = expires_at,
      study_id = normalized_study_id
    ),
    class = "sn_usage_consent"
  )
}

.sn_usage_consent_id_from_json <- function(value) {
  parsed <- tryCatch(
    jsonlite::fromJSON(value, simplifyVector = TRUE),
    error = function(error) list()
  )
  as.character(parsed$consent_id %||% NA_character_)
}

.sn_usage_pending_remote_rows <- function(store, limit, consent_id = NULL) {
  connection <- .sn_usage_connect(store$path)
  on.exit(DBI::dbDisconnect(connection), add = TRUE)
  sessions <- DBI::dbGetQuery(
    connection,
    paste(
      "SELECT sessions.* FROM sessions WHERE EXISTS (",
      "SELECT 1 FROM workflow_runs",
      "WHERE workflow_runs.session_id = sessions.session_id",
      "AND workflow_runs.remote_synced_at IS NULL",
      "AND workflow_runs.status <> 'running')"
    )
  )
  if (!is.null(consent_id) && nrow(sessions) > 0L) {
    recorded_ids <- vapply(
      sessions$consent_json,
      .sn_usage_consent_id_from_json,
      character(1)
    )
    sessions <- sessions[
      !is.na(recorded_ids) & recorded_ids == consent_id,
      ,
      drop = FALSE
    ]
  }
  if (nrow(sessions) == 0L) {
    runs <- DBI::dbGetQuery(
      connection,
      "SELECT * FROM workflow_runs WHERE 0 = 1"
    )
    return(list(runs = runs, sessions = sessions))
  }
  literals <- paste(
    as.character(DBI::dbQuoteString(connection, sessions$session_id)),
    collapse = ","
  )
  runs <- DBI::dbGetQuery(
    connection,
    paste0(
      "SELECT * FROM workflow_runs ",
      "WHERE remote_synced_at IS NULL AND status <> 'running' ",
      "AND session_id IN (", literals, ") ",
      "ORDER BY started_at, run_id LIMIT ?"
    ),
    params = list(as.integer(limit))
  )
  sessions <- sessions[
    match(unique(runs$session_id), sessions$session_id, nomatch = 0L),
    ,
    drop = FALSE
  ]
  list(runs = runs, sessions = sessions)
}

.sn_usage_remote_session_rows <- function(sessions, consent) {
  keep_versions <- "package_versions" %in% consent$data_categories
  keep_timing <- "timing" %in% consent$data_categories
  data.frame(
    session_id = as.character(sessions$session_id),
    started_at = if (keep_timing) as.character(sessions$started_at) else NA_character_,
    mode = as.character(sessions$mode),
    shennong_version = if (keep_versions) {
      as.character(sessions$shennong_version)
    } else {
      NA_character_
    },
    r_version = if (keep_versions) as.character(sessions$r_version) else NA_character_,
    packages_json = if (keep_versions) as.character(sessions$packages_json) else "{}",
    consent_id = consent$consent_id,
    consent_scope = consent$scope,
    policy_id = consent$policy_id,
    policy_version = consent$policy_version,
    study_id = consent$study_id,
    instrumentation_version = .sn_usage_instrumentation_version,
    stringsAsFactors = FALSE
  )
}

.sn_usage_remote_run_rows <- function(runs, consent) {
  keep_parameters <- "safe_parameters" %in% consent$data_categories
  keep_timing <- "timing" %in% consent$data_categories
  keep_errors <- "error_class" %in% consent$data_categories
  keep_acceleration <- "acceleration" %in% consent$data_categories
  data.frame(
    run_id = as.character(runs$run_id),
    session_id = as.character(runs$session_id),
    parent_run_id = as.character(runs$parent_run_id),
    root_run_id = as.character(runs$root_run_id),
    depth = as.integer(runs$depth),
    mode = as.character(runs$mode),
    workflow = as.character(runs$workflow),
    category = as.character(runs$category),
    call_scope = as.character(runs$call_scope),
    method = if (keep_parameters) as.character(runs$method) else NA_character_,
    backend = if (keep_parameters) as.character(runs$backend) else NA_character_,
    status = if (keep_errors) as.character(runs$status) else NA_character_,
    started_at = if (keep_timing) as.character(runs$started_at) else NA_character_,
    finished_at = if (keep_timing) as.character(runs$finished_at) else NA_character_,
    elapsed_ms = if (keep_timing) as.numeric(runs$elapsed_ms) else NA_real_,
    cpu_user_ms = if (keep_timing) as.numeric(runs$cpu_user_ms) else NA_real_,
    cpu_system_ms = if (keep_timing) as.numeric(runs$cpu_system_ms) else NA_real_,
    params_json = if (keep_parameters) as.character(runs$params_json) else "{}",
    params_sha256 = if (keep_parameters) {
      as.character(runs$params_sha256)
    } else {
      rep(
        digest::digest("{}", algo = "sha256", serialize = FALSE),
        nrow(runs)
      )
    },
    warning_count = if (keep_errors) as.integer(runs$warning_count) else NA_integer_,
    acceleration_json = if (keep_acceleration) {
      as.character(runs$acceleration_json)
    } else {
      "{}"
    },
    error_class = if (keep_errors) as.character(runs$error_class) else NA_character_,
    instrumentation_version = as.integer(runs$instrumentation_version),
    study_id = consent$study_id,
    stringsAsFactors = FALSE
  )
}

.sn_usage_mark_remote_sync <- function(store, run_ids, error = NULL) {
  connection <- .sn_usage_connect(store$path)
  on.exit(DBI::dbDisconnect(connection), add = TRUE)
  if (length(run_ids) == 0L) return(invisible(0L))
  for (run_id in run_ids) {
    if (is.null(error)) {
      DBI::dbExecute(
        connection,
        paste(
          "UPDATE workflow_runs SET remote_synced_at = ?,",
          "remote_sync_error = NULL WHERE run_id = ?"
        ),
        params = list(.sn_usage_now(), run_id)
      )
    } else {
      DBI::dbExecute(
        connection,
        "UPDATE workflow_runs SET remote_sync_error = ? WHERE run_id = ?",
        params = list(.sn_usage_redact_error(error), run_id)
      )
    }
  }
  invisible(length(run_ids))
}

#' Flush the local usage outbox to a managed remote DBI store
#'
#' Delivery is explicit, consent-gated, and idempotent by `session_id` and
#' `run_id`. Only sessions recorded under the exact consent receipt are sent;
#' other receipts remain pending. Optional remote fields follow the receipt's
#' data categories. Remote error messages and process identifiers are never
#' uploaded. A failed flush leaves rows in the local outbox for later retry.
#'
#' @param store A `backend = "dbi"` store from [sn_create_usage_store()].
#' @param consent Active `remote_research` consent from
#'   [sn_confirm_usage_consent()]. It must exactly match the recorded receipt.
#' @param limit Maximum completed outbox rows delivered in one call.
#' @param strict Stop on a remote error. `FALSE` warns and retains the outbox.
#'
#' @return A one-row delivery summary.
#' @export
sn_flush_usage_tracking <- function(store, consent, limit = 1000L, strict = TRUE) {
  .sn_usage_require_database()
  store <- .sn_usage_validate_store(store)
  if (!identical(store$backend, "dbi")) {
    stop("Only a `backend = \"dbi\"` usage store has a remote outbox.", call. = FALSE)
  }
  consent <- .sn_usage_validate_consent(consent, remote = TRUE)
  if (!"api_name" %in% consent$data_categories) {
    stop("Remote research consent must include the `api_name` data category.", call. = FALSE)
  }
  limit <- suppressWarnings(as.integer(limit))
  if (length(limit) != 1L || is.na(limit) || limit < 1L) {
    stop("`limit` must be one positive integer.", call. = FALSE)
  }
  .sn_usage_initialize_database(store$path)
  pending <- .sn_usage_pending_remote_rows(
    store,
    limit,
    consent_id = consent$consent_id
  )
  if (nrow(pending$runs) == 0L) {
    return(invisible(data.frame(
      pending = 0L, appended = 0L, already_present = 0L, synced = 0L
    )))
  }
  recorded_consent_ids <- vapply(
    pending$sessions$consent_json,
    .sn_usage_consent_id_from_json,
    character(1)
  )
  recorded_consent_ids <- unique(recorded_consent_ids[!is.na(recorded_consent_ids)])
  if (length(recorded_consent_ids) != 1L ||
      !identical(recorded_consent_ids, consent$consent_id)) {
    stop("The supplied consent does not match the pending local outbox session.", call. = FALSE)
  }
  result <- tryCatch({
    connection <- .sn_usage_remote_connection(store)
    on.exit(DBI::dbDisconnect(connection), add = TRUE)
    tables <- .sn_usage_initialize_remote_tables(connection, store)
    session_rows <- .sn_usage_remote_session_rows(pending$sessions, consent)
    run_rows <- .sn_usage_remote_run_rows(pending$runs, consent)
    existing_sessions <- .sn_usage_remote_existing_ids(
      connection,
      tables[["sessions"]],
      "session_id",
      session_rows$session_id
    )
    existing_runs <- .sn_usage_remote_existing_ids(
      connection,
      tables[["runs"]],
      "run_id",
      run_rows$run_id
    )
    new_sessions <- session_rows[!session_rows$session_id %in% existing_sessions, , drop = FALSE]
    new_runs <- run_rows[!run_rows$run_id %in% existing_runs, , drop = FALSE]
    DBI::dbWithTransaction(connection, {
      if (nrow(new_sessions) > 0L) {
        DBI::dbAppendTable(connection, tables[["sessions"]], new_sessions)
      }
      if (nrow(new_runs) > 0L) {
        DBI::dbAppendTable(connection, tables[["runs"]], new_runs)
      }
    })
    .sn_usage_mark_remote_sync(store, run_rows$run_id)
    data.frame(
      pending = nrow(run_rows),
      appended = nrow(new_runs),
      already_present = length(existing_runs),
      synced = nrow(run_rows),
      stringsAsFactors = FALSE
    )
  }, error = function(error) {
    .sn_usage_mark_remote_sync(store, pending$runs$run_id, error = error)
    if (isTRUE(strict)) stop(error)
    warning(
      "Remote usage flush failed; records remain in the local outbox: ",
      conditionMessage(error),
      call. = FALSE
    )
    data.frame(
      pending = nrow(pending$runs), appended = 0L,
      already_present = 0L, synced = 0L,
      stringsAsFactors = FALSE
    )
  })
  invisible(result)
}

#' Enable local workflow usage and timing records
#'
#' Instruments Shennong's high-level computational workflows for the current R
#' process and writes one local SQLite row per call. Tracking is always opt-in:
#' loading or attaching Shennong never creates a database, and no record is sent
#' over the network. The parameter recorder stores conservative scalar settings
#' but redacts objects, matrices, identifiers, paths, credentials, prompts, and
#' free text.
#'
#' @param path Local SQLite file path. It is retained for compatibility and
#'   cannot be combined with `store`.
#' @param store Optional object from [sn_create_usage_store()]. A remote DBI
#'   store still writes first to its local SQLite outbox.
#' @param consent Optional object from [sn_confirm_usage_consent()]. It is
#'   mandatory for a remote store.
#' @param remote_flush For a DBI store, flush manually or when tracking is
#'   disabled. Remote failures on disable warn and retain the local outbox.
#' @param mode Explicit execution context: `"development"`, `"production"`,
#'   `"test"`, or `"benchmark"`. Shennong never infers that a run is production.
#' @param display Show a compact elapsed-time message when each tracked call
#'   finishes.
#' @param nested Record nested high-level Shennong workflow calls. Summaries
#'   exclude nested calls by default to avoid double-counting elapsed time.
#' @param functions Optional exact character vector of exported functions to
#'   instrument. `NULL` records every public Shennong and rio adapter except the
#'   usage-tracking control/query functions themselves.
#' @param strict If `TRUE`, a database write failure stops the analysis. The
#'   default is fail-open: warn once and preserve the scientific call.
#'
#' @return Invisibly, the result of [sn_check_usage_tracking()].
#' @export
#'
#' @examples
#' \dontrun{
#' database <- file.path(tempdir(), "shennong-usage.sqlite")
#' sn_enable_usage_tracking(database, mode = "development")
#' object <- sn_run_cluster(object, integration_method = "unintegrated")
#' sn_disable_usage_tracking()
#' sn_summarize_usage(database)
#' }
sn_enable_usage_tracking <- function(
    path = NULL,
    store = NULL,
    consent = NULL,
    mode = c("development", "production", "test", "benchmark"),
    display = FALSE,
    nested = TRUE,
    functions = NULL,
    remote_flush = c("manual", "on_disable"),
    strict = FALSE) {
  if (isTRUE(.sn_usage_state$enabled)) {
    stop(
      "Usage tracking is already enabled; disable it before changing its configuration.",
      call. = FALSE
    )
  }
  .sn_usage_require_database()
  store <- .sn_usage_resolve_store(path = path, store = store)
  remote <- identical(store$backend, "dbi")
  consent <- .sn_usage_validate_consent(consent, remote = remote)
  path <- store$path
  mode <- match.arg(mode)
  remote_flush <- match.arg(remote_flush)
  display <- .sn_usage_validate_flag(display, "display")
  nested <- .sn_usage_validate_flag(nested, "nested")
  strict <- .sn_usage_validate_flag(strict, "strict")
  selected <- .sn_usage_registry(functions = functions)$workflow
  .sn_usage_initialize_database(path)

  config <- list(
    path = path,
    mode = mode,
    display = display,
    nested = nested,
    functions = selected,
    store = store,
    consent = consent,
    remote_flush = remote_flush,
    strict = strict,
    schema_version = .sn_usage_schema_version
  )
  session_id <- .sn_usage_create_session(
    path,
    mode,
    list(
      mode = mode,
      display = display,
      nested = nested,
      strict = strict,
      schema_version = .sn_usage_schema_version,
      instrumented_functions = length(selected),
      destination_kind = store$backend,
      remote_flush = remote_flush
    ),
    consent = consent,
    destination_kind = store$backend
  )
  .sn_usage_state$config <- config
  .sn_usage_state$session_id <- session_id
  .sn_usage_state$pid <- as.integer(Sys.getpid())
  .sn_usage_state$stack <- list()
  .sn_usage_state$warned <- character()
  .sn_usage_state$enabled <- TRUE
  instrumented <- tryCatch(
    .sn_usage_instrument(functions = selected),
    error = function(error) {
      .sn_usage_state$enabled <- FALSE
      .sn_usage_finish_session(path, session_id, strict = FALSE)
      .sn_usage_state$config <- NULL
      .sn_usage_state$session_id <- NULL
      stop("Could not instrument Shennong workflows: ", conditionMessage(error), call. = FALSE)
    }
  )
  if (length(instrumented) == 0L) {
    sn_disable_usage_tracking()
    stop("No Shennong workflow functions were available to instrument.", call. = FALSE)
  }
  invisible(sn_check_usage_tracking())
}

#' Disable local workflow usage tracking
#'
#' Restores the original Shennong function bindings and closes the logical
#' tracking session. There is no persistent database connection to close.
#'
#' @return Invisibly, `TRUE` when tracking was disabled and `FALSE` when it was
#'   already disabled.
#' @export
sn_disable_usage_tracking <- function() {
  if (!isTRUE(.sn_usage_state$enabled)) return(invisible(FALSE))
  if (length(.sn_usage_state$stack) > 0L) {
    stop(
      "Cannot disable usage tracking while a tracked function is still active.",
      call. = FALSE
    )
  }
  config <- .sn_usage_state$config
  session_id <- .sn_usage_state$session_id
  .sn_usage_restore_instrumentation()
  .sn_usage_state$enabled <- FALSE
  .sn_usage_finish_session(config$path, session_id, strict = config$strict)
  if (identical(config$store$backend, "dbi") &&
      identical(config$remote_flush, "on_disable")) {
    tryCatch(
      sn_flush_usage_tracking(
        store = config$store,
        consent = config$consent,
        strict = FALSE
      ),
      error = function(error) warning(
        "Remote usage flush failed; records remain in the local outbox: ",
        conditionMessage(error),
        call. = FALSE
      )
    )
  }
  .sn_usage_state$config <- NULL
  .sn_usage_state$session_id <- NULL
  .sn_usage_state$pid <- NULL
  .sn_usage_state$stack <- list()
  invisible(TRUE)
}

#' Check local workflow usage tracking state
#'
#' @return A one-row data frame describing the current process-local tracking
#'   configuration.
#' @export
sn_check_usage_tracking <- function() {
  config <- .sn_usage_state$config
  data.frame(
    enabled = isTRUE(.sn_usage_state$enabled),
    path = config$path %||% NA_character_,
    destination = config$store$backend %||% NA_character_,
    consent_scope = config$consent$scope %||% NA_character_,
    study_id = config$consent$study_id %||% NA_character_,
    mode = config$mode %||% NA_character_,
    session_id = .sn_usage_state$session_id %||% NA_character_,
    display = config$display %||% NA,
    nested = config$nested %||% NA,
    strict = config$strict %||% NA,
    instrumented_functions = length(.sn_usage_state$bindings),
    coverage_model = "namespace_bindings",
    covers_preexisting_function_references = FALSE,
    stringsAsFactors = FALSE
  )
}

#' Evaluate code with temporary local usage tracking
#'
#' @param expr Code to evaluate.
#' @inheritParams sn_enable_usage_tracking
#'
#' @return The value of `expr`, with its visibility preserved.
#' @export
sn_with_usage_tracking <- function(
    expr,
    path = NULL,
    store = NULL,
    consent = NULL,
    mode = c("development", "production", "test", "benchmark"),
    display = FALSE,
    nested = TRUE,
    functions = NULL,
    remote_flush = c("manual", "on_disable"),
    strict = FALSE) {
  if (isTRUE(.sn_usage_state$enabled)) {
    return(force(expr))
  }
  sn_enable_usage_tracking(
    path = path,
    store = store,
    consent = consent,
    mode = match.arg(mode),
    display = display,
    nested = nested,
    functions = functions,
    remote_flush = match.arg(remote_flush),
    strict = strict
  )
  on.exit(sn_disable_usage_tracking(), add = TRUE)
  force(expr)
}

#' Time an R expression
#'
#' Evaluates an expression once, optionally displays elapsed time, and records
#' it as an ad-hoc workflow when usage tracking is enabled. The expression's
#' value, visibility, warnings, and errors are preserved.
#'
#' @param expr Code to evaluate.
#' @param label Short, non-sensitive label used in timing output and the local
#'   usage database. It must use letters, numbers, dots, underscores, or dashes.
#' @param display Show the elapsed-time message.
#' @param record Record the call when usage tracking is enabled.
#'
#' @return The value of `expr`, with its visibility preserved.
#' @export
#'
#' @examples
#' result <- sn_time_call(sum(seq_len(1000)), label = "sum-example")
sn_time_call <- function(expr, label = "ad_hoc", display = TRUE, record = TRUE) {
  if (!is.character(label) || length(label) != 1L || is.na(label) ||
      !grepl("^[A-Za-z0-9_.-]+$", label)) {
    stop("`label` must be one short non-sensitive identifier.", call. = FALSE)
  }
  display <- .sn_usage_validate_flag(display, "display")
  record <- .sn_usage_validate_flag(record, "record")
  workflow <- paste0("timed:", label)
  token <- if (record && isTRUE(.sn_usage_state$enabled)) {
    .sn_usage_begin(workflow, call("sn_time_call", label = label), category = "ad_hoc")
  } else {
    token <- new.env(parent = emptyenv())
    token$started_proc <- proc.time()
    token$finished <- FALSE
    token$recorded <- FALSE
    token$tracking_active <- FALSE
    token$workflow <- workflow
    token$warning_count <- 0L
    token$acceleration_patches <- character()
    token
  }
  status <- "ok"
  captured <- NULL
  tryCatch(
    withCallingHandlers(
      force(expr),
      warning = function(condition) .sn_usage_count_warning(token)
    ),
    interrupt = function(condition) {
      status <<- "interrupted"
      captured <<- condition
      stop(condition)
    },
    error = function(condition) {
      status <<- "error"
      captured <<- condition
      stop(condition)
    },
    finally = {
      if (record && isTRUE(.sn_usage_state$enabled)) {
        previous <- .sn_usage_state$config$display
        .sn_usage_state$config$display <- display
        on.exit(.sn_usage_state$config$display <- previous, add = TRUE)
        .sn_usage_finish(token, status = status, condition = captured)
      } else if (display) {
        seconds <- unname((proc.time() - token$started_proc)[["elapsed"]])
        message(sprintf("[Shennong timing] %s: %.3f s [%s]", workflow, seconds, status))
      }
    }
  )
}

#' List local workflow usage records
#'
#' @param path Existing usage SQLite file. It may be omitted while tracking is
#'   enabled in the current process.
#' @param store Optional usage store. Use `source = "remote"` with a DBI store
#'   and reader-authorized connection factory.
#' @param source Read the local outbox or the managed remote table.
#' @param workflow Optional exact workflow name(s).
#' @param mode Optional execution mode(s).
#' @param status Optional run status values.
#' @param top_level If `TRUE`, keep only root calls; if `FALSE`, keep only
#'   nested calls; `NULL` keeps both.
#' @param limit Maximum number of newest rows returned.
#'
#' @return A data frame ordered from newest to oldest. `invocation_number` is
#'   the use number for the workflow; `parameter_set_invocation` is the use
#'   number for the same sanitized parameter fingerprint.
#' @export
sn_list_usage_runs <- function(
    path = NULL,
    store = NULL,
    source = c("local", "remote"),
    workflow = NULL,
    mode = NULL,
    status = NULL,
    top_level = NULL,
    limit = 1000L) {
  .sn_usage_require_database()
  store <- .sn_usage_resolve_query_store(path = path, store = store)
  source <- match.arg(source)
  if (!is.null(top_level)) .sn_usage_validate_flag(top_level, "top_level")
  limit <- suppressWarnings(as.integer(limit))
  if (length(limit) != 1L || is.na(limit) || limit < 1L) {
    stop("`limit` must be one positive integer.", call. = FALSE)
  }
  if (identical(source, "remote")) {
    if (!identical(store$backend, "dbi")) {
      stop("`source = \"remote\"` requires a DBI usage store.", call. = FALSE)
    }
    connection <- .sn_usage_remote_connection(store)
    on.exit(DBI::dbDisconnect(connection), add = TRUE)
    table <- .sn_usage_remote_table_names(store)[["runs"]]
    if (!DBI::dbExistsTable(connection, table)) {
      stop("The remote usage run table does not exist.", call. = FALSE)
    }
    clauses <- character()
    add_remote_filter <- function(column, values) {
      if (is.null(values)) return(invisible(NULL))
      if (!is.character(values) || anyNA(values) || any(!nzchar(values))) {
        stop("Usage query filters must be non-empty character vectors.", call. = FALSE)
      }
      column_sql <- as.character(DBI::dbQuoteIdentifier(connection, column))
      literals <- paste(
        as.character(DBI::dbQuoteString(connection, values)),
        collapse = ","
      )
      clauses <<- c(clauses, paste0(column_sql, " IN (", literals, ")"))
      invisible(NULL)
    }
    add_remote_filter("workflow", workflow)
    add_remote_filter("mode", mode)
    add_remote_filter("status", status)
    if (!is.null(top_level)) {
      clauses <- c(clauses, if (top_level) "depth = 0" else "depth > 0")
    }
    where <- if (length(clauses)) {
      paste("WHERE", paste(clauses, collapse = " AND "))
    } else {
      ""
    }
    query <- paste(
      "SELECT * FROM (SELECT usage_rows.* ,",
      "ROW_NUMBER() OVER (PARTITION BY workflow ORDER BY started_at, run_id)",
      "AS invocation_number,",
      "ROW_NUMBER() OVER (PARTITION BY workflow, params_sha256",
      "ORDER BY started_at, run_id) AS parameter_set_invocation",
      "FROM", as.character(DBI::dbQuoteIdentifier(connection, table)),
      "AS usage_rows) AS usage_ranked",
      where,
      "ORDER BY started_at DESC LIMIT", as.integer(limit)
    )
    rows <- DBI::dbGetQuery(connection, query)
    rows$elapsed_seconds <- as.numeric(rows$elapsed_ms) / 1000
    rows$cpu_user_seconds <- as.numeric(rows$cpu_user_ms) / 1000
    rows$cpu_system_seconds <- as.numeric(rows$cpu_system_ms) / 1000
    rows$invocation_number <- as.integer(rows$invocation_number)
    rows$parameter_set_invocation <- as.integer(rows$parameter_set_invocation)
    rows$remote_synced_at <- NA_character_
    rows$remote_sync_error <- NA_character_
    rownames(rows) <- NULL
    return(rows)
  }
  path <- .sn_usage_normalize_path(store$path, must_exist = TRUE)
  connection <- .sn_usage_connect(path)
  on.exit(DBI::dbDisconnect(connection), add = TRUE)
  clauses <- character()
  parameters <- list()
  add_filter <- function(column, values) {
    if (is.null(values)) return(invisible(NULL))
    if (!is.character(values) || anyNA(values) || any(!nzchar(values))) {
      stop("Usage query filters must be non-empty character vectors.", call. = FALSE)
    }
    clauses <<- c(
      clauses,
      paste0(column, " IN (", paste(rep("?", length(values)), collapse = ","), ")")
    )
    parameters <<- c(parameters, as.list(values))
    invisible(NULL)
  }
  add_filter("workflow", workflow)
  add_filter("mode", mode)
  add_filter("status", status)
  if (!is.null(top_level)) {
    clauses <- c(clauses, if (top_level) "depth = 0" else "depth > 0")
  }
  where <- if (length(clauses)) {
    paste("WHERE", paste(clauses, collapse = " AND "))
  } else {
    ""
  }
  parameters <- c(parameters, list(limit))
  rows <- DBI::dbGetQuery(
    connection,
    paste(
      "SELECT run_id, session_id, parent_run_id, root_run_id, depth, pid,",
      "mode, workflow, category, call_scope, method, backend, status, started_at,",
      "finished_at, elapsed_ms / 1000.0 AS elapsed_seconds,",
      "cpu_user_ms / 1000.0 AS cpu_user_seconds,",
      "cpu_system_ms / 1000.0 AS cpu_system_seconds, params_json,",
      "params_sha256, invocation_number, parameter_set_invocation,",
      "warning_count, acceleration_json, error_class, error_message_redacted,",
      "remote_synced_at, remote_sync_error",
      "FROM workflow_runs", where, "ORDER BY started_at DESC LIMIT ?"
    ),
    params = parameters
  )
  rownames(rows) <- NULL
  rows
}

#' Summarize local workflow usage and runtime
#'
#' @inheritParams sn_list_usage_runs
#' @param by_parameters Group separate sanitized parameter fingerprints. This
#'   is useful for comparing calls such as different `sn_run_cluster()` method
#'   settings.
#' @param sort_by Rank summaries by call count, cumulative elapsed time, or
#'   median elapsed time.
#'
#' @return A data frame ordered by call count and total elapsed time.
#' @export
sn_summarize_usage <- function(
    path = NULL,
    store = NULL,
    source = c("local", "remote"),
    workflow = NULL,
    mode = NULL,
    top_level = TRUE,
    by_parameters = FALSE,
    sort_by = c("calls", "total_seconds", "median_seconds")) {
  by_parameters <- .sn_usage_validate_flag(by_parameters, "by_parameters")
  sort_by <- match.arg(sort_by)
  source <- match.arg(source)
  rows <- sn_list_usage_runs(
    path = path,
    store = store,
    source = source,
    workflow = workflow,
    mode = mode,
    top_level = top_level,
    limit = .Machine$integer.max
  )
  if (nrow(rows) == 0L) {
    return(data.frame(
      workflow = character(), mode = character(), method = character(),
      calls = integer(), successful = integer(), errors = integer(),
      incomplete = integer(),
      median_seconds = numeric(), p95_seconds = numeric(),
      total_seconds = numeric(), first_used = character(), last_used = character(),
      stringsAsFactors = FALSE
    ))
  }
  group_names <- c("workflow", "mode")
  if (by_parameters) group_names <- c(group_names, "params_sha256", "method")
  grouping <- rows[group_names]
  grouping[] <- lapply(grouping, function(value) {
    value <- as.character(value)
    value[is.na(value)] <- "<none>"
    value
  })
  keys <- interaction(grouping, drop = TRUE, lex.order = TRUE)
  groups <- split(seq_len(nrow(rows)), keys)
  result <- lapply(groups, function(index) {
    selected <- rows[index, , drop = FALSE]
    elapsed <- selected$elapsed_seconds[is.finite(selected$elapsed_seconds)]
    output <- list(
      workflow = selected$workflow[[1L]],
      mode = selected$mode[[1L]]
    )
    if (by_parameters) {
      output$params_sha256 <- selected$params_sha256[[1L]]
      output$method <- selected$method[[1L]]
    } else {
      output$method <- NA_character_
    }
    output$calls <- nrow(selected)
    output$successful <- sum(selected$status == "ok")
    output$errors <- sum(selected$status %in% c("error", "interrupted"))
    output$incomplete <- sum(selected$status %in% c("running", "abandoned"))
    output$median_seconds <- if (length(elapsed)) stats::median(elapsed) else NA_real_
    output$p95_seconds <- if (length(elapsed)) {
      unname(stats::quantile(elapsed, 0.95, names = FALSE, type = 8))
    } else {
      NA_real_
    }
    output$total_seconds <- if (length(elapsed)) sum(elapsed) else 0
    output$first_used <- min(selected$started_at)
    output$last_used <- max(selected$started_at)
    as.data.frame(output, stringsAsFactors = FALSE)
  })
  result <- do.call(rbind, result)
  rownames(result) <- NULL
  primary <- result[[sort_by]]
  result[order(-primary, -result$total_seconds, result$workflow), , drop = FALSE]
}
