usage_test_metadata <- function() {
  data.frame(
    cluster = c("0", "0", "1", "1"),
    batch = c("a", "b", "a", "b"),
    stringsAsFactors = FALSE
  )
}

test_that("usage tracking is opt-in and restores workflow bindings", {
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")
  on.exit(Shennong::sn_disable_usage_tracking(), add = TRUE)

  path <- tempfile(fileext = ".sqlite")
  original <- Shennong::sn_calculate_cluster_entropy
  expect_false(Shennong::sn_check_usage_tracking()$enabled)
  result <- Shennong::sn_calculate_cluster_entropy(
    usage_test_metadata(),
    cluster_by = "cluster",
    label_by = "batch"
  )
  expect_equal(nrow(result), 2L)
  expect_false(file.exists(path))

  state <- Shennong::sn_enable_usage_tracking(
    path,
    mode = "test",
    display = FALSE,
    strict = TRUE
  )
  expect_true(state$enabled)
  expect_equal(
    state$instrumented_functions,
    nrow(Shennong:::.sn_usage_registry())
  )
  expect_false(identical(original, Shennong::sn_calculate_cluster_entropy))
  expect_true(Shennong::sn_disable_usage_tracking())
  expect_identical(Shennong::sn_calculate_cluster_entropy, original)
  expect_false(Shennong::sn_disable_usage_tracking())

  connection <- DBI::dbConnect(RSQLite::SQLite(), path)
  on.exit(DBI::dbDisconnect(connection), add = TRUE)
  sessions <- DBI::dbReadTable(connection, "sessions")
  expect_equal(nrow(sessions), 1L)
  expect_match(sessions$packages_json, "Shennong")
  expect_match(sessions$packages_json, "Matrix")
  expect_false(is.na(sessions$finished_at))
})

test_that("usage records preserve results, RNG, errors, and invocation numbers", {
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")
  on.exit(Shennong::sn_disable_usage_tracking(), add = TRUE)

  path <- tempfile(fileext = ".sqlite")
  set.seed(717)
  rng_before <- .Random.seed
  Shennong::sn_enable_usage_tracking(
    path,
    mode = "test",
    display = FALSE,
    strict = TRUE
  )
  expected <- Shennong::sn_disable_usage_tracking()
  expect_true(expected)
  untracked <- Shennong::sn_calculate_cluster_entropy(
    usage_test_metadata(),
    cluster_by = "cluster",
    label_by = "batch"
  )

  Shennong::sn_enable_usage_tracking(
    path,
    mode = "test",
    display = FALSE,
    strict = TRUE
  )
  tracked_one <- Shennong::sn_calculate_cluster_entropy(
    usage_test_metadata(),
    cluster_by = "cluster",
    label_by = "batch"
  )
  tracked_two <- Shennong::sn_calculate_cluster_entropy(
    usage_test_metadata(),
    cluster_by = "cluster",
    label_by = "batch"
  )
  captured <- expect_error(
    Shennong::sn_calculate_cluster_entropy(usage_test_metadata()),
    "must be supplied",
    class = "simpleError"
  )
  expect_equal(tracked_one, untracked)
  expect_equal(tracked_two, untracked)
  expect_identical(.Random.seed, rng_before)
  Shennong::sn_disable_usage_tracking()

  runs <- Shennong::sn_list_usage_runs(
    path,
    workflow = "sn_calculate_cluster_entropy"
  )
  runs <- runs[order(runs$invocation_number), , drop = FALSE]
  expect_identical(runs$invocation_number, 1:3)
  expect_identical(runs$parameter_set_invocation, c(1L, 2L, 1L))
  expect_identical(runs$status, c("ok", "ok", "error"))
  expect_true(all(is.finite(runs$elapsed_seconds)))
  expect_true(all(runs$elapsed_seconds >= 0))
  expect_match(runs$params_json[[1L]], '"cluster_by":"cluster"')
  expect_match(runs$error_class[[3L]], "simpleError")
  expect_match(runs$error_message_redacted[[3L]], "must be supplied")
})

test_that("tracking reports the pre-existing function-reference limitation", {
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")
  on.exit(Shennong::sn_disable_usage_tracking(), add = TRUE)

  saved <- Shennong::sn_calculate_cluster_entropy
  path <- tempfile(fileext = ".sqlite")
  state <- Shennong::sn_enable_usage_tracking(
    path,
    mode = "test",
    display = FALSE,
    functions = "sn_calculate_cluster_entropy"
  )
  expect_false(state$covers_preexisting_function_references)
  saved(usage_test_metadata(), cluster_by = "cluster", label_by = "batch")
  Shennong::sn_calculate_cluster_entropy(
    usage_test_metadata(),
    cluster_by = "cluster",
    label_by = "batch"
  )
  Shennong::sn_disable_usage_tracking()

  runs <- Shennong::sn_list_usage_runs(path)
  expect_equal(nrow(runs), 1L)
})

test_that("usage parameter records redact sensitive payloads", {
  call <- quote(sn_run_llm(
    provider = "openai",
    model = "gpt-5",
    prompt = "patient-123 secret",
    api_key = "sk-secret",
    path = "/private/patient-123.rds"
  ))
  parameters <- Shennong:::.sn_usage_call_parameters("sn_run_llm", call)
  encoded <- as.character(Shennong:::.sn_usage_json(parameters))
  expect_match(encoded, "openai")
  expect_match(encoded, "gpt-5")
  expect_false(grepl("patient-123|sk-secret|/private", encoded))
  expect_false(any(c("prompt", "api_key", "path") %in% names(parameters)))

  generic <- Shennong:::.sn_usage_call_parameters(
    "sn_run_cluster",
    quote(sn_run_cluster(
      object = patient_object,
      integration_method = "harmony",
      sample_ids = c("patient-a", "patient-b"),
      store_name = "patient-a-analysis",
      output_dir = "/private/output"
    ))
  )
  generic_json <- as.character(Shennong:::.sn_usage_json(generic))
  expect_match(generic_json, "harmony")
  expect_false(grepl("patient-a|patient-b|/private/output", generic_json))

  generic_input <- Shennong:::.sn_usage_call_parameters(
    "sn_enrich",
    quote(sn_enrich(
      x = c("BRCA1", "patient-gene-123"),
      genes = c("TP53", "private-gene"),
      database = "H"
    ))
  )
  generic_input_json <- as.character(Shennong:::.sn_usage_json(generic_input))
  expect_match(generic_input_json, '"database":"H"')
  expect_false(grepl("BRCA1|patient-gene|TP53|private-gene", generic_input_json))

  condition <- simpleError(
    "failed /private/patient-123.rds token=sk-secret user@example.org"
  )
  redacted <- Shennong:::.sn_usage_redact_error(condition)
  expect_false(grepl("/private|patient-123|sk-secret|user@example", redacted))
  expect_match(redacted, "<path>")
  expect_match(redacted, "<email>")
  expect_match(redacted, "<redacted>")
})

test_that("usage resolves only allowlisted selector promises and effective defaults", {
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")
  skip_if_not_installed("Seurat")
  on.exit(Shennong::sn_disable_usage_tracking(), add = TRUE)

  path <- tempfile(fileext = ".sqlite")
  Shennong::sn_enable_usage_tracking(
    path,
    mode = "test",
    display = FALSE,
    functions = "sn_run_cluster",
    strict = TRUE
  )
  selected_method <- "harmony"
  selected_modality <- "rna"
  expect_error(
    Shennong::sn_run_cluster(
      NULL,
      integration_method = selected_method,
      modality = selected_modality
    ),
    "Seurat object"
  )
  Shennong::sn_disable_usage_tracking()

  runs <- Shennong::sn_list_usage_runs(path, workflow = "sn_run_cluster")
  expect_equal(nrow(runs), 1L)
  expect_match(runs$params_json, '"integration_method":"harmony"')
  expect_match(runs$params_json, '"normalization_method":"seurat"')
  expect_match(runs$params_json, '"cluster_algorithm":"louvain"')
  expect_match(runs$params_json, '"modality":"rna"')
  expect_false(grepl("selected_method|selected_modality", runs$params_json))
  expect_match(runs$method, "integration_method=harmony")
})

test_that("consent receipts are stable and reject post-confirmation expansion", {
  consent_one <- Shennong::sn_confirm_usage_consent(
    scope = "remote_research",
    policy_id = "stable-policy",
    policy_version = "1.0",
    purposes = c("performance_research", "usage_research"),
    data_categories = c("api_name", "timing"),
    study_id = "stable-study"
  )
  consent_two <- Shennong::sn_confirm_usage_consent(
    scope = "remote_research",
    policy_id = "stable-policy",
    policy_version = "1.0",
    purposes = rev(c("performance_research", "usage_research")),
    data_categories = rev(c("api_name", "timing")),
    study_id = "stable-study"
  )
  expect_identical(consent_one$consent_id, consent_two$consent_id)

  expanded <- consent_one
  expanded$data_categories <- c(expanded$data_categories, "safe_parameters")
  expect_error(
    Shennong:::.sn_usage_validate_consent(expanded, remote = TRUE),
    "modified"
  )
})

test_that("usage records nested parentage and acceleration activation honestly", {
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")
  on.exit(Shennong::sn_disable_usage_tracking(), add = TRUE)

  path <- tempfile(fileext = ".sqlite")
  Shennong::sn_enable_usage_tracking(
    path,
    mode = "benchmark",
    display = FALSE,
    nested = TRUE,
    strict = TRUE
  )
  parent <- Shennong:::.sn_usage_begin(
    "parent-workflow",
    quote(parent_workflow(method = "one")),
    category = "test"
  )
  child <- Shennong:::.sn_usage_begin(
    "child-workflow",
    quote(child_workflow(method = "two")),
    category = "test"
  )
  Shennong:::.sn_usage_record_acceleration("lisi")
  Shennong:::.sn_usage_count_warning(parent)
  Shennong:::.sn_usage_count_warning(child)
  expect_equal(parent$warning_count, 0L)
  expect_equal(child$warning_count, 1L)
  Shennong:::.sn_usage_finish(child)
  Shennong:::.sn_usage_finish(parent)
  Shennong::sn_disable_usage_tracking()

  runs <- Shennong::sn_list_usage_runs(path, limit = 10L)
  parent_row <- runs[runs$workflow == "parent-workflow", , drop = FALSE]
  child_row <- runs[runs$workflow == "child-workflow", , drop = FALSE]
  expect_equal(parent_row$depth, 0L)
  expect_equal(child_row$depth, 1L)
  expect_identical(child_row$parent_run_id, parent_row$run_id)
  expect_identical(child_row$root_run_id, parent_row$run_id)
  expect_match(child_row$acceleration_json, "lisi")
  expect_match(child_row$acceleration_json, "scope_activation_only")
  expect_match(child_row$acceleration_json, '"fast_path_hit":null')
  expect_equal(parent_row$warning_count, 0L)
  expect_equal(child_row$warning_count, 1L)
})

test_that("nested tracking can be disabled without affecting the parent", {
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")
  on.exit(Shennong::sn_disable_usage_tracking(), add = TRUE)

  path <- tempfile(fileext = ".sqlite")
  Shennong::sn_enable_usage_tracking(
    path,
    mode = "test",
    display = FALSE,
    nested = FALSE,
    strict = TRUE
  )
  parent <- Shennong:::.sn_usage_begin(
    "parent-only",
    quote(parent_only()),
    category = "test"
  )
  child <- Shennong:::.sn_usage_begin(
    "unrecorded-child",
    quote(unrecorded_child()),
    category = "test"
  )
  expect_false(child$tracking_active)
  Shennong:::.sn_usage_finish(child)
  Shennong:::.sn_usage_finish(parent)
  Shennong::sn_disable_usage_tracking()

  runs <- Shennong::sn_list_usage_runs(path)
  expect_identical(runs$workflow, "parent-only")
  expect_equal(runs$depth, 0L)
})

test_that("finish-write failures never replace a scientific error", {
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")
  on.exit(Shennong::sn_disable_usage_tracking(), add = TRUE)

  path <- tempfile(fileext = ".sqlite")
  Shennong::sn_enable_usage_tracking(
    path,
    mode = "test",
    display = FALSE,
    strict = TRUE
  )
  token <- Shennong:::.sn_usage_begin(
    "scientific-error",
    quote(scientific_error()),
    category = "test"
  )
  expect_error(
    Shennong::sn_disable_usage_tracking(),
    "still active"
  )
  scientific_error <- simpleError("original scientific failure")
  testthat::local_mocked_bindings(
    .sn_usage_finish_run = function(...) stop("database finish failure"),
    .package = "Shennong"
  )
  expect_warning(
    expect_no_error(Shennong:::.sn_usage_finish(
      token,
      status = "error",
      condition = scientific_error
    )),
    "database finish failure"
  )
})

test_that("timing helper preserves visibility, warnings, and errors", {
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")
  on.exit(Shennong::sn_disable_usage_tracking(), add = TRUE)

  expect_false(withVisible(Shennong::sn_time_call(
    invisible(3L),
    label = "visibility",
    display = FALSE,
    record = FALSE
  ))$visible)

  path <- tempfile(fileext = ".sqlite")
  Shennong::sn_enable_usage_tracking(path, mode = "test", display = FALSE)
  expect_warning(
    value <- Shennong::sn_time_call(
      {
        warning("visible warning")
        7L
      },
      label = "warning-test",
      display = FALSE
    ),
    "visible warning"
  )
  expect_identical(value, 7L)
  expect_error(
    Shennong::sn_time_call(
      stop("timed error"),
      label = "error-test",
      display = FALSE
    ),
    "timed error"
  )
  Shennong::sn_disable_usage_tracking()
  runs <- Shennong::sn_list_usage_runs(path)
  warning_row <- runs[runs$workflow == "timed:warning-test", , drop = FALSE]
  error_row <- runs[runs$workflow == "timed:error-test", , drop = FALSE]
  expect_equal(warning_row$warning_count, 1L)
  expect_equal(warning_row$status, "ok")
  expect_equal(error_row$status, "error")
})

test_that("scientific timing excludes usage database writes", {
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")
  on.exit(Shennong::sn_disable_usage_tracking(), add = TRUE)

  path <- tempfile(fileext = ".sqlite")
  Shennong::sn_enable_usage_tracking(
    path,
    mode = "test",
    display = FALSE,
    functions = "sn_calculate_cluster_entropy",
    strict = TRUE
  )
  captured_timing <- NULL
  testthat::local_mocked_bindings(
    .sn_usage_insert_run = function(token) {
      Sys.sleep(0.2)
      list(
        run_id = "mock-run",
        root_run_id = "mock-run",
        invocation_number = 1L,
        parameter_set_invocation = 1L
      )
    },
    .sn_usage_finish_run = function(token, status, condition, timing) {
      captured_timing <<- timing
      Sys.sleep(0.2)
      invisible(NULL)
    },
    .package = "Shennong"
  )
  wall <- system.time({
    token <- Shennong:::.sn_usage_begin(
      "timing-boundary",
      quote(timing_boundary()),
      category = "test"
    )
    Sys.sleep(0.02)
    Shennong:::.sn_usage_finish(token)
  })[["elapsed"]]

  expect_true(is.numeric(captured_timing))
  expect_gte(unname(wall - captured_timing[["elapsed"]]), 0.3)
})

test_that("usage registry covers every safe exported function family", {
  registry <- Shennong:::.sn_usage_registry()
  exports <- getNamespaceExports("Shennong")
  expected <- setdiff(exports, Shennong:::.sn_usage_control_functions)
  expect_setequal(
    registry$workflow,
    expected
  )
  expect_true(all(c(
    "sn_run_cluster", "sn_plot_dim", "sn_get_result", "sn_list_methods",
    "sn_store_result", "sn_read", "sn_write", "sn_install_dependencies"
  ) %in% registry$workflow))
  expect_false(any(Shennong:::.sn_usage_control_functions %in% registry$workflow))
  expect_identical(
    Shennong:::.sn_usage_registry(c("sn_plot_dim", "sn_get_result"))$workflow,
    c("sn_get_result", "sn_plot_dim")
  )
  expect_error(
    Shennong:::.sn_usage_registry("sn_enable_usage_tracking"),
    "cannot instrument"
  )
})

test_that("plot, get, list, and store APIs create usage rows", {
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")
  skip_if_not_installed("SeuratObject")
  on.exit(Shennong::sn_disable_usage_tracking(), add = TRUE)

  path <- tempfile(fileext = ".sqlite")
  Shennong::sn_enable_usage_tracking(
    path,
    mode = "test",
    display = FALSE,
    nested = TRUE,
    strict = TRUE
  )
  expect_s3_class(Shennong::sn_list_palettes(display = "none"), "data.frame")
  expect_length(Shennong::sn_get_palette("Paired", n = 3L), 3L)
  expect_s3_class(
    Shennong::sn_plot_boxplot(
      data.frame(group = c("a", "b"), value = c(1, 2)),
      x = "group",
      y = "value"
    ),
    "ggplot"
  )
  object <- SeuratObject::CreateSeuratObject(
    counts = Matrix::sparseMatrix(
      i = c(1L, 2L),
      j = c(1L, 2L),
      x = c(1, 1),
      dims = c(2L, 2L),
      dimnames = list(c("gene1", "gene2"), c("cell1", "cell2"))
    )
  )
  result <- Shennong:::.sn_new_analysis_result(
    analysis_type = "usage_test",
    name = "usage_test",
    method = "test",
    tables = list(primary = data.frame(value = 1))
  )
  object <- Shennong::sn_store_result(
    object,
    type = "usage_test",
    name = "usage_test",
    result = result
  )
  expect_s4_class(object, "Seurat")
  Shennong::sn_disable_usage_tracking()

  workflows <- Shennong::sn_list_usage_runs(path)$workflow
  expect_true(all(c(
    "sn_list_palettes", "sn_get_palette", "sn_plot_boxplot", "sn_store_result"
  ) %in% workflows))
})

test_that("remote DBI delivery is consent-gated, idempotent, and queryable", {
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")
  on.exit(Shennong::sn_disable_usage_tracking(), add = TRUE)

  outbox <- tempfile(fileext = ".sqlite")
  remote <- tempfile(fileext = ".sqlite")
  connections <- 0L
  factory <- function() {
    connections <<- connections + 1L
    DBI::dbConnect(RSQLite::SQLite(), remote)
  }
  store <- Shennong::sn_create_usage_store(
    backend = "dbi",
    path = outbox,
    connect = factory,
    table_prefix = "usage_test",
    allow_remote_ddl = TRUE
  )
  expect_error(
    Shennong::sn_enable_usage_tracking(store = store, mode = "test"),
    "requires consent"
  )
  expect_equal(connections, 0L)
  expect_false(file.exists(remote))

  consent <- Shennong::sn_confirm_usage_consent(
    scope = "remote_research",
    policy_id = "usage-study",
    policy_version = "1.0",
    purposes = c("usage_research", "performance_research"),
    data_categories = c("api_name", "timing", "safe_parameters"),
    study_id = "test-cohort"
  )
  Shennong::sn_enable_usage_tracking(
    store = store,
    consent = consent,
    mode = "production",
    display = FALSE,
    functions = "sn_calculate_cluster_entropy",
    remote_flush = "manual",
    strict = TRUE
  )
  Shennong::sn_calculate_cluster_entropy(
    usage_test_metadata(),
    cluster_by = "cluster",
    label_by = "batch"
  )
  Shennong::sn_disable_usage_tracking()
  expect_false(file.exists(remote))

  first <- Shennong::sn_flush_usage_tracking(
    store,
    consent,
    strict = TRUE
  )
  expect_equal(first$pending, 1L)
  expect_equal(first$appended, 1L)
  expect_equal(first$synced, 1L)
  expect_true(file.exists(remote))
  second <- Shennong::sn_flush_usage_tracking(store, consent, strict = TRUE)
  expect_equal(second$pending, 0L)

  connection <- DBI::dbConnect(RSQLite::SQLite(), remote)
  on.exit(DBI::dbDisconnect(connection), add = TRUE)
  expect_true(DBI::dbExistsTable(connection, "usage_test_sessions"))
  expect_true(DBI::dbExistsTable(connection, "usage_test_runs"))
  remote_runs <- DBI::dbReadTable(connection, "usage_test_runs")
  expect_equal(nrow(remote_runs), 1L)
  expect_identical(remote_runs$workflow, "sn_calculate_cluster_entropy")
  expect_false(any(c("pid", "error_message_redacted") %in% names(remote_runs)))
  expect_match(remote_runs$params_json, '"cluster_by":"cluster"')
  expect_false(grepl("cell1|patient|/private", paste(remote_runs, collapse = " ")))

  queried <- Shennong::sn_list_usage_runs(
    store = store,
    source = "remote",
    workflow = "sn_calculate_cluster_entropy"
  )
  expect_equal(nrow(queried), 1L)
  expect_equal(queried$invocation_number, 1L)
  expect_equal(queried$parameter_set_invocation, 1L)
})

test_that("minimal remote consent uploads only its authorized projection", {
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")
  on.exit(Shennong::sn_disable_usage_tracking(), add = TRUE)

  outbox <- tempfile(fileext = ".sqlite")
  remote <- tempfile(fileext = ".sqlite")
  store <- Shennong::sn_create_usage_store(
    backend = "dbi",
    path = outbox,
    connect = function() DBI::dbConnect(RSQLite::SQLite(), remote),
    table_prefix = "usage_minimal",
    allow_remote_ddl = TRUE
  )
  consent <- Shennong::sn_confirm_usage_consent(
    scope = "remote_research",
    policy_id = "minimal-study",
    policy_version = "1.0",
    data_categories = "api_name"
  )
  Shennong::sn_enable_usage_tracking(
    store = store,
    consent = consent,
    mode = "test",
    display = FALSE,
    functions = "sn_calculate_cluster_entropy",
    strict = TRUE
  )
  Shennong::sn_calculate_cluster_entropy(
    usage_test_metadata(),
    cluster_by = "cluster",
    label_by = "batch"
  )
  Shennong::sn_disable_usage_tracking()
  expect_equal(
    Shennong::sn_flush_usage_tracking(store, consent)$synced,
    1L
  )

  connection <- DBI::dbConnect(RSQLite::SQLite(), remote)
  on.exit(DBI::dbDisconnect(connection), add = TRUE)
  sessions <- DBI::dbReadTable(connection, "usage_minimal_sessions")
  runs <- DBI::dbReadTable(connection, "usage_minimal_runs")
  expect_identical(runs$workflow, "sn_calculate_cluster_entropy")
  expect_true(all(is.na(sessions[, c(
    "started_at", "shennong_version", "r_version"
  )])))
  expect_identical(sessions$packages_json, "{}")
  expect_true(all(is.na(runs[, c(
    "method", "backend", "status", "started_at", "finished_at",
    "elapsed_ms", "cpu_user_ms", "cpu_system_ms", "warning_count",
    "error_class"
  )])))
  expect_identical(runs$params_json, "{}")
  expect_identical(
    runs$params_sha256,
    digest::digest("{}", algo = "sha256", serialize = FALSE)
  )
  expect_identical(runs$acceleration_json, "{}")

  local <- DBI::dbConnect(RSQLite::SQLite(), outbox)
  on.exit(DBI::dbDisconnect(local), add = TRUE)
  DBI::dbExecute(local, "UPDATE workflow_runs SET remote_synced_at = NULL")
  retry <- Shennong::sn_flush_usage_tracking(store, consent)
  expect_equal(retry$appended, 0L)
  expect_equal(retry$already_present, 1L)
  expect_equal(retry$synced, 1L)
  expect_equal(nrow(DBI::dbReadTable(connection, "usage_minimal_runs")), 1L)

  expect_error(
    DBI::dbAppendTable(connection, "usage_minimal_runs", runs),
    "UNIQUE|unique"
  )
})

test_that("one outbox flushes each consent receipt without crossing categories", {
  skip_if_not_installed("DBI")
  skip_if_not_installed("RSQLite")
  on.exit(Shennong::sn_disable_usage_tracking(), add = TRUE)

  outbox <- tempfile(fileext = ".sqlite")
  remote <- tempfile(fileext = ".sqlite")
  store <- Shennong::sn_create_usage_store(
    backend = "dbi",
    path = outbox,
    connect = function() DBI::dbConnect(RSQLite::SQLite(), remote),
    table_prefix = "usage_mixed",
    allow_remote_ddl = TRUE
  )
  consent_parameters <- Shennong::sn_confirm_usage_consent(
    scope = "remote_research",
    policy_id = "mixed-study",
    policy_version = "parameters",
    data_categories = c("api_name", "safe_parameters")
  )
  consent_timing <- Shennong::sn_confirm_usage_consent(
    scope = "remote_research",
    policy_id = "mixed-study",
    policy_version = "timing",
    data_categories = c("api_name", "timing")
  )
  for (consent in list(consent_parameters, consent_timing)) {
    Shennong::sn_enable_usage_tracking(
      store = store,
      consent = consent,
      mode = "test",
      display = FALSE,
      functions = "sn_calculate_cluster_entropy",
      strict = TRUE
    )
    Shennong::sn_calculate_cluster_entropy(
      usage_test_metadata(),
      cluster_by = "cluster",
      label_by = "batch"
    )
    Shennong::sn_disable_usage_tracking()
  }

  first <- Shennong::sn_flush_usage_tracking(store, consent_parameters)
  expect_equal(first$pending, 1L)
  local_rows <- Shennong::sn_list_usage_runs(outbox)
  expect_equal(sum(!is.na(local_rows$remote_synced_at)), 1L)
  expect_equal(sum(is.na(local_rows$remote_synced_at)), 1L)

  second <- Shennong::sn_flush_usage_tracking(store, consent_timing)
  expect_equal(second$pending, 1L)
  expect_true(all(!is.na(Shennong::sn_list_usage_runs(outbox)$remote_synced_at)))

  connection <- DBI::dbConnect(RSQLite::SQLite(), remote)
  on.exit(DBI::dbDisconnect(connection), add = TRUE)
  sessions <- DBI::dbReadTable(connection, "usage_mixed_sessions")
  runs <- DBI::dbReadTable(connection, "usage_mixed_runs")
  parameter_session <- sessions$session_id[
    sessions$consent_id == consent_parameters$consent_id
  ]
  parameter_run <- runs[runs$session_id == parameter_session, , drop = FALSE]
  timing_run <- runs[runs$session_id != parameter_session, , drop = FALSE]
  expect_false(identical(parameter_run$params_json, "{}"))
  expect_true(is.na(parameter_run$elapsed_ms))
  expect_identical(timing_run$params_json, "{}")
  expect_false(is.na(timing_run$elapsed_ms))
})
