#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", all_args, value = TRUE)
script_file <- if (length(script_arg)) {
  sub("^--file=", "", script_arg[[1]])
} else {
  file.path(getwd(), "scripts", "real-data", "benchmark-single-cell-autozyme.R")
}
script_file <- normalizePath(script_file, winslash = "/", mustWork = TRUE)
repo_root <- normalizePath(
  file.path(dirname(script_file), "..", ".."),
  winslash = "/",
  mustWork = TRUE
)

.option <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  inline <- args[startsWith(args, prefix)]
  if (length(inline)) {
    return(sub(prefix, "", inline[[length(inline)]], fixed = TRUE))
  }
  index <- match(paste0("--", name), args)
  if (!is.na(index) && index < length(args)) return(args[[index + 1L]])
  default
}

.flag <- function(name) paste0("--", name) %in% args

.absolute_path <- function(path) {
  path <- path.expand(path)
  if (!grepl("^/", path)) path <- file.path(repo_root, path)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

.integer_option <- function(name, default, minimum = 1L) {
  value <- suppressWarnings(as.integer(.option(name, as.character(default))))
  if (length(value) != 1L || is.na(value) || value < minimum) {
    stop("`--", name, "` must be one integer >= ", minimum, ".", call. = FALSE)
  }
  value
}

if (.flag("help")) {
  cat(
    "Usage: Rscript scripts/real-data/benchmark-single-cell-autozyme.R [options]\n",
    "\n",
    "Options:\n",
    "  --root PATH          Validated real-data root (default: environment or\n",
    "                       data-local/pkgdown-real)\n",
    "  --output PATH        JSON report path\n",
    "  --repetitions N      Fresh-process repetitions per condition (default: 3)\n",
    "  --seed N             Deterministic seed (default: 717)\n",
    "  --help               Show this message\n",
    sep = ""
  )
  quit(status = 0L)
}

root <- .absolute_path(.option(
  "root",
  Sys.getenv(
    "SHENNONG_REAL_DATA_DIR",
    unset = file.path(repo_root, "data-local", "pkgdown-real")
  )
))
output <- .absolute_path(.option(
  "output",
  file.path(root, "results", "single-cell-autozyme-benchmark.json")
))
repetitions <- .integer_option("repetitions", 3L)
seed <- .integer_option("seed", 717L, minimum = 0L)
worker <- .flag("worker")

required <- c(
  "Shennong", "autozyme", "digest", "jsonlite", "Matrix", "qs2",
  "scDblFinder", "Seurat", "SeuratObject", "SingleCellExperiment"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing benchmark package(s): ", paste(missing, collapse = ", "), ".", call. = FALSE)
}

.pbmc_path <- function() file.path(root, "single-cell", "kotliarov_pbmc.qs2")

.load_pbmc <- function() {
  object <- qs2::qs_read(.pbmc_path())
  if (!inherits(object, "Seurat") || ncol(object) != 2000L) {
    stop("Expected the validated 2,000-cell Kotliarov PBMC Seurat object.", call. = FALSE)
  }
  object
}

.prepare_sample_objects <- function(object) {
  sample_ids <- sort(unique(as.character(object$real_sample)))
  sample_ids <- sample_ids[!is.na(sample_ids) & nzchar(sample_ids)]
  objects <- lapply(sample_ids, function(sample_id) {
    cells <- colnames(object)[as.character(object$real_sample) == sample_id]
    counts <- SeuratObject::LayerData(
      object = object,
      assay = "RNA",
      layer = "counts"
    )[, cells, drop = FALSE]
    result <- SeuratObject::CreateSeuratObject(
      counts = counts,
      assay = "RNA",
      project = sample_id
    )
    result$real_sample <- sample_id
    result
  })
  names(objects) <- sample_ids
  objects
}

.assay_snapshot <- function(assay) {
  layers <- SeuratObject::Layers(assay)
  payload <- lapply(layers, function(layer) {
    matrix <- SeuratObject::LayerData(assay, layer = layer)
    list(
      class = class(matrix),
      dim = dim(matrix),
      dimnames = dimnames(matrix),
      value = digest::digest(matrix, algo = "sha256", serialize = TRUE)
    )
  })
  names(payload) <- layers
  list(
    class = class(assay),
    dim = dim(assay),
    dimnames = dimnames(assay),
    layers = payload,
    cells = digest::digest(methods::slot(assay, "cells"), algo = "sha256"),
    features = digest::digest(methods::slot(assay, "features"), algo = "sha256"),
    valid = isTRUE(methods::validObject(assay, test = TRUE))
  )
}

.result_snapshot <- function(operation, value) {
  if (identical(operation, "merge")) {
    list(
      dim = dim(value),
      cells = colnames(value),
      features = rownames(value),
      metadata = digest::digest(value[[]], algo = "sha256"),
      assay = .assay_snapshot(value[["RNA"]])
    )
  } else if (identical(operation, "joinlayers")) {
    .assay_snapshot(value)
  } else {
    data.frame(
      cell = colnames(value),
      class = as.character(value$scDblFinder.class),
      score = as.numeric(value$scDblFinder.score),
      stringsAsFactors = FALSE
    )
  }
}

.patch_for <- function(operation) {
  switch(
    operation,
    merge = "seurat_merge",
    joinlayers = "seurat_joinlayers",
    scdblfinder = "scdblfinder",
    stop("Unknown benchmark operation `", operation, "`.", call. = FALSE)
  )
}

.run_worker <- function() {
  operation <- .option("operation")
  if (!operation %in% c("merge", "joinlayers", "scdblfinder")) {
    stop("Worker requires a valid `--operation`.", call. = FALSE)
  }
  accelerated <- identical(tolower(.option("accelerated", "false")), "true")
  result_path <- .absolute_path(.option("result"))
  patch <- .patch_for(operation)

  suppressPackageStartupMessages(library(Shennong))
  old_options <- options(shennong.autozyme = FALSE)
  on.exit(options(old_options), add = TRUE)
  try(sn_disable_autozyme(patch), silent = TRUE)
  on.exit(try(sn_disable_autozyme(patch), silent = TRUE), add = TRUE)
  if (accelerated) sn_enable_autozyme(patch)

  set.seed(seed)
  pbmc <- .load_pbmc()
  sample_objects <- NULL
  if (operation %in% c("merge", "joinlayers")) {
    sample_objects <- .prepare_sample_objects(pbmc)
  }

  expression <- switch(
    operation,
    merge = quote(base::merge(
      x = sample_objects[[1L]],
      y = sample_objects[-1L],
      add.cell.ids = names(sample_objects),
      merge.data = FALSE
    )),
    joinlayers = quote({
      counts <- lapply(sample_objects, function(object) {
        SeuratObject::LayerData(object, assay = "RNA", layer = "counts")
      })
      assay <- SeuratObject::CreateAssay5Object(counts = counts)
      SeuratObject::JoinLayers(assay, layers = "counts", new = "counts")
    }),
    scdblfinder = quote(sn_find_doublets(
      pbmc,
      clusters = NULL,
      group_by = NULL,
      dbr_sd = NULL,
      ncores = 1L,
      assay = "RNA",
      layer = "counts",
      min_features = 200L
    ))
  )

  gc()
  value <- NULL
  timing <- system.time({
    value <- if (accelerated) {
      eval(expression)
    } else {
      autozyme::with_disabled(eval(expression))
    }
  })
  snapshot <- .result_snapshot(operation, value)
  fingerprint <- digest::digest(snapshot, algo = "sha256", serialize = TRUE)
  cells <- ncol(value)
  features <- nrow(value)
  sn_disable_autozyme(patch)
  active_after <- sn_check_autozyme(patch)$active

  result <- list(
    operation = operation,
    accelerated = accelerated,
    patch = patch,
    elapsed_seconds = unname(timing[["elapsed"]]),
    output_fingerprint = fingerprint,
    output_cells = cells,
    output_features = features,
    valid = if (operation == "scdblfinder") TRUE else snapshot$valid %||% snapshot$assay$valid,
    active_after = isTRUE(active_after)
  )
  dir.create(dirname(result_path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(result, result_path, auto_unbox = TRUE, pretty = TRUE)
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

if (worker) {
  .run_worker()
  quit(status = 0L)
}

validator <- file.path(repo_root, "scripts", "real-data", "validate-real-data.R")
validation_status <- system2(
  file.path(R.home("bin"), "Rscript"),
  c(shQuote(validator), "--bundle", "kotliarov_pbmc", "--root", shQuote(root))
)
if (!identical(as.integer(validation_status), 0L)) {
  stop("Kotliarov PBMC validation failed; benchmark was not run.", call. = FALSE)
}

time_bin <- Sys.which("time")
if (!nzchar(time_bin) || !identical(normalizePath(time_bin), "/usr/bin/time")) {
  time_bin <- "/usr/bin/time"
}
if (!file.exists(time_bin)) {
  stop("GNU `/usr/bin/time` is required to record process peak RSS.", call. = FALSE)
}

.parse_peak_rss <- function(path) {
  lines <- readLines(path, warn = FALSE)
  matched <- grep("Maximum resident set size \\(kbytes\\):", lines, value = TRUE)
  if (length(matched) != 1L) return(NA_real_)
  as.numeric(trimws(sub("^.*:", "", matched[[1L]]))) / 1024
}

.run_child <- function(operation, accelerated, repetition) {
  run_dir <- tempfile(paste0("shennong-autozyme-", operation, "-"))
  dir.create(run_dir, recursive = TRUE)
  result_path <- file.path(run_dir, "result.json")
  time_path <- file.path(run_dir, "time.txt")
  log_path <- file.path(run_dir, "worker.log")
  worker_args <- c(
    "-v", "-o", shQuote(time_path),
    shQuote(file.path(R.home("bin"), "Rscript")),
    shQuote(script_file),
    "--worker",
    "--root", shQuote(root),
    "--operation", operation,
    "--accelerated", if (accelerated) "true" else "false",
    "--result", shQuote(result_path),
    "--seed", as.character(seed)
  )
  status <- system2(
    time_bin,
    args = worker_args,
    stdout = log_path,
    stderr = log_path,
    wait = TRUE
  )
  if (!identical(as.integer(status), 0L) || !file.exists(result_path)) {
    log <- if (file.exists(log_path)) readLines(log_path, warn = FALSE) else character()
    stop(
      "Benchmark worker failed for ", operation, " (accelerated=", accelerated,
      ", repetition=", repetition, "):\n", paste(log, collapse = "\n"),
      call. = FALSE
    )
  }
  result <- jsonlite::read_json(result_path, simplifyVector = TRUE)
  data.frame(
    operation = operation,
    patch = result$patch,
    condition = if (accelerated) "accelerated" else "baseline",
    repetition = repetition,
    elapsed_seconds = as.numeric(result$elapsed_seconds),
    peak_rss_mb = .parse_peak_rss(time_path),
    output_fingerprint = result$output_fingerprint,
    output_cells = as.integer(result$output_cells),
    output_features = as.integer(result$output_features),
    valid = isTRUE(result$valid),
    active_after = isTRUE(result$active_after),
    stringsAsFactors = FALSE
  )
}

operations <- c("merge", "joinlayers", "scdblfinder")
runs <- list()
run_index <- 0L
for (repetition in seq_len(repetitions)) {
  conditions <- if (repetition %% 2L) c(FALSE, TRUE) else c(TRUE, FALSE)
  for (operation in operations) {
    for (accelerated in conditions) {
      run_index <- run_index + 1L
      message(
        "Running ", operation, " ", if (accelerated) "accelerated" else "baseline",
        " repetition ", repetition, "/", repetitions, "."
      )
      runs[[run_index]] <- .run_child(operation, accelerated, repetition)
    }
  }
}
runs <- do.call(rbind, runs)

if (any(!runs$valid) || any(runs$active_after) || any(!is.finite(runs$peak_rss_mb))) {
  stop("A benchmark worker failed validity, rollback, or RSS reporting.", call. = FALSE)
}

summaries <- lapply(operations, function(operation) {
  selected <- runs[runs$operation == operation, , drop = FALSE]
  if (length(unique(selected$output_fingerprint)) != 1L) {
    stop(operation, " baseline and accelerated outputs were not identical.", call. = FALSE)
  }
  baseline <- selected[selected$condition == "baseline", , drop = FALSE]
  accelerated <- selected[selected$condition == "accelerated", , drop = FALSE]
  baseline_seconds <- stats::median(baseline$elapsed_seconds)
  accelerated_seconds <- stats::median(accelerated$elapsed_seconds)
  baseline_rss <- stats::median(baseline$peak_rss_mb)
  accelerated_rss <- stats::median(accelerated$peak_rss_mb)
  data.frame(
    operation = operation,
    patch = unique(selected$patch),
    cells = unique(selected$output_cells),
    features = unique(selected$output_features),
    repetitions = repetitions,
    baseline_seconds = baseline_seconds,
    accelerated_seconds = accelerated_seconds,
    speedup = baseline_seconds / accelerated_seconds,
    baseline_peak_rss_mb = baseline_rss,
    accelerated_peak_rss_mb = accelerated_rss,
    peak_rss_reduction_percent = 100 * (baseline_rss - accelerated_rss) / baseline_rss,
    output_identical = TRUE,
    patch_inactive_after = TRUE,
    stringsAsFactors = FALSE
  )
})
summaries <- do.call(rbind, summaries)

manifest_path <- file.path(root, "manifests", "kotliarov_pbmc.json")
manifest <- jsonlite::read_json(manifest_path, simplifyVector = TRUE)
status <- Shennong::sn_check_autozyme(vapply(
  operations,
  .patch_for,
  character(1)
))
cpu_model <- if (file.exists("/proc/cpuinfo")) {
  lines <- readLines("/proc/cpuinfo", warn = FALSE)
  model <- sub("^model name[[:space:]]*:[[:space:]]*", "", grep("^model name", lines, value = TRUE))
  if (length(model)) model[[1L]] else NA_character_
} else {
  NA_character_
}

report <- list(
  schema_version = "1.0.0",
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  seed = seed,
  methodology = list(
    fresh_process_per_run = TRUE,
    repetitions_per_condition = repetitions,
    ordering = "alternating baseline and accelerated",
    timing = "R system.time elapsed around the analytical operation",
    memory = "GNU time maximum resident set size for the complete worker process",
    summary = "median"
  ),
  dataset = list(
    name = "Kotliarov PBMC RNA/ADT",
    source = manifest$sources,
    artifact_sha256 = manifest$artifact$sha256,
    cells = manifest$summary$cells,
    features = manifest$summary$features,
    samples = manifest$summary$samples
  ),
  environment = list(
    os = paste(Sys.info()[c("sysname", "release", "machine")], collapse = " "),
    cpu = cpu_model,
    logical_cores = parallel::detectCores(logical = TRUE),
    r = R.version.string,
    shennong = as.character(utils::packageVersion("Shennong")),
    autozyme = as.character(utils::packageVersion("autozyme")),
    seurat = as.character(utils::packageVersion("Seurat")),
    seurat_object = as.character(utils::packageVersion("SeuratObject")),
    scdblfinder = as.character(utils::packageVersion("scDblFinder"))
  ),
  patch_status = status,
  summary = summaries,
  runs = runs
)

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(
  report,
  output,
  dataframe = "rows",
  auto_unbox = TRUE,
  pretty = TRUE,
  null = "null",
  na = "null",
  digits = 8
)
utils::write.csv(
  summaries,
  sub("[.]json$", ".csv", output),
  row.names = FALSE,
  na = ""
)

print(summaries, row.names = FALSE)
cat("Wrote ", normalizePath(output, winslash = "/", mustWork = TRUE), "\n", sep = "")
