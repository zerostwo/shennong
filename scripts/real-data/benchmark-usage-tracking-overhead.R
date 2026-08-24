#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", all_args, value = TRUE)
script_file <- if (length(script_arg)) {
  sub("^--file=", "", script_arg[[1L]])
} else {
  file.path(getwd(), "scripts", "real-data", "benchmark-usage-tracking-overhead.R")
}
repo_root <- normalizePath(
  file.path(dirname(normalizePath(script_file, mustWork = TRUE)), "..", ".."),
  winslash = "/",
  mustWork = TRUE
)

.option <- function(name, default) {
  prefix <- paste0("--", name, "=")
  inline <- args[startsWith(args, prefix)]
  if (length(inline)) return(sub(prefix, "", inline[[length(inline)]], fixed = TRUE))
  index <- match(paste0("--", name), args)
  if (!is.na(index) && index < length(args)) return(args[[index + 1L]])
  default
}

.absolute_path <- function(path) {
  path <- path.expand(path)
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", path)) path <- file.path(repo_root, path)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

if ("--help" %in% args) {
  cat(
    "Usage: Rscript scripts/real-data/benchmark-usage-tracking-overhead.R [options]\n",
    "\n",
    "Options:\n",
    "  --iterations N  Calls per condition (default: 100)\n",
    "  --output PATH   JSON output path\n",
    "  --help          Show this message\n",
    sep = ""
  )
  quit(status = 0L)
}

iterations <- suppressWarnings(as.integer(.option("iterations", "100")))
if (length(iterations) != 1L || is.na(iterations) || iterations < 10L) {
  stop("`--iterations` must be one integer >= 10.", call. = FALSE)
}
output <- .absolute_path(.option(
  "output",
  file.path("data-local", "autozyme-benchmark", "usage-tracking-overhead.json")
))
required <- c("Shennong", "DBI", "RSQLite", "jsonlite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing benchmark package(s): ", paste(missing, collapse = ", "), ".", call. = FALSE)
}

suppressPackageStartupMessages(library(Shennong))
on.exit(Shennong::sn_disable_usage_tracking(), add = TRUE)
metadata <- data.frame(
  cluster = rep(sprintf("cluster_%02d", 1:10), each = 100),
  batch = rep(sprintf("batch_%02d", 1:5), length.out = 1000),
  stringsAsFactors = FALSE
)
call_metric <- function() {
  invisible(Shennong::sn_calculate_cluster_entropy(
    metadata,
    cluster_by = "cluster",
    label_by = "batch"
  ))
}
for (index in seq_len(10L)) call_metric()

set.seed(717)
rng_before <- .Random.seed
baseline <- system.time(for (index in seq_len(iterations)) call_metric())[["elapsed"]]

database <- sub("[.]json$", ".sqlite", output)
if (file.exists(database)) {
  stop("Refusing to overwrite existing usage database: ", database, call. = FALSE)
}
Shennong::sn_enable_usage_tracking(
  database,
  mode = "benchmark",
  display = FALSE,
  nested = TRUE,
  strict = TRUE
)
tracked <- system.time(for (index in seq_len(iterations)) call_metric())[["elapsed"]]
Shennong::sn_disable_usage_tracking()
runs <- Shennong::sn_list_usage_runs(
  database,
  workflow = "sn_calculate_cluster_entropy",
  top_level = TRUE,
  limit = iterations + 10L
)

report <- list(
  schema_version = "shennong.usage-overhead/v1",
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  iterations = iterations,
  workload = list(
    workflow = "sn_calculate_cluster_entropy",
    rows = nrow(metadata),
    clusters = length(unique(metadata$cluster)),
    batches = length(unique(metadata$batch))
  ),
  results = list(
    baseline_total_seconds = unname(baseline),
    tracked_total_seconds = unname(tracked),
    baseline_wall_ms_per_call = unname(1000 * baseline / iterations),
    tracked_wall_ms_per_call = unname(1000 * tracked / iterations),
    incremental_wall_ms_per_call = unname(1000 * (tracked - baseline) / iterations),
    stored_median_analysis_ms = 1000 * stats::median(runs$elapsed_seconds),
    stored_p95_analysis_ms = 1000 * unname(stats::quantile(
      runs$elapsed_seconds,
      0.95,
      names = FALSE,
      type = 8
    )),
    recorded_rows = nrow(runs),
    invocation_numbers_exact = identical(sort(runs$invocation_number), seq_len(iterations)),
    rng_unchanged = identical(rng_before, .Random.seed)
  ),
  environment = list(
    r = R.version.string,
    shennong = as.character(utils::packageVersion("Shennong")),
    DBI = as.character(utils::packageVersion("DBI")),
    RSQLite = as.character(utils::packageVersion("RSQLite"))
  ),
  database = normalizePath(database, winslash = "/", mustWork = TRUE)
)
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(
  report,
  output,
  auto_unbox = TRUE,
  pretty = TRUE,
  null = "null",
  digits = 10
)
print(as.data.frame(report$results), row.names = FALSE)
cat("Wrote ", normalizePath(output, winslash = "/", mustWork = TRUE), "\n", sep = "")
