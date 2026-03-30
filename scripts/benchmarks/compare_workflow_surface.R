args <- commandArgs(trailingOnly = TRUE)

root_dir <- if (length(args) == 0L) {
  getwd()
} else {
  args[[1]]
}

benchmark_dir <- file.path(root_dir, "scripts", "benchmarks")
script_dir <- file.path(benchmark_dir, "scripts")
result_dir <- file.path(benchmark_dir, "results")

dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

workflow_files <- c(
  shennong = file.path(script_dir, "pbmc_shennong_workflow.R"),
  seurat_baseline = file.path(script_dir, "pbmc_seurat_baseline.R")
)

missing_files <- workflow_files[!file.exists(workflow_files)]
if (length(missing_files) > 0L) {
  stop(
    "Missing benchmark script(s):\n",
    paste(unname(missing_files), collapse = "\n"),
    call. = FALSE
  )
}

read_executable_lines <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- trimws(lines)
  lines <- lines[nzchar(lines)]
  lines[!grepl("^#", lines)]
}

extract_function_calls <- function(lines) {
  matches <- gregexpr("\\b[A-Za-z][A-Za-z0-9_.:]*\\s*\\(", lines, perl = TRUE)
  tokens <- regmatches(lines, matches)
  calls <- unlist(tokens, use.names = FALSE)
  calls <- trimws(sub("\\($", "", calls))
  calls[nzchar(calls)]
}

count_regex <- function(lines, pattern) {
  sum(grepl(pattern, lines, perl = TRUE))
}

measure_script <- function(name, path) {
  lines <- read_executable_lines(path)
  calls <- extract_function_calls(lines)

  data.frame(
    workflow = name,
    file = normalizePath(path, winslash = "/", mustWork = TRUE),
    executable_lines = length(lines),
    assignment_steps = count_regex(lines, "<-"),
    unique_function_calls = length(unique(calls)),
    shennong_calls = count_regex(lines, "\\bsn_[A-Za-z0-9_]+\\s*\\("),
    seurat_calls = count_regex(lines, "(\\bSeurat::|\\bNormalizeData\\s*\\(|\\bFindVariableFeatures\\s*\\(|\\bScaleData\\s*\\(|\\bRunPCA\\s*\\(|\\bFindNeighbors\\s*\\(|\\bFindClusters\\s*\\(|\\bRunUMAP\\s*\\(|\\bFindAllMarkers\\s*\\()"),
    persistence_calls = count_regex(lines, "(\\bsaveRDS\\s*\\(|\\bwrite\\.[A-Za-z0-9_]+\\s*\\(|\\bsn_store_[A-Za-z0-9_]+\\s*\\(|\\bsn_list_results\\s*\\(|\\bsn_get_[A-Za-z0-9_]+\\s*\\()"),
    conversion_calls = count_regex(lines, "(\\bas\\.Seurat\\s*\\(|\\bas\\.SingleCellExperiment\\s*\\(|\\bCreateSeuratObject\\s*\\(|\\bSingleCellExperiment::|\\banndataR::|\\bsn_add_data_from_anndata\\s*\\()"),
    stringsAsFactors = FALSE
  )
}

metrics <- do.call(
  rbind,
  Map(measure_script, names(workflow_files), unname(workflow_files))
)

metrics_csv <- file.path(result_dir, "workflow_surface_metrics.csv")
write.csv(metrics, metrics_csv, row.names = FALSE)

render_table <- function(df) {
  header <- paste(names(df), collapse = " | ")
  divider <- paste(rep("---", ncol(df)), collapse = " | ")
  rows <- apply(df, 1, function(row) paste(row, collapse = " | "))
  c(
    paste0("| ", header, " |"),
    paste0("| ", divider, " |"),
    paste0("| ", rows, " |")
  )
}

summary_lines <- c(
  "# Workflow Surface Metrics",
  "",
  paste("Generated:", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  "",
  "This phase-0 benchmark is a static script audit intended to quantify early",
  "workflow-surface differences between a Shennong-first and a Seurat-first",
  "implementation of a routine clustering plus marker-discovery workflow.",
  "",
  render_table(metrics)
)

summary_md <- file.path(result_dir, "workflow_surface_metrics.md")
writeLines(summary_lines, summary_md)

print(metrics)
message("Wrote: ", metrics_csv)
message("Wrote: ", summary_md)
