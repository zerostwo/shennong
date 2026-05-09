#!/usr/bin/env Rscript

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- sub("^--", "", args[[i]])
    value <- if (i < length(args)) args[[i + 1L]] else NA_character_
    out[[key]] <- value
    i <- i + 2L
  }
  out
}

parse_vector <- function(x, default, type = c("character", "integer")) {
  type <- match.arg(type)
  x <- x %||% default
  values <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  values <- values[nzchar(values)]
  if (identical(type, "integer")) {
    values <- as.integer(values)
  }
  values
}

parse_bool <- function(x, default = FALSE) {
  if (is.null(x)) {
    return(default)
  }
  tolower(x) %in% c("true", "t", "1", "yes", "y")
}

parse_elapsed <- function(x) {
  if (length(x) == 0L || is.na(x)) {
    return(NA_real_)
  }
  x <- sub(".*: ", "", x)
  parts <- strsplit(x, ":", fixed = TRUE)[[1]]
  parts <- as.numeric(parts)
  if (length(parts) == 3L) {
    return(parts[[1]] * 3600 + parts[[2]] * 60 + parts[[3]])
  }
  if (length(parts) == 2L) {
    return(parts[[1]] * 60 + parts[[2]])
  }
  parts[[1]]
}

parse_time_log <- function(path) {
  if (!file.exists(path)) {
    return(list(max_rss_kb = NA_real_, elapsed_sec_time = NA_real_))
  }
  lines <- readLines(path, warn = FALSE)
  rss <- grep("Maximum resident set size", lines, value = TRUE)
  elapsed <- grep("Elapsed \\(wall clock\\) time", lines, value = TRUE)
  list(
    max_rss_kb = if (length(rss)) as.numeric(sub(".*: ", "", rss[[1]])) else NA_real_,
    elapsed_sec_time = parse_elapsed(elapsed[[1]] %||% NA_character_)
  )
}

write_csv_row <- function(path, row) {
  row <- lapply(row, function(x) {
    if (is.null(x) || length(x) == 0L) {
      return(NA_character_)
    }
    if (length(x) > 1L) {
      return(paste(x, collapse = ";"))
    }
    x
  })
  df <- as.data.frame(row, stringsAsFactors = FALSE)
  if (!file.exists(path)) {
    utils::write.table(df, file = path, sep = ",", row.names = FALSE, col.names = TRUE, quote = TRUE)
  } else {
    utils::write.table(df, file = path, sep = ",", row.names = FALSE, col.names = FALSE, quote = TRUE, append = TRUE)
  }
}

prepare_input <- function(base, input_path, n_cells, n_batches, assay, layer, seed) {
  if (file.exists(input_path)) {
    return(invisible(input_path))
  }

  set.seed(seed + n_cells + n_batches)
  total_cells <- ncol(base)
  replace <- n_cells > total_cells
  idx <- sample(seq_len(total_cells), size = n_cells, replace = replace)
  cell_names <- sprintf("bench_cell_%07d", seq_len(n_cells))

  counts <- SeuratObject::LayerData(base, assay = assay, layer = "counts")[, idx, drop = FALSE]
  decontam <- SeuratObject::LayerData(base, assay = assay, layer = layer)[, idx, drop = FALSE]
  colnames(counts) <- cell_names
  colnames(decontam) <- cell_names

  metadata <- base[[]][idx, , drop = FALSE]
  rownames(metadata) <- cell_names
  metadata$sample <- rep(sprintf("bench_batch_%02d", seq_len(n_batches)), length.out = n_cells)
  metadata$sample <- sample(metadata$sample, size = n_cells, replace = FALSE)

  object <- Seurat::CreateSeuratObject(counts = counts, assay = assay, meta.data = metadata, project = "coralysis_capacity")
  SeuratObject::LayerData(object, assay = assay, layer = layer) <- decontam
  SeuratObject::DefaultLayer(object[[assay]]) <- layer

  dir.create(dirname(input_path), recursive = TRUE, showWarnings = FALSE)
  qs::qsave(object, input_path, preset = "fast")
  rm(counts, decontam, metadata, object)
  invisible(gc(verbose = FALSE))
  invisible(input_path)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

base_path <- args$base %||% "/home/sduan/projects/immune-atlas/data/processed/pbmc.qs"
outdir <- args$outdir %||% file.path(getwd(), "benchmarks", "coralysis_capacity", "results")
shennong_path <- args$shennong_path %||% getwd()
assay <- args$assay %||% "RNA"
layer <- args$layer %||% "decontaminated_counts"
cell_counts <- parse_vector(args$cell_counts, "3000,10000", type = "integer")
batch_counts <- parse_vector(args$batch_counts, "3,10", type = "integer")
methods <- parse_vector(args$methods, "coralysis,coralysis2", type = "character")
threads <- parse_vector(args$threads, "1", type = "integer")
nfeatures <- as.integer(args$nfeatures %||% "2000")
npcs <- as.integer(args$npcs %||% "10")
dims <- as.integer(args$dims %||% as.character(npcs))
L <- as.integer(args$L %||% "3")
k <- as.integer(args$k %||% "8")
seed <- as.integer(args$seed %||% "717")
store_sce <- args$store_sce %||% "false"
train_with_bnn <- args$train_with_bnn %||% "true"
build_train_set <- args$build_train_set %||% "true"
pca_chunk_size <- as.integer(args$pca_chunk_size %||% "20000")
predict_chunk_size <- as.integer(args$predict_chunk_size %||% "20000")
skip_existing <- parse_bool(args$skip_existing, default = TRUE)
max_runs <- as.integer(args$max_runs %||% "0")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
input_dir <- file.path(outdir, "inputs")
run_dir <- file.path(outdir, "runs")
dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
summary_csv <- file.path(outdir, "summary.csv")

suppressPackageStartupMessages({
  library(qs)
  library(jsonlite)
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
})

message("Loading base object: ", base_path)
base <- qs::qread(base_path)
message("Base dimensions: ", nrow(base), " features x ", ncol(base), " cells")

cmdline <- commandArgs(FALSE)
file_arg <- grep("^--file=", cmdline, value = TRUE)
script_dir <- if (length(file_arg) > 0L) {
  dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = TRUE))
} else {
  file.path(getwd(), "benchmarks", "coralysis_capacity")
}
bench_one <- file.path(script_dir, "bench_one.R")
if (!file.exists(bench_one)) {
  bench_one <- file.path(getwd(), "benchmarks", "coralysis_capacity", "bench_one.R")
}
if (!file.exists(bench_one)) {
  stop("Could not locate bench_one.R", call. = FALSE)
}

run_count <- 0L
for (n_cells in cell_counts) {
  for (n_batches in batch_counts) {
    input_path <- file.path(input_dir, sprintf("pbmc_n%06d_b%02d.qs", n_cells, n_batches))
    prepare_input(
      base = base,
      input_path = input_path,
      n_cells = n_cells,
      n_batches = n_batches,
      assay = assay,
      layer = layer,
      seed = seed
    )

    for (method in methods) {
      for (thread_count in threads) {
        run_id <- sprintf(
          "%s_n%06d_b%02d_t%02d_L%02d_k%02d",
          method, n_cells, n_batches, thread_count, L, k
        )
        run_prefix <- file.path(run_dir, run_id)
        result_json <- paste0(run_prefix, ".json")
        time_log <- paste0(run_prefix, ".time.txt")
        stdout_log <- paste0(run_prefix, ".stdout.txt")
        stderr_log <- paste0(run_prefix, ".stderr.txt")
        if (skip_existing && file.exists(result_json) && file.exists(time_log)) {
          message("Skipping existing run: ", run_id)
          next
        }

        run_count <- run_count + 1L
        if (max_runs > 0L && run_count > max_runs) {
          message("Reached --max-runs = ", max_runs)
          quit(status = 0L)
        }

        message("Running ", run_id)
        cmd_args <- c(
          "-v",
          "-o", time_log,
          file.path(R.home("bin"), "Rscript"),
          bench_one,
          "--input", input_path,
          "--output", result_json,
          "--method", method,
          "--batch", "sample",
          "--assay", assay,
          "--layer", layer,
          "--threads", as.character(thread_count),
          "--nfeatures", as.character(nfeatures),
          "--npcs", as.character(npcs),
          "--dims", as.character(dims),
          "--L", as.character(L),
          "--k", as.character(k),
          "--store_sce", store_sce,
          "--train_with_bnn", train_with_bnn,
          "--build_train_set", build_train_set,
          "--pca_chunk_size", as.character(pca_chunk_size),
          "--predict_chunk_size", as.character(predict_chunk_size),
          "--shennong_path", shennong_path
        )
        status <- system2("/usr/bin/time", args = cmd_args, stdout = stdout_log, stderr = stderr_log)
        if (is.null(status)) {
          status <- 0L
        }
        time_info <- parse_time_log(time_log)
        result <- if (file.exists(result_json)) {
          jsonlite::read_json(result_json, simplifyVector = TRUE)
        } else {
          list(status = "missing_result", error = NA_character_)
        }
        row <- c(
          list(
            run_id = run_id,
            exit_status = status,
            max_rss_kb = time_info$max_rss_kb,
            max_rss_gb = time_info$max_rss_kb / 1024^2,
            elapsed_sec_time = time_info$elapsed_sec_time
          ),
          result
        )
        write_csv_row(summary_csv, row)
      }
    }
  }
}

message("Benchmark summary: ", summary_csv)
