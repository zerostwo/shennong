#!/usr/bin/env Rscript

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

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L || (length(x) == 1L && is.na(x))) y else x
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x)) {
    return(default)
  }
  tolower(x) %in% c("true", "t", "1", "yes", "y")
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
input <- args$input
output <- args$output
method <- args$method %||% "coralysis2"
batch <- args$batch %||% "sample"
assay <- args$assay %||% "RNA"
layer <- args$layer %||% "decontaminated_counts"
threads <- as.integer(args$threads %||% "1")
nfeatures <- as.integer(args$nfeatures %||% "2000")
npcs <- as.integer(args$npcs %||% "10")
dims <- seq_len(as.integer(args$dims %||% as.character(npcs)))
L <- as.integer(args$L %||% "3")
k <- as.integer(args$k %||% "8")
resolution <- as.numeric(args$resolution %||% "0.2")
store_sce <- as_bool(args$store_sce, default = FALSE)
train_with_bnn <- as_bool(args$train_with_bnn, default = TRUE)
build_train_set <- as_bool(args$build_train_set, default = TRUE)
store_joint_probability <- as_bool(args$store_joint_probability, default = identical(method, "coralysis"))
shennong_path <- args$shennong_path %||% normalizePath(file.path(getwd()), mustWork = FALSE)
pca_chunk_size <- as.integer(args$pca_chunk_size %||% "20000")
predict_chunk_size <- as.integer(args$predict_chunk_size %||% "20000")

if (is.null(input) || is.null(output)) {
  stop("Usage: bench_one.R --input <input.qs> --output <result.json> [options]", call. = FALSE)
}

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(qs)
  library(jsonlite)
  library(SeuratObject)
  library(Matrix)
})

if (requireNamespace("devtools", quietly = TRUE) && file.exists(file.path(shennong_path, "DESCRIPTION"))) {
  suppressPackageStartupMessages(devtools::load_all(shennong_path, quiet = TRUE))
} else {
  suppressPackageStartupMessages(library(Shennong))
}

options(
  Coralysis.pca.chunk.size = pca_chunk_size,
  Coralysis.predict.chunk.size = predict_chunk_size
)

result <- list(
  input = input,
  method = method,
  batch = batch,
  assay = assay,
  layer = layer,
  threads = threads,
  nfeatures = nfeatures,
  npcs = npcs,
  dims = paste(dims, collapse = ","),
  L = L,
  k = k,
  resolution = resolution,
  store_sce = store_sce,
  train_with_bnn = train_with_bnn,
  build_train_set = build_train_set,
  store_joint_probability = store_joint_probability,
  pca_chunk_size = pca_chunk_size,
  predict_chunk_size = predict_chunk_size,
  status = "started",
  error = NA_character_
)

started <- Sys.time()
obj <- qs::qread(input)
result$cells <- ncol(obj)
result$features <- nrow(obj)
result$batches <- length(unique(as.character(obj[[batch, drop = TRUE]])))
result$input_object_gb <- as.numeric(object.size(obj)) / 1024^3

icp_args <- list(
  k = k,
  L = L,
  threads = threads,
  RNGseed = 717,
  train.with.bnn = train_with_bnn,
  build.train.set = build_train_set
)
if (identical(method, "coralysis2")) {
  icp_args$store.joint.probability <- store_joint_probability
}

pca_args <- list(
  p = npcs,
  pca.method = "stats",
  return.model = TRUE
)
if (identical(method, "coralysis2")) {
  pca_args$pca.chunk.size <- pca_chunk_size
}

control <- list(
  store_sce = store_sce,
  icp_args = icp_args,
  pca_args = pca_args
)

run_result <- tryCatch(
  {
    out <- sn_run_cluster(
      object = obj,
      integration_method = method,
      batch_by = batch,
      assay = assay,
      layer = layer,
      normalization_method = "seurat",
      nfeatures = nfeatures,
      npcs = npcs,
      dims = dims,
      resolution = resolution,
      block_genes = NULL,
      reuse = FALSE,
      integration_control = control,
      verbose = FALSE
    )
    list(status = "ok", object = out, error = NA_character_)
  },
  error = function(e) {
    list(status = "error", object = NULL, error = conditionMessage(e))
  }
)

finished <- Sys.time()
result$status <- run_result$status
result$error <- run_result$error
result$elapsed_sec_r <- as.numeric(difftime(finished, started, units = "secs"))
if (!is.null(run_result$object)) {
  result$output_object_gb <- as.numeric(object.size(run_result$object)) / 1024^3
  result$reductions <- paste(names(run_result$object@reductions), collapse = ",")
  result$integration_reduction <- run_result$object@misc$integration$reduction %||% NA_character_
  result$integration_features <- length(run_result$object@misc$integration$input_features %||% character(0))
} else {
  result$output_object_gb <- NA_real_
  result$reductions <- NA_character_
  result$integration_reduction <- NA_character_
  result$integration_features <- NA_integer_
}

jsonlite::write_json(result, output, auto_unbox = TRUE, pretty = TRUE, null = "null")
if (!identical(run_result$status, "ok")) {
  quit(status = 2L)
}
