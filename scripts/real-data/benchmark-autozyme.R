#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", all_args, value = TRUE)
script_file <- if (length(script_arg)) {
  sub("^--file=", "", script_arg[[1]])
} else {
  file.path(getwd(), "scripts", "real-data", "benchmark-autozyme.R")
}
repo_root <- normalizePath(
  file.path(dirname(script_file), "..", ".."),
  winslash = "/",
  mustWork = TRUE
)

.option <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  inline <- args[startsWith(args, prefix)]
  if (length(inline)) return(sub(prefix, "", inline[[length(inline)]], fixed = TRUE))
  index <- match(paste0("--", name), args)
  if (!is.na(index) && index < length(args)) return(args[[index + 1L]])
  default
}

if ("--help" %in% args) {
  cat(
    "Usage: Rscript scripts/real-data/benchmark-autozyme.R [options]\n",
    "\n",
    "Options:\n",
    "  --root PATH    Validated local real-data root (default: environment or\n",
    "                 data-local/pkgdown-real)\n",
    "  --output PATH  Local JSON result (default: ROOT/results/autozyme-benchmark.json)\n",
    "  --seed INTEGER Deterministic seed (default: 717)\n",
    "  --help         Show this message\n",
    sep = ""
  )
  quit(status = 0L)
}

.absolute_path <- function(path) {
  path <- path.expand(path)
  if (!grepl("^/", path)) path <- file.path(repo_root, path)
  normalizePath(path, winslash = "/", mustWork = FALSE)
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
  file.path(root, "results", "autozyme-benchmark.json")
))
seed <- suppressWarnings(as.integer(.option("seed", "717")))
if (length(seed) != 1L || is.na(seed)) stop("`--seed` must be one integer.", call. = FALSE)

required <- c(
  "Shennong", "autozyme", "CellChat", "Seurat", "SeuratObject", "WGCNA",
  "qs2", "jsonlite"
)
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing benchmark package(s): ", paste(missing, collapse = ", "), ".", call. = FALSE)
}

validator <- file.path(repo_root, "scripts", "real-data", "validate-real-data.R")
status <- system2(
  file.path(R.home("bin"), "Rscript"),
  c(shQuote(validator), "--bundle", "all", "--root", shQuote(root))
)
if (!identical(as.integer(status), 0L)) {
  stop("Real-data validation failed; AutoZyme was not benchmarked.", call. = FALSE)
}

suppressPackageStartupMessages({
  library(Shennong)
  library(Seurat)
  library(CellChat)
})

on.exit({
  try(sn_disable_autozyme(c("cellchat", "wgcna")), silent = TRUE)
}, add = TRUE)

.elapsed <- function(expression) {
  value <- NULL
  timing <- system.time(value <- force(expression))
  list(value = value, seconds = unname(timing[["elapsed"]]))
}

.cellchat_run <- function(object, database, accelerated) {
  if (accelerated) {
    sn_enable_autozyme("cellchat")
  } else {
    sn_disable_autozyme("cellchat")
  }
  on.exit(sn_disable_autozyme("cellchat"), add = TRUE)
  gc()
  .elapsed(sn_run_cell_communication(
    object,
    method = "cellchat",
    group_by = "seurat_clusters",
    assay = "RNA",
    layer = "data",
    species = "human",
    cellchat_db = database,
    min_cells = 10,
    raw_use = TRUE,
    store_name = if (accelerated) "accelerated" else "baseline",
    return_object = FALSE,
    nboot = 10,
    seed.use = seed
  ))
}

set.seed(seed)
pbmc <- qs2::qs_read(file.path(root, "single-cell", "kotliarov_pbmc.qs2"))
selected_cells <- unlist(lapply(
  split(colnames(pbmc), pbmc$real_response),
  function(cells) sample(cells, min(150L, length(cells)))
), use.names = FALSE)
pbmc <- subset(pbmc, cells = selected_cells)
pbmc <- sn_run_cluster(
  pbmc,
  normalization_method = "seurat",
  nfeatures = 1000L,
  npcs = 15L,
  dims = 1:10,
  resolution = 0.3,
  species = "human",
  verbose = FALSE
)
data(CellChatDB.human, package = "CellChat")
cellchat_db <- CellChat::subsetDB(
  CellChatDB.human,
  search = "Secreted Signaling",
  key = "annotation"
)
cellchat_baseline <- .cellchat_run(pbmc, cellchat_db, accelerated = FALSE)
cellchat_accelerated <- .cellchat_run(pbmc, cellchat_db, accelerated = TRUE)

cellchat_columns <- c("source", "target", "ligand", "receptor", "score")
.ordered_cellchat <- function(result) {
  table <- result[["tables"]][["primary"]][, cellchat_columns, drop = FALSE]
  table[do.call(order, table[seq_len(4L)]), , drop = FALSE]
}
cellchat_a <- .ordered_cellchat(cellchat_baseline$value)
cellchat_b <- .ordered_cellchat(cellchat_accelerated$value)
cellchat_keys_equal <- identical(cellchat_a[seq_len(4L)], cellchat_b[seq_len(4L)])
cellchat_max_diff <- if (nrow(cellchat_a)) {
  max(abs(cellchat_a$score - cellchat_b$score), na.rm = TRUE)
} else {
  0
}
if (!cellchat_keys_equal || !is.finite(cellchat_max_diff) || cellchat_max_diff > 1e-12) {
  stop("CellChat baseline and accelerated results are not exactly equivalent.", call. = FALSE)
}

.wgcna_run <- function(expression, metadata, accelerated) {
  if (accelerated) {
    sn_enable_autozyme("wgcna")
  } else {
    sn_disable_autozyme("wgcna")
  }
  on.exit(sn_disable_autozyme("wgcna"), add = TRUE)
  gc()
  .elapsed(sn_run_wgcna(
    expression,
    metadata = metadata,
    traits = c("event", "sample_type"),
    power = 6,
    min_module_size = 20L,
    merge_cut_height = 0.25,
    backend_control = list(
      blockwise = list(maxBlockSize = 400L, randomSeed = seed, nThreads = 1L)
    )
  ))
}

tcga <- qs2::qs_read(file.path(root, "bulk", "tcga_skcm.qs2"))
bulk_expression <- tcga$log2_uq_rsem_tpm
feature_variance <- apply(bulk_expression, 1L, stats::var)
selected_features <- names(sort(feature_variance, decreasing = TRUE))[seq_len(400L)]
bulk_expression <- bulk_expression[selected_features, , drop = FALSE]
bulk_metadata <- tcga$sample_data[colnames(bulk_expression), , drop = FALSE]
wgcna_baseline <- .wgcna_run(bulk_expression, bulk_metadata, accelerated = FALSE)
wgcna_accelerated <- .wgcna_run(bulk_expression, bulk_metadata, accelerated = TRUE)

.ordered_modules <- function(result) {
  table <- result[["tables"]][["modules"]]
  table[order(table$gene), , drop = FALSE]
}
wgcna_modules_equal <- identical(
  .ordered_modules(wgcna_baseline$value),
  .ordered_modules(wgcna_accelerated$value)
)

eigengene_a <- wgcna_baseline$value[["tables"]][["eigengenes"]]
eigengene_b <- wgcna_accelerated$value[["tables"]][["eigengenes"]]
eigengene_names <- intersect(
  setdiff(names(eigengene_a), "sample"),
  setdiff(names(eigengene_b), "sample")
)
eigengene_a <- as.matrix(eigengene_a[, eigengene_names, drop = FALSE])
eigengene_b <- as.matrix(eigengene_b[, eigengene_names, drop = FALSE])
for (index in seq_len(ncol(eigengene_a))) {
  if (stats::cor(eigengene_a[, index], eigengene_b[, index], use = "complete.obs") < 0) {
    eigengene_b[, index] <- -eigengene_b[, index]
  }
}
wgcna_eigengene_diff <- max(abs(eigengene_a - eigengene_b), na.rm = TRUE)

.ordered_traits <- function(result) {
  table <- result[["tables"]][["trait_associations"]]
  table[order(table$module, table$trait), , drop = FALSE]
}
trait_a <- .ordered_traits(wgcna_baseline$value)
trait_b <- .ordered_traits(wgcna_accelerated$value)
wgcna_trait_diff <- max(abs(trait_a$correlation - trait_b$correlation), na.rm = TRUE)
if (!wgcna_modules_equal || wgcna_eigengene_diff > 1e-5 || wgcna_trait_diff > 1e-5) {
  stop("WGCNA baseline and accelerated results exceed the declared tolerance.", call. = FALSE)
}

active_after <- sn_check_autozyme(c("cellchat", "wgcna"))$active
if (any(active_after)) stop("AutoZyme patches remained active after scoped validation.", call. = FALSE)

autozyme_status <- sn_check_autozyme(c("cellchat", "wgcna"))
report <- list(
  schema_version = "1.0.0",
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  seed = seed,
  inputs = list(
    kotliarov_cells = ncol(pbmc),
    kotliarov_clusters = length(unique(pbmc$seurat_clusters)),
    tcga_samples = ncol(bulk_expression),
    tcga_features = nrow(bulk_expression)
  ),
  autozyme = list(
    version = unique(autozyme_status$autozyme_version),
    remote_sha = unique(autozyme_status$autozyme_remote_sha),
    source_match = all(autozyme_status$autozyme_source_match),
    upstream_versions_match = all(autozyme_status$version_match),
    active_after = any(active_after)
  ),
  cellchat = list(
    baseline_seconds = cellchat_baseline$seconds,
    accelerated_seconds = cellchat_accelerated$seconds,
    speedup = cellchat_baseline$seconds / cellchat_accelerated$seconds,
    interactions = nrow(cellchat_a),
    keys_identical = cellchat_keys_equal,
    max_score_absolute_difference = cellchat_max_diff,
    equivalence = "exact_scoped",
    acceleration_provenance = cellchat_accelerated$value$provenance$acceleration
  ),
  wgcna = list(
    baseline_seconds = wgcna_baseline$seconds,
    accelerated_seconds = wgcna_accelerated$seconds,
    speedup = wgcna_baseline$seconds / wgcna_accelerated$seconds,
    modules_identical = wgcna_modules_equal,
    max_eigengene_absolute_difference = wgcna_eigengene_diff,
    max_trait_correlation_absolute_difference = wgcna_trait_diff,
    tolerance = 1e-5,
    equivalence = "numeric_tolerance_scoped",
    acceleration_provenance = wgcna_accelerated$value$provenance$acceleration
  )
)

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
jsonlite::write_json(report, output, auto_unbox = TRUE, pretty = TRUE, null = "null")
cat(
  sprintf(
    "CellChat %.2fx (%d interactions; max difference %.3g)\n",
    report$cellchat$speedup,
    report$cellchat$interactions,
    report$cellchat$max_score_absolute_difference
  ),
  sprintf(
    "WGCNA %.2fx (modules identical; eigengene max difference %.3g)\n",
    report$wgcna$speedup,
    report$wgcna$max_eigengene_absolute_difference
  ),
  "Wrote ", normalizePath(output, winslash = "/", mustWork = TRUE), "\n",
  sep = ""
)
