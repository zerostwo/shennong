#!/usr/bin/env Rscript

# Build a deterministic 30,000-cell, three-capture PBMC benchmark fixture from
# independently generated public 10x datasets. Raw source objects are kept
# outside the repository; the derived fixture and its manifest belong under an
# ignored data-local directory.

args <- commandArgs(trailingOnly = TRUE)
all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", all_args, value = TRUE)
script_file <- if (length(script_arg)) {
  sub("^--file=", "", script_arg[[1L]])
} else {
  file.path(getwd(), "scripts", "real-data", "prepare-autozyme-pbmc30k.R")
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

.absolute_path <- function(path, must_work = FALSE) {
  path <- path.expand(path)
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", path)) path <- file.path(repo_root, path)
  normalizePath(path, winslash = "/", mustWork = must_work)
}

if ("--help" %in% args) {
  cat(
    "Usage: Rscript scripts/real-data/prepare-autozyme-pbmc30k.R [options]\n",
    "\n",
    "Options:\n",
    "  --source-dir PATH  Directory containing pbmc10k_seurat.rds,\n",
    "                     pbmc33k_seurat.rds, and pbmc68k_seurat.rds\n",
    "  --output PATH      Derived .qs2 fixture path\n",
    "  --cells-per-capture N  Cells sampled from each capture (default: 10000)\n",
    "  --seed N           Sampling seed (default: 717)\n",
    "  --help             Show this message\n",
    sep = ""
  )
  quit(status = 0L)
}

source_dir <- .absolute_path(.option(
  "source-dir",
  "/home/pi/dev/optimizing/UCell_AddModuleScore_UCell/data"
), must_work = TRUE)
output <- .absolute_path(.option(
  "output",
  file.path("data-local", "autozyme-benchmark", "pbmc-captures-30k.qs2")
))
cells_per_capture <- suppressWarnings(as.integer(.option("cells-per-capture", "10000")))
seed <- suppressWarnings(as.integer(.option("seed", "717")))
if (length(cells_per_capture) != 1L || is.na(cells_per_capture) || cells_per_capture < 1L) {
  stop("`--cells-per-capture` must be one positive integer.", call. = FALSE)
}
if (length(seed) != 1L || is.na(seed) || seed < 0L) {
  stop("`--seed` must be one non-negative integer.", call. = FALSE)
}

required <- c("digest", "jsonlite", "Matrix", "qs2", "SeuratObject")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing fixture package(s): ", paste(missing, collapse = ", "), ".", call. = FALSE)
}

sources <- list(
  pbmc10k_v3 = list(
    file = "pbmc10k_seurat.rds",
    sha256 = "c30f8ba7d4af1e136f5100cafd21d16f27090560646948a57e39044ba3b1cadb",
    source = "10x Genomics PBMC 10k, 3' v3, GRCh38",
    upstream_url = paste0(
      "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/",
      "pbmc_10k_v3_filtered_feature_bc_matrix.tar.gz"
    )
  ),
  pbmc33k = list(
    file = "pbmc33k_seurat.rds",
    sha256 = "263b81043cfd489861e72d17da9b62f283306ffe271d3391c462e37c5e513966",
    source = "TENxPBMCData::TENxPBMCData('pbmc33k', as.sparse = TRUE)",
    upstream_url = "https://doi.org/10.18129/B9.bioc.TENxPBMCData"
  ),
  pbmc68k_v1 = list(
    file = "pbmc68k_seurat.rds",
    sha256 = "fd0b4cfb7737d9d5e48d02f8ce5a64979a0e7226d78bf5c9eaaeee532c03551d",
    source = "10x Genomics Fresh PBMC 68k donor A, 3' v1, hg19",
    upstream_url = paste0(
      "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/fresh_68k_pbmc_donor_a/",
      "fresh_68k_pbmc_donor_a_filtered_gene_bc_matrices.tar.gz"
    )
  )
)

paths <- vapply(sources, function(spec) {
  path <- file.path(source_dir, spec$file)
  if (!file.exists(path)) stop("Missing source fixture: ", path, call. = FALSE)
  actual <- digest::digest(file = path, algo = "sha256", serialize = FALSE)
  if (!identical(actual, spec$sha256)) {
    stop("SHA-256 mismatch for source fixture: ", path, call. = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}, character(1))

feature_sets <- lapply(paths, function(path) rownames(readRDS(path)))
common_features <- Reduce(intersect, feature_sets)
common_features <- feature_sets[[1L]][feature_sets[[1L]] %in% common_features]
if (length(common_features) < 10000L) {
  stop("Fewer than 10,000 common gene symbols remain across captures.", call. = FALSE)
}

count_matrices <- vector("list", length(paths))
metadata <- vector("list", length(paths))
names(count_matrices) <- names(paths)
names(metadata) <- names(paths)
source_dimensions <- list()
selected_cells <- list()
for (index in seq_along(paths)) {
  capture <- names(paths)[[index]]
  object <- readRDS(paths[[index]])
  if (!inherits(object, "Seurat")) {
    stop("Source fixture is not a Seurat object: ", paths[[index]], call. = FALSE)
  }
  if (ncol(object) < cells_per_capture) {
    stop("Capture `", capture, "` has fewer requested cells.", call. = FALSE)
  }
  source_dimensions[[capture]] <- list(features = nrow(object), cells = ncol(object))
  set.seed(seed + index - 1L)
  cells <- sort(sample(colnames(object), cells_per_capture, replace = FALSE))
  selected_cells[[capture]] <- digest::digest(cells, algo = "sha256", serialize = TRUE)
  counts <- SeuratObject::LayerData(object, assay = "RNA", layer = "counts")[
    common_features,
    cells,
    drop = FALSE
  ]
  colnames(counts) <- paste(capture, colnames(counts), sep = "_")
  count_matrices[[capture]] <- methods::as(counts, "dgCMatrix")
  current_metadata <- object[[]][cells, , drop = FALSE]
  rownames(current_metadata) <- colnames(counts)
  current_metadata$capture <- capture
  current_metadata$source_cell_hash <- vapply(
    cells,
    function(cell) digest::digest(cell, algo = "sha256", serialize = FALSE),
    character(1)
  )
  current_metadata$barcode <- NULL
  metadata[[capture]] <- current_metadata
  rm(object, counts, current_metadata)
  invisible(gc())
}

counts <- do.call(cbind, count_matrices)
metadata <- do.call(rbind, unname(metadata))
metadata <- metadata[colnames(counts), , drop = FALSE]
object <- SeuratObject::CreateSeuratObject(
  counts = counts,
  assay = "RNA",
  project = "autozyme_pbmc_captures_30k",
  meta.data = metadata,
  min.cells = 0,
  min.features = 0
)
object@misc$benchmark_fixture <- list(
  schema_version = "shennong.autozyme-pbmc30k/v1",
  seed = seed,
  cells_per_capture = cells_per_capture,
  common_features = length(common_features),
  captures = names(paths),
  source_sha256 = vapply(sources, `[[`, character(1), "sha256"),
  selected_cell_sha256 = unlist(selected_cells, use.names = TRUE)
)

dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
qs2::qs_save(object, output)
artifact_sha256 <- digest::digest(file = output, algo = "sha256", serialize = FALSE)
capture_values <- as.character(object[["capture", drop = TRUE]])
manifest <- list(
  schema_version = "shennong.autozyme-pbmc30k/v1",
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  output = normalizePath(output, winslash = "/", mustWork = TRUE),
  artifact_sha256 = artifact_sha256,
  seed = seed,
  cells_per_capture = cells_per_capture,
  summary = list(
    features = nrow(object),
    cells = ncol(object),
    captures = length(unique(capture_values)),
    capture_counts = as.list(table(capture_values))
  ),
  common_feature_sha256 = digest::digest(common_features, algo = "sha256", serialize = TRUE),
  sources = Map(function(spec, path, dimensions, selected_hash) {
    c(
      spec,
      list(
        local_path = path,
        source_dimensions = dimensions,
        selected_cell_sha256 = selected_hash
      )
    )
  }, sources, paths, source_dimensions, selected_cells),
  package_versions = list(
    R = R.version.string,
    SeuratObject = as.character(utils::packageVersion("SeuratObject")),
    Matrix = as.character(utils::packageVersion("Matrix")),
    qs2 = as.character(utils::packageVersion("qs2"))
  )
)
manifest_path <- paste0(sub("[.]qs2$", "", output, ignore.case = TRUE), ".manifest.json")
jsonlite::write_json(
  manifest,
  manifest_path,
  auto_unbox = TRUE,
  pretty = TRUE,
  null = "null",
  digits = NA
)

cat(
  "Prepared ", nrow(object), " features x ", ncol(object), " cells across ",
  length(unique(capture_values)), " captures.\n",
  "Fixture: ", normalizePath(output, winslash = "/", mustWork = TRUE), "\n",
  "SHA-256: ", artifact_sha256, "\n",
  "Manifest: ", normalizePath(manifest_path, winslash = "/", mustWork = TRUE), "\n",
  sep = ""
)
