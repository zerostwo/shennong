.get_info <- function(file) {
  rio:::.check_file(file, single_only = TRUE)
  if (tolower(file) == "clipboard") {
    return(rio:::.query_format(input = "clipboard", file = "clipboard"))
  }
  if (isFALSE(R.utils::isUrl(file))) {
    ext <- tolower(tools::file_ext(file))
  } else {
    parsed <- strsplit(strsplit(file, "?", fixed = TRUE)[[1]][1], "/", fixed = TRUE)[[1]]
    url_file <- parsed[length(parsed)]
    ext <- tolower(tools::file_ext(url_file))
  }
  # if (ext == "") {
  #   print("This is a dir")
  # }
  return(rio:::.query_format(input = ext, file = file))
}

#' @export
sn_read <- function(path, from = NULL, to = NULL, row_names = NULL, ...) {
  f <- rio:::find_compress(f = path)
  if (!is.na(f$compress)) {
    file <- f$file
  } else {
    file <- path
  }

  info <- .get_info(file = file)
  input <- info$input
  format <- from %||% info$format
  if (input == "h5ad") {
    to <- to %||% "InMemoryAnnData"
  } else {
    to <- to %||% "data.frame"
  }

  if (R.utils::isUrl(path)) {
    temp_file <- tempfile(fileext = paste0(".", input))
    u <- curl::curl_fetch_memory(url = path)
    writeBin(object = u$content, con = temp_file)
    path <- temp_file
  }

  rlang::arg_match(arg = to, values = c(
    "data.frame", "tbl_df", "tbl", "tibble", "arrow", "arrow_table", "data.table",
    "Seurat", "SeuratObject", "SingleCellExperiment", "sce", "InMemoryAnnData", "anndata", "BPCells"
  ))

  if (input == "") {
    files <- list.files(path = path)
    if (all(any(grepl("matrix.mtx.gz$", files)) && any(grepl("barcodes.tsv.gz$", files)) && any(grepl("features.tsv.gz$", files)))) {
      logger::log_info("Detected 10x format, reading with `Seurat::Read10X()`")
      x <- Seurat::Read10X(data.dir = path)
    } else if (all(any(grepl("matrix.mtx$", files)) && any(grepl("barcodes.tsv$", files)) && any(grepl("features.tsv$", files)))) {
      x <- Seurat::ReadSTARsolo(data.dir = path, ...)
    } else if (all(grepl("val$", files))) {
      x <- BPCells::open_matrix_dir(dir = path, ...)
    } else if (any(grepl("filtered_feature_bc_matrix.h5$", files))) {
      x <- Seurat::Load10X_Spatial(data.dir = path, ...)
    }
  } else {
    if (is.na(format)) {
      if (input == "h5") {
        x <- Seurat::Read10X_h5(filename = path, ...)
      } else if (input == "h5ad") {
        x <- anndataR::read_h5ad(path = path, to = to, ...)
      } else if (input %in% c("gtf", "gff", "gff1", "gff2", "gff3", "wig", "bed", "bed15", "bedGraph", "bw")) {
        x <- rtracklayer::import(con = path, format = format, ...)
      } else if (input == "gmt") {
        x <- clusterProfiler::read.gmt(gmtfile = path)
      }
    } else {
      # With extension file and fio support format
      x <- rio::import(file = path, format = format, setclass = to, ...)
    }
  }

  if (!is.null(x = row_names) && inherits(x = x, what = "data.frame")) {
    rownames(x = x) <- x[, row_names]
    x <- x[, -row_names]
  }

  return(x)
}

#' @export
sn_add_data_from_anndata <- function(object, umap_path = NULL, metadata_path = NULL) {
  if (!is_null(x = umap_path)) {
    umap <- read.csv(file = umap_path, row.names = 1) |>
      Matrix::as.matrix()
    object <- object[, rownames(umap)]
    object[["umap"]] <- Seurat::CreateDimReducObject(umap, key = "umap")
  }
  if (!is_null(x = metadata_path)) {
    metadata <- read.csv(file = metadata_path, row.names = 1)
    object@meta.data <- metadata
  }
  object
}

#' @export
sn_write_h5ad <- function(object, path, mode = "w") {
  dir_path <- dirname(path = path)
  if (!dir.exists(paths = dir_path)) {
    dir.create(path = dir_path, recursive = TRUE)
  }

  if (!inherits(x = object, what = "SingleCellExperiment")) {
    object <- Seurat::as.SingleCellExperiment(x = object)
  }

  anndataR::write_h5ad(object = object, path = path, mode = mode)

  rlang::inform(paste("Successfully wrote h5ad file to:", path))
}
