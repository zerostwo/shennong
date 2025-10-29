#' @export
sn_read <- function(path,
                    format,
                    to = NULL,
                    which,
                    row_names = NULL,
                    ...) {
  # Copy from rio::import(), modified to support some  Bioinformatics formats
  check_installed(pkg = "rio", reason = "to use `sn_read()` function.")

  rio:::.check_file(path, single_only = TRUE)
  if (R.utils::isUrl(path)) {
    path <- rio:::remote_to_local(path, format = format)
  }
  if ((path != "clipboard") && !file.exists(path)) {
    stop("No such path: ", path, call. = FALSE)
  }
  # Compressed path, f is a pretty bad name; but export() uses it.
  f <- rio:::find_compress(path)
  if (!is.na(f$compress)) {
    cpath <- path
    path <- f$path
    which <- ifelse(missing(which), 1, which)
    path <- rio:::parse_archive(cpath, which = which, path_type = f$compress)
    format <- rio:::.get_compressed_format(cpath, path, f$compress, format)
    # Reset which if `path` is zip or tar. #412
    which <- rio:::.reset_which(file_type = f$compress, which = which)
  }
  if (missing(format)) {
    format <- .get_info(path)$format
  } else {
    # Format such as "|"
    format <- rio:::.standardize_format(format)
  }
  class(path) <- c(paste0("rio_", format), class(path))
  if (missing(which)) {
    x <- rio:::.import(file = path, ...)
  } else {
    x <- rio:::.import(file = path, which = which, ...)
  }

  # Return non data.frame objects as is (if not setting class),
  # for example, dgCMatrix, Seurat, etc.
  ext_classes <- c(
    "dgCMatrix", "dgTMatrix", "dgRMatrix", "dsCMatrix", "dsTMatrix",
    "dsRMatrix", "CsparseMatrix", "RsparseMatrix", "sparseMatrix",
    "Seurat", "SingleCellExperiment", "InMemoryAnnData", "anndata",
    "BPCells"
  )
  if (inherits(x, ext_classes)) {
    return(x)
  }

  if (!is_null(x = row_names) && inherits(x = x, what = "data.frame")) {
    rownames(x = x) <- x[, row_names]
    x <- x[, -row_names]
  }

  # f R serialized object, just return it without setting object class
  if (inherits(path, c("rio_rdata", "rio_rds", "rio_json", "rio_qs", "rio_qs2")) && !inherits(x, "data.frame")) {
    return(x)
  }
  # Otherwise, make sure it's a data frame (or requested class)
  return(rio:::set_class(x, class = to))
}

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
  info <- rio:::.query_format(input = ext, file = file)
  if (info$format == "") {
    # Guess the format of the directory
    dir_format <- .guess_dir_format(path = file)
    if (!is_null(dir_format)) {
      info$format <- dir_format
      info$input <- dir_format
    }
  }
  return(info)
}

# If .get_info return format is "", that means it's a directory of files,
# So, we need to use a function to guess the format based on the files inside the directory.
.guess_dir_format <- function(path) {
  files <- list.files(path = path)
  if (all(any(grepl("matrix.mtx.gz$", files)) && any(grepl("barcodes.tsv.gz$", files)) && any(grepl("features.tsv.gz$", files)))) {
    log_info("Detected 10x format in directory: ", path)
    return("10x")
  } else if (all(any(grepl("matrix.mtx$", files)) && any(grepl("barcodes.tsv$", files)) && any(grepl("features.tsv$", files)))) {
    log_info("Detected STARsolo format in directory: ", path)
    return("starsolo")
  } else if (any(grepl("val$", files))) {
    log_info("Detected BPCells format in directory: ", path)
    return("bpcells")
  } else if (any(grepl("filtered_feature_bc_matrix.h5$", files))) {
    log_info("Detected 10x spatial format in directory: ", path)
    return("10x_spatial")
  } else {
    log_info("Cannot guess the format of the directory: ", path)
    return(NULL)
  }
}

#' @export
.import.rio_bpcells <- function(file, ...) {
  check_installed(pkg = "BPCells", reason = "to read BPCells format files.")
  mat <- BPCells::open_matrix_dir(dir = file, ...)

  return(mat)
}

#' @export
.import.rio_10x <- function(file, ...) {
  check_installed(pkg = "Seurat", reason = "to read 10x format files.")
  mat <- Seurat::Read10X(data.dir = file, ...)

  return(mat)
}

#' @export
.import.rio_starsolo <- function(file, ...) {
  check_installed(pkg = "Seurat", reason = "to read STARsolo format files.")
  mat <- Seurat::ReadSTARsolo(data.dir = file, ...)

  return(mat)
}

#' @export
.import.rio_h5ad <- function(file, ...) {
  check_installed(pkg = "anndataR", reason = "to read h5ad files.")
  mat <- anndataR::read_h5ad(path = file, ...)

  return(mat)
}

#' @export
.import.rio_h5 <- function(file, ...) {
  check_installed(pkg = c("Seurat", "hdf5r"), reason = "to read 10x h5 files.")
  mat <- Seurat::Read10X_h5(file)

  return(mat)
}

#' @title Import GMT files (from clusterProfiler)
#' @description
#' This function imports Gene Matrix Transposed (GMT) files into a two-column
#' data.frame of term–gene pairs.
#'
#' Adapted from `clusterProfiler::read.gmt()` (Y. Yu et al., 2012–2024),
#' which is licensed under the Artistic License 2.0.
#'
#' @details
#' The original implementation is available in the clusterProfiler package:
#' https://github.com/YuLab-SMU/clusterProfiler
#'
#' @return A data.frame with two columns: `term` and `gene`.
#' @seealso [clusterProfiler::read.gmt()]
#' @author Yu Guangchuang (original implementation)
#' @export
.import.rio_gmt <- function(file, ...) {
  x <- readLines(gmtfile)
  res <- strsplit(x, "\t")
  names(res) <- vapply(res, function(y) y[1], character(1))
  res <- lapply(res, "[", -c(1:2))
  ont2gene <- stack(res)
  ont2gene <- ont2gene[, c("ind", "values")]
  colnames(ont2gene) <- c("term", "gene")

  return(ont2gene)
}

#' Internal shared importer for genome annotation / track files
#' @keywords internal
.import_genomic_tracks <- function(file, ext = NULL, ...) {
  check_installed(
    pkg = "rtracklayer",
    reason = "to import genomic annotation/track files (gtf/gff/bed/bw/etc.)"
  )

  if (is.null(ext)) {
    ext <- tolower(tools::file_ext(file))
  }

  fmt_map <- list(
    gtf = "gtf",
    gff = "gff",
    gff1 = "gff1",
    gff2 = "gff2",
    gff3 = "gff3",
    wig = "wig",
    bed = "bed",
    bed15 = "bed15",
    bedgraph = "bedGraph",
    bw = "bigWig"
  )

  format <- fmt_map[[ext]]
  if (is.null(format)) {
    # fallback: let rtracklayer guess
    format <- ext
  }

  gr <- rtracklayer::import(con = file, format = format, ...)

  if (ext %in% c("gtf", "gff", "gff1", "gff2", "gff3", "bed", "bed15")) {
    return(as.data.frame(gr))
  } else {
    return(gr)
  }
}

# Vector of extensions you care about
.genomic_exts <- c("gtf", "gff", "gff1", "gff2", "gff3", "wig", "bed", "bed15", "bedGraph", "bw")
for (ext in .genomic_exts) {
  fn <- function(file, ...) {
    .import_genomic_tracks(file, ext = tolower(ext), ...)
  }
  assign(paste0(".import.rio_", ext), fn, envir = globalenv())
}

#' @export
sn_write <- function(x, path = NULL, to = NULL, ...) {
  check_installed(pkg = "rio", reason = "to use `sn_write()` function.")
  # Copy from rio::export(), modified to support some  Bioinformatics formats
  rio:::.check_file(path, single_only = TRUE)
  if (is.null(path) && is.null(to)) {
    stop("Must specify 'path' and/or 'to'")
  }

  if (is.null(path)) {
    format <- rio:::.standardize_format(to)
    path <- paste0(as.character(substitute(x)), ".", format)
    compress <- NA_character_
  } else {
    cfile <- path
    f <- rio:::find_compress(path)
    path <- f$file
    compress <- f$compress
    format <- ifelse(is.null(to), .get_info(path)$input, tolower(to))
  }

  rio:::.check_tar_support(compress, getRversion())
  format <- rio:::.standardize_format(format)
  outfile <- path
  if (is.matrix(x) || inherits(x, "ArrowTabular")) {
    x <- as.data.frame(x)
  }
  if (!is.data.frame(x) && is.list(x) && length(x) == 1 && is.data.frame(x[[1]]) &&
    !format %in% c("xlsx", "html", "rdata", "rds", "json", "qs", "qs2", "fods", "ods")) {
    x <- x[[1]] ## fix 385
  }

  if (!is.data.frame(x) && !format %in% c(
    "xlsx", "html", "rdata", "rds", "json", "qs", "qs2", "fods", "ods",
    "bpcells", "h5ad", "h5"
  )) {
    stop("'x' is not a data.frame or matrix", call. = FALSE)
  }
  if (format == "gz") {
    format <- .get_info(tools::file_path_sans_ext(path, compression = FALSE))$format
    if (format != "csv") {
      stop("gz is only supported for csv (for now).", call. = FALSE)
    }
  }
  rio:::.create_directory_if_not_exists(file = path) ## fix 347
  class(path) <- c(paste0("rio_", format), class(path))
  rio:::.export(file = path, x = x, ...)

  if (!is.na(compress)) {
    cfile <- rio:::compress_out(cfile = cfile, filename = path, type = compress)
    unlink(path)
    return(invisible(cfile))
  }
  invisible(unclass(outfile))
}

#' @export
.export.rio_bpcells <- function(file, x, overwrite = FALSE, ...) {
  check_installed_github(
    pkg = "BPCells",
    repo = "bnprks/BPCells/r",
    reason = "to write BPCells format files."
  )
  BPCells::write_matrix_dir(mat = x, dir = file, overwrite = overwrite, ...)
}

#' @export
.export.rio_h5ad <- function(file, x, mode = "w", ...) {
  check_installed_github(
    pkg = "anndataR", repo = "scverse/anndataR",
    reason = "to write h5ad files."
  )
  check_installed(pkg = c("SingleCellExperiment", "rhdf5"), reason = "to write h5ad files.")
  if (!inherits(x = x, what = "SingleCellExperiment")) {
    sce <- Seurat::as.SingleCellExperiment(x = x)
  }
  anndataR::write_h5ad(object = sce, path = file, mode = mode, ...)
}

#' @export
.export.rio_h5 <- function(file, x, ...) {
  check_installed_github(
    pkg = "BPCells",
    repo = "bnprks/BPCells/r",
    reason = "to write 10x h5 files."
  )
  BPCells::write_matrix_10x_hdf5(
    mat = x, path = file, ...
  )
}

#' @export
sn_add_data_from_anndata <- function(object, metadata_path = NULL, umap_path = NULL, key = "umap") {
  if (!is_null(x = umap_path)) {
    umap <- sn_read(path = umap_path, row_names = 1) |>
      Matrix::as.matrix()
    object <- object[, rownames(umap)]
    object[[key]] <- Seurat::CreateDimReducObject(umap, key = key)
  }
  if (!is_null(x = metadata_path)) {
    metadata <- sn_read(path = metadata_path, row_names = 1)
    # object@meta.data <- metadata
    object <- Seurat::AddMetaData(object = object, metadata = metadata)
  }
  object
}
