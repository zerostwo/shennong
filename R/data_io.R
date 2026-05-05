#' Read tabular and bioinformatics file formats
#'
#' A Shennong wrapper around `rio::import()` with support for common
#' single-cell and omics formats such as 10x directories, `.h5`, `.h5ad`,
#' BPCells directories, and GMT files.
#'
#' @param path Path, URL, or supported directory to import.
#' @param format Optional explicit format override.
#' @param to Optional target class forwarded to `rio`.
#' @param which Optional archive member index for compressed archives.
#' @param row_names Optional column index or name to promote to row names when
#'   importing a data frame.
#' @param file Input path used by the exported `rio` adapter methods.
#' @param ... Additional arguments forwarded to the underlying importer.
#'
#' @return The imported object. Tabular inputs are typically returned as data
#'   frames, while supported matrix-like formats are returned in their native
#'   classes.
#'
#' @examples
#' tmp <- tempfile(fileext = ".csv")
#' write.csv(mtcars[1:3, 1:3], tmp, row.names = FALSE)
#' sn_read(tmp)
#'
#' @seealso [sn_write()]
#' @export
sn_read <- function(path,
                    format,
                    to = NULL,
                    which,
                    row_names = NULL,
                    ...) {
  check_installed(pkg = "rio", reason = "to use `sn_read()` function.")

  format_supplied <- !missing(format)
  format <- if (format_supplied) .sn_standardize_format(format) else .sn_infer_format(path)

  if ((path != "clipboard") && !.sn_is_url(path) && !file.exists(path)) {
    stop("No such path: ", path, call. = FALSE)
  }

  if (.sn_is_custom_format(format)) {
    local_path <- .sn_download_custom_source(path, format)
    on.exit(if (!identical(local_path, path)) unlink(local_path), add = TRUE)
    return(.sn_dispatch_custom_reader(local_path, format = format, ...))
  }

  args <- c(
    list(file = path, setclass = to %||% getOption("rio.import.class", "data.frame")),
    if (format_supplied) list(format = format) else list(),
    if (!missing(which)) list(which = which) else list(),
    list(...)
  )
  x <- do.call(rio::import, args = args)

  if (!is_null(x = row_names) && inherits(x, "data.frame")) {
    if (length(row_names) != 1L || !(is.character(row_names) || is.numeric(row_names))) {
      stop("`row_names` must be a single column name or column index.", call. = FALSE)
    }
    row_index <- if (is.character(row_names)) {
      match(row_names, colnames(x))
    } else {
      as.integer(row_names)
    }
    if (is.na(row_index) || row_index < 1L || row_index > ncol(x)) {
      stop("`row_names` must identify a column in the imported data frame.", call. = FALSE)
    }
    rownames(x) <- x[[row_index]]
    x <- x[, -row_index, drop = FALSE]
  }

  x
}

#' List detected 10x Genomics output paths under a root folder
#'
#' This helper scans a root directory for typical 10x Genomics \code{outs/}
#' layouts and returns one selected path per detected sample. By default it
#' returns the \code{outs/} directory itself so the output can be passed
#' directly to [sn_initialize_seurat_object()]. The \code{path_type} argument
#' can instead return the filtered matrix path, raw matrix path, H5 files, or
#' \code{metrics_summary.csv} file.
#'
#' @param path Root directory to scan.
#' @param recursive Logical; if \code{TRUE}, scan subdirectories recursively.
#'   Defaults to \code{TRUE}.
#' @param include_root Logical; if \code{TRUE}, also test \code{path} itself.
#'   Defaults to \code{TRUE}.
#' @param path_type Which detected 10x path to return. Defaults to
#'   \code{"outs"}. Supported values are \code{"outs"}, \code{"filtered"},
#'   \code{"raw"}, \code{"filtered_h5"}, \code{"raw_h5"}, and
#'   \code{"metrics"}.
#'
#' @return A named character vector of matching paths, where names are inferred
#'   sample identifiers.
#'
#' @examples
#' root <- tempfile("tenx-root-")
#' dir.create(root)
#' sample_dir <- file.path(root, "sample1", "outs")
#' dir.create(file.path(sample_dir, "filtered_feature_bc_matrix"), recursive = TRUE)
#' dir.create(file.path(sample_dir, "raw_feature_bc_matrix"), recursive = TRUE)
#' file.create(file.path(sample_dir, "metrics_summary.csv"))
#' sn_list_10x_paths(root)
#' sn_list_10x_paths(root, path_type = "filtered")
#'
#' @export
sn_list_10x_paths <- function(path,
                              recursive = TRUE,
                              include_root = TRUE,
                              path_type = c("outs", "filtered", "raw", "filtered_h5", "raw_h5", "metrics")) {
  if (!dir.exists(path)) {
    stop("`path` must be an existing directory.", call. = FALSE)
  }

  path_type <- match.arg(path_type)

  candidates <- .sn_find_10x_candidates(
    path = path,
    recursive = recursive,
    include_root = include_root
  )

  resolved <- lapply(candidates, .sn_locate_10x_outs_paths)
  resolved <- resolved[!vapply(resolved, is.null, logical(1))]
  if (length(resolved) == 0) {
    return(character(0))
  }

  outs_paths <- vapply(resolved, `[[`, character(1), "outs_path")
  resolved <- resolved[!duplicated(outs_paths)]

  selected <- vapply(
    resolved,
    FUN = function(info) .sn_select_10x_path(info, path_type = path_type) %||% NA_character_,
    FUN.VALUE = character(1),
    USE.NAMES = FALSE
  )
  keep <- !is.na(selected)
  selected <- selected[keep]
  sample_names <- vapply(
    resolved[keep],
    FUN = function(info) info$sample_name %||% basename(info$outs_path),
    FUN.VALUE = character(1),
    USE.NAMES = FALSE
  )
  names(selected) <- sample_names
  selected <- selected[order(names(selected), selected)]
  selected
}

.sn_find_10x_candidates <- function(path, recursive = TRUE, include_root = TRUE) {
  normalized_root <- normalizePath(path, winslash = "/", mustWork = TRUE)

  candidates <- if (isTRUE(recursive)) {
    .sn_find_10x_candidates_recursive(normalized_root)
  } else {
    list.dirs(path = normalized_root, recursive = FALSE, full.names = TRUE)
  }

  if (isTRUE(include_root)) {
    candidates <- c(normalized_root, candidates)
  }

  candidates <- unique(candidates[dir.exists(candidates)])
  if (!isTRUE(include_root)) {
    candidates <- setdiff(candidates, normalized_root)
  }

  candidates
}

.sn_find_10x_candidates_recursive <- function(path) {
  find_bin <- Sys.which("find")[["find"]]
  if (.Platform$OS.type == "unix" && nzchar(find_bin)) {
    find_args <- c(
      path,
      "-type", "d",
      "-name", "outs",
      "-print"
    )
    found <- tryCatch(
      suppressWarnings(system2(find_bin, args = find_args, stdout = TRUE, stderr = FALSE)),
      error = function(e) character(0)
    )
    found <- found[nzchar(found)]
    if (length(found) > 0) {
      return(found)
    }
  }

  list.dirs(path = path, recursive = TRUE, full.names = TRUE)
}

.sn_select_10x_path <- function(info, path_type = c("outs", "filtered", "raw", "filtered_h5", "raw_h5", "metrics")) {
  path_type <- match.arg(path_type)

  switch(path_type,
    outs = info$outs_path %||% NULL,
    filtered = info$filtered_path %||% NULL,
    raw = info$raw_path %||% NULL,
    filtered_h5 = if (!is_null(info$filtered_path) && file.exists(info$filtered_path) && !dir.exists(info$filtered_path)) info$filtered_path else NULL,
    raw_h5 = if (!is_null(info$raw_path) && file.exists(info$raw_path) && !dir.exists(info$raw_path)) info$raw_path else NULL,
    metrics = info$metrics_path %||% NULL
  )
}

.sn_is_url <- function(path) {
  is.character(path) &&
    length(path) == 1 &&
    grepl("^(https?|ftp)://", path, ignore.case = TRUE)
}

.sn_standardize_format <- function(format) {
  if (format %in% c(",", ";", "|")) {
    return(format)
  }
  tolower(format)
}

.sn_infer_format <- function(path) {
  if (tolower(path) == "clipboard") {
    return("clipboard")
  }

  if (.sn_is_url(path)) {
    path <- strsplit(path, "?", fixed = TRUE)[[1]][1]
  }

  if (dir.exists(path)) {
    return(.guess_dir_format(path))
  }

  ext <- tolower(tools::file_ext(path))
  if (ext %in% c("gz", "bz2", "xz")) {
    path <- tools::file_path_sans_ext(path)
    ext <- tolower(tools::file_ext(path))
  }

  if (nzchar(ext)) ext else NULL
}

.sn_is_custom_format <- function(format) {
  !is_null(format) && format %in% c(
    "10x", "10x_spatial", "starsolo", "bpcells", "h5ad", "h5", "gmt", "qs", "qs2", .genomic_exts
  )
}

.sn_download_custom_source <- function(path, format) {
  if (!.sn_is_url(path)) {
    return(path)
  }

  ext <- tools::file_ext(strsplit(path, "?", fixed = TRUE)[[1]][1])
  if (!nzchar(ext)) {
    ext <- format
  }
  tmp <- tempfile(fileext = paste0(".", ext))
  curl::curl_download(url = path, destfile = tmp, quiet = TRUE)
  tmp
}

.sn_dispatch_custom_reader <- function(path, format, ...) {
  if (format %in% .genomic_exts) {
    return(.import_genomic_tracks(file = path, ext = format, ...))
  }

  reader <- switch(format,
    bpcells = .import.rio_bpcells,
    `10x` = .import.rio_10x,
    `10x_spatial` = .import.rio_10x_spatial,
    starsolo = .import.rio_starsolo,
    h5ad = .import.rio_h5ad,
    h5 = .import.rio_h5,
    gmt = .import.rio_gmt,
    qs = .import.rio_qs,
    qs2 = .import.rio_qs2,
    NULL
  )

  if (is_null(reader)) {
    stop("Unsupported import format: ", format, call. = FALSE)
  }

  reader(file = path, ...)
}

# If .get_info return format is "", that means it's a directory of files,
# So, we need to use a function to guess the format based on the files inside the directory.
.guess_dir_format <- function(path) {
  files <- list.files(path = path)
  if (all(any(grepl("matrix.mtx.gz$", files)) && any(grepl("barcodes.tsv.gz$", files)) && any(grepl("features.tsv.gz$", files)))) {
    .sn_log_info("Detected 10x format in directory: {path}.")
    return("10x")
  } else if (all(any(grepl("matrix.mtx$", files)) && any(grepl("barcodes.tsv$", files)) && any(grepl("features.tsv$", files)))) {
    .sn_log_info("Detected STARsolo format in directory: {path}.")
    return("starsolo")
  } else if (any(grepl("val$", files))) {
    .sn_log_info("Detected BPCells format in directory: {path}.")
    return("bpcells")
  } else if (any(grepl("filtered_feature_bc_matrix.h5$", files))) {
    .sn_log_info("Detected 10x spatial format in directory: {path}.")
    return("10x_spatial")
  } else {
    .sn_log_info("Could not infer a directory format for: {path}.")
    return(NULL)
  }
}

#' Import a BPCells matrix directory through `rio`
#'
#' This adapter is exported so `rio` can dispatch to it.
#'
#' @rdname sn_read
#' @export
.import.rio_bpcells <- function(file, ...) {
  check_installed(pkg = "BPCells", reason = "to read BPCells format files.")
  mat <- BPCells::open_matrix_dir(dir = file, ...)

  return(mat)
}

#' @rdname sn_read
#' @export
.import.rio_10x <- function(file, ...) {
  check_installed(pkg = "Seurat", reason = "to read 10x format files.")
  mat <- Seurat::Read10X(data.dir = file, ...)

  return(mat)
}

#' @rdname sn_read
#' @export
.import.rio_10x_spatial <- function(file, ...) {
  check_installed(pkg = "Seurat", reason = "to read 10x spatial format files.")
  Seurat::Load10X_Spatial(data.dir = file, ...)
}

#' @rdname sn_read
#' @export
.import.rio_starsolo <- function(file, ...) {
  check_installed(pkg = "Seurat", reason = "to read STARsolo format files.")
  mat <- Seurat::ReadSTARsolo(data.dir = file, ...)

  return(mat)
}

#' @rdname sn_read
#' @export
.import.rio_h5ad <- function(file, ...) {
  check_installed(pkg = "anndataR", reason = "to read h5ad files.")
  mat <- anndataR::read_h5ad(path = file, ...)

  return(mat)
}

#' @rdname sn_read
#' @export
.import.rio_h5 <- function(file, ...) {
  check_installed(pkg = c("Seurat", "hdf5r"), reason = "to read 10x h5 files.")
  mat <- Seurat::Read10X_h5(file)

  return(mat)
}

#' @rdname sn_read
#' @export
.import.rio_qs <- function(file, ...) {
  check_installed(pkg = "qs", reason = "to read qs files.")
  getExportedValue("qs", "qread")(file = file, ...)
}

#' @rdname sn_read
#' @export
.import.rio_qs2 <- function(file, ...) {
  check_installed(pkg = "qs2", reason = "to read qs2 files.")
  qs2::qs_read(file = file, ...)
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
#' @seealso The original `read.gmt()` implementation in the `clusterProfiler`
#'   package.
#' @author Yu Guangchuang (original implementation)
#' @rdname sn_read
#' @export
.import.rio_gmt <- function(file, ...) {
  x <- readLines(file)
  res <- strsplit(x, "\t")
  names(res) <- vapply(res, function(y) y[1], character(1))
  res <- lapply(res, "[", -c(1:2))
  ont2gene <- utils::stack(res)
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

.genomic_exts <- c("gtf", "gff", "gff1", "gff2", "gff3", "wig", "bed", "bed15", "bedGraph", "bw")

#' Write tabular and bioinformatics file formats
#'
#' A Shennong wrapper around `rio::export()` with support for common
#' single-cell formats such as BPCells, `.h5ad`, and 10x `.h5`.
#'
#' @param x Object to write.
#' @param path Output path. Missing parent directories are created
#'   automatically before dispatching the selected writer.
#' @param to Optional format override when it cannot be inferred from `path`.
#' @param file Output path used by the exported `rio` adapter methods.
#' @param overwrite Logical; overwrite an existing BPCells directory.
#' @param mode File mode passed through to `anndataR::write_h5ad()`.
#' @param ... Additional arguments forwarded to the underlying exporter.
#'
#' @return Invisibly returns the output path.
#'
#' @examples
#' tmp <- tempfile(fileext = ".csv")
#' sn_write(mtcars[1:3, 1:3], tmp)
#'
#' @seealso [sn_read()]
#' @export
sn_write <- function(x, path = NULL, to = NULL, ...) {
  check_installed(pkg = "rio", reason = "to use `sn_write()` function.")

  if (is.null(path) && is.null(to)) {
    stop("Must specify 'path' and/or 'to'")
  }

  format <- if (!is_null(to)) {
    .sn_standardize_format(to)
  } else {
    .sn_infer_format(path)
  }

  if (is.null(path)) {
    path <- paste0(as.character(substitute(x)), ".", format)
  }
  .sn_ensure_output_parent(path)

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

  if (.sn_is_custom_format(format)) {
    .sn_dispatch_custom_writer(path = path, x = x, format = format, ...)
    return(invisible(path))
  }

  args <- c(
    list(x = x, file = path),
    if (!is_null(to)) list(format = format) else list(),
    list(...)
  )
  do.call(rio::export, args = args)
  invisible(path)
}

.sn_ensure_output_parent <- function(path) {
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop("`path` must be a non-empty character scalar.", call. = FALSE)
  }

  parent <- dirname(path)
  if (is.na(parent) || !nzchar(parent) || identical(parent, ".")) {
    return(invisible(path))
  }

  if (!dir.exists(parent)) {
    dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  }
  if (!dir.exists(parent)) {
    stop("Failed to create output directory: ", parent, call. = FALSE)
  }

  invisible(path)
}

.sn_dispatch_custom_writer <- function(path, x, format, ...) {
  writer <- switch(format,
    bpcells = .export.rio_bpcells,
    h5ad = .export.rio_h5ad,
    h5 = .export.rio_h5,
    qs = .export.rio_qs,
    qs2 = .export.rio_qs2,
    NULL
  )

  if (is_null(writer)) {
    stop("Unsupported export format: ", format, call. = FALSE)
  }

  writer(file = path, x = x, ...)
}

#' Export a matrix to a BPCells directory through `rio`
#'
#' This adapter is exported so `rio` can dispatch to it.
#'
#' @rdname sn_write
#' @export
.export.rio_bpcells <- function(file, x, overwrite = FALSE, ...) {
  check_installed_github(
    pkg = "BPCells",
    repo = "bnprks/BPCells/r",
    reason = "to write BPCells format files."
  )
  BPCells::write_matrix_dir(mat = x, dir = file, overwrite = overwrite, ...)
}

#' @rdname sn_write
#' @export
.export.rio_h5ad <- function(file, x, mode = "w", ...) {
  check_installed_github(
    pkg = "anndataR", repo = "scverse/anndataR",
    reason = "to write h5ad files."
  )
  check_installed(pkg = c("SingleCellExperiment", "rhdf5"), reason = "to write h5ad files.")
  if (inherits(x = x, what = "SingleCellExperiment")) {
    sce <- x
  } else {
    check_installed(pkg = "Seurat", reason = "to convert Seurat objects before writing h5ad files.")
    sce <- Seurat::as.SingleCellExperiment(x = x)
  }
  anndataR::write_h5ad(object = sce, path = file, mode = mode, ...)
}

#' @rdname sn_write
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

#' @rdname sn_write
#' @export
.export.rio_qs <- function(file, x, ...) {
  check_installed(pkg = "qs", reason = "to write qs files.")
  getExportedValue("qs", "qsave")(x = x, file = file, ...)
}

#' @rdname sn_write
#' @export
.export.rio_qs2 <- function(file, x, ...) {
  check_installed(pkg = "qs2", reason = "to write qs2 files.")
  qs2::qs_save(object = x, file = file, ...)
}

#' Add metadata and embeddings exported from AnnData
#'
#' Reads metadata and/or embedding tables from disk and adds them back to a
#' Seurat object.
#'
#' @param object A Seurat object.
#' @param metadata_path Optional path to a metadata table whose first column
#'   contains cell names.
#' @param umap_path Optional path to an embedding table whose first column
#'   contains cell names.
#' @param key Reduction key used when storing the imported embedding.
#'
#' @return The updated Seurat object.
#'
#' @examples
#' \dontrun{
#' object <- sn_add_data_from_anndata(
#'   object,
#'   metadata_path = "obs.csv",
#'   umap_path = "umap.csv"
#' )
#' }
#'
#' @export
sn_add_data_from_anndata <- function(object, metadata_path = NULL, umap_path = NULL, key = "umap") {
  read_anndata_table <- function(path) {
    if (identical(.sn_infer_format(path), "csv")) {
      return(utils::read.csv(
        file = path,
        row.names = 1,
        check.names = FALSE
      ))
    }

    sn_read(path = path, row_names = 1)
  }

  if (!is_null(x = umap_path)) {
    umap <- read_anndata_table(umap_path) |>
      Matrix::as.matrix()
    object <- object[, rownames(umap)]
    reduction_key <- if (grepl("_$", key)) key else paste0(key, "_")
    object[[key]] <- Seurat::CreateDimReducObject(
      embeddings = umap,
      key = reduction_key,
      assay = SeuratObject::DefaultAssay(object)
    )
  }
  if (!is_null(x = metadata_path)) {
    metadata <- read_anndata_table(metadata_path)
    object <- Seurat::AddMetaData(object = object, metadata = metadata)
  }
  .sn_log_seurat_command(object = object, name = "sn_add_data_from_anndata")
}
