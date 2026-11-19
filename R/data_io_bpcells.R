# A QS2 snapshot retains BPCells' lazy operations. Record the backing inputs
# once per directory, even if several assays or layers reference the same one.
.sn_map_bpcells_dirs <- function(object, visit) {
  map_matrix <- function(matrix) {
    if (!inherits(matrix, "IterableMatrix")) return(matrix)
    inputs <- BPCells::all_matrix_inputs(matrix)
    changed <- FALSE
    for (i in seq_along(inputs)) {
      if (!inherits(inputs[[i]], "MatrixDir")) next
      path <- visit(inputs[[i]])
      if (!identical(path, inputs[[i]]@dir)) {
        inputs[[i]]@dir <- path
        changed <- TRUE
      }
    }
    if (changed) BPCells::all_matrix_inputs(matrix) <- inputs
    matrix
  }
  if (!inherits(object, "Seurat")) return(map_matrix(object))
  for (i in seq_along(object@assays)) {
    assay <- object@assays[[i]]
    if ("layers" %in% methods::slotNames(assay)) {
      # Work with the stored layer orientation and feature/cell maps.
      assay@layers <- lapply(assay@layers, map_matrix)
      object@assays[[i]] <- assay
    }
  }
  object
}

.sn_relative_bpcells_dir <- function(path, snapshot_dir) {
  # Normalize first so a moved directory tree resolves relative to its new
  # parent even when the working directory or the spelling of a path changes.
  parent <- normalizePath(snapshot_dir, winslash = "/", mustWork = TRUE)
  up <- character()
  repeat {
    prefix <- paste0(sub("/$", "", parent), "/")
    if (identical(parent, path)) return(paste(c(up, "."), collapse = "/"))
    if (startsWith(path, prefix)) {
      return(paste(c(up, substring(path, nchar(prefix) + 1L)), collapse = "/"))
    }
    next_parent <- dirname(parent)
    if (identical(next_parent, parent)) return(NULL) # different Windows volumes
    parent <- next_parent
    up <- c(up, "..")
  }
}

.sn_bpcells_files <- function(matrix) {
  index <- if (matrix@compressed) c("index_data", "index_idx", "index_idx_offsets", "index_starts") else "index"
  values <- if (matrix@compressed && matrix@type == "uint32_t") c("val_data", "val_idx", "val_idx_offsets") else "val"
  c("version", "shape", "storage_order", "col_names", "row_names", "idxptr", index, values)
}

.sn_check_bpcells_overwrite <- function(path) {
  matrix <- tryCatch(BPCells::open_matrix_dir(path), error = function(e) NULL)
  if (is.null(matrix)) stop("Refusing to overwrite a directory that is not a BPCells matrix: ", path, call. = FALSE)
  files <- list.files(path, all.files = TRUE, no.. = TRUE)
  links <- Sys.readlink(c(path, file.path(path, files)))
  if (length(setdiff(files, .sn_bpcells_files(matrix))) ||
      any(nzchar(links[!is.na(links)])) || any(dir.exists(file.path(path, files)))) {
    stop("BPCells directory contains unrecognized files or links: ", path, call. = FALSE)
  }
}

.sn_bpcells_fingerprint <- function(path, matrix) {
  # Hash bytes on disk without materializing the gene-by-cell matrix. Include
  # file names to distinguish shape, ordering, names, indices and values.
  files <- list.files(path, all.files = TRUE, no.. = TRUE)
  files <- files[!dir.exists(file.path(path, files))]
  if (!length(files) || !"version" %in% files) {
    stop("Missing BPCells matrix directory or version file: ", path, call. = FALSE)
  }
  version <- readLines(file.path(path, "version"), n = 1L, warn = FALSE)
  if (length(version) == 1L && grepl("^(packed|unpacked)-(uint|float|double)-matrix-v2$", version)) {
    expected <- .sn_bpcells_files(matrix)
    required <- setdiff(expected, c("row_names", "col_names"))
    if (!all(required %in% files)) stop("Incomplete BPCells matrix identity: ", path, call. = FALSE)
    if (!all(c("row_names", "col_names") %in% files)) {
      # Let the installed reader decide whether this layout permits absent
      # name files. Unnamed v2 matrices generally include empty files.
      BPCells::open_matrix_dir(path)
    }
    files <- intersect(expected, files)
  }
  # For future layout versions, conservatively include all files until the
  # format is known. The ordering must be independent of the reader's locale.
  files <- sort(files, method = "radix")
  checksums <- vapply(files, function(name) {
    digest::digest(file = file.path(path, name), algo = "sha256")
  }, character(1))
  digest::digest(checksums, algo = "sha256", serializeVersion = 2L)
}

.sn_prepare_qs2_bpcells <- function(object, file) {
  records <- list()
  remember <- function(matrix) {
    source <- matrix@dir
    if (is.null(records[[source]])) {
      path <- normalizePath(source, winslash = "/", mustWork = TRUE)
      records[[source]] <<- list(
        relative = .sn_relative_bpcells_dir(path, dirname(file)),
        sha256 = .sn_bpcells_fingerprint(path, matrix)
      )
    }
    source
  }
  object <- .sn_map_bpcells_dirs(object, remember)
  if (length(records)) {
    attr(object, "shennong_qs2_bpcells") <- list(version = 1L, directories = records)
  }
  object
}

.sn_restore_qs2_bpcells <- function(object, file, bpcells_paths) {
  if (!inherits(object, c("Seurat", "IterableMatrix")) && is.null(bpcells_paths)) return(object)
  manifest <- attr(object, "shennong_qs2_bpcells", exact = TRUE)
  if (!is.null(manifest) &&
      (!identical(manifest$version, 1L) || !length(manifest$directories))) {
    stop("Unsupported or incomplete QS2 BPCells directory manifest.", call. = FALSE)
  }
  if (!is.null(bpcells_paths) &&
      (!is.character(bpcells_paths) || is.null(names(bpcells_paths)) ||
       anyNA(bpcells_paths) || anyNA(names(bpcells_paths)) || anyDuplicated(names(bpcells_paths)) ||
       any(!nzchar(bpcells_paths)) || any(!nzchar(names(bpcells_paths))))) {
    stop("`bpcells_paths` must map original absolute paths to replacement directories.", call. = FALSE)
  }
  resolved <- list()
  restore <- function(matrix) {
    source <- matrix@dir
    if (!is.null(resolved[[source]])) return(resolved[[source]])
    record <- manifest$directories[[source]]
    if (!is.null(manifest) && is.null(record)) stop("Missing BPCells identity for ", source, call. = FALSE)
    target <- source
    if (source %in% names(bpcells_paths)) {
      target <- bpcells_paths[[source]]
    } else if (length(record$relative)) {
      relative <- file.path(dirname(file), record$relative)
      if (dir.exists(relative)) target <- relative
    }
    if (!dir.exists(target)) {
      stop("Missing BPCells directory: ", source,
           ". Provide verified locations with `bpcells_paths`.", call. = FALSE)
    }
    target <- normalizePath(target, winslash = "/", mustWork = TRUE)
    if (!is.null(record)) {
      if (!identical(.sn_bpcells_fingerprint(target, matrix), record$sha256)) {
        stop("BPCells matrix identity differs from the QS2 snapshot: ", target, call. = FALSE)
      }
    } else {
      backing <- BPCells::open_matrix_dir(target)
      expected <- if (identical(backing@transpose, matrix@transpose)) dim(matrix) else rev(dim(matrix))
      if (!identical(dim(backing), expected) || !identical(backing@type, matrix@type) ||
          !identical(backing@compressed, matrix@compressed)) {
        stop("Legacy BPCells directory does not match the stored matrix layout: ", target, call. = FALSE)
      }
    }
    resolved[[source]] <<- target
    target
  }
  object <- .sn_map_bpcells_dirs(object, restore)
  if (length(setdiff(names(bpcells_paths), names(resolved)))) {
    stop("`bpcells_paths` names a source not present in the QS2 object.", call. = FALSE)
  }
  if (length(setdiff(names(manifest$directories), names(resolved)))) {
    stop("QS2 BPCells manifest contains directories absent from the object.", call. = FALSE)
  }
  if (is.null(manifest) && length(setdiff(names(resolved), names(bpcells_paths)))) {
    warning("Legacy QS2 snapshot has no BPCells content fingerprints. Verify the original directories before resaving with sn_write().", call. = FALSE)
  }
  if (inherits(object, "Seurat") && length(object@misc$bpcells_layers)) {
    object@misc$bpcells_layers <- lapply(object@misc$bpcells_layers, function(path) resolved[[path]] %||% path)
  }
  attr(object, "shennong_qs2_bpcells") <- NULL
  object
}
