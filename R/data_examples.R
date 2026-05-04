# Internal registry for example datasets
.sn_example_data_catalog <- function() {
  data.frame(
    dataset = c("pbmc1k", "pbmc3k", "pbmc4k", "pbmc8k"),
    zenodo_record = c("14884845", "14884845", "14884845", "14884845"),
    species = c("human", "human", "human", "human"),
    stringsAsFactors = FALSE
  )
}

#' Load example datasets from Zenodo
#'
#' This function downloads (if not already cached) and loads processed
#' single-cell RNA-seq example datasets from a Zenodo-backed registry.
#' The current registry contains PBMC example datasets (`pbmc1k`, `pbmc3k`,
#' `pbmc4k`, and `pbmc8k`), and the interface is intentionally generalized so
#' additional example datasets can be added later without introducing another
#' top-level loader function.
#'
#' Instead of storing serialized R objects, the function caches the original
#' Zenodo \code{.h5} files locally for reproducibility and version consistency.
#' When requested, it dynamically constructs and returns either:
#' \itemize{
#'   \item a Seurat object (for filtered data), or
#'   \item a sparse count matrix (for raw data).
#' }
#' Users can also choose to only download/cache the data without loading it
#' into memory.
#'
#' @param dataset Character vector. Which example dataset(s) to load.
#'   Currently one of \code{"pbmc1k"}, \code{"pbmc3k"}, \code{"pbmc4k"},
#'   or \code{"pbmc8k"}. Default: \code{"pbmc3k"}.
#'
#' @param matrix_type Character scalar. Which matrix type to load.
#'   One of:
#'   \itemize{
#'     \item \code{"filtered"}: the \code{filtered_feature_bc_matrix.h5}
#'           (Cell Ranger filtered barcodes).
#'     \item \code{"raw"}: the \code{raw_feature_bc_matrix.h5}
#'           (unfiltered barcodes; useful for ambient RNA correction tools
#'           such as SoupX).
#'   }
#'   Default: \code{"filtered"}.
#'
#' @param save_dir Character scalar. Local cache directory where the downloaded
#'   Zenodo files will be stored. The directory will be created if it does not
#'   exist. Default: \code{"~/.shennong/data"}.
#'
#' @param return_object Logical. If \code{TRUE} (default), the function returns
#'   an in-memory object. If \code{FALSE}, the function only ensures that the
#'   file is downloaded locally and then returns (invisibly) the local file
#'   path.
#'
#' @param species Character scalar. Species label passed to
#'   \code{sn_initialize_seurat_object()} when constructing Seurat objects
#'   (i.e. when \code{matrix_type == "filtered"}). Ignored if
#'   \code{matrix_type == "raw"}. If \code{NULL}, use the dataset default from
#'   the example-data registry. When loading multiple datasets, supply either
#'   one species label for all datasets or one label per dataset.
#'
#' @param token Optional Zenodo access token. Public example datasets do not
#'   require a token. Supply this only for restricted/private records.
#'
#' @param overwrite Logical. If \code{TRUE}, re-download cached files.
#'
#' @param quiet Logical. If \code{TRUE}, suppress Zenodo download progress
#'   messages.
#'
#' @details
#' All datasets were re-aligned using Cell Ranger v9.0.1 with a custom reference
#' based on GENCODE v48 (GRCh38.p14).
#'
#' \strong{Caching strategy}
#'
#' The function caches the original files from Zenodo on first use:
#'
#' \preformatted{
#' {dataset}_{matrix_type}_feature_bc_matrix.h5
#' }
#'
#' For example, after loading \code{pbmc3k} you may see:
#'
#' \preformatted{
#' ~/.shennong/data/
#' |- pbmc1k_filtered_feature_bc_matrix.h5
#' |- pbmc1k_raw_feature_bc_matrix.h5
#' |- pbmc3k_filtered_feature_bc_matrix.h5
#' |- pbmc3k_raw_feature_bc_matrix.h5
#' |- pbmc4k_filtered_feature_bc_matrix.h5
#' |- pbmc4k_raw_feature_bc_matrix.h5
#' |- pbmc8k_filtered_feature_bc_matrix.h5
#' \\- pbmc8k_raw_feature_bc_matrix.h5
#' }
#'
#' When \code{return_object = TRUE}, the cached \code{.h5} file is read on the
#' fly:
#' \itemize{
#'   \item If \code{matrix_type == "filtered"}:
#'         a Seurat object is constructed via
#'         \code{sn_initialize_seurat_object()}. When multiple datasets are
#'         requested, the per-dataset Seurat objects are merged and each source
#'         is recorded in the \code{sample} metadata column.
#'   \item If \code{matrix_type == "raw"}:
#'         the function returns a sparse count matrix
#'         (typically a \code{dgCMatrix}), suitable for ambient RNA
#'         estimation / SoupX workflows. When multiple raw datasets are
#'         requested, the matrices are returned as a named list.
#' }
#'
#' When \code{return_object = FALSE}, no object is constructed;
#' only the file is ensured to exist locally. Multiple datasets return a named
#' character vector of cached paths.
#'
#' @return
#' One of:
#'
#' \itemize{
#'   \item If \code{matrix_type == "filtered"} and
#'         \code{return_object = TRUE}:
#'         a Seurat object. Multiple datasets are returned as one merged Seurat
#'         object.
#'
#'   \item If \code{matrix_type == "raw"} and
#'         \code{return_object = TRUE}:
#'         a sparse count matrix (e.g. \code{dgCMatrix}) or, for multiple
#'         datasets, a named list of sparse count matrices.
#'
#'   \item If \code{return_object = FALSE}:
#'         the local file path to the cached \code{.h5} file, or a named
#'         character vector of paths for multiple datasets, returned invisibly.
#' }
#'
#' @examples
#' \dontrun{
#' # 1. Load filtered PBMC3k as a Seurat object:
#' pbmc <- sn_load_data()
#'
#' # 2. Load raw PBMC3k counts as a sparse matrix (for SoupX etc.):
#' pbmc_raw <- sn_load_data(matrix_type = "raw")
#'
#' # 3. Only download/cache PBMC8k, don't construct anything in-memory:
#' sn_load_data(dataset = "pbmc8k", return_object = FALSE)
#'
#' # 4. Load and merge multiple filtered PBMC examples:
#' pbmc_merged <- sn_load_data(dataset = c("pbmc1k", "pbmc3k"))
#'
#' # 5. Use a custom cache directory:
#' pbmc4k <- sn_load_data(
#'   dataset  = "pbmc4k",
#'   save_dir = "~/datasets/pbmc_cache"
#' )
#' }
#'
#' @references
#' Duan S. (2025).
#' Processed PBMC datasets re-aligned with Cell Ranger v9.0.1
#' (GENCODE v48, GRCh38.p14). Zenodo.
#' DOI: 10.5281/zenodo.14884845
#'
#' @seealso
#' \code{\link{sn_initialize_seurat_object}},
#' \code{\link{sn_download_zenodo}},
#' \code{\link{sn_write}},
#' \code{\link{sn_read}}
#'
#' @export
sn_load_data <- function(dataset = "pbmc3k",
                         matrix_type = c("filtered", "raw"),
                         save_dir = "~/.shennong/data",
                         return_object = TRUE,
                         species = NULL,
                         token = NULL,
                         overwrite = FALSE,
                         quiet = FALSE) {
  catalog <- .sn_example_data_catalog()
  dataset <- .sn_match_example_datasets(dataset = dataset, catalog = catalog)
  matrix_type <- match.arg(matrix_type)
  species <- .sn_resolve_example_species(species = species, dataset = dataset, catalog = catalog)
  save_dir <- sn_set_path(save_dir)

  local_h5 <- vapply(dataset, function(dataset_name) {
    dataset_info <- catalog[catalog$dataset == dataset_name, , drop = FALSE]
    file_name <- paste0(dataset_name, "_", matrix_type, "_feature_bc_matrix.h5")
    paths <- sn_download_zenodo(
      record_id = dataset_info$zenodo_record[[1]],
      files = file_name,
      save_dir = save_dir,
      token = token,
      overwrite = overwrite,
      quiet = quiet
    )
    unname(paths[[file_name]])
  }, character(1), USE.NAMES = TRUE)

  if (!return_object) {
    if (length(local_h5) == 1L) {
      return(invisible(unname(local_h5)))
    }
    return(invisible(local_h5))
  }

  counts <- lapply(local_h5, sn_read)

  if (matrix_type == "raw") {
    if (length(counts) == 1L) {
      return(counts[[1]])
    }
    return(counts)
  }

  objects <- lapply(seq_along(counts), function(i) {
    sn_initialize_seurat_object(
      x = counts[[i]],
      species = species[[i]],
      sample_name = names(local_h5)[[i]],
      project = names(local_h5)[[i]]
    )
  })
  names(objects) <- names(local_h5)

  .sn_merge_example_objects(objects)
}

.sn_match_example_datasets <- function(dataset, catalog) {
  if (!is.character(dataset) || length(dataset) == 0L || any(!nzchar(dataset))) {
    stop("`dataset` must be a non-empty character vector.", call. = FALSE)
  }
  invalid <- setdiff(dataset, catalog$dataset)
  if (length(invalid) > 0L) {
    stop(
      "`dataset` must contain only known example datasets: ",
      paste(catalog$dataset, collapse = ", "),
      ". Invalid value(s): ",
      paste(invalid, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  duplicated_dataset <- unique(dataset[duplicated(dataset)])
  if (length(duplicated_dataset) > 0L) {
    stop(
      "`dataset` cannot contain duplicated values: ",
      paste(duplicated_dataset, collapse = ", "),
      ".",
      call. = FALSE
    )
  }
  dataset
}

.sn_resolve_example_species <- function(species = NULL, dataset, catalog) {
  dataset_species <- catalog$species[match(dataset, catalog$dataset)]
  if (is.null(species)) {
    return(dataset_species)
  }
  if (!is.character(species) || length(species) == 0L || any(!nzchar(species))) {
    stop("`species` must be `NULL`, a character scalar, or a character vector parallel to `dataset`.", call. = FALSE)
  }
  if (length(species) == 1L) {
    return(rep(species, length(dataset)))
  }
  if (length(species) != length(dataset)) {
    stop("`species` must have length 1 or the same length as `dataset`.", call. = FALSE)
  }
  species
}

.sn_merge_example_objects <- function(objects) {
  if (length(objects) == 1L) {
    return(objects[[1]])
  }
  merged <- merge(
    x = objects[[1]],
    y = objects[-1],
    add.cell.ids = names(objects),
    project = "Shennong"
  )
  Seurat::Misc(merged, "loaded_datasets") <- list(
    datasets = names(objects),
    loader = "sn_load_data",
    merged = TRUE
  )
  merged
}
