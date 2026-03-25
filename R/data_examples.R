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
#' @param dataset Character scalar. Which example dataset to load.
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
#'   the example-data registry.
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
#'         \code{sn_initialize_seurat_object()}.
#'   \item If \code{matrix_type == "raw"}:
#'         the function returns a sparse count matrix
#'         (typically a \code{dgCMatrix}), suitable for ambient RNA
#'         estimation / SoupX workflows.
#' }
#'
#' When \code{return_object = FALSE}, no object is constructed;
#' only the file is ensured to exist locally.
#'
#' @return
#' One of:
#'
#' \itemize{
#'   \item If \code{matrix_type == "filtered"} and
#'         \code{return_object = TRUE}:
#'         a Seurat object.
#'
#'   \item If \code{matrix_type == "raw"} and
#'         \code{return_object = TRUE}:
#'         a sparse count matrix (e.g. \code{dgCMatrix}).
#'
#'   \item If \code{return_object = FALSE}:
#'         the local file path to the cached \code{.h5} file,
#'         returned invisibly.
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
#' # 4. Use a custom cache directory:
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
#' \code{\link{sn_write}},
#' \code{\link{sn_read}}
#'
#' @export
sn_load_data <- function(dataset = "pbmc3k",
                         matrix_type = c("filtered", "raw"),
                         save_dir = "~/.shennong/data",
                         return_object = TRUE,
                         species = NULL) {
  catalog <- .sn_example_data_catalog()
  dataset <- arg_match(dataset, values = catalog$dataset)
  matrix_type <- match.arg(matrix_type)
  dataset_info <- catalog[catalog$dataset == dataset, , drop = FALSE]
  species <- species %||% dataset_info$species[[1]]
  zenodo_record <- dataset_info$zenodo_record[[1]]

  save_dir <- sn_set_path(save_dir)

  local_h5 <- glue("{save_dir}/{dataset}_{matrix_type}_feature_bc_matrix.h5")

  if (!file.exists(local_h5)) {
    remote_url <- glue::glue(
      "https://zenodo.org/records/{zenodo_record}/files/{dataset}_{matrix_type}_feature_bc_matrix.h5"
    )
    cli::cli_inform("Fetching {dataset} ({matrix_type}) from Zenodo cache...")
    curl::curl_download(url = remote_url, destfile = local_h5)
  }

  if (!return_object) {
    return(invisible(local_h5))
  }

  counts <- sn_read(local_h5)

  if (matrix_type == "filtered") {
    return(sn_initialize_seurat_object(
      x = counts,
      species = species
    ))
  }

  counts
}
