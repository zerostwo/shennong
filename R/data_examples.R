# Internal registry for example datasets
.sn_example_data_catalog <- function() {
  data.frame(
    dataset = c("pbmc1k", "pbmc3k", "pbmc4k", "pbmc8k"),
    zenodo_record = c("14884845", "14884845", "14884845", "14884845"),
    species = c("human", "human", "human", "human"),
    stringsAsFactors = FALSE
  )
}

.sn_public_data_record_id <- function(record_id = NULL) {
  as.character(record_id %||% getOption("shennong.public_data_record", "20044788"))
}

.sn_public_data_cache_dir <- function(save_dir, record_id) {
  file.path(save_dir, paste0("zenodo_", record_id))
}

.sn_public_data_file_type <- function(matrix_type) {
  switch(
    matrix_type,
    filtered = "filtered_h5",
    raw = "raw_h5",
    metrics = "metrics_summary"
  )
}

.sn_read_public_data_index <- function(record_id = NULL,
                                       save_dir = "~/.shennong/data",
                                       token = NULL,
                                       overwrite = FALSE,
                                       quiet = TRUE) {
  record_id <- .sn_public_data_record_id(record_id)
  cache_dir <- .sn_public_data_cache_dir(sn_set_path(save_dir), record_id)
  index_path <- sn_download_zenodo(
    record_id = record_id,
    files = "shennong_index.json",
    save_dir = cache_dir,
    token = token,
    overwrite = overwrite,
    quiet = quiet
  )
  jsonlite::fromJSON(unname(index_path[["shennong_index.json"]]), simplifyVector = FALSE)
}

.sn_public_accessions <- function(accessions, name) {
  values <- accessions[[name]] %||% character(0)
  paste(as.character(values), collapse = ";")
}

.sn_public_data_catalog <- function(record_id = NULL,
                                    save_dir = "~/.shennong/data",
                                    token = NULL,
                                    overwrite = FALSE,
                                    quiet = TRUE) {
  index <- .sn_read_public_data_index(
    record_id = record_id,
    save_dir = save_dir,
    token = token,
    overwrite = overwrite,
    quiet = quiet
  )
  record_id <- as.character(index$zenodo$record_id %||% .sn_public_data_record_id(record_id))
  doi <- as.character(index$zenodo$doi %||% NA_character_)
  project_id <- as.character(index$project$project_id %||% NA_character_)
  rows <- unlist(lapply(index$studies %||% list(), function(study) {
    samples <- study$samples %||% list()
    lapply(samples, function(sample) {
      processing <- sample$processing %||% list()
      accessions <- sample$source_accessions %||% list()
      files <- sample$files %||% list()
      data.frame(
        source = "public",
        dataset = as.character(sample$sample_id %||% NA_character_),
        sample_id = as.character(sample$sample_id %||% NA_character_),
        sample_display_name = as.character(sample$display_name %||% NA_character_),
        study_id = as.character(study$study_id %||% NA_character_),
        study_display_name = as.character(study$display_name %||% NA_character_),
        organism = as.character(sample$organism %||% NA_character_),
        taxon_id = as.character(sample$taxon_id %||% NA_character_),
        assay = as.character(sample$assay %||% NA_character_),
        technology = as.character(sample$technology %||% NA_character_),
        pipeline = as.character(processing$pipeline %||% NA_character_),
        pipeline_version = as.character(processing$pipeline_version %||% NA_character_),
        reference = as.character(processing$reference %||% NA_character_),
        zenodo_record = record_id,
        doi = doi,
        project_id = project_id,
        zip_file = as.character(study$zip_file %||% NA_character_),
        filtered_h5 = as.character(files$filtered_h5 %||% NA_character_),
        raw_h5 = as.character(files$raw_h5 %||% NA_character_),
        metrics_summary = as.character(files$metrics_summary %||% NA_character_),
        geo_accession = .sn_public_accessions(accessions, "geo"),
        bioproject_accession = .sn_public_accessions(accessions, "bioproject"),
        arrayexpress_accession = .sn_public_accessions(accessions, "arrayexpress"),
        stringsAsFactors = FALSE
      )
    })
  }), recursive = FALSE)
  if (length(rows) == 0L) {
    return(data.frame())
  }
  do.call(rbind, rows)
}

.sn_example_catalog_for_listing <- function() {
  catalog <- .sn_example_data_catalog()
  data.frame(
    source = "example",
    dataset = catalog$dataset,
    sample_id = catalog$dataset,
    sample_display_name = catalog$dataset,
    study_id = "pbmc_examples",
    study_display_name = "PBMC example datasets",
    organism = ifelse(catalog$species == "human", "Homo sapiens", catalog$species),
    taxon_id = ifelse(catalog$species == "human", "9606", NA_character_),
    assay = "single_cell_transcriptomics",
    technology = "10x_genomics",
    pipeline = NA_character_,
    pipeline_version = NA_character_,
    reference = NA_character_,
    zenodo_record = catalog$zenodo_record,
    doi = "10.5281/zenodo.14884845",
    project_id = "pbmc_examples",
    zip_file = NA_character_,
    filtered_h5 = paste0(catalog$dataset, "_filtered_feature_bc_matrix.h5"),
    raw_h5 = paste0(catalog$dataset, "_raw_feature_bc_matrix.h5"),
    metrics_summary = NA_character_,
    geo_accession = NA_character_,
    bioproject_accession = NA_character_,
    arrayexpress_accession = NA_character_,
    stringsAsFactors = FALSE
  )
}

#' List datasets available through Shennong
#'
#' `sn_list_datasets()` returns a sample-level registry for datasets that can be
#' loaded with \code{\link{sn_load_data}}. By default it reads the current
#' Shennong public Zenodo collection, whose machine-readable index follows the
#' `UPLOAD_RULES.md` layout stored in the record.
#'
#' @param source Which registry to list. \code{"public"} lists the current
#'   Shennong public Zenodo collection; \code{"examples"} lists the legacy PBMC
#'   examples; \code{"all"} returns both.
#' @param record_id Zenodo record ID for the public Shennong collection.
#'   Defaults to option \code{shennong.public_data_record} or \code{"20044788"}.
#' @param save_dir Local cache directory used to store \code{shennong_index.json}.
#' @param token Optional Zenodo access token for restricted/private records.
#' @param overwrite Logical; if \code{TRUE}, re-download the public index.
#' @param quiet Logical; if \code{TRUE}, suppress Zenodo download progress.
#'
#' @return A data frame with one row per loadable sample-level dataset.
#'
#' @examples
#' \dontrun{
#' datasets <- sn_list_datasets()
#' head(datasets[, c("dataset", "study_id", "organism", "technology")])
#' }
#'
#' @export
sn_list_datasets <- function(source = c("public", "examples", "all"),
                             record_id = NULL,
                             save_dir = "~/.shennong/data",
                             token = NULL,
                             overwrite = FALSE,
                             quiet = TRUE) {
  source <- match.arg(source)
  public <- if (source %in% c("public", "all")) {
    .sn_public_data_catalog(
      record_id = record_id,
      save_dir = save_dir,
      token = token,
      overwrite = overwrite,
      quiet = quiet
    )
  } else {
    NULL
  }
  examples <- if (source %in% c("examples", "all")) .sn_example_catalog_for_listing() else NULL
  rows <- c(list(public), list(examples))
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0L) {
    return(data.frame())
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
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
#' @param dataset Character vector. Which dataset(s) to load. Legacy example
#'   values include \code{"pbmc1k"}, \code{"pbmc3k"}, \code{"pbmc4k"}, and
#'   \code{"pbmc8k"}. Public Shennong collection values are sample IDs returned
#'   by \code{sn_list_datasets()}, or study IDs when \code{sample_id} is also
#'   supplied. Default: \code{"pbmc3k"}.
#'
#' @param sample_id Optional sample ID(s) for public Shennong collection
#'   records. Use this when \code{dataset} names a study ID such as
#'   \code{"AndrewDHildreth2021"}. If \code{dataset} already names sample IDs,
#'   leave this as \code{NULL}.
#'
#' @param matrix_type Character scalar. Which matrix type to load.
#'   One of:
#'   \itemize{
#'     \item \code{"filtered"}: the \code{filtered_feature_bc_matrix.h5}
#'           (Cell Ranger filtered barcodes).
#'     \item \code{"raw"}: the \code{raw_feature_bc_matrix.h5}
#'           (unfiltered barcodes; useful for ambient RNA correction tools
#'           such as SoupX).
#'     \item \code{"metrics"}: the Cell Ranger \code{metrics_summary.csv}
#'           for public Shennong collection samples.
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
#' @param record_id Zenodo record ID for public Shennong collection samples.
#'   Defaults to option \code{shennong.public_data_record} or \code{"20044788"}.
#'
#' @param validate Logical; if \code{TRUE}, validate extracted public
#'   collection files against \code{manifest.tsv} MD5 checksums.
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
#'
#' # 6. Load a sample from the Shennong public Zenodo collection:
#' data_index <- sn_list_datasets()
#' sample_obj <- sn_load_data(dataset = data_index$dataset[[1]])
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
                         sample_id = NULL,
                         matrix_type = c("filtered", "raw", "metrics"),
                         save_dir = "~/.shennong/data",
                         return_object = TRUE,
                         species = NULL,
                         token = NULL,
                         overwrite = FALSE,
                         quiet = FALSE,
                         record_id = NULL,
                         validate = FALSE) {
  catalog <- .sn_example_data_catalog()
  matrix_type <- match.arg(matrix_type)
  dataset <- .sn_normalize_dataset_argument(dataset)
  if (is.null(sample_id) && all(dataset %in% catalog$dataset)) {
    return(.sn_load_example_data(
      dataset = dataset,
      matrix_type = matrix_type,
      save_dir = save_dir,
      return_object = return_object,
      species = species,
      token = token,
      overwrite = overwrite,
      quiet = quiet
    ))
  }
  if (is.null(sample_id) && any(dataset %in% catalog$dataset)) {
    stop("Legacy example datasets and public collection datasets cannot be mixed in one `sn_load_data()` call.", call. = FALSE)
  }

  .sn_load_public_data(
    dataset = dataset,
    sample_id = sample_id,
    matrix_type = matrix_type,
    save_dir = save_dir,
    return_object = return_object,
    species = species,
    token = token,
    overwrite = overwrite,
    quiet = quiet,
    record_id = record_id,
    validate = validate
  )
}

.sn_load_example_data <- function(dataset,
                                  matrix_type,
                                  save_dir,
                                  return_object,
                                  species,
                                  token,
                                  overwrite,
                                  quiet) {
  catalog <- .sn_example_data_catalog()
  dataset <- .sn_match_example_datasets(dataset = dataset, catalog = catalog)
  if (identical(matrix_type, "metrics")) {
    stop("`matrix_type = \"metrics\"` is available only for public Shennong collection samples.", call. = FALSE)
  }
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

.sn_normalize_dataset_argument <- function(dataset) {
  if (!is.character(dataset) || length(dataset) == 0L || any(!nzchar(dataset))) {
    stop("`dataset` must be a non-empty character vector.", call. = FALSE)
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

.sn_match_public_data_rows <- function(dataset, sample_id = NULL, catalog) {
  if (is.null(sample_id)) {
    sample_rows <- catalog[catalog$sample_id %in% dataset, , drop = FALSE]
    if (nrow(sample_rows) == length(dataset)) {
      return(sample_rows[match(dataset, sample_rows$sample_id), , drop = FALSE])
    }
    study_rows <- catalog[catalog$study_id %in% dataset, , drop = FALSE]
    if (nrow(study_rows) > 0L) {
      stop(
        "`sample_id` must be supplied when `dataset` names public studies. ",
        "Use `sn_list_datasets()` to inspect samples for study_id(s): ",
        paste(unique(study_rows$study_id), collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    invalid <- setdiff(dataset, catalog$sample_id)
    stop(
      "`dataset` must contain known sample IDs from `sn_list_datasets()` or legacy examples. Invalid value(s): ",
      paste(invalid, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  if (!is.character(sample_id) || length(sample_id) == 0L || any(!nzchar(sample_id))) {
    stop("`sample_id` must be `NULL` or a non-empty character vector.", call. = FALSE)
  }
  if (length(dataset) == 1L && length(sample_id) > 1L) {
    dataset <- rep(dataset, length(sample_id))
  }
  if (length(sample_id) == 1L && length(dataset) > 1L) {
    sample_id <- rep(sample_id, length(dataset))
  }
  if (length(dataset) != length(sample_id)) {
    stop("`sample_id` must have length 1 or the same length as `dataset`.", call. = FALSE)
  }
  key <- paste(dataset, sample_id, sep = "\r")
  row_key <- paste(catalog$study_id, catalog$sample_id, sep = "\r")
  matched <- match(key, row_key)
  if (any(is.na(matched))) {
    invalid <- key[is.na(matched)]
    stop(
      "Unknown public study/sample combination(s): ",
      paste(gsub("\r", "/", invalid, fixed = TRUE), collapse = ", "),
      ". Use `sn_list_datasets()` to inspect available values.",
      call. = FALSE
    )
  }
  catalog[matched, , drop = FALSE]
}

.sn_public_species <- function(organism) {
  organism <- tolower(as.character(organism))
  ifelse(
    organism %in% c("homo sapiens", "human"),
    "human",
    ifelse(organism %in% c("mus musculus", "mouse"), "mouse", organism)
  )
}

.sn_public_manifest <- function(record_id,
                                save_dir,
                                token,
                                overwrite,
                                quiet) {
  cache_dir <- .sn_public_data_cache_dir(save_dir, record_id)
  manifest_path <- sn_download_zenodo(
    record_id = record_id,
    files = "manifest.tsv",
    save_dir = cache_dir,
    token = token,
    overwrite = overwrite,
    quiet = quiet
  )
  utils::read.delim(unname(manifest_path[["manifest.tsv"]]), check.names = FALSE, stringsAsFactors = FALSE)
}

.sn_download_public_study_zip <- function(row,
                                          save_dir,
                                          token,
                                          overwrite,
                                          quiet) {
  cache_dir <- .sn_public_data_cache_dir(save_dir, row$zenodo_record[[1]])
  zip_path <- sn_download_zenodo(
    record_id = row$zenodo_record[[1]],
    files = row$zip_file[[1]],
    save_dir = cache_dir,
    token = token,
    overwrite = overwrite,
    quiet = quiet
  )
  unname(zip_path[[row$zip_file[[1]]]])
}

.sn_validate_public_file <- function(path, row, file_type, manifest) {
  if (is.null(manifest)) {
    return(invisible(TRUE))
  }
  manifest_row <- manifest[
    manifest$study_id == row$study_id[[1]] &
      manifest$sample_id == row$sample_id[[1]] &
      manifest$file_type == file_type,
    ,
    drop = FALSE
  ]
  if (nrow(manifest_row) != 1L) {
    stop("Could not find a unique manifest row for `", row$study_id[[1]], "/", row$sample_id[[1]], "`.", call. = FALSE)
  }
  md5 <- unname(tools::md5sum(path))
  if (!identical(tolower(md5), tolower(manifest_row$md5[[1]]))) {
    stop("Checksum validation failed for extracted file: ", path, call. = FALSE)
  }
  invisible(TRUE)
}

.sn_prepare_public_sample_file <- function(row,
                                           matrix_type,
                                           save_dir,
                                           token,
                                           overwrite,
                                           quiet,
                                           validate,
                                           manifest = NULL) {
  file_type <- .sn_public_data_file_type(matrix_type)
  internal_path <- row[[file_type]][[1]]
  if (is.na(internal_path) || !nzchar(internal_path)) {
    stop("The selected public sample does not provide file type `", file_type, "`.", call. = FALSE)
  }
  cache_dir <- .sn_public_data_cache_dir(save_dir, row$zenodo_record[[1]])
  local_path <- file.path(cache_dir, internal_path)
  if (!file.exists(local_path) || isTRUE(overwrite)) {
    zip_path <- .sn_download_public_study_zip(
      row = row,
      save_dir = save_dir,
      token = token,
      overwrite = overwrite,
      quiet = quiet
    )
    utils::unzip(zipfile = zip_path, files = internal_path, exdir = cache_dir, overwrite = TRUE)
  }
  if (!file.exists(local_path)) {
    stop("Failed to extract `", internal_path, "` from `", row$zip_file[[1]], "`.", call. = FALSE)
  }
  local_path <- normalizePath(local_path, winslash = "/", mustWork = TRUE)
  if (isTRUE(validate)) {
    .sn_validate_public_file(path = local_path, row = row, file_type = file_type, manifest = manifest)
  }
  local_path
}

.sn_load_public_data <- function(dataset,
                                 sample_id,
                                 matrix_type,
                                 save_dir,
                                 return_object,
                                 species,
                                 token,
                                 overwrite,
                                 quiet,
                                 record_id,
                                 validate) {
  save_dir <- sn_set_path(save_dir)
  catalog <- .sn_public_data_catalog(
    record_id = record_id,
    save_dir = save_dir,
    token = token,
    overwrite = overwrite,
    quiet = quiet
  )
  rows <- .sn_match_public_data_rows(dataset = dataset, sample_id = sample_id, catalog = catalog)
  if (any(duplicated(rows$sample_id))) {
    stop("Public sample IDs must be unique in one `sn_load_data()` call.", call. = FALSE)
  }
  species <- .sn_resolve_example_species(
    species = species,
    dataset = rows$sample_id,
    catalog = data.frame(dataset = rows$sample_id, species = .sn_public_species(rows$organism), stringsAsFactors = FALSE)
  )
  record_id <- unique(rows$zenodo_record)
  manifest <- if (isTRUE(validate)) {
    .sn_public_manifest(
      record_id = record_id[[1]],
      save_dir = save_dir,
      token = token,
      overwrite = overwrite,
      quiet = quiet
    )
  } else {
    NULL
  }
  local_files <- vapply(seq_len(nrow(rows)), function(i) {
    .sn_prepare_public_sample_file(
      row = rows[i, , drop = FALSE],
      matrix_type = matrix_type,
      save_dir = save_dir,
      token = token,
      overwrite = overwrite,
      quiet = quiet,
      validate = validate,
      manifest = manifest
    )
  }, character(1))
  names(local_files) <- rows$sample_id

  if (!return_object) {
    if (length(local_files) == 1L) {
      return(invisible(unname(local_files)))
    }
    return(invisible(local_files))
  }

  if (identical(matrix_type, "metrics")) {
    metrics <- lapply(local_files, utils::read.csv, check.names = FALSE)
    if (length(metrics) == 1L) {
      return(metrics[[1]])
    }
    return(metrics)
  }

  counts <- lapply(local_files, sn_read)
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
      sample_name = names(local_files)[[i]],
      project = rows$study_id[[i]]
    )
  })
  names(objects) <- names(local_files)
  result <- .sn_merge_example_objects(objects)
  if (inherits(result, "Seurat")) {
    Seurat::Misc(result, "loaded_datasets") <- list(
      datasets = rows$sample_id,
      studies = rows$study_id,
      loader = "sn_load_data",
      zenodo_record = rows$zenodo_record[[1]],
      merged = length(objects) > 1L
    )
  }
  result
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
