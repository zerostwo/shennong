#' Upload reusable data files to Zenodo
#'
#' `sn_upload_zenodo()` creates or updates a Zenodo draft through
#' \pkg{zen4R}, uploads user-selected files, and writes a Shennong manifest
#' alongside those files so downstream users can verify the exact dataset
#' version and file checksums they reused.
#'
#' @param files Character vector of local files to upload.
#' @param title Record title.
#' @param description Record description. When \code{NULL}, Shennong writes a
#'   short reusable-data description.
#' @param creators Creator names. A character vector such as
#'   \code{c("Duan, Songqi")} is the simplest form. A data frame or list with
#'   \code{name}, \code{firstname}, \code{lastname}, \code{orcid}, and
#'   \code{affiliation}/\code{affiliations} fields is also accepted.
#' @param version Dataset version string recorded in Zenodo metadata and the
#'   Shennong manifest.
#' @param record_id Optional Zenodo record ID. When supplied,
#'   \code{new_version = FALSE} updates that draft; \code{new_version = TRUE}
#'   creates a new version from the published record.
#' @param new_version If \code{TRUE}, create a new Zenodo version for
#'   \code{record_id}.
#' @param publish If \code{TRUE}, publish the record after upload. Defaults to
#'   \code{FALSE} so users can inspect the draft first.
#' @param sandbox If \code{TRUE}, use the Zenodo sandbox service.
#' @param token Zenodo access token. When \code{NULL}, Shennong checks
#'   \code{ZENODO_SANDBOX_TOKEN} for sandbox uploads and then
#'   \code{ZENODO_TOKEN} / \code{ZENODO_PAT}.
#' @param license License identifier stored in Zenodo rights metadata.
#' @param resource_type Zenodo resource type ID. Defaults to \code{"dataset"}.
#' @param access File and record access policy. Defaults to \code{"public"}.
#' @param keywords Optional subject keywords.
#' @param manifest If \code{TRUE}, upload a
#'   \code{shennong_zenodo_manifest.json} file with checksums and provenance.
#' @param manifest_dir Directory used to write the temporary manifest file.
#' @param reserve_doi Whether to reserve a DOI for new drafts.
#' @param dry_run If \code{TRUE}, build and return the planned record and
#'   manifest without contacting Zenodo.
#' @param logger Logger level passed to \code{zen4R::ZenodoManager$new()}.
#'
#' @return A list containing the Zenodo identifiers, uploaded-file manifest,
#'   manifest path, publication status, and the \code{ZenodoRecord} object.
#'
#' @examples
#' \dontrun{
#' result <- sn_upload_zenodo(
#'   files = c("data/reference_atlas.qs", "data/reference_markers.csv"),
#'   title = "Reusable PBMC reference atlas",
#'   creators = "Duan, Songqi",
#'   version = "2026.05.04",
#'   sandbox = TRUE
#' )
#' result$doi
#' result$manifest
#' }
#'
#' @export
sn_upload_zenodo <- function(files,
                             title,
                             description = NULL,
                             creators = NULL,
                             version = NULL,
                             record_id = NULL,
                             new_version = FALSE,
                             publish = FALSE,
                             sandbox = FALSE,
                             token = NULL,
                             license = "cc-by-4.0",
                             resource_type = "dataset",
                             access = c("public", "restricted"),
                             keywords = NULL,
                             manifest = TRUE,
                             manifest_dir = tempdir(),
                             reserve_doi = TRUE,
                             dry_run = FALSE,
                             logger = NULL) {
  check_installed("zen4R", reason = "to upload reusable data to Zenodo.")
  access <- rlang::arg_match(access)

  if (missing(files) || length(files) == 0L) {
    stop("`files` must contain at least one local file.", call. = FALSE)
  }
  if (missing(title) || !is.character(title) || length(title) != 1L || !nzchar(title)) {
    stop("`title` must be a non-empty string.", call. = FALSE)
  }
  if (!is.logical(new_version) || length(new_version) != 1L || is.na(new_version)) {
    stop("`new_version` must be TRUE or FALSE.", call. = FALSE)
  }
  if (isTRUE(new_version) && (is.null(record_id) || !nzchar(as.character(record_id)))) {
    stop("`record_id` is required when `new_version = TRUE`.", call. = FALSE)
  }

  files <- .sn_prepare_zenodo_files(files)
  record <- .sn_build_zenodo_record(
    title = title,
    description = description,
    creators = creators,
    version = version,
    license = license,
    resource_type = resource_type,
    access = access,
    keywords = keywords
  )

  file_manifest <- .sn_zenodo_file_manifest(files)
  manifest_path <- NULL
  upload_files <- files
  if (isTRUE(manifest)) {
    manifest_path <- .sn_write_zenodo_manifest(
      file_manifest = file_manifest,
      title = title,
      version = version,
      sandbox = sandbox,
      record_id = record_id,
      manifest_dir = manifest_dir
    )
    upload_files <- c(upload_files, manifest_path)
  }

  if (isTRUE(dry_run)) {
    return(.sn_zenodo_upload_summary(
      record = record,
      files = file_manifest,
      manifest_path = manifest_path,
      sandbox = sandbox,
      published = FALSE,
      dry_run = TRUE
    ))
  }

  manager <- .sn_zenodo_manager(token = token, sandbox = sandbox, logger = logger)

  if (isTRUE(new_version)) {
    base_record <- manager$getRecordById(as.character(record_id))
    if (inherits(base_record, "ZenodoException") || is.null(base_record)) {
      stop("Could not fetch published Zenodo record `", record_id, "` for versioning.", call. = FALSE)
    }
    base_record$metadata <- record$metadata
    base_record$access <- record$access
    deposited <- manager$depositRecordVersion(
      record = base_record,
      delete_latest_files = TRUE,
      files = upload_files,
      publish = publish
    )
  } else {
    if (!is.null(record_id)) {
      record$id <- as.character(record_id)
    }
    deposited <- manager$depositRecord(record = record, reserveDOI = reserve_doi, publish = FALSE)
    if (inherits(deposited, "ZenodoException") || is.null(deposited)) {
      stop("Zenodo record deposition failed.", call. = FALSE)
    }
    for (path in upload_files) {
      manager$uploadFile(path, record = deposited)
    }
    if (isTRUE(publish)) {
      deposited <- manager$publishRecord(deposited$id)
    }
  }

  if (inherits(deposited, "ZenodoException") || is.null(deposited)) {
    stop("Zenodo upload failed.", call. = FALSE)
  }

  .sn_zenodo_upload_summary(
    record = deposited,
    files = file_manifest,
    manifest_path = manifest_path,
    sandbox = sandbox,
    published = isTRUE(publish),
    dry_run = FALSE
  )
}

#' Download files from a Zenodo record
#'
#' `sn_download_zenodo()` downloads one or more files from a Zenodo record into
#' a local cache directory. Public Zenodo records do not require a token. Supply
#' `token` only for restricted files or private records that your account can
#' access.
#'
#' @param record_id Zenodo record ID.
#' @param files Character vector of file names to download. When `NULL`, the
#'   function queries the Zenodo record metadata and downloads all listed files.
#' @param save_dir Local directory where files are cached.
#' @param token Optional Zenodo access token for restricted/private downloads.
#' @param sandbox If `TRUE`, use the Zenodo sandbox service.
#' @param overwrite If `TRUE`, re-download files that already exist locally.
#' @param quiet If `TRUE`, suppress download progress messages.
#'
#' @return A named character vector of local file paths.
#'
#' @examples
#' \dontrun{
#' paths <- sn_download_zenodo(
#'   record_id = "14884845",
#'   files = "pbmc3k_filtered_feature_bc_matrix.h5"
#' )
#' paths
#' }
#'
#' @export
sn_download_zenodo <- function(record_id,
                               files = NULL,
                               save_dir = "~/.shennong/data",
                               token = NULL,
                               sandbox = FALSE,
                               overwrite = FALSE,
                               quiet = FALSE) {
  if (missing(record_id) || length(record_id) != 1L || !nzchar(as.character(record_id))) {
    stop("`record_id` must be a non-empty Zenodo record ID.", call. = FALSE)
  }
  if (!is.null(files) && (!is.character(files) || length(files) == 0L || any(!nzchar(files)))) {
    stop("`files` must be `NULL` or a non-empty character vector.", call. = FALSE)
  }
  if (!is.logical(overwrite) || length(overwrite) != 1L || is.na(overwrite)) {
    stop("`overwrite` must be TRUE or FALSE.", call. = FALSE)
  }

  record_id <- as.character(record_id)
  save_dir <- sn_set_path(save_dir)
  token <- if (is.null(token)) NULL else as.character(token[[1]])
  if (!is.null(token) && !nzchar(token)) {
    token <- NULL
  }

  if (is.null(files)) {
    files <- .sn_list_zenodo_record_files(record_id = record_id, token = token, sandbox = sandbox)
  }

  paths <- vapply(files, function(file_name) {
    destfile <- file.path(save_dir, basename(file_name))
    if (file.exists(destfile) && !isTRUE(overwrite)) {
      return(normalizePath(destfile, winslash = "/", mustWork = TRUE))
    }

    url <- .sn_zenodo_file_download_url(
      record_id = record_id,
      file_name = file_name,
      sandbox = sandbox
    )
    .sn_download_zenodo_file(
      url = url,
      destfile = destfile,
      token = token,
      quiet = quiet
    )
  }, character(1), USE.NAMES = TRUE)

  names(paths) <- basename(files)
  paths
}

.sn_prepare_zenodo_files <- function(files) {
  files <- path.expand(as.character(files))
  if (any(!nzchar(files))) {
    stop("`files` cannot contain empty paths.", call. = FALSE)
  }
  missing_files <- files[!file.exists(files)]
  if (length(missing_files) > 0L) {
    stop("Zenodo upload file(s) not found: ", paste(missing_files, collapse = ", "), call. = FALSE)
  }
  directories <- files[dir.exists(files)]
  if (length(directories) > 0L) {
    stop("Zenodo uploads require files, not directories: ", paste(directories, collapse = ", "), call. = FALSE)
  }
  normalizePath(files, winslash = "/", mustWork = TRUE)
}

.sn_build_zenodo_record <- function(title,
                                    description = NULL,
                                    creators = NULL,
                                    version = NULL,
                                    license = "cc-by-4.0",
                                    resource_type = "dataset",
                                    access = "public",
                                    keywords = NULL) {
  record <- zen4R::ZenodoRecord$new(logger = NULL)
  record$metadata$title <- title
  record$metadata$description <- description %||% paste(
    "Reusable data uploaded with Shennong",
    as.character(utils::packageVersion("Shennong"))
  )
  record$metadata$resource_type <- list(id = resource_type)
  record$metadata$rights <- list(list(id = license))
  record$metadata$creators <- .sn_normalize_zenodo_creators(creators)
  if (!is.null(version)) {
    record$metadata$version <- as.character(version)
  }
  if (!is.null(keywords)) {
    record$metadata$subjects <- lapply(as.character(keywords), function(keyword) {
      list(subject = keyword)
    })
  }
  record$access$record <- access
  record$access$files <- access
  record
}

.sn_normalize_zenodo_creators <- function(creators = NULL) {
  if (is.null(creators)) {
    creators <- "Duan, Songqi"
  }
  if (is.data.frame(creators)) {
    creators <- lapply(seq_len(nrow(creators)), function(i) as.list(creators[i, , drop = FALSE]))
  }
  if (is.character(creators)) {
    creators <- as.list(creators)
  }
  if (!is.list(creators) || length(creators) == 0L) {
    stop("`creators` must be a character vector, data frame, or list.", call. = FALSE)
  }
  lapply(creators, .sn_normalize_zenodo_creator)
}

.sn_normalize_zenodo_creator <- function(creator) {
  if (is.character(creator) && length(creator) == 1L) {
    creator <- list(name = creator)
  }
  if (!is.list(creator)) {
    stop("Each creator must be a string or a named list/data-frame row.", call. = FALSE)
  }
  creator <- lapply(creator, function(x) {
    if (length(x) == 0L || all(is.na(x))) NULL else as.character(x[[1]])
  })
  name <- creator$name %||% paste(stats::na.omit(c(creator$lastname, creator$firstname)), collapse = ", ")
  if (!nzchar(name)) {
    stop("Each Zenodo creator must include `name` or `firstname`/`lastname`.", call. = FALSE)
  }
  person <- list(
    type = if (!is.null(creator$firstname) || !is.null(creator$lastname)) "personal" else "organizational",
    name = name,
    identifiers = list()
  )
  if (!is.null(creator$firstname)) {
    person$given_name <- creator$firstname
  }
  if (!is.null(creator$lastname)) {
    person$family_name <- creator$lastname
  }
  if (!is.null(creator$orcid)) {
    person$identifiers[[length(person$identifiers) + 1L]] <- list(scheme = "orcid", identifier = creator$orcid)
  }
  affiliations <- creator$affiliations %||% creator$affiliation %||% NULL
  affiliations <- if (is.null(affiliations)) {
    list()
  } else {
    lapply(strsplit(affiliations, "\\s*;\\s*")[[1]], function(affiliation) list(name = affiliation))
  }
  list(person_or_org = person, affiliations = affiliations)
}

.sn_zenodo_file_manifest <- function(files) {
  info <- file.info(files)
  tibble::tibble(
    file = basename(files),
    path = files,
    size = as.numeric(info$size),
    md5 = unname(tools::md5sum(files)),
    sha256 = unname(tools::sha256sum(files))
  )
}

.sn_write_zenodo_manifest <- function(file_manifest,
                                      title,
                                      version = NULL,
                                      sandbox = FALSE,
                                      record_id = NULL,
                                      manifest_dir = tempdir()) {
  dir.create(manifest_dir, recursive = TRUE, showWarnings = FALSE)
  manifest_path <- file.path(manifest_dir, "shennong_zenodo_manifest.json")
  payload <- list(
    schema_version = "1.0.0",
    title = title,
    dataset_version = version %||% NA_character_,
    package = "Shennong",
    package_version = as.character(utils::packageVersion("Shennong")),
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    sandbox = isTRUE(sandbox),
    record_id = record_id %||% NA_character_,
    files = file_manifest
  )
  jsonlite::write_json(payload, path = manifest_path, auto_unbox = TRUE, pretty = TRUE, dataframe = "rows")
  normalizePath(manifest_path, winslash = "/", mustWork = TRUE)
}

.sn_zenodo_base_url <- function(sandbox = FALSE) {
  if (isTRUE(sandbox)) "https://sandbox.zenodo.org" else "https://zenodo.org"
}

.sn_zenodo_api_url <- function(record_id, sandbox = FALSE) {
  paste0(.sn_zenodo_base_url(sandbox), "/api/records/", record_id)
}

.sn_zenodo_file_download_url <- function(record_id, file_name, sandbox = FALSE) {
  paste0(
    .sn_zenodo_base_url(sandbox),
    "/records/",
    record_id,
    "/files/",
    utils::URLencode(file_name, reserved = TRUE),
    "?download=1"
  )
}

.sn_zenodo_handle <- function(token = NULL) {
  handle <- curl::new_handle()
  if (!is.null(token) && nzchar(token)) {
    curl::handle_setheaders(handle, Authorization = paste("Bearer", token))
  }
  handle
}

.sn_download_zenodo_file <- function(url, destfile, token = NULL, quiet = FALSE) {
  dir.create(dirname(destfile), recursive = TRUE, showWarnings = FALSE)
  if (!isTRUE(quiet)) {
    cli::cli_inform("Fetching {basename(destfile)} from Zenodo...")
  }
  curl::curl_download(
    url = url,
    destfile = destfile,
    quiet = quiet,
    mode = "wb",
    handle = .sn_zenodo_handle(token)
  )
  normalizePath(destfile, winslash = "/", mustWork = TRUE)
}

.sn_list_zenodo_record_files <- function(record_id, token = NULL, sandbox = FALSE) {
  response <- curl::curl_fetch_memory(
    url = .sn_zenodo_api_url(record_id = record_id, sandbox = sandbox),
    handle = .sn_zenodo_handle(token)
  )
  if (response$status_code >= 400L) {
    stop("Could not read Zenodo record `", record_id, "` (HTTP ", response$status_code, ").", call. = FALSE)
  }
  record <- jsonlite::fromJSON(rawToChar(response$content), simplifyVector = FALSE)
  files <- record$files %||% list()
  file_names <- vapply(files, function(file) {
    file$key %||% file$filename %||% file$id %||% NA_character_
  }, character(1))
  file_names <- file_names[!is.na(file_names) & nzchar(file_names)]
  if (length(file_names) == 0L) {
    stop("Zenodo record `", record_id, "` does not list downloadable files.", call. = FALSE)
  }
  file_names
}

.sn_zenodo_token <- function(token = NULL, sandbox = FALSE) {
  candidates <- c(
    token,
    if (isTRUE(sandbox)) Sys.getenv("ZENODO_SANDBOX_TOKEN", unset = NA_character_) else NA_character_,
    Sys.getenv("ZENODO_TOKEN", unset = NA_character_),
    Sys.getenv("ZENODO_PAT", unset = NA_character_)
  )
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  if (length(candidates) == 0L) {
    return(NULL)
  }
  candidates[[1]]
}

.sn_zenodo_manager <- function(token = NULL, sandbox = FALSE, logger = NULL) {
  token <- .sn_zenodo_token(token = token, sandbox = sandbox)
  if (is.null(token)) {
    token_env <- if (isTRUE(sandbox)) "`ZENODO_SANDBOX_TOKEN` or `ZENODO_TOKEN`" else "`ZENODO_TOKEN`"
    stop("A Zenodo token is required. Supply `token` or set ", token_env, ".", call. = FALSE)
  }
  zen4R::ZenodoManager$new(token = token, sandbox = sandbox, logger = logger, keyring_backend = "env")
}

.sn_zenodo_upload_summary <- function(record,
                                      files,
                                      manifest_path = NULL,
                                      sandbox = FALSE,
                                      published = FALSE,
                                      dry_run = FALSE) {
  pids <- record$pids %||% list()
  doi <- pids$doi$identifier %||% record$metadata$doi %||% NA_character_
  concept_doi <- pids$conceptdoi$identifier %||% NA_character_
  links <- record$links %||% list()
  list(
    record_id = record$id %||% NA_character_,
    concept_id = record$conceptid %||% record$parent$id %||% NA_character_,
    doi = doi,
    concept_doi = concept_doi,
    url = links$record_html %||% links$self_html %||% links$latest_html %||% NA_character_,
    version = record$metadata$version %||% NA_character_,
    sandbox = isTRUE(sandbox),
    published = isTRUE(published),
    dry_run = isTRUE(dry_run),
    files = files,
    manifest_path = manifest_path %||% NA_character_,
    record = record
  )
}
