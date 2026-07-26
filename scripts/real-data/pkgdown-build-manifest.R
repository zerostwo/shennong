# Helpers shared by the real-data pkgdown builder and publisher.
#
# The manifest lives next to the ignored local runtime reports, not in the
# published site. This lets it cover every site file without a self-hash
# exception and keeps local build provenance out of gh-pages.

.pkgdown_manifest_filename <- "pkgdown-real-build-manifest.json"

.pkgdown_sha256 <- function(paths) {
  if (!length(paths)) return(character())
  paths <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  hashes <- unname(tools::sha256sum(paths))
  if (length(hashes) != length(paths) || anyNA(hashes) || any(!nzchar(hashes))) {
    stop("Unable to calculate SHA-256 for all build inputs.", call. = FALSE)
  }
  hashes
}

.pkgdown_relative_paths <- function(paths, root) {
  root <- sub("/+$", "", normalizePath(root, winslash = "/", mustWork = TRUE))
  paths <- normalizePath(paths, winslash = "/", mustWork = TRUE)
  prefix <- paste0(root, "/")
  if (any(!startsWith(paths, prefix))) {
    stop("A manifest file is outside its declared root.", call. = FALSE)
  }
  substring(paths, nchar(prefix) + 1L)
}

.pkgdown_site_files <- function(site_dir) {
  site_dir <- normalizePath(site_dir, winslash = "/", mustWork = TRUE)
  files <- list.files(
    site_dir,
    recursive = TRUE,
    all.files = TRUE,
    full.names = TRUE,
    include.dirs = FALSE,
    no.. = TRUE
  )
  files <- files[file.exists(files) & !dir.exists(files)]
  files[order(.pkgdown_relative_paths(files, site_dir), method = "radix")]
}

.pkgdown_forbidden_data_pattern <- paste0(
  "[.]",
  "(log|out|err|qs|qs2|rds|rda|rdata|rds2|h5|hdf5|h5ad|h5mu|loom|",
  "mtx|csv|tsv|tab|parquet|feather|arrow|xls|xlsx|ods|sqlite|sqlite3|",
  "db|duckdb|sql|sas7bdat|dta|sav|pkl|pickle|npy|npz|mat|",
  "fa|fna|fasta|fq|fastq|sam|bam|cram|bed|bedgraph|bigwig|bw|",
  "vcf|bcf|gtf|gff|gff3|zarr|cel|idat|mzml|mzxml|",
  "zip|tar|tgz|gz|bgz|bz2|xz|7z|rar|zst)$"
)

.pkgdown_allowed_static_extensions <- c(
  "html", "htm", "css", "js", "map", "md", "svg", "png", "jpg",
  "jpeg", "gif", "webp", "avif", "ico", "pdf", "woff", "woff2",
  "ttf", "otf", "eot", "webmanifest", "json", "yml", "yaml", "xml",
  "txt"
)

.pkgdown_assert_static_allowlist <- function(site_dir) {
  files <- .pkgdown_site_files(site_dir)
  if (!length(files)) stop("Rendered site contains no static files.", call. = FALSE)

  relative <- .pkgdown_relative_paths(files, site_dir)
  base <- basename(relative)
  forbidden <- grepl(
    .pkgdown_forbidden_data_pattern,
    base,
    ignore.case = TRUE,
    perl = TRUE
  )
  if (any(forbidden)) {
    stop(
      "Rendered site contains forbidden log/raw/data/table/database/archive files: ",
      paste(relative[forbidden], collapse = ", "),
      call. = FALSE
    )
  }

  allowed_extensionless <- base %in% c("CNAME", ".nojekyll")
  has_extension <- grepl("[.][^.]+$", base, perl = TRUE)
  suspicious_extensionless <- !has_extension & !allowed_extensionless
  if (any(suspicious_extensionless)) {
    stop(
      "Rendered site contains suspicious extensionless files: ",
      paste(relative[suspicious_extensionless], collapse = ", "),
      call. = FALSE
    )
  }

  extension <- tolower(sub("^.*[.]", "", base))
  allowed <- allowed_extensionless |
    (has_extension & extension %in% .pkgdown_allowed_static_extensions)
  if (any(!allowed)) {
    stop(
      "Rendered site contains files outside the conservative static-site allowlist: ",
      paste(relative[!allowed], collapse = ", "),
      call. = FALSE
    )
  }

  # JSON/YAML/XML/TXT are useful pkgdown metadata, but arbitrary files with
  # these extensions are also common data-leak vectors. Keep only known static
  # metadata names produced by pkgdown or explicitly intended for crawlers/LLMs.
  structured <- extension %in% c("json", "yml", "yaml", "xml", "txt")
  known_structured <- relative %in% c(
    "search.json", "pkgdown.yml", "sitemap.xml", "robots.txt", "llms.txt",
    "deps/data-deps.txt"
  )
  if (any(structured & !known_structured)) {
    stop(
      "Rendered site contains unrecognized structured/text artifacts: ",
      paste(relative[structured & !known_structured], collapse = ", "),
      call. = FALSE
    )
  }

  invisible(relative)
}

.pkgdown_remove_generated_logs <- function(site_dir) {
  if (!dir.exists(site_dir)) return(character())
  files <- .pkgdown_site_files(site_dir)
  logs <- files[grepl("[.]log$", basename(files), ignore.case = TRUE, perl = TRUE)]
  if (!length(logs)) return(character())
  relative <- .pkgdown_relative_paths(logs, site_dir)
  unlink(logs, force = TRUE)
  remaining <- logs[file.exists(logs)]
  if (length(remaining)) {
    stop(
      "Unable to remove generated log files from the pkgdown site: ",
      paste(remaining, collapse = ", "),
      call. = FALSE
    )
  }
  relative
}

.pkgdown_expected_reference_pages <- function(repo_root) {
  man_dir <- file.path(repo_root, "man")
  rd_files <- list.files(
    man_dir,
    pattern = "[.]Rd$",
    full.names = TRUE
  )
  if (!length(rd_files)) {
    stop("No generated Rd files were found for reference-page validation.", call. = FALSE)
  }

  canonical <- tools::file_path_sans_ext(basename(rd_files))
  public_aliases <- unlist(lapply(rd_files, function(path) {
    lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
    aliases <- grep("^\\\\alias\\{.*\\}$", lines, value = TRUE)
    aliases <- sub("^\\\\alias\\{(.*)\\}$", "\\1", aliases)
    aliases[!startsWith(aliases, ".")]
  }), use.names = FALSE)

  sort(unique(file.path(
    "reference",
    paste0(c("index", canonical, public_aliases), ".html")
  )), method = "radix")
}

.pkgdown_actual_reference_pages <- function(site_dir) {
  reference_dir <- file.path(site_dir, "reference")
  if (!dir.exists(reference_dir)) return(character())
  pages <- list.files(
    reference_dir,
    pattern = "[.]html$",
    recursive = FALSE,
    full.names = TRUE
  )
  sort(file.path("reference", basename(pages)), method = "radix")
}

.pkgdown_assert_reference_pages <- function(site_dir, repo_root) {
  expected <- .pkgdown_expected_reference_pages(repo_root)
  actual <- .pkgdown_actual_reference_pages(site_dir)
  missing <- setdiff(expected, actual)
  extra <- setdiff(actual, expected)
  if (length(missing) || length(extra)) {
    details <- c(
      if (length(missing)) paste0("missing: ", paste(missing, collapse = ", ")),
      if (length(extra)) paste0("stale/extra: ", paste(extra, collapse = ", "))
    )
    stop(
      "Rendered reference pages do not match the current Rd/export contract (",
      paste(details, collapse = "; "), ").",
      call. = FALSE
    )
  }
  invisible(expected)
}

.pkgdown_git_output <- function(repo_root, arguments, fail = TRUE) {
  output <- suppressWarnings(system2(
    "git",
    c("-C", shQuote(repo_root), arguments),
    stdout = TRUE,
    stderr = TRUE
  ))
  status <- attr(output, "status")
  status <- if (is.null(status)) 0L else as.integer(status)
  if (fail && status != 0L) {
    stop(
      "Git command failed while creating pkgdown provenance: ",
      paste(output, collapse = "\n"),
      call. = FALSE
    )
  }
  list(status = status, output = output)
}

.pkgdown_source_files <- function(repo_root) {
  listed <- .pkgdown_git_output(
    repo_root,
    c("ls-files", "--cached", "--others", "--exclude-standard")
  )$output
  listed <- sort(unique(listed[nzchar(listed)]), method = "radix")
  paths <- file.path(repo_root, listed)
  keep <- file.exists(paths) & !dir.exists(paths)
  list(relative = listed[keep], absolute = paths[keep])
}

.pkgdown_hash_records <- function(paths, root) {
  relative <- .pkgdown_relative_paths(paths, root)
  hashes <- .pkgdown_sha256(paths)
  sizes <- unname(file.info(paths)$size)
  order_index <- order(relative, method = "radix")
  lapply(order_index, function(index) {
    list(
      path = relative[[index]],
      sha256 = hashes[[index]],
      size = as.numeric(sizes[[index]])
    )
  })
}

.pkgdown_records_fingerprint <- function(records) {
  payload <- paste(vapply(records, function(record) {
    paste(record$path, record$sha256, format(record$size, scientific = FALSE), sep = "\t")
  }, character(1)), collapse = "\n")
  temp <- tempfile("shennong-pkgdown-fingerprint-")
  on.exit(unlink(temp, force = TRUE), add = TRUE)
  writeBin(charToRaw(payload), temp)
  .pkgdown_sha256(temp)[[1]]
}

.pkgdown_source_identity <- function(repo_root) {
  source_files <- .pkgdown_source_files(repo_root)
  records <- .pkgdown_hash_records(source_files$absolute, repo_root)
  status <- .pkgdown_git_output(
    repo_root,
    c("status", "--porcelain", "--untracked-files=normal")
  )$output
  sha <- .pkgdown_git_output(repo_root, c("rev-parse", "HEAD"))$output
  if (length(sha) != 1L || !grepl("^[0-9a-f]{40}$", sha[[1]])) {
    stop("Unable to resolve the source Git SHA.", call. = FALSE)
  }
  dirty <- length(status) > 0L
  fingerprint <- .pkgdown_records_fingerprint(records)
  list(
    kind = if (dirty) "dirty_fingerprint" else "git_sha",
    value = if (dirty) fingerprint else sha[[1]],
    git_sha = sha[[1]],
    dirty = dirty,
    fingerprint = fingerprint,
    file_count = length(records)
  )
}

.pkgdown_write_build_manifest <- function(
    site_dir,
    repo_root,
    runtime_report,
    manifest_path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Install `jsonlite` to write the pkgdown build manifest.", call. = FALSE)
  }
  site_dir <- normalizePath(site_dir, winslash = "/", mustWork = TRUE)
  repo_root <- normalizePath(repo_root, winslash = "/", mustWork = TRUE)
  runtime_report <- normalizePath(runtime_report, winslash = "/", mustWork = TRUE)
  .pkgdown_assert_static_allowlist(site_dir)
  reference_pages <- .pkgdown_assert_reference_pages(site_dir, repo_root)
  site_files <- .pkgdown_site_files(site_dir)
  site_records <- .pkgdown_hash_records(site_files, site_dir)

  manifest <- list(
    schema_version = "1.0.0",
    build_time_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    source = .pkgdown_source_identity(repo_root),
    runtime_report = list(
      name = "runtime-coverage-core-summary.json",
      sha256 = .pkgdown_sha256(runtime_report)[[1]]
    ),
    reference_pages = reference_pages,
    site = list(
      file_count = length(site_records),
      files = site_records
    )
  )
  dir.create(dirname(manifest_path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(
    manifest,
    manifest_path,
    auto_unbox = TRUE,
    pretty = TRUE,
    null = "null"
  )
  invisible(normalizePath(manifest_path, winslash = "/", mustWork = TRUE))
}

.pkgdown_manifest_scalar <- function(value, label) {
  if (length(value) != 1L || !is.character(value) || is.na(value) || !nzchar(value)) {
    stop("Pkgdown build manifest has invalid `", label, "`.", call. = FALSE)
  }
  value
}

.pkgdown_verify_build_manifest <- function(
    site_dir,
    repo_root,
    runtime_report,
    manifest_path) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Install `jsonlite` to verify the pkgdown build manifest.", call. = FALSE)
  }
  if (!file.exists(manifest_path)) {
    stop(
      "The same-build pkgdown manifest is missing; rebuild with ",
      "scripts/build-pkgdown.R --real.",
      call. = FALSE
    )
  }
  manifest <- jsonlite::read_json(manifest_path, simplifyVector = FALSE)
  if (!identical(manifest$schema_version, "1.0.0")) {
    stop("Pkgdown build manifest schema_version must be 1.0.0.", call. = FALSE)
  }
  build_time <- .pkgdown_manifest_scalar(manifest$build_time_utc, "build_time_utc")
  parsed_build_time <- suppressWarnings(as.POSIXct(build_time, tz = "UTC"))
  if (is.na(parsed_build_time)) {
    stop("Pkgdown build manifest has an invalid UTC build time.", call. = FALSE)
  }

  current_source <- .pkgdown_source_identity(repo_root)
  recorded_fingerprint <- .pkgdown_manifest_scalar(
    manifest$source$fingerprint,
    "source.fingerprint"
  )
  if (!identical(recorded_fingerprint, current_source$fingerprint)) {
    stop(
      "Pkgdown build manifest source fingerprint is stale; rebuild the real-data site ",
      "(recorded ", substr(recorded_fingerprint, 1L, 12L),
      ", current ", substr(current_source$fingerprint, 1L, 12L), ").",
      call. = FALSE
    )
  }
  recorded_kind <- .pkgdown_manifest_scalar(manifest$source$kind, "source.kind")
  recorded_value <- .pkgdown_manifest_scalar(manifest$source$value, "source.value")
  if (!recorded_kind %in% c("git_sha", "dirty_fingerprint")) {
    stop("Pkgdown build manifest source.kind is unsupported.", call. = FALSE)
  }
  expected_value <- if (identical(recorded_kind, "git_sha")) {
    current_source$git_sha
  } else {
    current_source$fingerprint
  }
  if (!identical(recorded_value, expected_value)) {
    stop("Pkgdown build manifest source identity does not match the checkout.", call. = FALSE)
  }

  runtime_hash <- .pkgdown_sha256(runtime_report)[[1]]
  runtime_name <- .pkgdown_manifest_scalar(
    manifest$runtime_report$name,
    "runtime_report.name"
  )
  if (!identical(runtime_name, "runtime-coverage-core-summary.json")) {
    stop("Pkgdown build manifest is bound to an unexpected runtime report.", call. = FALSE)
  }
  recorded_runtime_hash <- .pkgdown_manifest_scalar(
    manifest$runtime_report$sha256,
    "runtime_report.sha256"
  )
  if (!identical(recorded_runtime_hash, runtime_hash)) {
    stop(
      "Pkgdown build manifest runtime-report SHA-256 does not match; ",
      "rerun runtime coverage and the site build together.",
      call. = FALSE
    )
  }

  expected_references <- .pkgdown_expected_reference_pages(repo_root)
  recorded_references <- sort(
    unique(unlist(manifest$reference_pages, use.names = FALSE)),
    method = "radix"
  )
  if (!identical(recorded_references, expected_references)) {
    stop(
      "Pkgdown build manifest reference-page contract is stale.",
      call. = FALSE
    )
  }
  .pkgdown_assert_reference_pages(site_dir, repo_root)

  recorded_files <- manifest$site$files
  if (!is.list(recorded_files) || !length(recorded_files)) {
    stop("Pkgdown build manifest contains no site file hashes.", call. = FALSE)
  }
  recorded_file_count <- manifest$site$file_count
  if (length(recorded_file_count) != 1L || !is.numeric(recorded_file_count) ||
      is.na(recorded_file_count) || recorded_file_count != length(recorded_files)) {
    stop("Pkgdown build manifest site.file_count is invalid.", call. = FALSE)
  }
  recorded_paths <- vapply(recorded_files, function(record) {
    .pkgdown_manifest_scalar(record$path, "site.files[].path")
  }, character(1))
  if (anyDuplicated(recorded_paths) || any(startsWith(recorded_paths, "/")) ||
      any(grepl("(^|/)\\.\\.(/|$)", recorded_paths, perl = TRUE))) {
    stop("Pkgdown build manifest contains unsafe or duplicate site paths.", call. = FALSE)
  }

  site_files <- .pkgdown_site_files(site_dir)
  current_paths <- .pkgdown_relative_paths(site_files, site_dir)
  missing <- setdiff(recorded_paths, current_paths)
  extra <- setdiff(current_paths, recorded_paths)
  if (length(missing) || length(extra)) {
    details <- c(
      if (length(missing)) paste0("missing: ", paste(missing, collapse = ", ")),
      if (length(extra)) paste0("extra: ", paste(extra, collapse = ", "))
    )
    stop(
      "Pkgdown site file set differs from the same-build manifest (",
      paste(details, collapse = "; "), ").",
      call. = FALSE
    )
  }

  index <- match(recorded_paths, current_paths)
  hashes <- .pkgdown_sha256(site_files[index])
  sizes <- unname(file.info(site_files[index])$size)
  recorded_hashes <- vapply(recorded_files, function(record) {
    .pkgdown_manifest_scalar(record$sha256, "site.files[].sha256")
  }, character(1))
  recorded_sizes <- vapply(recorded_files, function(record) {
    value <- record$size
    if (length(value) != 1L || !is.numeric(value) || is.na(value) || value < 0) {
      stop("Pkgdown build manifest has invalid site file size.", call. = FALSE)
    }
    as.numeric(value)
  }, numeric(1))
  changed <- hashes != recorded_hashes | sizes != recorded_sizes
  if (any(changed)) {
    stop(
      "Pkgdown site files changed after the manifest was created: ",
      paste(recorded_paths[changed], collapse = ", "),
      call. = FALSE
    )
  }

  invisible(manifest)
}
