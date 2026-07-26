#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", all_args, value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else file.path(getwd(), "scripts", "real-data", "validate-real-data.R")
repo_root <- normalizePath(file.path(dirname(script_file), "..", ".."), winslash = "/", mustWork = TRUE)
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

.option <- function(name, default = NULL) {
  prefix <- paste0("--", name, "=")
  inline <- args[startsWith(args, prefix)]
  if (length(inline)) return(sub(prefix, "", inline[[length(inline)]], fixed = TRUE))
  index <- match(paste0("--", name), args)
  if (!is.na(index) && index < length(args)) return(args[[index + 1L]])
  default
}

.flag <- function(name) paste0("--", name) %in% args

if (.flag("help")) {
  cat(
    "Usage: Rscript scripts/real-data/validate-real-data.R [options]\n",
    "\n",
    "Options:\n",
    "  --bundle NAME       Validate one bundle or all bundles (default: all)\n",
    "  --root PATH         Local data root (default: SHENNONG_REAL_DATA_DIR or\n",
    "                      data-local/pkgdown-real)\n",
    "  --allow-missing     Report missing local artifacts without failing\n",
    "  --coverage-only     Validate sources.yml and coverage.csv only\n",
    "  --help              Show this message\n",
    sep = ""
  )
  quit(status = 0L)
}

for (package in c("yaml", "jsonlite")) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop("Install `", package, "` to validate real-data infrastructure.", call. = FALSE)
  }
}

source_file <- file.path(repo_root, "scripts", "real-data", "sources.yml")
coverage_file <- file.path(repo_root, "scripts", "real-data", "coverage.csv")
coverage_exclusion_file <- file.path(
  repo_root, "scripts", "real-data", "coverage-exclusions.csv"
)
sources <- yaml::read_yaml(source_file)
bundle_names <- names(sources$bundles)
bundle <- .option("bundle", "all")
if (!bundle %in% c("all", bundle_names)) {
  stop("Unknown bundle `", bundle, "`.", call. = FALSE)
}

.absolute_path <- function(path) {
  path <- path.expand(path)
  if (!grepl("^/", path)) path <- file.path(repo_root, path)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

default_root <- Sys.getenv(
  sources$root_environment_variable %||% "SHENNONG_REAL_DATA_DIR",
  unset = file.path(repo_root, sources$default_root %||% "data-local/pkgdown-real")
)
root <- .absolute_path(.option("root", default_root))
allow_missing <- .flag("allow-missing")
coverage_only <- .flag("coverage-only")
errors <- character()
warnings <- character()

.record_error <- function(...) errors <<- c(errors, paste0(...))
.record_warning <- function(...) warnings <<- c(warnings, paste0(...))
.expect <- function(condition, ...) {
  if (!isTRUE(condition)) .record_error(...)
  invisible(condition)
}

.sha256 <- function(path) {
  command <- Sys.which("sha256sum")
  if (!nzchar(command)) stop("`sha256sum` is required for artifact validation.", call. = FALSE)
  output <- system2(command, path, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status") %||% 0L
  if (!identical(as.integer(status), 0L) || !length(output)) return(NA_character_)
  strsplit(output[[1]], "[[:space:]]+")[[1]][[1]]
}

.validate_storage_boundaries <- function() {
  gitignore <- readLines(file.path(repo_root, ".gitignore"), warn = FALSE)
  rbuildignore <- readLines(file.path(repo_root, ".Rbuildignore"), warn = FALSE)
  .expect(any(grepl("(^|/)data-local/?$", gitignore)), ".gitignore must exclude data-local/.")
  .expect(any(grepl("data-local", rbuildignore, fixed = TRUE)), ".Rbuildignore must exclude data-local.")
  repo <- paste0(normalizePath(repo_root, winslash = "/", mustWork = TRUE), "/")
  candidate <- paste0(sub("/+$", "", root), "/")
  local_root <- paste0(normalizePath(file.path(repo_root, "data-local"), winslash = "/", mustWork = FALSE), "/")
  .expect(!startsWith(candidate, repo) || startsWith(candidate, local_root), "Local data root is inside the repository but outside data-local/: ", root)
  tracked <- tryCatch(
    system2("git", c("-C", repo_root, "ls-files", "--", "data-local"), stdout = TRUE, stderr = TRUE),
    error = function(e) character()
  )
  .expect(length(tracked) == 0L, "Real-data files are tracked by Git: ", paste(tracked, collapse = ", "))
}

.validate_sources <- function() {
  .expect(identical(as.character(sources$schema_version), "1.0.0"), "sources.yml schema_version must be 1.0.0.")
  .expect(setequal(bundle_names, c("kotliarov_pbmc", "melanoma_bridge", "hermann_spermatogenesis", "visium_lymph_node")), "sources.yml must define exactly the four planned bundles.")
  for (current in bundle_names) {
    specification <- sources$bundles[[current]]
    .expect(length(specification$output) >= 1L, "Bundle ", current, " has no output path.")
    .expect(length(specification$sources) >= 1L, "Bundle ", current, " has no public source.")
    .expect(length(specification$intended_coverage) >= 1L, "Bundle ", current, " has no intended coverage.")
    .expect(!is.null(specification$subset$strategy), "Bundle ", current, " lacks a deterministic subset strategy.")
    for (index in seq_along(specification$sources)) {
      source <- specification$sources[[index]]
      .expect(!is.null(source$provider), "Bundle ", current, " source ", index, " lacks provider.")
      .expect(!is.null(source$accession), "Bundle ", current, " source ", index, " lacks accession.")
      .expect(!is.null(source$landing_page), "Bundle ", current, " source ", index, " lacks landing_page.")
      .expect(!is.null(source$license), "Bundle ", current, " source ", index, " lacks license/data-use terms.")
      .expect(!is.null(source$citation), "Bundle ", current, " source ", index, " lacks citation.")
    }
  }
}

.namespace_exports <- function() {
  namespace <- readLines(file.path(repo_root, "NAMESPACE"), warn = FALSE)
  sub("^export\\(([^)]+)\\)$", "\\1", grep("^export\\(sn_", namespace, value = TRUE))
}

.validate_coverage <- function() {
  .expect(file.exists(coverage_file), "Missing scripts/real-data/coverage.csv.")
  if (!file.exists(coverage_file)) return(invisible(NULL))
  coverage <- utils::read.csv(coverage_file, stringsAsFactors = FALSE, check.names = FALSE)
  required <- c(
    "function", "category", "bundle", "level", "backend", "notes", "article"
  )
  .expect(all(required %in% names(coverage)), "coverage.csv must contain: ", paste(required, collapse = ", "), ".")
  if (!all(required %in% names(coverage))) return(invisible(NULL))
  .expect(!anyDuplicated(coverage[["function"]]), "coverage.csv contains duplicate functions.")
  .expect(all(coverage$category %in% c("analysis", "visualization")), "coverage.csv category must be analysis or visualization.")
  .expect(all(coverage$bundle %in% bundle_names), "coverage.csv contains an unknown bundle.")
  .expect(all(coverage$level %in% c("core", "extended")), "coverage.csv level must be core or extended.")
  coverage[["article"]] <- trimws(coverage[["article"]])
  mapped <- nzchar(coverage[["article"]])
  valid_article_path <- !mapped | grepl(
    "^vignettes/[A-Za-z0-9._-]+[.]Rmd$",
    coverage[["article"]]
  )
  .expect(
    all(valid_article_path),
    "coverage.csv article paths must name repository vignettes: ",
    paste(unique(coverage[["article"]][!valid_article_path]), collapse = ", "),
    "."
  )
  for (index in which(mapped & valid_article_path)) {
    article <- coverage[["article"]][[index]]
    article_path <- file.path(repo_root, article)
    .expect(
      file.exists(article_path),
      "coverage.csv article does not exist for ",
      coverage[["function"]][[index]],
      ": ",
      article,
      "."
    )
    if (!file.exists(article_path)) next
    article_source <- readLines(article_path, warn = FALSE)
    function_pattern <- paste0(
      "\\b", coverage[["function"]][[index]], "\\s*\\("
    )
    .expect(
      any(grepl(function_pattern, article_source, perl = TRUE)),
      "coverage.csv maps ",
      coverage[["function"]][[index]],
      " to ",
      article,
      ", but the function does not appear in that article's source."
    )
  }
  exports <- .namespace_exports()
  unknown <- setdiff(coverage[["function"]], exports)
  .expect(length(unknown) == 0L, "coverage.csv contains non-exported functions: ", paste(unknown, collapse = ", "), ".")
  .expect(
    file.exists(coverage_exclusion_file),
    "Missing scripts/real-data/coverage-exclusions.csv."
  )
  if (file.exists(coverage_exclusion_file)) {
    exclusions <- utils::read.csv(
      coverage_exclusion_file,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    exclusion_columns <- c("function", "category", "reason")
    .expect(
      all(exclusion_columns %in% names(exclusions)),
      "coverage-exclusions.csv must contain: ",
      paste(exclusion_columns, collapse = ", "), "."
    )
    if (all(exclusion_columns %in% names(exclusions))) {
      allowed_exclusion_categories <- c(
        "acceleration_control", "backend_adapter", "interpretation_reporting",
        "io_project", "result_management", "runtime_admin",
        "signature_management"
      )
      .expect(
        !anyDuplicated(exclusions[["function"]]),
        "coverage-exclusions.csv contains duplicate functions."
      )
      .expect(
        all(nzchar(trimws(exclusions[["reason"]]))),
        "Every coverage exclusion must have a reason."
      )
      .expect(
        all(exclusions[["category"]] %in% allowed_exclusion_categories),
        "coverage-exclusions.csv contains an unknown category."
      )
      .expect(
        length(intersect(coverage[["function"]], exclusions[["function"]])) == 0L,
        "Functions cannot appear in both coverage.csv and coverage-exclusions.csv."
      )
      missing_classification <- setdiff(
        exports,
        c(coverage[["function"]], exclusions[["function"]])
      )
      stale_exclusions <- setdiff(exclusions[["function"]], exports)
      .expect(
        length(missing_classification) == 0L,
        "Exported sn_* functions lack a runtime-coverage classification: ",
        paste(missing_classification, collapse = ", "), "."
      )
      .expect(
        length(stale_exclusions) == 0L,
        "coverage-exclusions.csv contains non-exported functions: ",
        paste(stale_exclusions, collapse = ", "), "."
      )
    }
  }
  plot_engine <- c(
    "sn_apply_figure_profile", "sn_export_figure", "sn_export_figure_bundle",
    "sn_figure_spec", "sn_recommend_figure_size", "sn_save_figure", "sn_validate_figure"
  )
  expected_visuals <- sort(unique(c(exports[grepl("^sn_plot_", exports)], plot_engine)))
  covered_visuals <- sort(coverage[["function"]][coverage$category == "visualization"])
  missing_visuals <- setdiff(expected_visuals, covered_visuals)
  extra_visuals <- setdiff(covered_visuals, expected_visuals)
  .expect(length(missing_visuals) == 0L, "coverage.csv misses visualization exports: ", paste(missing_visuals, collapse = ", "), ".")
  .expect(length(extra_visuals) == 0L, "coverage.csv misclassifies visualization exports: ", paste(extra_visuals, collapse = ", "), ".")
  .expect(sum(coverage$category == "analysis") >= 79L, "coverage.csv must declare at least 79 analysis APIs.")
  .expect(sum(coverage$category == "visualization") == length(expected_visuals), "coverage.csv visualization count differs from NAMESPACE.")
  message(
    "Coverage declarations: ", sum(coverage$category == "analysis"), " analysis; ",
    sum(coverage$category == "visualization"), " visualization; ",
    sum(coverage$level == "core"), " core; ", sum(coverage$level == "extended"),
    " extended; ", sum(mapped), " article mappings; ", sum(!mapped), " unmapped."
  )
}

.manifest_path <- function(artifact_path) {
  file.path(root, "manifests", paste0(tools::file_path_sans_ext(basename(artifact_path)), ".json"))
}

.validate_manifest <- function(artifact_path, bundle) {
  manifest_path <- .manifest_path(artifact_path)
  .expect(file.exists(manifest_path), "Missing manifest for ", artifact_path, ".")
  if (!file.exists(manifest_path)) return(NULL)
  manifest <- jsonlite::read_json(manifest_path, simplifyVector = TRUE)
  .expect(identical(as.character(manifest$schema_version), "1.0.0"), "Manifest schema mismatch: ", manifest_path)
  .expect(identical(as.character(manifest$bundle), bundle), "Manifest bundle mismatch: ", manifest_path)
  actual_bytes <- unname(file.info(artifact_path)$size)
  .expect(identical(as.numeric(manifest$artifact$bytes), as.numeric(actual_bytes)), "Artifact size mismatch: ", artifact_path)
  actual_sha <- .sha256(artifact_path)
  .expect(identical(as.character(manifest$artifact$sha256), actual_sha), "Artifact SHA-256 mismatch: ", artifact_path)
  .expect(length(manifest$sources) >= 1L, "Manifest source IDs are empty: ", manifest_path)
  .expect(!is.null(manifest$selection), "Manifest selection metadata are empty: ", manifest_path)
  manifest_packages <- manifest[["packages"]]
  .expect(
    length(manifest_packages) > 0L &&
      !is.null(names(manifest_packages)) &&
      all(nzchar(names(manifest_packages))) &&
      all(nzchar(as.character(manifest_packages))),
    "Manifest package versions must be a named package-to-version map: ",
    manifest_path
  )
  manifest
}

.read_qs2 <- function(path) {
  if (!requireNamespace("qs2", quietly = TRUE)) stop("Install `qs2` to validate local artifacts.", call. = FALSE)
  qs2::qs_read(path)
}

.expect_seurat <- function(object, label) {
  .expect(inherits(object, "Seurat"), label, " is not a Seurat object.")
  if (!inherits(object, "Seurat")) return(FALSE)
  .expect(nrow(object) > 0L && ncol(object) > 0L, label, " is empty.")
  TRUE
}

.validate_kotliarov <- function(path) {
  object <- .read_qs2(path)
  if (!.expect_seurat(object, "Kotliarov artifact")) return()
  .expect(all(c("RNA", "ADT") %in% names(object@assays)), "Kotliarov artifact must contain RNA and ADT assays.")
  .expect(ncol(object) >= 1000L, "Kotliarov artifact has too few cells.")
  .expect(nrow(object) >= 10000L, "Kotliarov artifact must retain the broad RNA feature space for QC.")
  metadata <- object[[]]
  .expect("real_sample" %in% names(metadata), "Kotliarov artifact lacks real_sample.")
  .expect("real_response" %in% names(metadata), "Kotliarov artifact lacks real_response.")
  if ("real_sample" %in% names(metadata)) .expect(length(unique(metadata$real_sample)) >= 10L, "Kotliarov artifact retains too few biological samples.")
  .expect(identical(object@misc$real_data$source, "scRNAseq/kotliarov-pbmc"), "Kotliarov provenance source mismatch.")
}

.validate_melanoma <- function(path) {
  object <- .read_qs2(path)
  if (!.expect_seurat(object, "GSE72056 artifact")) return()
  layers <- SeuratObject::Layers(object[["RNA"]])
  .expect("data" %in% layers, "GSE72056 normalized expression must be stored in the data layer.")
  .expect(!"counts" %in% layers, "GSE72056 normalized expression must not be mislabeled as raw counts.")
  .expect(all(c("tumor", "malignant_call", "cell_type") %in% names(object[[]])), "GSE72056 artifact lacks author annotations.")
  .expect(length(unique(object$tumor)) >= 10L, "GSE72056 artifact retains too few tumors.")
  .expect(all(c("malignant", "non_malignant") %in% unique(object$malignant_call)), "GSE72056 artifact must contain malignant and non-malignant cells.")
}

.validate_tcga <- function(path) {
  artifact <- .read_qs2(path)
  .expect(is.list(artifact), "TCGA-SKCM artifact must be a list.")
  required <- c(
    "rsem_tpm", "log2_uq_rsem_tpm", "sample_data",
    "feature_data", "provenance"
  )
  .expect(all(required %in% names(artifact)), "TCGA-SKCM artifact lacks: ", paste(setdiff(required, names(artifact)), collapse = ", "), ".")
  if (!all(required %in% names(artifact))) return()
  .expect(
    identical(dim(artifact$rsem_tpm), dim(artifact$log2_uq_rsem_tpm)),
    "TCGA expression assay dimensions differ."
  )
  .expect(
    identical(colnames(artifact$rsem_tpm), rownames(artifact$sample_data)) &&
      identical(colnames(artifact$log2_uq_rsem_tpm), rownames(artifact$sample_data)),
    "TCGA sample metadata are not aligned to both expression assays."
  )
  .expect(
    all(is.finite(artifact$rsem_tpm)),
    "TCGA RSEM TPM contains missing or non-finite values."
  )
  .expect(
    all(is.finite(artifact$log2_uq_rsem_tpm)),
    "TCGA log2 upper-quartile RSEM TPM contains missing or non-finite values."
  )
  .expect(all(c("time", "event", "sample_type") %in% names(artifact$sample_data)), "TCGA sample_data lacks survival or sample-type fields.")
  .expect(all(stats::na.omit(unique(artifact$sample_data$event)) %in% c(0, 1)), "TCGA event must be binary 0/1.")
  .expect(sum(is.finite(artifact$sample_data$time) & artifact$sample_data$time > 0) >= 50L, "TCGA artifact has too few complete survival samples.")
  .expect(
    identical(as.character(artifact$provenance$data_version), "2.1.1"),
    "TCGA artifact must record curatedTCGAData version 2.1.1."
  )
}

.validate_hermann <- function(path) {
  object <- .read_qs2(path)
  if (!.expect_seurat(object, "Hermann artifact")) return()
  layers <- SeuratObject::Layers(object[["RNA"]])
  .expect(all(c("counts", "spliced", "unspliced") %in% layers), "Hermann artifact must contain counts, spliced, and unspliced layers.")
  .expect("real_cell_type" %in% names(object[[]]), "Hermann artifact lacks real_cell_type.")
  if ("real_cell_type" %in% names(object[[]])) .expect(length(unique(object$real_cell_type)) >= 4L, "Hermann artifact retains too few trajectory states.")
  .expect(ncol(object) >= 500L && nrow(object) >= 1000L, "Hermann artifact is too small for trajectory and velocity tests.")
}

.validate_visium <- function(path) {
  artifact <- .read_qs2(path)
  .expect(is.list(artifact) && all(c("sections", "provenance") %in% names(artifact)), "Visium artifact must contain sections and provenance.")
  if (!is.list(artifact) || !"sections" %in% names(artifact)) return()
  .expect(length(artifact$sections) == 2L, "Visium artifact must contain two real sections.")
  for (name in names(artifact$sections)) {
    object <- artifact$sections[[name]]
    if (!.expect_seurat(object, paste0("Visium section ", name))) next
    .expect(ncol(object) >= 200L, "Visium section ", name, " has too few spots.")
    .expect(length(object@images) >= 1L, "Visium section ", name, " lacks spatial image/coordinates.")
    .expect(all(c("real_section", "real_platform") %in% names(object[[]])), "Visium section ", name, " lacks standardized metadata.")
  }
}

.validate_artifact <- function(bundle, relative_path) {
  path <- file.path(root, relative_path)
  if (!file.exists(path)) {
    message("MISSING ", relative_path)
    if (!allow_missing) .record_error("Missing local artifact: ", relative_path)
    return(invisible(NULL))
  }
  .validate_manifest(path, bundle)
  if (identical(relative_path, "single-cell/kotliarov_pbmc.qs2")) .validate_kotliarov(path)
  if (identical(relative_path, "cancer/melanoma_gse72056.qs2")) .validate_melanoma(path)
  if (identical(relative_path, "bulk/tcga_skcm.qs2")) .validate_tcga(path)
  if (identical(relative_path, "trajectory/hermann_spermatogenesis.qs2")) .validate_hermann(path)
  if (identical(relative_path, "spatial/visium_lymph_node.qs2")) .validate_visium(path)
  message("VALID ", relative_path)
}

.validate_storage_boundaries()
.validate_sources()
.validate_coverage()

if (!coverage_only) {
  selected <- if (identical(bundle, "all")) bundle_names else bundle
  for (current in selected) {
    for (relative_path in sources$bundles[[current]]$output) {
      .validate_artifact(current, relative_path)
    }
  }
}

if (length(warnings)) {
  message("Warnings:\n- ", paste(unique(warnings), collapse = "\n- "))
}
if (length(errors)) {
  stop("Real-data validation failed:\n- ", paste(unique(errors), collapse = "\n- "), call. = FALSE)
}
message("Real-data validation passed for ", if (coverage_only) "infrastructure only" else bundle, ".")
