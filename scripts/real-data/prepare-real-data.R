#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
all_args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", all_args, value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else file.path(getwd(), "scripts", "real-data", "prepare-real-data.R")
repo_root <- normalizePath(file.path(dirname(script_file), "..", ".."), winslash = "/", mustWork = TRUE)

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
    "Usage: Rscript scripts/real-data/prepare-real-data.R [options]\n",
    "\n",
    "Options:\n",
    "  --bundle NAME       kotliarov_pbmc, melanoma_bridge,\n",
    "                      hermann_spermatogenesis, visium_lymph_node, or all\n",
    "  --root PATH         Local data root (default: SHENNONG_REAL_DATA_DIR or\n",
    "                      data-local/pkgdown-real)\n",
    "  --seed INTEGER      Deterministic subsetting seed (default: 717)\n",
    "  --force             Replace existing raw downloads and derived artifacts\n",
    "  --list              Print bundle names without preparing data\n",
    "  --help              Show this message\n",
    sep = ""
  )
  quit(status = 0L)
}

source_file <- file.path(repo_root, "scripts", "real-data", "sources.yml")
if (!requireNamespace("yaml", quietly = TRUE)) {
  stop("Install `yaml` to read scripts/real-data/sources.yml.", call. = FALSE)
}
sources <- yaml::read_yaml(source_file)
bundle_names <- names(sources$bundles)
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

if (.flag("list")) {
  cat(paste(bundle_names, collapse = "\n"), "\n")
  quit(status = 0L)
}

bundle <- .option("bundle", "all")
if (!bundle %in% c("all", bundle_names)) {
  stop("Unknown bundle `", bundle, "`. Use --list to inspect valid names.", call. = FALSE)
}

seed <- suppressWarnings(as.integer(.option("seed", sources$rules$default_seed %||% 717L)))
if (length(seed) != 1L || is.na(seed)) stop("`--seed` must be one integer.", call. = FALSE)
replace_existing <- .flag("force")

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

.assert_safe_root <- function(path) {
  repo <- paste0(normalizePath(repo_root, winslash = "/", mustWork = TRUE), "/")
  local_root <- paste0(normalizePath(file.path(repo_root, "data-local"), winslash = "/", mustWork = FALSE), "/")
  candidate <- paste0(sub("/+$", "", path), "/")
  if (startsWith(candidate, repo) && !startsWith(candidate, local_root)) {
    stop(
      "Refusing to write real data inside the repository outside `data-local/`: ",
      path,
      call. = FALSE
    )
  }
  invisible(path)
}

.assert_safe_root(root)
for (directory in c(
  "raw", "single-cell", "cancer", "bulk", "trajectory", "spatial",
  "manifests", "cache", "results", "logs"
)) {
  dir.create(file.path(root, directory), recursive = TRUE, showWarnings = FALSE)
}

.require <- function(packages, reason) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing package(s) ", paste(missing, collapse = ", "), " required ", reason, ".",
      call. = FALSE
    )
  }
}

.download <- function(url, destination) {
  if (file.exists(destination) && !replace_existing && file.info(destination)$size > 0) {
    message("Reusing ", destination)
    return(destination)
  }
  .require("curl", "to download public real-data sources")
  dir.create(dirname(destination), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(destination, ".partial")
  unlink(temporary)
  message("Downloading ", url)
  curl::curl_download(url, temporary, quiet = FALSE, mode = "wb")
  if (!file.exists(temporary) || file.info(temporary)$size <= 0) {
    stop("Download produced an empty file: ", url, call. = FALSE)
  }
  if (!file.rename(temporary, destination)) {
    stop("Could not move completed download to ", destination, call. = FALSE)
  }
  destination
}

.sha256 <- function(path) {
  command <- Sys.which("sha256sum")
  if (!nzchar(command)) stop("`sha256sum` is required to validate local artifacts.", call. = FALSE)
  output <- system2(command, path, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status") %||% 0L
  if (!identical(as.integer(status), 0L) || !length(output)) {
    stop("Could not calculate SHA-256 for ", path, call. = FALSE)
  }
  strsplit(output[[1]], "[[:space:]]+")[[1]][[1]]
}

.save_qs2 <- function(object, path) {
  .require("qs2", "to save local real-data artifacts")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temporary <- paste0(path, ".partial")
  unlink(temporary)
  qs2::qs_save(object = object, file = temporary)
  if (!file.rename(temporary, path)) stop("Could not finalize ", path, call. = FALSE)
  path
}

.package_versions <- function(packages) {
  installed <- packages[vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  stats::setNames(
    vapply(installed, function(package) as.character(utils::packageVersion(package)), character(1)),
    installed
  )
}

.write_manifest <- function(path, bundle, component, source_ids, summary, selection) {
  .require("jsonlite", "to write real-data manifests")
  info <- file.info(path)
  manifest <- list(
    schema_version = "1.0.0",
    bundle = bundle,
    component = component,
    artifact = list(
      relative_path = substring(normalizePath(path, winslash = "/", mustWork = TRUE), nchar(root) + 2L),
      bytes = unname(info$size),
      sha256 = .sha256(path)
    ),
    sources = source_ids,
    selection = selection,
    summary = summary,
    prepared_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    seed = seed,
    packages = as.list(.package_versions(c(
      "scRNAseq", "SingleCellExperiment", "SummarizedExperiment", "Seurat",
      "SeuratObject", "Matrix", "data.table", "curatedTCGAData",
      "MultiAssayExperiment", "qs2"
    ))),
    script = "scripts/real-data/prepare-real-data.R"
  )
  manifest_path <- file.path(
    root, "manifests", paste0(tools::file_path_sans_ext(basename(path)), ".json")
  )
  jsonlite::write_json(manifest, manifest_path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  message("Wrote ", path, " (", manifest$artifact$sha256, ")")
  invisible(manifest)
}

.find_column <- function(data, candidates, required = TRUE) {
  lowered <- tolower(names(data))
  exact <- match(tolower(candidates), lowered)
  exact <- exact[!is.na(exact)]
  if (length(exact)) return(names(data)[exact[[1]]])
  patterns <- paste(candidates, collapse = "|")
  fuzzy <- grep(patterns, names(data), ignore.case = TRUE, value = TRUE)
  if (length(fuzzy)) return(fuzzy[[1]])
  if (isTRUE(required)) {
    stop("Could not find a required metadata column matching: ", paste(candidates, collapse = ", "), call. = FALSE)
  }
  NULL
}

.stratified_target <- function(indices, strata, target) {
  indices <- as.integer(indices)
  strata <- as.character(strata)
  strata[is.na(strata) | !nzchar(strata)] <- "missing"
  groups <- split(indices, strata)
  cap <- max(1L, ceiling(target / length(groups)))
  selected <- unlist(lapply(groups, function(current) {
    if (length(current) <= cap) current else sample(current, cap)
  }), use.names = FALSE)
  if (length(selected) < target) {
    remaining <- setdiff(indices, selected)
    selected <- c(selected, sample(remaining, min(length(remaining), target - length(selected))))
  }
  if (length(selected) > target) selected <- sample(selected, target)
  sort(unique(selected))
}

.as_sparse <- function(matrix) {
  .require("Matrix", "to create sparse local fixtures")
  if (inherits(matrix, "sparseMatrix")) return(methods::as(matrix, "dgCMatrix"))
  methods::as(Matrix::Matrix(as.matrix(matrix), sparse = TRUE), "dgCMatrix")
}

.row_variance <- function(matrix) {
  means <- Matrix::rowMeans(matrix)
  pmax(Matrix::rowMeans(matrix ^ 2) - means ^ 2, 0)
}

.existing_or_prepare <- function(path, expression) {
  manifest_path <- file.path(
    root, "manifests", paste0(tools::file_path_sans_ext(basename(path)), ".json")
  )
  if (file.exists(path) && file.exists(manifest_path) && !replace_existing) {
    message("Reusing derived artifact ", path)
    return(invisible(path))
  }
  if (file.exists(path) && !file.exists(manifest_path) && !replace_existing) {
    message("Rebuilding derived artifact because its manifest is missing: ", path)
  }
  base::force(expression)
}

.prepare_kotliarov <- function() {
  output <- file.path(root, "single-cell", "kotliarov_pbmc.qs2")
  .existing_or_prepare(output, {
    .require(
      c("scRNAseq", "SingleCellExperiment", "SummarizedExperiment", "SeuratObject", "Matrix"),
      "to prepare the Kotliarov PBMC CITE-seq fixture"
    )
    set.seed(seed)
    message("Retrieving Kotliarov PBMC CITE-seq through scRNAseq/gypsum.")
    rna_sce <- scRNAseq::KotliarovPBMCData(
      mode = "rna", ensembl = FALSE, location = FALSE, legacy = FALSE
    )
    adt_sce <- scRNAseq::KotliarovPBMCData(
      mode = "adt", ensembl = FALSE, location = FALSE, legacy = FALSE
    )
    common_cells <- intersect(colnames(rna_sce), colnames(adt_sce))
    if (!length(common_cells)) stop("Kotliarov RNA and ADT objects have no shared cell barcodes.", call. = FALSE)
    rna_sce <- rna_sce[, common_cells]
    adt_sce <- adt_sce[, common_cells]
    if (!identical(colnames(rna_sce), colnames(adt_sce))) {
      stop("Kotliarov RNA and ADT cell barcodes could not be exactly aligned.", call. = FALSE)
    }
    full_dimensions <- list(
      rna = c(features = nrow(rna_sce), cells = ncol(rna_sce)),
      adt = c(features = nrow(adt_sce), cells = ncol(adt_sce))
    )
    metadata <- as.data.frame(SummarizedExperiment::colData(rna_sce))
    sample_column <- .find_column(metadata, c("sample", "sample_id", "sampleid"))
    eligible <- seq_len(ncol(rna_sce))
    singlet_column <- .find_column(metadata, c("joint_classification_global", "joint.*classification"), required = FALSE)
    if (!is.null(singlet_column)) {
      singlets <- grepl("singlet", as.character(metadata[[singlet_column]]), ignore.case = TRUE)
      if (sum(singlets, na.rm = TRUE) >= 1000L) eligible <- which(singlets)
    }
    selected <- .stratified_target(
      eligible,
      metadata[[sample_column]][eligible],
      target = min(2000L, length(eligible))
    )
    rna_sce <- rna_sce[, selected]
    adt_sce <- adt_sce[, selected]
    metadata <- metadata[selected, , drop = FALSE]
    rownames(metadata) <- colnames(rna_sce)

    assay_names <- SummarizedExperiment::assayNames(rna_sce)
    rna_assay <- if ("counts" %in% assay_names) "counts" else assay_names[[1]]
    rna <- .as_sparse(SummarizedExperiment::assay(rna_sce, rna_assay))
    object <- SeuratObject::CreateSeuratObject(counts = rna, meta.data = metadata, project = "KotliarovPBMC")

    adt_names <- SummarizedExperiment::assayNames(adt_sce)
    adt_assay <- if ("counts" %in% adt_names) "counts" else adt_names[[1]]
    adt <- .as_sparse(SummarizedExperiment::assay(adt_sce, adt_assay))
    adt <- adt[, colnames(object), drop = FALSE]
    object[["ADT"]] <- SeuratObject::CreateAssay5Object(counts = adt)

    response_column <- .find_column(metadata, c("adjmfc.time", "response", "responder"), required = FALSE)
    batch_column <- .find_column(metadata, c("batch"), required = FALSE)
    lane_column <- .find_column(metadata, c("tenx_lane", "lane"), required = FALSE)
    object$real_sample <- as.character(metadata[[sample_column]])
    if (!is.null(response_column)) object$real_response <- as.character(metadata[[response_column]])
    if (!is.null(batch_column)) object$real_batch <- as.character(metadata[[batch_column]])
    if (!is.null(lane_column)) object$real_lane <- as.character(metadata[[lane_column]])
    object@misc$real_data <- list(
      bundle = "kotliarov_pbmc",
      source = "scRNAseq/kotliarov-pbmc",
      source_expression = "RNA UMI counts and ADT counts",
      source_dimensions = full_dimensions,
      subset_seed = seed,
      subset_strategy = "up to 100 eligible cells per biological sample",
      license = "CC0",
      citation = sources$bundles$kotliarov_pbmc$sources[[1]]$citation
    )
    .save_qs2(object, output)
    .write_manifest(
      output, "kotliarov_pbmc", "rna_adt", "scRNAseq/kotliarov-pbmc",
      summary = list(
        features = nrow(object), cells = ncol(object), assays = names(object@assays),
        samples = length(unique(object$real_sample)),
        responses = if ("real_response" %in% colnames(object[[]])) sort(unique(object$real_response)) else NULL
      ),
      selection = list(target_cells = 2000L, retained_cells = ncol(object), sample_column = sample_column)
    )
  })
  invisible(output)
}

.prepare_hermann <- function() {
  output <- file.path(root, "trajectory", "hermann_spermatogenesis.qs2")
  .existing_or_prepare(output, {
    .require(
      c("scRNAseq", "SummarizedExperiment", "SeuratObject", "Matrix"),
      "to prepare the Hermann spermatogenesis fixture"
    )
    set.seed(seed)
    message("Retrieving Hermann spermatogenesis through scRNAseq/gypsum.")
    sce <- scRNAseq::HermannSpermatogenesisData(strip = TRUE, location = TRUE, legacy = FALSE)
    full_dimensions <- c(features = nrow(sce), cells = ncol(sce))
    assay_names <- SummarizedExperiment::assayNames(sce)
    spliced_candidates <- assay_names[grepl("spliced", assay_names, ignore.case = TRUE) &
      !grepl("unspliced", assay_names, ignore.case = TRUE)]
    unspliced_candidates <- assay_names[grepl("unspliced", assay_names, ignore.case = TRUE)]
    if (!length(spliced_candidates) || !length(unspliced_candidates)) {
      stop("Hermann source does not expose both spliced and unspliced assays.", call. = FALSE)
    }
    spliced_name <- if ("spliced" %in% assay_names) "spliced" else spliced_candidates[[1]]
    unspliced_name <- if ("unspliced" %in% assay_names) "unspliced" else unspliced_candidates[[1]]
    metadata <- as.data.frame(SummarizedExperiment::colData(sce))
    type_column <- .find_column(metadata, c("cell_type", "celltype", "cell.type", "celltypes", "cluster"))
    selected_cells <- .stratified_target(
      seq_len(ncol(sce)), metadata[[type_column]], target = min(1000L, ncol(sce))
    )
    sce <- sce[, selected_cells]
    metadata <- metadata[selected_cells, , drop = FALSE]
    rownames(metadata) <- colnames(sce)
    spliced <- .as_sparse(SummarizedExperiment::assay(sce, spliced_name))
    unspliced <- .as_sparse(SummarizedExperiment::assay(sce, unspliced_name))
    total <- spliced + unspliced
    variances <- .row_variance(total)
    feature_order <- order(variances, decreasing = TRUE, na.last = NA)
    selected_features <- rownames(total)[utils::head(feature_order, min(2000L, length(feature_order)))]
    original_features <- selected_features
    unique_features <- make.unique(selected_features)
    spliced <- spliced[selected_features, , drop = FALSE]
    unspliced <- unspliced[selected_features, , drop = FALSE]
    total <- total[selected_features, , drop = FALSE]
    rownames(spliced) <- unique_features
    rownames(unspliced) <- unique_features
    rownames(total) <- unique_features

    object <- SeuratObject::CreateSeuratObject(
      counts = total, meta.data = metadata, project = "HermannSpermatogenesis"
    )
    SeuratObject::LayerData(object, assay = "RNA", layer = "spliced") <- spliced
    SeuratObject::LayerData(object, assay = "RNA", layer = "unspliced") <- unspliced
    object$real_cell_type <- as.character(metadata[[type_column]])
    object@misc$real_data <- list(
      bundle = "hermann_spermatogenesis",
      source = "scRNAseq/hermann-spermatogenesis",
      source_dimensions = full_dimensions,
      subset_seed = seed,
      subset_strategy = "cell-type-stratified sampling followed by total-count variance",
      original_features = stats::setNames(original_features, unique_features),
      license = "CC0",
      citation = sources$bundles$hermann_spermatogenesis$sources[[1]]$citation
    )
    .save_qs2(object, output)
    .write_manifest(
      output, "hermann_spermatogenesis", "spliced_unspliced",
      "scRNAseq/hermann-spermatogenesis",
      summary = list(
        features = nrow(object), cells = ncol(object),
        layers = SeuratObject::Layers(object[["RNA"]]),
        cell_types = sort(unique(object$real_cell_type))
      ),
      selection = list(
        target_cells = 1000L, retained_cells = ncol(object),
        target_features = 2000L, retained_features = nrow(object), cell_type_column = type_column
      )
    )
  })
  invisible(output)
}

.prepare_melanoma_single_cell <- function() {
  output <- file.path(root, "cancer", "melanoma_gse72056.qs2")
  .existing_or_prepare(output, {
    .require(c("data.table", "SeuratObject", "Matrix"), "to prepare GSE72056")
    set.seed(seed)
    source_spec <- sources$bundles$melanoma_bridge$sources[[1]]
    raw_file <- .download(
      source_spec$download_url,
      file.path(root, "raw", "gse72056", basename(source_spec$download_url))
    )
    header <- data.table::fread(raw_file, nrows = 3L, check.names = FALSE, showProgress = FALSE)
    if (nrow(header) != 3L || ncol(header) < 10L) stop("Unexpected GSE72056 table structure.", call. = FALSE)
    cell_names <- names(header)[-1L]
    tumor <- as.character(unlist(header[1L, -1L, with = FALSE], use.names = FALSE))
    malignant_code <- suppressWarnings(as.integer(unlist(header[2L, -1L, with = FALSE], use.names = FALSE)))
    nonmalignant_code <- suppressWarnings(as.integer(unlist(header[3L, -1L, with = FALSE], use.names = FALSE)))
    nonmalignant_labels <- c("1" = "T", "2" = "B", "3" = "Macrophage", "4" = "Endothelial", "5" = "CAF", "6" = "NK")
    malignant_call <- ifelse(malignant_code == 2L, "malignant", ifelse(malignant_code == 1L, "non_malignant", "unresolved"))
    cell_type <- ifelse(
      malignant_code == 2L, "Malignant",
      ifelse(malignant_code == 1L, unname(nonmalignant_labels[as.character(nonmalignant_code)]), "Unresolved")
    )
    cell_type[is.na(cell_type)] <- "Non_malignant_other"
    selected <- .stratified_target(
      seq_along(cell_names), paste(tumor, cell_type, sep = "::"),
      target = min(1500L, length(cell_names))
    )
    selected_names <- cell_names[selected]
    table <- data.table::fread(
      raw_file, select = c(names(header)[[1]], selected_names),
      check.names = FALSE, showProgress = TRUE
    )
    features <- as.character(table[[1]])
    metadata_rows <- features %in% as.character(header[[1]])
    expression_table <- table[!metadata_rows, -1L, with = FALSE]
    features <- features[!metadata_rows]
    expression <- as.matrix(expression_table)
    storage.mode(expression) <- "double"
    rownames(expression) <- features
    colnames(expression) <- selected_names
    variances <- apply(expression, 1L, stats::var)
    feature_order <- order(variances, decreasing = TRUE, na.last = NA)
    selected_feature_indices <- utils::head(feature_order, min(6000L, length(feature_order)))
    expression <- .as_sparse(expression[selected_feature_indices, , drop = FALSE])
    cell_data <- data.frame(
      tumor = tumor[selected],
      malignant_code = malignant_code[selected],
      malignant_call = malignant_call[selected],
      cell_type = cell_type[selected],
      row.names = selected_names,
      check.names = FALSE
    )
    assay <- SeuratObject::CreateAssay5Object(data = expression)
    object <- SeuratObject::CreateSeuratObject(counts = assay, meta.data = cell_data, project = "GSE72056")
    object@misc$real_data <- list(
      bundle = "melanoma_bridge",
      source = "GSE72056",
      source_expression = "author-provided normalized Smart-seq2 expression; stored only in the data layer",
      source_dimensions = c(features = length(features), cells = length(cell_names)),
      subset_seed = seed,
      license = source_spec$license,
      citation = source_spec$citation
    )
    .save_qs2(object, output)
    .write_manifest(
      output, "melanoma_bridge", "gse72056_single_cell", "GSE72056",
      summary = list(
        features = nrow(object), cells = ncol(object), tumors = length(unique(object$tumor)),
        cell_types = sort(unique(object$cell_type)), layers = SeuratObject::Layers(object[["RNA"]])
      ),
      selection = list(target_cells = 1500L, target_features = 6000L, stratification = "tumor x author cell class")
    )
  })
  invisible(output)
}

.clinical_survival <- function(clinical) {
  vital_column <- .find_column(clinical, c("vital_status"))
  death_column <- .find_column(clinical, c("days_to_death"), required = FALSE)
  follow_column <- .find_column(clinical, c("days_to_last_follow_up", "days_to_last_followup"), required = FALSE)
  death <- if (is.null(death_column)) rep(NA_real_, nrow(clinical)) else suppressWarnings(as.numeric(clinical[[death_column]]))
  follow <- if (is.null(follow_column)) rep(NA_real_, nrow(clinical)) else suppressWarnings(as.numeric(clinical[[follow_column]]))
  vital <- clinical[[vital_column]]
  if (is.numeric(vital) || is.integer(vital) || is.logical(vital)) {
    event <- as.integer(vital == 1)
  } else {
    event <- as.integer(tolower(as.character(vital)) %in% c("1", "dead", "deceased"))
  }
  time <- ifelse(event == 1L & is.finite(death), death, follow)
  list(time = time, event = event)
}

.prepare_tcga_skcm <- function() {
  output <- file.path(root, "bulk", "tcga_skcm.qs2")
  .existing_or_prepare(output, {
    .require(
      c("curatedTCGAData", "MultiAssayExperiment", "SummarizedExperiment", "Matrix"),
      "to prepare the curated TCGA-SKCM bridge"
    )
    set.seed(seed)
    message("Retrieving curatedTCGAData SKCM RNASeq2Gene assays (data version 2.1.1).")
    mae <- curatedTCGAData::curatedTCGAData(
      "SKCM", assays = c("RNASeq2Gene", "RNASeq2GeneNorm"),
      version = "2.1.1", dry.run = FALSE
    )
    experiments <- MultiAssayExperiment::experiments(mae)
    experiment_names <- names(experiments)
    raw_names <- experiment_names[
      grepl("RNASeq2Gene-", experiment_names, fixed = TRUE) &
        !grepl("RNASeq2GeneNorm", experiment_names, fixed = TRUE)
    ]
    normalized_names <- experiment_names[
      grepl("RNASeq2GeneNorm", experiment_names, fixed = TRUE)
    ]
    if (length(raw_names) != 1L || length(normalized_names) != 1L) {
      stop("Could not uniquely resolve RNASeq2Gene and RNASeq2GeneNorm experiments.", call. = FALSE)
    }
    raw_name <- raw_names[[1]]
    normalized_name <- normalized_names[[1]]
    raw_expression <- SummarizedExperiment::assay(experiments[[raw_name]], 1L)
    normalized_expression <- SummarizedExperiment::assay(experiments[[normalized_name]], 1L)

    sample_map <- as.data.frame(MultiAssayExperiment::sampleMap(mae), stringsAsFactors = FALSE)
    raw_map <- sample_map[sample_map$assay == raw_name, c("primary", "colname"), drop = FALSE]
    normalized_map <- sample_map[
      sample_map$assay == normalized_name, c("primary", "colname"), drop = FALSE
    ]
    names(raw_map)[names(raw_map) == "colname"] <- "raw_colname"
    names(normalized_map)[names(normalized_map) == "colname"] <- "normalized_colname"
    raw_map$sample_barcode <- substr(raw_map$raw_colname, 1L, 15L)
    normalized_map$sample_barcode <- substr(normalized_map$normalized_colname, 1L, 15L)
    paired <- merge(
      raw_map, normalized_map, by = c("primary", "sample_barcode"),
      all = FALSE, sort = TRUE
    )
    paired <- paired[
      order(paired$primary, paired$sample_barcode, paired$raw_colname, paired$normalized_colname),
      , drop = FALSE
    ]
    paired <- paired[!duplicated(paired$primary), , drop = FALSE]

    clinical <- as.data.frame(MultiAssayExperiment::colData(mae), stringsAsFactors = FALSE)
    survival <- .clinical_survival(clinical)
    clinical$real_time <- survival$time
    clinical$real_event <- survival$event
    group_column <- .find_column(clinical, c("ALL_PRIMARY_VS_METASTATIC"))
    paired$curated_group <- as.character(clinical[paired$primary, group_column])
    paired$time <- as.numeric(clinical[paired$primary, "real_time"])
    paired$event <- as.integer(clinical[paired$primary, "real_event"])
    paired <- paired[
      paired$curated_group %in% c("All_Primaries", "All_Metastases") &
        is.finite(paired$time) & paired$time > 0,
      , drop = FALSE
    ]
    selected_rows <- unlist(
      lapply(split(seq_len(nrow(paired)), paired$curated_group), function(current) {
        if (length(current) <= 90L) current else sample(current, 90L)
      }),
      use.names = FALSE
    )
    paired <- paired[sort(unique(selected_rows)), , drop = FALSE]
    if (nrow(paired) < 50L) {
      stop("Too few survival-complete primary/metastatic TCGA-SKCM patients.", call. = FALSE)
    }

    common_features <- intersect(rownames(raw_expression), rownames(normalized_expression))
    raw_for_selection <- raw_expression[
      common_features, paired$raw_colname, drop = FALSE
    ]
    normalized_for_selection <- normalized_expression[
      common_features, paired$normalized_colname, drop = FALSE
    ]
    complete_features <- rowSums(!is.finite(raw_for_selection)) == 0L &
      rowSums(!is.finite(normalized_for_selection)) == 0L
    common_features <- common_features[complete_features]
    normalized_for_selection <- normalized_for_selection[
      complete_features, , drop = FALSE
    ]
    if (length(common_features) < 3000L) {
      stop("Fewer than 3,000 TCGA-SKCM genes are finite in both selected assays.", call. = FALSE)
    }
    variances <- apply(normalized_for_selection, 1L, stats::var)
    feature_order <- names(sort(variances, decreasing = TRUE, na.last = NA))
    feature_order <- feature_order[!duplicated(feature_order)]
    melanoma_path <- file.path(root, "cancer", "melanoma_gse72056.qs2")
    if (file.exists(melanoma_path) && requireNamespace("qs2", quietly = TRUE)) {
      melanoma_features <- rownames(qs2::qs_read(melanoma_path))
      feature_order <- c(
        feature_order[feature_order %in% melanoma_features],
        feature_order[!feature_order %in% melanoma_features]
      )
    }
    selected_features <- utils::head(feature_order, min(3000L, length(feature_order)))
    rsem_tpm <- .as_sparse(raw_expression[
      selected_features, paired$raw_colname, drop = FALSE
    ])
    log2_uq_rsem_tpm <- normalized_expression[
      selected_features, paired$normalized_colname, drop = FALSE
    ]
    colnames(rsem_tpm) <- paired$primary
    colnames(log2_uq_rsem_tpm) <- paired$primary

    optional_clinical <- function(candidates) {
      column <- .find_column(clinical, candidates, required = FALSE)
      if (is.null(column)) rep(NA, nrow(paired)) else clinical[paired$primary, column]
    }
    sample_data <- data.frame(
      patient = paired$primary,
      sample = paired$sample_barcode,
      sample_type = unname(c(
        All_Primaries = "Primary Tumor", All_Metastases = "Metastatic"
      )[paired$curated_group]),
      curated_group = paired$curated_group,
      time = paired$time,
      event = paired$event,
      vital_status = clinical[paired$primary, "vital_status"],
      days_to_death = clinical[paired$primary, "days_to_death"],
      days_to_last_followup = clinical[paired$primary, "days_to_last_followup"],
      age_at_diagnosis = optional_clinical(c(
        "age_at_initial_pathologic_diagnosis", "patient.age_at_initial_pathologic_diagnosis"
      )),
      sex = optional_clinical(c("gender", "patient.gender")),
      stage = optional_clinical(c("pathologic_stage", "clinical_stage", "patient.stage_event.pathologic_stage")),
      raw_colname = paired$raw_colname,
      normalized_colname = paired$normalized_colname,
      row.names = paired$primary,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    artifact <- list(
      rsem_tpm = rsem_tpm,
      log2_uq_rsem_tpm = log2_uq_rsem_tpm,
      sample_data = sample_data,
      feature_data = data.frame(
        feature = rownames(rsem_tpm), row.names = rownames(rsem_tpm)
      ),
      provenance = list(
        bundle = "melanoma_bridge", source = "TCGA-SKCM",
        retrieval = "curatedTCGAData", data_version = "2.1.1",
        assays = c(rsem_tpm = raw_name, log2_uq_rsem_tpm = normalized_name),
        expression_semantics = c(
          rsem_tpm = "continuous RNASeq2Gene RSEM TPM; not raw counts",
          log2_uq_rsem_tpm = "continuous RNASeq2GeneNorm upper-quartile log2-normalized RSEM TPM; not raw counts"
        ),
        sample_alignment = "MultiAssayExperiment sampleMap by patient and 15-character sample barcode",
        open_access = TRUE, subset_seed = seed,
        citation = sources$bundles$melanoma_bridge$sources[[2]]$citation
      )
    )
    .save_qs2(artifact, output)
    .write_manifest(
      output, "melanoma_bridge", "tcga_skcm_bulk", "TCGA-SKCM",
      summary = list(
        features = nrow(rsem_tpm), samples = ncol(rsem_tpm),
        events = sum(sample_data$event == 1L, na.rm = TRUE),
        sample_types = sort(unique(sample_data$sample_type))
      ),
      selection = list(
        target_samples_per_type = 90L, retained_samples = nrow(sample_data),
        target_features = 3000L, retained_features = nrow(rsem_tpm),
        finite_shared_features = length(common_features),
        endpoint = "overall survival"
      )
    )
  })
  invisible(output)
}

.prepare_melanoma_bridge <- function() {
  .prepare_melanoma_single_cell()
  .prepare_tcga_skcm()
  invisible(c(
    file.path(root, "cancer", "melanoma_gse72056.qs2"),
    file.path(root, "bulk", "tcga_skcm.qs2")
  ))
}

.prepare_visium_section <- function(source_spec, sample_name, platform) {
  .require(c("Seurat", "SeuratObject"), "to prepare 10x Visium fixtures")
  raw_directory <- file.path(root, "raw", "visium", sample_name)
  dir.create(raw_directory, recursive = TRUE, showWarnings = FALSE)
  matrix_path <- .download(source_spec$matrix_url, file.path(raw_directory, basename(source_spec$matrix_url)))
  spatial_archive <- .download(source_spec$spatial_url, file.path(raw_directory, basename(source_spec$spatial_url)))
  spatial_directory <- file.path(raw_directory, "spatial")
  if (!dir.exists(spatial_directory) || replace_existing) {
    if (replace_existing && dir.exists(spatial_directory)) unlink(spatial_directory, recursive = TRUE)
    utils::untar(spatial_archive, exdir = raw_directory)
  }
  object <- Seurat::Load10X_Spatial(
    data.dir = raw_directory, filename = basename(matrix_path), assay = "Spatial", slice = sample_name
  )
  coordinates <- Seurat::GetTissueCoordinates(object)
  coordinate_columns <- intersect(c("x", "y", "imagerow", "imagecol"), names(coordinates))
  if (length(coordinate_columns) < 2L) {
    numeric_columns <- names(coordinates)[vapply(coordinates, is.numeric, logical(1))]
    coordinate_columns <- utils::head(numeric_columns, 2L)
  } else {
    coordinate_columns <- utils::tail(coordinate_columns, 2L)
  }
  if (length(coordinate_columns) != 2L) stop("Could not resolve Visium coordinates.", call. = FALSE)
  center <- vapply(coordinates[, coordinate_columns, drop = FALSE], stats::median, numeric(1), na.rm = TRUE)
  distance <- rowSums((sweep(as.matrix(coordinates[, coordinate_columns, drop = FALSE]), 2L, center, "-")) ^ 2)
  selected <- rownames(coordinates)[order(distance)][seq_len(min(400L, nrow(coordinates)))]
  object <- subset(object, cells = selected)
  object$real_section <- sample_name
  object$real_platform <- platform
  object@misc$real_data <- list(
    bundle = "visium_lymph_node", accession = source_spec$accession,
    subset_strategy = "400 tissue spots nearest spatial median", subset_seed = seed,
    license = source_spec$license, citation = source_spec$citation
  )
  object
}

.prepare_visium <- function() {
  output <- file.path(root, "spatial", "visium_lymph_node.qs2")
  .existing_or_prepare(output, {
    specs <- sources$bundles$visium_lymph_node$sources
    fresh <- .prepare_visium_section(specs[[1]], "fresh_frozen", "Visium whole transcriptome")
    ffpe <- .prepare_visium_section(specs[[2]], "ffpe_cytassist", "Visium CytAssist FFPE")
    artifact <- list(
      sections = list(fresh_frozen = fresh, ffpe_cytassist = ffpe),
      provenance = list(
        bundle = "visium_lymph_node", sources = vapply(specs, `[[`, character(1), "accession"),
        licenses = unique(vapply(specs, `[[`, character(1), "license")), subset_seed = seed
      )
    )
    .save_qs2(artifact, output)
    .write_manifest(
      output, "visium_lymph_node", "two_sections",
      vapply(specs, `[[`, character(1), "accession"),
      summary = list(
        sections = names(artifact$sections),
        cells = vapply(artifact$sections, ncol, numeric(1)),
        features = vapply(artifact$sections, nrow, numeric(1))
      ),
      selection = list(target_spots_per_section = 400L, strategy = "nearest to spatial median")
    )
  })
  invisible(output)
}

preparers <- list(
  kotliarov_pbmc = .prepare_kotliarov,
  melanoma_bridge = .prepare_melanoma_bridge,
  hermann_spermatogenesis = .prepare_hermann,
  visium_lymph_node = .prepare_visium
)

selected_bundles <- if (identical(bundle, "all")) bundle_names else bundle
message("Local real-data root: ", root)
for (current in selected_bundles) {
  message("\nPreparing bundle: ", current)
  preparers[[current]]()
}
message("\nReal-data preparation completed for: ", paste(selected_bundles, collapse = ", "))
