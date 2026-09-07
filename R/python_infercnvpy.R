# inferCNVpy Python backend helpers.
#
# This module owns inferCNVpy control validation, input and configuration
# preparation, managed execution, and result import.

.sn_validate_infercnvpy_controls <- function(object,
                                              reference_key,
                                              reference_cat,
                                              run_pca,
                                              run_neighbors,
                                              run_leiden,
                                              run_umap,
                                              score,
                                              cnv_score_group_by) {
  controls <- list(
    run_pca = run_pca,
    run_neighbors = run_neighbors,
    run_leiden = run_leiden,
    run_umap = run_umap,
    score = score
  )
  invalid <- names(controls)[!vapply(
    controls,
    function(value) is.logical(value) && length(value) == 1L && !is.na(value),
    logical(1)
  )]
  if (length(invalid) > 0L) {
    stop("infercnvpy control(s) must be TRUE or FALSE: ", paste(invalid, collapse = ", "), call. = FALSE)
  }
  if (isTRUE(run_neighbors) && !isTRUE(run_pca)) {
    stop("infercnvpy `run_neighbors = TRUE` requires `run_pca = TRUE`.", call. = FALSE)
  }
  if (isTRUE(run_leiden) && !isTRUE(run_neighbors)) {
    stop("infercnvpy `run_leiden = TRUE` requires `run_neighbors = TRUE`.", call. = FALSE)
  }
  if (isTRUE(run_umap) && !isTRUE(run_neighbors)) {
    stop("infercnvpy `run_umap = TRUE` requires `run_neighbors = TRUE`.", call. = FALSE)
  }
  if (isTRUE(score) && is.null(cnv_score_group_by) && !isTRUE(run_leiden)) {
    stop(
      "infercnvpy `score = TRUE` requires `run_leiden = TRUE` unless `cnv_score_group_by` is supplied.",
      call. = FALSE
    )
  }
  if (!is.null(cnv_score_group_by) &&
      (!is.character(cnv_score_group_by) || length(cnv_score_group_by) != 1L ||
        is.na(cnv_score_group_by) || !nzchar(cnv_score_group_by) ||
        !cnv_score_group_by %in% colnames(object[[]]))) {
    stop("`cnv_score_group_by` must name one metadata column in `object`.", call. = FALSE)
  }
  if (!is.null(reference_cat) && is.null(reference_key)) {
    stop("`reference_cat` requires `reference_key`/`reference_by`.", call. = FALSE)
  }
  invisible(TRUE)
}

.sn_select_infercnvpy_layer <- function(object, assay) {
  layers <- SeuratObject::Layers(object[[assay]])
  if ("data" %in% layers || any(grepl("^data\\.", layers))) {
    return("data")
  }
  stop(
    "infercnvpy requires normalized, log-transformed expression, but assay '",
    assay, "' has no `data`/`data.*` layer. Run normalization first or select ",
    "an assay with a normalized `data` layer; raw counts are not used as a fallback.",
    call. = FALSE
  )
}

.sn_infercnvpy_script_path <- function(script = NULL) {
  if (!is.null(script) && nzchar(script)) {
    script <- path.expand(script)
    if (!file.exists(script)) {
      stop("infercnvpy runner script does not exist: ", script, call. = FALSE)
    }
    return(normalizePath(script, winslash = "/", mustWork = TRUE))
  }

  installed <- system.file("pixi", "infercnvpy", "scripts", "infercnvpy_run.py", package = "Shennong")
  if (nzchar(installed) && file.exists(installed)) {
    return(normalizePath(installed, winslash = "/", mustWork = TRUE))
  }
  source_path <- file.path(getwd(), "inst", "pixi", "infercnvpy", "scripts", "infercnvpy_run.py")
  if (file.exists(source_path)) {
    return(normalizePath(source_path, winslash = "/", mustWork = TRUE))
  }
  stop("Could not locate Shennong's infercnvpy Python runner.", call. = FALSE)
}

.sn_write_infercnvpy_input <- function(object,
                                       input_dir,
                                       assay,
                                       layer,
                                       species = NULL,
                                       gene_order = NULL,
                                       gtf_file = NULL,
                                       metadata_columns = NULL) {
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  expr <- .sn_get_seurat_layer_data(object = object, assay = assay, layer = layer)
  expr <- .sn_as_sparse_matrix(expr)
  expr <- expr[, colnames(object), drop = FALSE]

  gene_positions <- .sn_prepare_infercnvpy_gene_positions(
    features = rownames(expr),
    object = object,
    species = species,
    gene_order = gene_order,
    gtf_file = gtf_file
  )
  keep_features <- gene_positions$feature_id
  expr <- expr[keep_features, , drop = FALSE]
  gene_positions <- gene_positions[match(rownames(expr), gene_positions$feature_id), , drop = FALSE]

  matrix_path <- file.path(input_dir, "matrix.mtx")
  obs_path <- file.path(input_dir, "obs.csv")
  var_path <- file.path(input_dir, "var.csv")
  Matrix::writeMM(obj = expr, file = matrix_path)

  full_metadata <- object[[]][colnames(expr), , drop = FALSE]
  metadata_columns <- unique(as.character(unlist(metadata_columns, use.names = FALSE)))
  metadata_columns <- metadata_columns[!is.na(metadata_columns) & nzchar(metadata_columns)]
  missing_metadata <- setdiff(metadata_columns, colnames(full_metadata))
  if (length(missing_metadata) > 0L) {
    stop(
      "Required infercnvpy metadata column(s) were not found: ",
      paste(missing_metadata, collapse = ", "),
      call. = FALSE
    )
  }
  obs <- full_metadata[, metadata_columns, drop = FALSE]
  obs <- data.frame(cell_id = rownames(full_metadata), obs, check.names = FALSE)
  utils::write.csv(obs, file = obs_path, row.names = FALSE, quote = TRUE)

  utils::write.csv(gene_positions, file = var_path, row.names = FALSE, quote = TRUE)
  list(
    input_dir = normalizePath(input_dir, winslash = "/", mustWork = TRUE),
    matrix_path = normalizePath(matrix_path, winslash = "/", mustWork = TRUE),
    obs_path = normalizePath(obs_path, winslash = "/", mustWork = TRUE),
    var_path = normalizePath(var_path, winslash = "/", mustWork = TRUE),
    assay = assay,
    layer = layer,
    features = rownames(expr),
    cells = colnames(expr),
    metadata_columns = metadata_columns,
    n_features = nrow(expr),
    n_cells = ncol(expr)
  )
}

.sn_prepare_infercnvpy_gene_positions <- function(features,
                                                  object,
                                                  species = NULL,
                                                  gene_order = NULL,
                                                  gtf_file = NULL) {
  features <- as.character(features)
  if (!is.null(gene_order)) {
    positions <- .sn_match_user_gene_order(features = features, gene_order = gene_order)
  } else if (!is.null(gtf_file)) {
    if (!file.exists(path.expand(gtf_file))) {
      stop("`gtf_file` does not exist: ", gtf_file, call. = FALSE)
    }
    positions <- data.frame(feature_id = features, stringsAsFactors = FALSE)
  } else {
    species <- sn_get_species(object = object, species = species)
    annotations <- .sn_get_gene_annotation_table(species = species)
    positions <- .sn_match_annotation_gene_order(features = features, annotations = annotations)
  }

  if (!all(c("feature_id") %in% colnames(positions))) {
    stop("Internal error: infercnvpy gene positions must include `feature_id`.", call. = FALSE)
  }
  if (is.null(gtf_file) && !all(c("chromosome", "start", "end") %in% colnames(positions))) {
    stop("Gene positions must include `chromosome`, `start`, and `end`.", call. = FALSE)
  }
  positions <- positions[positions$feature_id %in% features, , drop = FALSE]
  positions <- positions[!duplicated(positions$feature_id), , drop = FALSE]
  positions <- positions[match(intersect(features, positions$feature_id), positions$feature_id), , drop = FALSE]
  if (nrow(positions) == 0L) {
    stop("No input genes could be matched to genomic positions for infercnvpy.", call. = FALSE)
  }
  missing <- setdiff(features, positions$feature_id)
  if (length(missing) > 0L) {
    .sn_log_warn("Dropping {length(missing)} gene(s) without genomic positions before infercnvpy.")
  }
  positions
}

.sn_match_annotation_gene_order <- function(features, annotations) {
  feature_base <- sub("\\..*$", "", features)
  by_name <- annotations[!is.na(annotations$gene_name) & nzchar(annotations$gene_name), , drop = FALSE]
  by_name <- by_name[!duplicated(by_name$gene_name), , drop = FALSE]
  rownames(by_name) <- by_name$gene_name
  by_id <- annotations[!is.na(annotations$gene_id) & nzchar(annotations$gene_id), , drop = FALSE]
  by_id <- by_id[!duplicated(by_id$gene_id), , drop = FALSE]
  rownames(by_id) <- by_id$gene_id
  by_base <- annotations[!is.na(annotations$gene_id_base) & nzchar(annotations$gene_id_base), , drop = FALSE]
  by_base <- by_base[!duplicated(by_base$gene_id_base), , drop = FALSE]
  rownames(by_base) <- by_base$gene_id_base

  matched <- vector("list", length(features))
  for (i in seq_along(features)) {
    feature <- features[[i]]
    row <- NULL
    matched_by <- NA_character_
    if (feature %in% rownames(by_name)) {
      row <- by_name[feature, , drop = FALSE]
      matched_by <- "gene_name"
    } else if (feature %in% rownames(by_id)) {
      row <- by_id[feature, , drop = FALSE]
      matched_by <- "gene_id"
    } else if (feature_base[[i]] %in% rownames(by_base)) {
      row <- by_base[feature_base[[i]], , drop = FALSE]
      matched_by <- "gene_id_base"
    }
    if (!is.null(row)) {
      matched[[i]] <- data.frame(
        feature_id = feature,
        chromosome = row$seqname,
        start = row$start,
        end = row$end,
        gene_id = row$gene_id,
        gene_id_base = row$gene_id_base,
        gene_name = row$gene_name,
        matched_by = matched_by,
        stringsAsFactors = FALSE
      )
    }
  }
  matched <- matched[!vapply(matched, is.null, logical(1))]
  if (length(matched) == 0L) {
    return(data.frame(feature_id = character(0), stringsAsFactors = FALSE))
  }
  do.call(rbind, matched)
}

.sn_match_user_gene_order <- function(features, gene_order) {
  gene_order <- as.data.frame(gene_order, stringsAsFactors = FALSE)
  id_col <- intersect(c("feature_id", "feature", "gene", "gene_name", "gene_id"), colnames(gene_order))
  if (length(id_col) == 0L) {
    stop(
      "`gene_order` must contain one of: feature_id, feature, gene, gene_name, gene_id.",
      call. = FALSE
    )
  }
  id_col <- id_col[[1]]
  chr_col <- intersect(c("chromosome", "seqname", "chr"), colnames(gene_order))
  if (length(chr_col) == 0L || !"start" %in% colnames(gene_order) || !"end" %in% colnames(gene_order)) {
    stop("`gene_order` must contain chromosome/seqname, start, and end columns.", call. = FALSE)
  }
  chr_col <- chr_col[[1]]
  gene_order <- gene_order[!duplicated(gene_order[[id_col]]), , drop = FALSE]
  rownames(gene_order) <- as.character(gene_order[[id_col]])
  keep <- features[features %in% rownames(gene_order)]
  out <- gene_order[keep, , drop = FALSE]
  data.frame(
    feature_id = keep,
    chromosome = out[[chr_col]],
    start = out[["start"]],
    end = out[["end"]],
    stringsAsFactors = FALSE
  )
}

.sn_execute_infercnvpy_pixi <- function(script,
                                        input_dir,
                                        output_dir,
                                        config_path,
                                        ...) {
  sn_call_pixi_environment(
    "infercnvpy",
    command = "python",
    args = c(
      shQuote(script),
      "--input-dir", shQuote(input_dir),
      "--output-dir", shQuote(output_dir),
      "--config", shQuote(config_path)
    ),
    ...
  )
}

.sn_import_infercnvpy_results <- function(object,
                                          output_dir,
                                          run_dir,
                                          assay,
                                          input,
                                          config,
                                          metadata_prefix = "infercnvpy_",
                                          result_name = "infercnvpy",
                                          return_object = TRUE,
                                          retain_run_dir = TRUE,
                                          max_artifact_import_gb = 0.5) {
  obs_path <- file.path(output_dir, "obs.csv")
  manifest_path <- file.path(output_dir, "manifest.json")
  if (!file.exists(obs_path)) {
    stop("infercnvpy output is missing `obs.csv`: ", obs_path, call. = FALSE)
  }
  if (!file.exists(manifest_path)) {
    stop("infercnvpy output is missing required `manifest.json`: ", manifest_path, call. = FALSE)
  }

  .sn_assert_python_artifact_budget(obs_path, max_artifact_import_gb, "infercnvpy metadata output")
  .sn_assert_python_artifact_budget(manifest_path, max_artifact_import_gb, "infercnvpy manifest")
  metadata <- utils::read.csv(obs_path, row.names = 1, check.names = FALSE)
  .sn_validate_exact_python_ids(
    rownames(metadata), input$cells,
    "infercnvpy metadata output cell"
  )
  if (ncol(metadata) == 0L) {
    stop("infercnvpy metadata output contains no result columns.", call. = FALSE)
  }
  metadata <- metadata[colnames(object), , drop = FALSE]
  colnames(metadata) <- .sn_prefix_metadata_columns(colnames(metadata), metadata_prefix)
  object <- Seurat::AddMetaData(object = object, metadata = metadata)

  pca_path <- file.path(output_dir, "cnv_pca.csv")
  umap_path <- file.path(output_dir, "cnv_umap.csv")
  reductions <- character(0)
  if (file.exists(pca_path)) {
    .sn_assert_python_artifact_budget(pca_path, max_artifact_import_gb, "infercnvpy PCA output")
    pca <- .sn_read_embedding_csv(pca_path, cells = colnames(object))
    reduction_name <- paste0(metadata_prefix, "cnv_pca")
    colnames(pca) <- paste0("CNVPCA_", seq_len(ncol(pca)))
    object[[reduction_name]] <- Seurat::CreateDimReducObject(
      embeddings = pca,
      key = "CNVPCA_",
      assay = assay
    )
    reductions <- c(reductions, reduction_name)
  }
  if (file.exists(umap_path)) {
    .sn_assert_python_artifact_budget(umap_path, max_artifact_import_gb, "infercnvpy UMAP output")
    umap <- .sn_read_embedding_csv(umap_path, cells = colnames(object))
    reduction_name <- paste0(metadata_prefix, "cnv_umap")
    colnames(umap) <- paste0("CNVUMAP_", seq_len(ncol(umap)))
    object[[reduction_name]] <- Seurat::CreateDimReducObject(
      embeddings = umap,
      key = "CNVUMAP_",
      assay = assay
    )
    reductions <- c(reductions, reduction_name)
  }

  chromosome_path <- file.path(output_dir, "cnv_chromosome.csv")
  if (!file.exists(chromosome_path)) {
    stop("infercnvpy output is missing `cnv_chromosome.csv`: ", chromosome_path, call. = FALSE)
  }
  .sn_assert_python_artifact_budget(
    chromosome_path,
    max_artifact_import_gb,
    "infercnvpy chromosome output"
  )
  chromosome <- utils::read.csv(chromosome_path, row.names = 1, check.names = FALSE)
  .sn_validate_exact_python_ids(
    rownames(chromosome), input$cells,
    "infercnvpy chromosome output cell"
  )
  chromosome <- as.matrix(chromosome[input$cells, , drop = FALSE])
  storage.mode(chromosome) <- "numeric"
  if (ncol(chromosome) < 1L || any(!is.finite(chromosome))) {
    stop("infercnvpy chromosome output must contain finite numeric values.", call. = FALSE)
  }

  backend_manifest <- jsonlite::read_json(manifest_path, simplifyVector = TRUE)
  if (!is.list(backend_manifest) ||
      !identical(tolower(as.character(backend_manifest$method %||% "")), "infercnvpy")) {
    stop("infercnvpy manifest must declare `method: infercnvpy`.", call. = FALSE)
  }
  if (is.null(names(backend_manifest)) || anyNA(names(backend_manifest)) ||
      any(!nzchar(names(backend_manifest))) || anyDuplicated(names(backend_manifest))) {
    stop("infercnvpy manifest fields must have unique, non-empty names.", call. = FALSE)
  }
  manifest_cells <- suppressWarnings(as.integer(backend_manifest$n_cells %||% NA_integer_))
  if (length(manifest_cells) != 1L || is.na(manifest_cells) || manifest_cells != input$n_cells) {
    stop(
      "infercnvpy manifest must report `n_cells = ", input$n_cells,
      "` for the current input.",
      call. = FALSE
    )
  }
  manifest_input_features <- suppressWarnings(as.integer(backend_manifest$input_n_features %||% NA_integer_))
  if (length(manifest_input_features) != 1L || is.na(manifest_input_features) ||
      manifest_input_features != input$n_features) {
    stop(
      "infercnvpy manifest must report `input_n_features = ", input$n_features,
      "` for the current input.",
      call. = FALSE
    )
  }
  manifest_features <- suppressWarnings(as.integer(backend_manifest$n_features %||% NA_integer_))
  if (length(manifest_features) != 1L || is.na(manifest_features) ||
      manifest_features < 1L || manifest_features > input$n_features) {
    stop("infercnvpy manifest reports an invalid retained `n_features` value.", call. = FALSE)
  }
  run_manifest <- list(
    method = "infercnvpy",
    result_name = result_name,
    assay = assay,
    source_layer = input$layer,
    run_dir = if (isTRUE(retain_run_dir)) normalizePath(run_dir, winslash = "/", mustWork = TRUE) else NULL,
    output_dir = if (isTRUE(retain_run_dir)) normalizePath(output_dir, winslash = "/", mustWork = TRUE) else NULL,
    run_dir_retained = isTRUE(retain_run_dir),
    input_features = input$features,
    n_features = input$n_features,
    n_retained_features = manifest_features,
    n_cells = input$n_cells,
    imported_metadata = colnames(metadata),
    imported_reductions = reductions,
    config = .sn_safe_python_config(config)
  )
  safe_backend_manifest <- .sn_sanitize_python_backend_manifest(
    backend_manifest,
    retain_paths = isTRUE(retain_run_dir)
  )
  safe_backend_manifest <- safe_backend_manifest[
    !names(safe_backend_manifest) %in% names(run_manifest)
  ]
  run_manifest <- c(run_manifest, safe_backend_manifest)
  if (anyDuplicated(names(run_manifest))) {
    stop("Internal error: infercnvpy run manifest contains duplicate field names.", call. = FALSE)
  }

  object@misc$infercnvpy <- object@misc$infercnvpy %||% list()
  object@misc$infercnvpy[[result_name]] <- run_manifest
  object <- .sn_log_seurat_command(object = object, name = "sn_run_infercnvpy")
  if (isTRUE(return_object)) {
    return(object)
  }
  run_manifest
}
