# Seurat-to-Python bridge.
#
# Extracted from package_tools.R: helpers that serialize Seurat objects for pixi Python backends, execute runner scripts, and import results (including infercnvpy).
.sn_run_python_object_method <- function(object,
                                         environment,
                                         script_name,
                                         method,
                                         assay = NULL,
                                         layer = NULL,
                                         reference_object = NULL,
                                         reference_assay = NULL,
                                         reference_layer = NULL,
                                         spatial_cols = NULL,
                                         output_dir = NULL,
                                         runtime_dir = NULL,
                                         metadata_prefix = paste0(method, "_"),
                                         result_name = method,
                                         return_object = TRUE,
                                         config = list(),
                                         ...) {
  check_installed(pkg = "Seurat", reason = glue::glue("to run {method} on a Seurat object."))
  assay <- assay %||% SeuratObject::DefaultAssay(object = object)
  layer <- layer %||% .sn_select_python_object_layer(object = object, assay = assay)
  output_dir <- output_dir %||% .sn_default_python_run_dir(method, runtime_dir = runtime_dir)
  input_dir <- file.path(output_dir, "input")
  result_dir <- file.path(output_dir, "output")
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  spatial_cols <- .sn_resolve_spatial_cols(object = object, spatial_cols = spatial_cols, required = environment %in% c("cell2location", "tangram", "squidpy", "spatialdata", "stlearn"))
  query <- .sn_write_python_object_input(
    object = object,
    input_dir = file.path(input_dir, "query"),
    assay = assay,
    layer = layer,
    spatial_cols = spatial_cols
  )

  reference <- NULL
  if (!is.null(reference_object)) {
    reference_assay <- reference_assay %||% SeuratObject::DefaultAssay(object = reference_object)
    reference_layer <- reference_layer %||% .sn_select_python_object_layer(object = reference_object, assay = reference_assay)
    reference <- .sn_write_python_object_input(
      object = reference_object,
      input_dir = file.path(input_dir, "reference"),
      assay = reference_assay,
      layer = reference_layer,
      spatial_cols = NULL
    )
  }

  config <- .sn_prepare_python_object_config(config = config, input_dir = input_dir)
  config <- c(
    list(
      method = method,
      assay = assay,
      layer = layer,
      reference_assay = reference_assay,
      reference_layer = reference_layer
    ),
    config
  )
  config_path <- .sn_write_json_file(config, file.path(output_dir, paste0(method, "_config.json")))
  script <- .sn_pixi_script_path(environment = environment, script_name = script_name)
  .sn_execute_python_object_pixi(
    environment = environment,
    script = script,
    input_dir = input_dir,
    output_dir = result_dir,
    config_path = config_path,
    ...
  )
  .sn_import_python_object_results(
    object = object,
    method = method,
    result_name = result_name,
    output_dir = result_dir,
    run_dir = output_dir,
    assay = assay,
    query = query,
    reference = reference,
    config = config,
    metadata_prefix = metadata_prefix,
    return_object = return_object
  )
}

.sn_select_python_object_layer <- function(object, assay) {
  layers <- SeuratObject::Layers(object[[assay]])
  if ("data" %in% layers || any(grepl("^data\\.", layers))) {
    return("data")
  }
  "counts"
}

.sn_pixi_script_path <- function(environment, script_name) {
  installed <- system.file("pixi", environment, "scripts", script_name, package = "Shennong")
  if (nzchar(installed) && file.exists(installed)) {
    return(normalizePath(installed, winslash = "/", mustWork = TRUE))
  }
  source_path <- file.path(getwd(), "inst", "pixi", environment, "scripts", script_name)
  if (file.exists(source_path)) {
    return(normalizePath(source_path, winslash = "/", mustWork = TRUE))
  }
  stop("Could not locate Shennong pixi runner script: ", file.path(environment, "scripts", script_name), call. = FALSE)
}

.sn_write_python_object_input <- function(object,
                                          input_dir,
                                          assay,
                                          layer,
                                          spatial_cols = NULL) {
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  expr <- .sn_get_seurat_layer_data(object = object, assay = assay, layer = layer)
  expr <- .sn_as_sparse_matrix(expr)
  expr <- expr[, colnames(object), drop = FALSE]

  matrix_path <- file.path(input_dir, "matrix.mtx")
  obs_path <- file.path(input_dir, "obs.csv")
  var_path <- file.path(input_dir, "var.csv")
  Matrix::writeMM(obj = expr, file = matrix_path)

  obs <- object[[]][colnames(expr), , drop = FALSE]
  obs <- data.frame(cell_id = rownames(obs), obs, check.names = FALSE)
  utils::write.csv(obs, file = obs_path, row.names = FALSE, quote = TRUE)
  var <- data.frame(feature_id = rownames(expr), stringsAsFactors = FALSE)
  utils::write.csv(var, file = var_path, row.names = FALSE, quote = TRUE)

  spatial_path <- NULL
  if (!is.null(spatial_cols)) {
    spatial <- obs[, spatial_cols, drop = FALSE]
    rownames(spatial) <- obs$cell_id
    spatial_path <- file.path(input_dir, "spatial.csv")
    utils::write.csv(spatial, file = spatial_path, row.names = TRUE, quote = TRUE)
  }

  list(
    input_dir = normalizePath(input_dir, winslash = "/", mustWork = TRUE),
    matrix_path = normalizePath(matrix_path, winslash = "/", mustWork = TRUE),
    obs_path = normalizePath(obs_path, winslash = "/", mustWork = TRUE),
    var_path = normalizePath(var_path, winslash = "/", mustWork = TRUE),
    spatial_path = if (!is.null(spatial_path)) normalizePath(spatial_path, winslash = "/", mustWork = TRUE) else NULL,
    assay = assay,
    layer = layer,
    features = rownames(expr),
    n_features = nrow(expr),
    n_cells = ncol(expr)
  )
}

.sn_resolve_spatial_cols <- function(object, spatial_cols = NULL, required = FALSE) {
  metadata_cols <- colnames(object[[]])
  if (!is.null(spatial_cols)) {
    if (length(spatial_cols) != 2L || !all(spatial_cols %in% metadata_cols)) {
      stop("`spatial_cols` must contain two metadata columns present in `object`.", call. = FALSE)
    }
    return(spatial_cols)
  }
  candidates <- list(
    c("x", "y"),
    c("spatial_x", "spatial_y"),
    c("imagecol", "imagerow"),
    c("pxl_col_in_fullres", "pxl_row_in_fullres"),
    c("array_col", "array_row")
  )
  for (candidate in candidates) {
    if (all(candidate %in% metadata_cols)) {
      return(candidate)
    }
  }
  if (isTRUE(required)) {
    stop("Spatial object workflow requires `spatial_cols` or recognized coordinate metadata columns.", call. = FALSE)
  }
  NULL
}

.sn_prepare_python_object_config <- function(config, input_dir) {
  if (!is.null(config$reference_signatures)) {
    signatures <- config$reference_signatures
    if (is.data.frame(signatures) || is.matrix(signatures)) {
      signatures_path <- file.path(input_dir, "reference_signatures.csv")
      utils::write.csv(as.data.frame(signatures), file = signatures_path, row.names = TRUE, quote = TRUE)
      config$reference_signatures <- normalizePath(signatures_path, winslash = "/", mustWork = TRUE)
    } else if (is.character(signatures) && length(signatures) == 1L && file.exists(path.expand(signatures))) {
      config$reference_signatures <- normalizePath(path.expand(signatures), winslash = "/", mustWork = TRUE)
    }
  }
  config
}

.sn_execute_python_object_pixi <- function(environment,
                                          script,
                                          input_dir,
                                          output_dir,
                                          config_path,
                                          ...) {
  sn_call_pixi_environment(
    environment = environment,
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

.sn_import_python_object_results <- function(object,
                                            method,
                                            result_name,
                                            output_dir,
                                            run_dir,
                                            assay,
                                            query,
                                            reference = NULL,
                                            config = list(),
                                            metadata_prefix = paste0(method, "_"),
                                            return_object = TRUE) {
  obs_path <- file.path(output_dir, "obs.csv")
  imported_metadata <- character(0)
  if (file.exists(obs_path)) {
    metadata <- utils::read.csv(obs_path, row.names = 1, check.names = FALSE)
    shared_cells <- intersect(colnames(object), rownames(metadata))
    if (length(shared_cells) > 0L && ncol(metadata) > 0L) {
      metadata <- metadata[shared_cells, , drop = FALSE]
      colnames(metadata) <- .sn_prefix_metadata_columns(colnames(metadata), metadata_prefix)
      imported_metadata <- colnames(metadata)
      object <- Seurat::AddMetaData(object = object, metadata = metadata)
    }
  }

  reductions <- character(0)
  embedding_files <- list.files(output_dir, pattern = "\\.csv$", full.names = TRUE)
  embedding_files <- embedding_files[grepl("(latent|pca|umap|embedding)", basename(embedding_files), ignore.case = TRUE)]
  for (embedding_file in embedding_files) {
    embedding <- tryCatch(.sn_read_embedding_csv(embedding_file, cells = colnames(object)), error = function(e) NULL)
    if (is.null(embedding)) {
      next
    }
    reduction_name <- paste0(metadata_prefix, tools::file_path_sans_ext(basename(embedding_file)))
    reduction_key <- paste0(gsub("[^A-Za-z0-9]", "", toupper(reduction_name)), "_")
    colnames(embedding) <- paste0(reduction_key, seq_len(ncol(embedding)))
    object[[reduction_name]] <- Seurat::CreateDimReducObject(
      embeddings = embedding,
      key = reduction_key,
      assay = assay
    )
    reductions <- c(reductions, reduction_name)
  }

  manifest_path <- file.path(output_dir, "manifest.json")
  backend_manifest <- if (file.exists(manifest_path)) {
    jsonlite::read_json(manifest_path, simplifyVector = TRUE)
  } else {
    list()
  }
  run_manifest <- c(
    list(
      method = method,
      result_name = result_name,
      assay = assay,
      source_layer = query$layer,
      run_dir = normalizePath(run_dir, winslash = "/", mustWork = TRUE),
      output_dir = normalizePath(output_dir, winslash = "/", mustWork = TRUE),
      input_features = query$features,
      n_features = query$n_features,
      n_cells = query$n_cells,
      reference = reference,
      imported_metadata = imported_metadata,
      imported_reductions = reductions,
      config = config
    ),
    backend_manifest
  )
  object@misc[[method]] <- object@misc[[method]] %||% list()
  object@misc[[method]][[result_name]] <- run_manifest
  object <- .sn_log_seurat_command(object = object, name = paste0("sn_run_", method))
  if (isTRUE(return_object)) {
    return(object)
  }
  run_manifest
}

.sn_run_infercnvpy_object <- function(object,
                                      assay = NULL,
                                      layer = NULL,
                                      species = NULL,
                                      reference_key = NULL,
                                      reference_cat = NULL,
                                      gene_order = NULL,
                                      gtf_file = NULL,
                                      gtf_gene_id = "gene_name",
                                      adata_gene_id = NULL,
                                      output_dir = NULL,
                                      runtime_dir = NULL,
                                      key_added = "cnv",
                                      window_size = 100,
                                      step = 10,
                                      dynamic_threshold = 1.5,
                                      exclude_chromosomes = c("chrX", "chrY"),
                                      chunksize = 5000,
                                      n_jobs = NULL,
                                      calculate_gene_values = FALSE,
                                      lfc_clip = 3,
                                      run_pca = TRUE,
                                      run_neighbors = TRUE,
                                      run_leiden = TRUE,
                                      run_umap = FALSE,
                                      score = TRUE,
                                      leiden_resolution = 1,
                                      cnv_score_group_by = NULL,
                                      metadata_prefix = "infercnvpy_",
                                      result_name = "infercnvpy",
                                      return_object = TRUE,
                                      ...) {
  check_installed(pkg = "Seurat", reason = "to run infercnvpy on a Seurat object.")
  assay <- assay %||% SeuratObject::DefaultAssay(object = object)
  layer <- layer %||% .sn_select_infercnvpy_layer(object = object, assay = assay)
  if (!is.null(reference_key) && (!nzchar(reference_key) || !reference_key %in% colnames(object[[]]))) {
    stop("`reference_key` must name a metadata column in `object`.", call. = FALSE)
  }

  output_dir <- output_dir %||% .sn_default_python_run_dir("infercnvpy", runtime_dir = runtime_dir)
  input_dir <- file.path(output_dir, "input")
  result_dir <- file.path(output_dir, "output")
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

  input <- .sn_write_infercnvpy_input(
    object = object,
    input_dir = input_dir,
    assay = assay,
    layer = layer,
    species = species,
    gene_order = gene_order,
    gtf_file = gtf_file
  )
  script <- .sn_infercnvpy_script_path()
  config <- list(
    method = "infercnvpy",
    assay = assay,
    layer = NULL,
    source_layer = layer,
    reference_key = reference_key,
    reference_cat = reference_cat,
    gtf_file = if (!is.null(gtf_file)) normalizePath(path.expand(gtf_file), winslash = "/", mustWork = TRUE) else NULL,
    gtf_gene_id = gtf_gene_id,
    adata_gene_id = adata_gene_id,
    key_added = key_added,
    window_size = window_size,
    step = step,
    dynamic_threshold = dynamic_threshold,
    exclude_chromosomes = exclude_chromosomes,
    chunksize = chunksize,
    n_jobs = n_jobs,
    calculate_gene_values = isTRUE(calculate_gene_values),
    lfc_clip = lfc_clip,
    run_pca = isTRUE(run_pca),
    run_neighbors = isTRUE(run_neighbors),
    run_leiden = isTRUE(run_leiden),
    run_umap = isTRUE(run_umap),
    score = isTRUE(score),
    leiden_resolution = leiden_resolution,
    cnv_score_groupby = cnv_score_group_by,
    write_h5ad = TRUE
  )
  config_path <- .sn_write_json_file(config, file.path(output_dir, "infercnvpy_config.json"))

  .sn_execute_infercnvpy_pixi(
    script = script,
    input_dir = input$input_dir,
    output_dir = result_dir,
    config_path = config_path,
    ...
  )

  manifest <- .sn_import_infercnvpy_results(
    object = object,
    output_dir = result_dir,
    run_dir = output_dir,
    assay = assay,
    input = input,
    config = config,
    metadata_prefix = metadata_prefix,
    result_name = result_name,
    return_object = return_object
  )
  manifest
}

.sn_select_infercnvpy_layer <- function(object, assay) {
  layers <- SeuratObject::Layers(object[[assay]])
  if ("data" %in% layers || any(grepl("^data\\.", layers))) {
    return("data")
  }
  .sn_log_warn(
    "No normalized `data` layer was found for assay '{assay}'; using `counts`. ",
    "infercnvpy expects normalized log-transformed expression."
  )
  "counts"
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
                                       gtf_file = NULL) {
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

  obs <- object[[]][colnames(expr), , drop = FALSE]
  obs <- data.frame(cell_id = rownames(obs), obs, check.names = FALSE)
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
  sn_call_infercnvpy(
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
                                          return_object = TRUE) {
  obs_path <- file.path(output_dir, "obs.csv")
  manifest_path <- file.path(output_dir, "manifest.json")
  if (!file.exists(obs_path)) {
    stop("infercnvpy output is missing `obs.csv`: ", obs_path, call. = FALSE)
  }

  metadata <- utils::read.csv(obs_path, row.names = 1, check.names = FALSE)
  shared_cells <- intersect(colnames(object), rownames(metadata))
  if (length(shared_cells) == 0L) {
    stop("infercnvpy metadata output does not contain cells from the input object.", call. = FALSE)
  }
  metadata <- metadata[shared_cells, , drop = FALSE]
  colnames(metadata) <- .sn_prefix_metadata_columns(colnames(metadata), metadata_prefix)
  object <- Seurat::AddMetaData(object = object, metadata = metadata)

  pca_path <- file.path(output_dir, "cnv_pca.csv")
  umap_path <- file.path(output_dir, "cnv_umap.csv")
  reductions <- character(0)
  if (file.exists(pca_path)) {
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

  backend_manifest <- if (file.exists(manifest_path)) {
    jsonlite::read_json(manifest_path, simplifyVector = TRUE)
  } else {
    list()
  }
  run_manifest <- c(
    list(
      method = "infercnvpy",
      result_name = result_name,
      assay = assay,
      source_layer = input$layer,
      run_dir = normalizePath(run_dir, winslash = "/", mustWork = TRUE),
      output_dir = normalizePath(output_dir, winslash = "/", mustWork = TRUE),
      input_features = input$features,
      n_features = input$n_features,
      n_cells = input$n_cells,
      imported_metadata = colnames(metadata),
      imported_reductions = reductions,
      config = config
    ),
    backend_manifest
  )

  object@misc$infercnvpy <- object@misc$infercnvpy %||% list()
  object@misc$infercnvpy[[result_name]] <- run_manifest
  object <- .sn_log_seurat_command(object = object, name = "sn_run_infercnvpy")
  if (isTRUE(return_object)) {
    return(object)
  }
  run_manifest
}

.sn_prefix_metadata_columns <- function(columns, prefix) {
  if (is.null(prefix) || !nzchar(prefix)) {
    return(columns)
  }
  ifelse(startsWith(columns, prefix), columns, paste0(prefix, columns))
}

.sn_read_embedding_csv <- function(path, cells) {
  embedding <- utils::read.csv(path, row.names = 1, check.names = FALSE)
  missing_cells <- setdiff(cells, rownames(embedding))
  if (length(missing_cells) > 0L) {
    stop("Embedding output is missing cells from the input object: ", path, call. = FALSE)
  }
  embedding <- as.matrix(embedding[cells, , drop = FALSE])
  storage.mode(embedding) <- "numeric"
  embedding
}
