# Internal helper to normalize the clustering dimension arguments.
.sn_resolve_cluster_dims <- function(dims = NULL, npcs = 50) {
  dims <- dims %||% seq_len(min(20, npcs))
  dims <- as.integer(dims)
  dims[dims > 0]
}

.sn_supported_integration_methods <- function() {
  c(
    "unintegrated", "harmony", "coralysis", "seurat_cca", "seurat_rpca",
    "scvi", "scanvi", "scpoli", "bbknn", "totalvi", "mmochi"
  )
}

.sn_cluster_performance <- function(code) {
  invisible(gc(reset = TRUE))
  started_at <- Sys.time()
  elapsed_start <- proc.time()[["elapsed"]]
  value <- force(code)
  elapsed_seconds <- unname(proc.time()[["elapsed"]] - elapsed_start)
  memory <- gc()
  peak_r_memory_mb <- if (ncol(memory) >= 6L) {
    sum(memory[, 6L], na.rm = TRUE)
  } else {
    NA_real_
  }
  list(
    value = value,
    performance = list(
      elapsed_seconds = elapsed_seconds,
      peak_memory_mb = peak_r_memory_mb,
      peak_r_memory_mb = peak_r_memory_mb,
      backend_peak_rss_mb = NA_real_,
      memory_scope = "R_heap",
      started_at = as.character(started_at),
      completed_at = as.character(Sys.time())
    )
  )
}

.sn_complete_integration_performance <- function(performance, object) {
  backend <- object@misc$integration$backend_performance %||% NULL
  backend_peak <- backend$backend_peak_rss_mb %||% NA_real_
  performance$backend_peak_rss_mb <- backend_peak
  peaks <- c(performance$peak_r_memory_mb, backend_peak)
  peaks <- peaks[is.finite(peaks)]
  performance$peak_memory_mb <- if (length(peaks) > 0L) max(peaks) else NA_real_
  performance$memory_scope <- if (is.finite(backend_peak)) {
    "max_of_R_heap_and_backend_process_tree_rss"
  } else {
    "R_heap"
  }
  performance$backend <- backend
  performance
}

.sn_external_integration_control <- function() {
  list(
    runtime_dir = NULL,
    pixi_project = NULL,
    pixi_project_dir = NULL,
    pixi_home = NULL,
    run_dir = NULL,
    pixi = NULL,
    manifest_path = NULL,
    manifest_lines = NULL,
    overwrite_manifest = FALSE,
    platforms = NULL,
    install_pixi = TRUE,
    pixi_version = "latest",
    pixi_download_url = NULL,
    mirror = "default",
    mirror_append_original = TRUE,
    script = NULL,
    environment = NULL
  )
}

.sn_integration_control_templates <- function() {
  external <- .sn_external_integration_control()
  accelerator <- c(external, list(accelerator = "auto", cuda_version = NULL))
  list(
    unintegrated = list(),
    harmony = list(theta = 2, group_by_vars = NULL),
    coralysis = list(
      icp_args = list(
        threads = 1L,
        verbose = TRUE,
        RNGseed = 717L,
        build.train.params = list(nhvg = 2000L, p = 30L)
      ),
      pca_args = list(
        assay.name = "joint.probability",
        p = 30L,
        dimred.name = "Coralysis",
        return.model = TRUE
      ),
      store_sce = TRUE
    ),
    seurat_cca = list(
      orig.reduction = "pca", assay = "RNA", features = NULL,
      dims = 1:30, new.reduction = "integrated.cca", verbose = TRUE
    ),
    seurat_rpca = list(
      orig.reduction = "pca", assay = "RNA", features = NULL,
      dims = 1:30, new.reduction = "integrated.rpca", verbose = TRUE
    ),
    scvi = c(accelerator, list(
      reduction = "scvi", label_by = NULL, unlabeled_category = "Unknown",
      n_latent = 30L, seed = 717L, max_epochs = NULL,
      model_args = list(), train_args = list(), write_h5ad = TRUE
    )),
    scanvi = c(accelerator, list(
      reduction = "scanvi", label_by = NULL, unlabeled_category = "Unknown",
      n_latent = 30L, seed = 717L, max_epochs = NULL,
      scanvi_max_epochs = NULL, model_args = list(), train_args = list(),
      scanvi_model_args = list(), scanvi_train_args = list(), write_h5ad = TRUE
    )),
    scpoli = c(accelerator, list(
      reduction = "scpoli", label_by = NULL, n_latent = 10L,
      embedding_dims = 5L, latent_batch_size = 2048L, seed = 717L,
      n_epochs = 100L, max_epochs = NULL, pretraining_epochs = 90L,
      model_args = list(), train_args = list(), write_h5ad = TRUE,
      save_model = TRUE
    )),
    bbknn = c(external, list(
      graph_name = "bbknn_snn", umap_reduction = "umap", seed = 717L,
      bbknn_args = list(), umap_args = list()
    )),
    totalvi = c(accelerator, list(
      reduction = "totalvi", label_by = NULL, n_latent = 30L, seed = 717L,
      max_epochs = NULL, model_args = list(), train_args = list(),
      totalvi_model_args = list(), totalvi_train_args = list(),
      protein_assay = NULL, protein_layer = "counts", protein_features = NULL,
      adt_assay = "ADT", adt_layer = "counts", adt_features = NULL,
      protein_obsm_key = "protein_expression", write_h5ad = TRUE
    )),
    mmochi = c(external, list(
      protein_layer = "data", data_key = "protein",
      key_added = "landmark_protein", single_peaks = list(),
      marker_bandwidths = list(), peak_overrides = list(),
      inclusion_mask = NULL, landmark_args = list(), show = FALSE,
      reduction = "mmochi", corrected_layer = "mmochi.data",
      store_corrected_layer = TRUE,
      single_sample_batch_key = ".sn_mmochi_single_sample",
      keep_single_sample_batch = FALSE
    ))
  )
}

#' Inspect complete integration-control templates
#'
#' Returns the Shennong-supported `integration_control` fields and their
#' defaults for each integration backend. Values that depend on the input data,
#' such as Coralysis PCA rank or Seurat integration features, are illustrative
#' defaults and are resolved against the object by [sn_run_cluster()]. Extra
#' fields supplied for Seurat CCA/RPCA are forwarded to
#' `Seurat::IntegrateLayers()`.
#'
#' @details The returned templates enumerate every field consumed directly by
#'   Shennong. The main method-specific controls are:
#'
#'   - `unintegrated`: no backend controls.
#'   - `harmony`: `theta`, `group_by_vars`.
#'   - `coralysis`: `icp_args`, `pca_args`, `store_sce`.
#'   - `seurat_cca`, `seurat_rpca`: `orig.reduction`, `assay`, `features`,
#'     `dims`, `new.reduction`, `verbose`; additional fields are forwarded to
#'     `Seurat::IntegrateLayers()`.
#'   - `scvi`: runtime/pixi controls plus `accelerator`, `cuda_version`,
#'     `reduction`, `label_by`, `unlabeled_category`, `n_latent`, `seed`,
#'     `max_epochs`, `model_args`, `train_args`, `write_h5ad`.
#'   - `scanvi`: all scVI controls plus `scanvi_max_epochs`,
#'     `scanvi_model_args`, and `scanvi_train_args`; `label_by` is required.
#'   - `scpoli`: runtime/pixi/accelerator controls plus `reduction`, `label_by`,
#'     `n_latent`, `embedding_dims`, `latent_batch_size`, `seed`, `n_epochs`,
#'     `max_epochs`, `pretraining_epochs`, `model_args`, `train_args`,
#'     `write_h5ad`, and `save_model`.
#'   - `bbknn`: runtime/pixi controls plus `graph_name`, `umap_reduction`,
#'     `seed`, `bbknn_args`, and `umap_args`.
#'   - `totalvi`: all scVI runtime/accelerator controls plus
#'     `totalvi_model_args`, `totalvi_train_args`, `protein_assay`,
#'     `protein_layer`, `protein_features`, `adt_assay`, `adt_layer`,
#'     `adt_features`, and `protein_obsm_key`.
#'   - `mmochi`: runtime/pixi controls plus `protein_layer`, `data_key`,
#'     `key_added`, `single_peaks`, `marker_bandwidths`, `peak_overrides`,
#'     `inclusion_mask`, `landmark_args`, `show`, `reduction`,
#'     `corrected_layer`, `store_corrected_layer`, `single_sample_batch_key`,
#'     and `keep_single_sample_batch`.
#'
#'   Runtime/pixi fields are `runtime_dir`, `pixi_project`,
#'   `pixi_project_dir`, `pixi_home`, `run_dir`, `pixi`, `manifest_path`,
#'   `manifest_lines`, `overwrite_manifest`, `platforms`, `install_pixi`,
#'   `pixi_version`, `pixi_download_url`, `mirror`,
#'   `mirror_append_original`, `script`, and `environment`.
#'
#' @param method Optional integration method name or vector. When `NULL`, return
#'   templates for every supported method.
#'
#' @return A named list keyed by integration method, or one named control list
#'   when a single `method` is requested.
#' @examples
#' sn_get_integration_control_template("harmony")
#' sn_get_integration_control_template(c("scvi", "scanvi"))
#' @export
sn_get_integration_control_template <- function(method = NULL) {
  templates <- .sn_integration_control_templates()
  if (is.null(method)) {
    return(templates)
  }
  method <- .sn_normalize_integration_methods(method)
  unknown <- setdiff(method, names(templates))
  if (length(unknown) > 0L) {
    stop("Unsupported integration method(s): ", paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  selected <- templates[method]
  if (length(selected) == 1L) selected[[1L]] else selected
}

.sn_normalize_integration_methods <- function(integration_method) {
  if (!is.character(integration_method) || length(integration_method) == 0L ||
      anyNA(integration_method) || any(!nzchar(integration_method))) {
    stop("`integration_method` must contain one or more method names.", call. = FALSE)
  }
  if (any(integration_method == "unintergrated")) {
    warning(
      "`integration_method = \"unintergrated\"` is deprecated; use \"unintegrated\".",
      call. = FALSE
    )
    integration_method[integration_method == "unintergrated"] <- "unintegrated"
  }
  unknown <- setdiff(integration_method, .sn_supported_integration_methods())
  if (length(unknown) > 0L) {
    stop(
      "Unsupported `integration_method`: ", paste(unique(unknown), collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (anyDuplicated(integration_method)) {
    stop("`integration_method` must not contain duplicate methods.", call. = FALSE)
  }
  integration_method
}

.sn_reduction_key <- function(prefix, method) {
  paste0(toupper(gsub("[^A-Za-z0-9]", "", paste0(prefix, method))), "_")
}

.sn_resolve_multimodal_method <- function(multimodal_method = NULL,
                                          integration_method = "harmony",
                                          integration_method_supplied = FALSE) {
  supported <- c("wnn", "coralysis", "totalvi", "mmochi")
  if (is.null(multimodal_method)) {
    if (isTRUE(integration_method_supplied) && integration_method %in% supported) {
      return(integration_method)
    }
    return("wnn")
  }
  match.arg(multimodal_method, supported)
}

.sn_cluster_requires_rna_workflow <- function(modality, multimodal_method = NULL) {
  !(identical(modality, "cite_seq") && multimodal_method %in% c("coralysis", "mmochi"))
}

.sn_cluster_requires_rna_pca <- function(modality,
                                         multimodal_method = NULL,
                                         batch = NULL,
                                         integration_method = "harmony") {
  if (identical(modality, "cite_seq")) {
    return(identical(multimodal_method, "wnn"))
  }
  is.null(batch) || integration_method %in% c("unintegrated", "harmony", "seurat_cca", "seurat_rpca", "bbknn")
}

.sn_cluster_requires_adt_data <- function(modality,
                                          multimodal_method = NULL,
                                          integration_control = list()) {
  if (!identical(modality, "cite_seq")) {
    return(FALSE)
  }
  if (multimodal_method %in% c("wnn", "coralysis")) {
    return(TRUE)
  }
  if (identical(multimodal_method, "mmochi")) {
    protein_layer <- integration_control$protein_layer %||% integration_control$adt_layer %||% "data"
    return(identical(protein_layer, "data"))
  }
  FALSE
}

.sn_cluster_requires_adt_pca <- function(modality, multimodal_method = NULL) {
  identical(modality, "cite_seq") && identical(multimodal_method, "wnn")
}

.sn_coralysis_store_sce <- function(integration_control = list()) {
  !identical(integration_control$store_sce, FALSE)
}

.sn_resolve_find_clusters_algorithm <- function(cluster_algorithm = c("louvain", "louvain_multilevel", "slm", "leiden")) {
  if (is.numeric(cluster_algorithm) && length(cluster_algorithm) == 1L) {
    algorithm <- as.integer(cluster_algorithm)
    if (!algorithm %in% 1:4) {
      stop("`cluster_algorithm` must be one of 1, 2, 3, 4 or a supported algorithm name.", call. = FALSE)
    }
    return(algorithm)
  }

  cluster_algorithm <- match.arg(cluster_algorithm)
  switch(
    cluster_algorithm,
    louvain = 1L,
    louvain_multilevel = 2L,
    slm = 3L,
    leiden = 4L
  )
}

.sn_cluster_stage_order <- c(
  normalize = 1L,
  cell_cycle = 2L,
  hvg = 3L,
  pca = 4L,
  adt = 5L,
  integration = 6L,
  neighbors = 7L,
  clusters = 8L,
  umap = 9L,
  tsne = 10L
)

.sn_cluster_matrix_block_payload <- function(x) {
  sparse <- tryCatch(
    {
      if (inherits(x, "sparseMatrix")) {
        methods::as(x, "dgCMatrix")
      } else if (inherits(x, "Matrix") || is.matrix(x)) {
        methods::as(Matrix::Matrix(x, sparse = TRUE), "dgCMatrix")
      } else {
        methods::as(x, "dgCMatrix")
      }
    },
    error = function(error) NULL
  )
  if (!is.null(sparse)) {
    sparse <- Matrix::drop0(sparse)
    return(list(
      storage = "dgCMatrix",
      dim = dim(sparse),
      p = sparse@p,
      i = sparse@i,
      x = sparse@x
    ))
  }

  dense <- as.matrix(x)
  list(
    storage = "matrix",
    dim = dim(dense),
    values = unname(dense)
  )
}

.sn_cluster_matrix_signature <- function(x, block_ncol = 1024L) {
  dimensions <- dim(x)
  if (length(dimensions) != 2L) {
    stop("Clustering analysis inputs must be two-dimensional.", call. = FALSE)
  }
  block_ncol <- max(1L, as.integer(block_ncol))
  starts <- if (dimensions[[2L]] > 0L) {
    seq.int(1L, dimensions[[2L]], by = block_ncol)
  } else {
    integer(0)
  }
  block_hashes <- vapply(starts, function(start) {
    finish <- min(start + block_ncol - 1L, dimensions[[2L]])
    block <- x[, seq.int(start, finish), drop = FALSE]
    digest::digest(
      .sn_cluster_matrix_block_payload(block),
      algo = "sha256",
      serialize = TRUE
    )
  }, character(1))

  digest::digest(
    list(
      schema_version = 1L,
      dim = dimensions,
      rownames = rownames(x),
      colnames = colnames(x),
      block_ncol = block_ncol,
      block_hashes = block_hashes
    ),
    algo = "sha256",
    serialize = TRUE
  )
}

.sn_cluster_analysis_input_signature <- function(object,
                                                 assay = "RNA",
                                                 layer = "counts") {
  input <- .sn_get_seurat_layer_data(
    object = object,
    assay = assay,
    layer = layer
  )
  digest::digest(
    list(
      schema_version = 1L,
      assay = assay,
      layer = layer,
      values = .sn_cluster_matrix_signature(input)
    ),
    algo = "sha256",
    serialize = TRUE
  )
}

.sn_cluster_metadata_signature <- function(object, columns = character(0)) {
  columns <- unique(as.character(columns %||% character(0)))
  columns <- columns[!is.na(columns) & nzchar(columns)]
  if (length(columns) == 0L) {
    return(NULL)
  }

  metadata <- object[[]]
  values <- lapply(columns, function(column) {
    if (!column %in% colnames(metadata)) {
      return(structure(list(), class = "shennong_missing_metadata_column"))
    }
    metadata[[column]]
  })
  names(values) <- columns
  digest::digest(
    list(
      schema_version = 1L,
      cells = rownames(metadata),
      columns = columns,
      values = values
    ),
    algo = "sha256",
    serialize = TRUE
  )
}

.sn_cluster_control_metadata_columns <- function(control) {
  if (!is.list(control)) {
    return(character(0))
  }
  metadata_fields <- c(
    "group_by_vars", "label_by", "batch_by", "batch_key",
    "categorical_covariate_keys", "continuous_covariate_keys"
  )
  current <- unlist(
    control[intersect(metadata_fields, names(control))],
    use.names = FALSE
  )
  nested <- unlist(
    lapply(Filter(is.list, control), .sn_cluster_control_metadata_columns),
    use.names = FALSE
  )
  unique(c(as.character(current), nested))
}

.sn_resolve_cluster_rerun_from <- function(rerun_from = NULL) {
  if (is.null(rerun_from)) {
    return(NULL)
  }
  rerun_from <- match.arg(rerun_from, names(.sn_cluster_stage_order))
  rerun_from
}

.sn_can_reuse_cluster_stage <- function(object,
                                        stage,
                                        signature,
                                        reuse = TRUE,
                                        rerun_from = NULL,
                                        required = TRUE) {
  if (!isTRUE(reuse)) {
    return(FALSE)
  }
  if (!is.null(rerun_from) && .sn_cluster_stage_order[[stage]] >= .sn_cluster_stage_order[[rerun_from]]) {
    return(FALSE)
  }
  stages <- object@misc$sn_run_cluster$stages %||% list()
  stage_info <- stages[[stage]] %||% NULL
  if (is.null(stage_info) || !identical(stage_info$signature, signature)) {
    return(FALSE)
  }
  if (is.function(required)) {
    return(isTRUE(required(object, stage_info)))
  }
  isTRUE(required)
}

.sn_record_cluster_stage <- function(object, stage, signature, ...) {
  metadata <- list(...)
  object@misc$sn_run_cluster <- object@misc$sn_run_cluster %||% list()
  object@misc$sn_run_cluster$stages <- object@misc$sn_run_cluster$stages %||% list()
  object@misc$sn_run_cluster$stages[[stage]] <- c(
    list(
      signature = signature,
      created_at = as.character(Sys.time())
    ),
    metadata
  )
  object
}

.sn_has_seurat_layer <- function(object, assay = "RNA", layer = "data") {
  assay %in% names(object@assays) &&
    length(.sn_match_seurat_layers(object = object, assay = assay, layer = layer)) > 0L
}

.sn_ensure_cluster_algorithm_dependencies <- function(cluster_algorithm_value,
                                                      leiden_method = "leidenbase",
                                                      auto_install = TRUE,
                                                      repos = getOption("repos"),
                                                      ask = FALSE) {
  if (!identical(as.integer(cluster_algorithm_value), 4L) || !identical(leiden_method, "leidenbase")) {
    return(invisible(TRUE))
  }
  missing <- .sn_find_missing_packages("leidenbase")
  if (length(missing) == 0L) {
    return(invisible(TRUE))
  }
  if (!isTRUE(auto_install)) {
    check_installed("leidenbase", reason = "to run Leiden clustering with `cluster_algorithm = \"leiden\"`.")
  }

  .sn_log_info("Installing missing clustering package for Leiden: leidenbase.")
  sn_install_dependencies(
    packages = "leidenbase",
    missing_only = TRUE,
    repos = repos,
    ask = ask,
    github_dependencies = NA
  )
  if (length(.sn_find_missing_packages("leidenbase")) > 0L) {
    check_installed("leidenbase", reason = "to run Leiden clustering with `cluster_algorithm = \"leiden\"`.")
  }
  invisible(TRUE)
}

.sn_valid_reduction_dims <- function(object, reduction, dims) {
  embeddings <- Seurat::Embeddings(object = object[[reduction]])
  max_dim <- ncol(embeddings)
  valid_dims <- dims[dims <= max_dim]
  if (length(valid_dims) == 0L) {
    stop(glue("No requested `dims` are available in reduction '{reduction}'."), call. = FALSE)
  }
  valid_dims
}

.sn_select_variable_features <- function(object,
                                         nfeatures = 3000,
                                         split_by = NULL,
                                         assay = NULL,
                                         layer = NULL,
                                         verbose = TRUE) {
  if (is_null(split_by) || identical(split_by, "global")) {
    object <- .sn_with_default_seurat_acceleration(
      Seurat::FindVariableFeatures(object, nfeatures = nfeatures, verbose = verbose),
      object = object,
      assay = assay
    )
    return(list(
      object = object,
      features = Seurat::VariableFeatures(object = object)
    ))
  }

  if (!split_by %in% colnames(object[[]])) {
    stop(glue("`hvg_group_by` must be NULL or a metadata column name. '{split_by}' was not found."))
  }

  metadata <- object[[]]
  split_values <- as.character(metadata[[split_by]])
  names(split_values) <- rownames(metadata)
  analysis_cells <- colnames(object)
  assay <- assay %||% SeuratObject::DefaultAssay(object = object)
  if (!is.null(layer) && assay %in% names(object@assays)) {
    matched_layers <- .sn_match_seurat_layers(object = object, assay = assay, layer = layer)
    if (length(matched_layers) > 0L) {
      layer_cells <- unique(unlist(lapply(matched_layers, function(current_layer) {
        colnames(SeuratObject::LayerData(object = object, assay = assay, layer = current_layer))
      }), use.names = FALSE))
      analysis_cells <- intersect(analysis_cells, layer_cells)
    }
  }
  analysis_cells <- intersect(analysis_cells, names(split_values))
  split_values <- split_values[analysis_cells]
  valid_group <- !is.na(split_values) & nzchar(split_values)
  if (verbose && any(!valid_group)) {
    .sn_log_warn(
      "[hvg] Skipping {sum(!valid_group)} cell(s) with missing or empty `{split_by}` values during grouped HVG selection."
    )
  }
  analysis_cells <- analysis_cells[valid_group]
  split_values <- split_values[valid_group]
  if (length(analysis_cells) == 0L) {
    stop(glue("No cells with non-missing `{split_by}` values are available for grouped HVG selection."), call. = FALSE)
  }

  hvg_object <- object
  SeuratObject::DefaultAssay(object = hvg_object) <- assay
  extra_assays <- setdiff(names(hvg_object@assays), assay)
  if (length(extra_assays) > 0L) {
    for (extra_assay in extra_assays) {
      hvg_object[[extra_assay]] <- NULL
    }
  }

  split_levels <- unique(split_values)
  skipped_groups <- character(0)
  feature_lists <- lapply(split_levels, function(current_group) {
    current_cells <- analysis_cells[split_values == current_group]
    if (length(current_cells) < 2L) {
      skipped_groups <<- c(skipped_groups, current_group)
      return(character(0))
    }
    current_object <- hvg_object[, current_cells]
    current_object <- .sn_with_default_seurat_acceleration(
      Seurat::FindVariableFeatures(
        current_object,
        nfeatures = nfeatures,
        verbose = FALSE
      ),
      object = current_object,
      assay = assay
    )
    Seurat::VariableFeatures(current_object)
  })
  feature_lists <- Filter(length, feature_lists)
  if (verbose && length(skipped_groups) > 0L) {
    .sn_log_warn(
      "[hvg] Skipping {length(skipped_groups)} `{split_by}` group(s) with fewer than 2 analyzable cells: ",
      "{paste(utils::head(skipped_groups, 5), collapse = ', ')}",
      "{if (length(skipped_groups) > 5) '...' else ''}."
    )
  }
  if (length(feature_lists) == 0L) {
    stop(glue("No `{split_by}` groups have enough analyzable cells for grouped HVG selection."), call. = FALSE)
  }

  feature_frequency <- sort(table(unlist(feature_lists, use.names = FALSE)), decreasing = TRUE)
  all_features <- names(feature_frequency)
  feature_ranks <- lapply(feature_lists, function(features) {
    stats::setNames(seq_along(features), features)
  })
  mean_rank <- vapply(all_features, function(feature) {
    mean(vapply(feature_ranks, function(rank_map) {
      if (feature %in% names(rank_map)) {
        rank_map[[feature]]
      } else {
        nfeatures + 1
      }
    }, numeric(1)))
  }, numeric(1))
  selected <- all_features[order(-as.numeric(feature_frequency), mean_rank, all_features)]
  selected <- utils::head(selected, nfeatures)

  Seurat::VariableFeatures(object = object) <- selected
  list(object = object, features = selected)
}

.sn_resolve_user_hvg_features <- function(features,
                                          object,
                                          arg_name = "hvg_features",
                                          verbose = TRUE) {
  if (is_null(features)) {
    return(list(features = character(0), missing = character(0)))
  }
  if (!is.character(features)) {
    stop("`", arg_name, "` must be NULL or a character vector of feature names.", call. = FALSE)
  }

  features <- unique(features[!is.na(features) & nzchar(features)])
  if (length(features) == 0L) {
    return(list(features = character(0), missing = character(0)))
  }

  present <- intersect(features, rownames(object))
  missing <- setdiff(features, present)
  if (length(missing) > 0L && isTRUE(verbose)) {
    .sn_log_warn(
      "`{arg_name}` contains {length(missing)} feature(s) not present in the object; ",
      "ignoring examples: {paste(utils::head(missing, 5), collapse = ', ')}."
    )
  }

  list(features = present, missing = missing)
}

.sn_normalize_block_signature_query <- function(query) {
  aliases <- c(
    g1s = "Programs/cellCycle.G1S",
    g2m = "Programs/cellCycle.G2M"
  )
  normalized_query <- .sn_signature_normalize_key(query)
  alias <- unname(aliases[normalized_query])
  if (!is.na(alias)) alias else query
}

.sn_match_signature_query <- function(query, leaf_table) {
  query <- .sn_normalize_block_signature_query(query)

  current_matches <- which(leaf_table$path == query)
  if (length(current_matches) == 0L) {
    current_matches <- which(leaf_table$name == query)
  }
  if (length(current_matches) == 0L) {
    normalized_query <- .sn_signature_normalize_key(query)
    normalized_path <- .sn_signature_normalize_key(leaf_table$path)
    normalized_name <- .sn_signature_normalize_key(leaf_table$name)
    current_matches <- which(normalized_path == normalized_query)
    if (length(current_matches) == 0L) {
      current_matches <- which(normalized_name == normalized_query)
    }
  }

  if (length(current_matches) > 1L) {
    stop(
      glue(
        "Signature query '{query}' is ambiguous. ",
        "Use one of the full paths instead: ",
        "{paste(leaf_table$path[current_matches], collapse = ', ')}."
      ),
      call. = FALSE
    )
  }

  current_matches
}

.sn_resolve_block_genes <- function(block_genes,
                                    species = c("human", "mouse"),
                                    verbose = TRUE) {
  if (is_null(block_genes)) {
    return(NULL)
  }
  if (!is.character(block_genes)) {
    stop("`block_genes` must be NULL or a character vector.", call. = FALSE)
  }

  block_genes <- unique(trimws(block_genes[!is.na(block_genes) & nzchar(block_genes)]))
  if (length(block_genes) == 0L) {
    return(character(0))
  }

  species <- rlang::arg_match(species, c("human", "mouse"))
  leaf_table <- .sn_signature_leaf_table(species = species)
  signature_matches <- lapply(block_genes, .sn_match_signature_query, leaf_table = leaf_table)
  matched_signature <- lengths(signature_matches) > 0L
  signature_queries <- block_genes[matched_signature]
  custom_queries <- block_genes[!matched_signature]
  signature_rows <- unique(unlist(signature_matches[matched_signature], use.names = FALSE))
  signature_genes <- unique(unlist(leaf_table$genes[signature_rows], use.names = FALSE))

  custom_genes <- character(0)
  if (length(custom_queries) > 0L) {
    check_installed(
      "HGNChelper",
      reason = "to validate custom `block_genes` symbols."
    )
    checked <- HGNChelper::checkGeneSymbols(custom_queries, species = species)
    suggested <- checked$Suggested.Symbol
    valid <- !is.na(suggested) & nzchar(suggested)
    custom_genes <- unique(suggested[valid])
    invalid <- custom_queries[!valid]

    if (length(invalid) > 0L) {
      .sn_log_warn(
        "Removed {length(invalid)} invalid genes from block list ",
        "(e.g. {paste(utils::head(invalid, 3), collapse=', ')})."
      )
    }
  }

  resolved <- unique(c(signature_genes, custom_genes))
  if (verbose) {
    if (length(signature_queries) > 0L) {
      .sn_log_info(
        "Loaded {length(signature_genes)} blocked genes from ",
        "{length(signature_queries)} built-in signature query(s)."
      )
    }
    if (length(custom_genes) > 0L) {
      .sn_log_info("Using {length(custom_genes)} custom blocked genes.")
    }
    .sn_log_info("Resolved {length(resolved)} blocked genes total.")
  }

  resolved
}

.sn_resolve_assay_features <- function(object,
                                       assay,
                                       features = NULL,
                                       arg_name = "features",
                                       verbose = TRUE) {
  available <- rownames(object[[assay]])
  if (is.null(features)) {
    return(list(features = available, missing = character(0)))
  }
  if (!is.character(features)) {
    stop("`", arg_name, "` must be NULL or a character vector of feature names.", call. = FALSE)
  }

  features <- unique(features[!is.na(features) & nzchar(features)])
  present <- intersect(features, available)
  missing <- setdiff(features, present)
  if (length(missing) > 0L && isTRUE(verbose)) {
    .sn_log_warn(
      "`{arg_name}` contains {length(missing)} feature(s) not present in assay '{assay}'; ",
      "ignoring examples: {paste(utils::head(missing, 5), collapse = ', ')}."
    )
  }
  if (length(present) == 0L) {
    stop("`", arg_name, "` did not contain any features present in assay '", assay, "'.", call. = FALSE)
  }

  list(features = present, missing = missing)
}

.sn_cluster_tail_defaults <- function() {
  list(
    cluster_control = list(),
    reuse = TRUE,
    rerun_from = NULL,
    auto_install = TRUE,
    install_repos = getOption("repos"),
    install_ask = FALSE,
    hvg_group_by = NULL,
    rare_feature_method = "none",
    rare_feature_group_by = NULL,
    rare_feature_n = 200,
    rare_feature_control = list(),
    block_genes = c("heatshock", "ribo", "mito", "tcr", "immunoglobulins", "pseudogenes"),
    theta = 2,
    group_by_vars = NULL,
    npcs = 50,
    dims = NULL,
    species = NULL,
    assay = "RNA",
    layer = "counts",
    modality = c("rna", "cite_seq"),
    multimodal_method = NULL,
    adt_assay = "ADT",
    adt_layer = "counts",
    adt_features = NULL,
    adt_npcs = 30,
    adt_dims = NULL,
    wnn_control = list(),
    umap_control = list(),
    run_tsne = FALSE,
    tsne_control = list(),
    checkpoint_dir = NULL,
    resume = TRUE,
    checkpoint_compress = FALSE,
    return_cluster = FALSE,
    verbose = TRUE
  )
}

.sn_resolve_cluster_tail_args <- function(dots) {
  defaults <- .sn_cluster_tail_defaults()
  if (length(dots) == 0L) {
    return(list(values = defaults, supplied = character()))
  }

  dot_names <- names(dots)
  if (is.null(dot_names)) {
    dot_names <- rep("", length(dots))
  }
  named <- nzchar(dot_names)
  unknown <- setdiff(unique(dot_names[named]), names(defaults))
  if (length(unknown) > 0L) {
    stop("Unused clustering argument(s): ", paste(unknown, collapse = ", "), ".", call. = FALSE)
  }
  duplicated_names <- unique(dot_names[named][duplicated(dot_names[named])])
  if (length(duplicated_names) > 0L) {
    stop("Clustering argument(s) supplied more than once: ", paste(duplicated_names, collapse = ", "), ".", call. = FALSE)
  }

  unnamed <- which(!named)
  available <- setdiff(names(defaults), dot_names[named])
  if (length(unnamed) > length(available)) {
    stop("Too many positional clustering arguments were supplied after `leiden_objective_function`.", call. = FALSE)
  }
  if (length(unnamed) > 0L) {
    dot_names[unnamed] <- utils::head(available, length(unnamed))
  }
  if (anyDuplicated(dot_names)) {
    duplicated_names <- unique(dot_names[duplicated(dot_names)])
    stop("Clustering argument(s) supplied more than once: ", paste(duplicated_names, collapse = ", "), ".", call. = FALSE)
  }

  for (index in seq_along(dots)) {
    defaults[dot_names[[index]]] <- list(dots[[index]])
  }
  list(values = defaults, supplied = dot_names)
}

#' Run clustering for a single dataset or batch integration workflow
#'
#' This function is the main clustering entry point in `Shennong`.
#' When `batch = NULL`, it performs single-dataset clustering with either the
#' standard Seurat workflow, an SCTransform workflow, or a single-sample
#' CITE-seq workflow. When `batch` is supplied, it performs batch integration
#' followed by clustering and UMAP.
#'
#' @param object A \code{Seurat} object.
#' @param batch A column name in \code{object@meta.data} specifying the batch
#'   labels used for integration. If \code{NULL}, no RNA batch integration is
#'   performed. CITE-seq MMoCHi runs in single-sample mode by passing an
#'   internal constant batch key to the Python backend.
#' @param normalization_method One of \code{"seurat"}, \code{"scran"}, or
#'   \code{"sctransform"}. The \code{"seurat"} and \code{"scran"} workflows can
#'   be followed by any supported \code{integration_method} when \code{batch}
#'   is supplied. The SCTransform workflow can currently be
#'   combined with \code{integration_method = "harmony"} by supplying
#'   \code{batch}.
#' @param integration_method One or more batch-analysis methods used when
#'   \code{batch} is supplied. \code{"unintegrated"} keeps the PCA baseline;
#'   multiple values run against the same normalized/HVG/PCA preparation and
#'   retain method-specific reductions, graphs, cluster columns, UMAP, and
#'   optional t-SNE results in one object. Scalar analysis parameters supplied
#'   as vectors are expanded as a conditional Cartesian grid; parameters that
#'   are naturally vector-valued, such as \code{dims}, \code{hvg_features},
#'   \code{vars_to_regress}, and \code{block_genes}, remain intact. Supported
#'   integration backends are
#'   \code{"harmony"},
#'   \code{"coralysis"}, \code{"seurat_cca"}, \code{"seurat_rpca"},
#'   \code{"scvi"}, \code{"scanvi"}, \code{"scpoli"}, \code{"bbknn"}, and
#'   \code{"totalvi"}. \code{"mmochi"} is
#'   accepted as a CITE-seq convenience alias and requires
#'   \code{modality = "cite_seq"}.
#'   \code{"harmony"} preserves the historical Shennong behavior.
#'   \code{"coralysis"} runs native Coralysis on the selected log-normalized
#'   feature set and stores the integrated embedding as the \code{"coralysis"}
#'   reduction. \code{"scvi"}, \code{"scanvi"}, and \code{"scpoli"} export
#'   the selected \code{assay}/\code{layer} count matrix to a pixi-managed
#'   environment under \code{~/.shennong/pixi/}, run the Python backend, and
#'   import the latent representation as a Seurat reduction. \code{"bbknn"}
#'   computes a batch-balanced graph from the selected-layer PCA and uses that
#'   graph directly for clustering and UMAP. \code{"totalvi"} is used for
#'   RNA+ADT CITE-seq workflows and is
#'   usually selected through \code{modality = "cite_seq"} and
#'   \code{multimodal_method = "totalvi"}.
#'   Python expression/protein inputs remain sparse; learned backends may create
#'   bounded dense minibatch tensors, while imported latent/PCA/UMAP results are
#'   low-dimensional dense outputs.
#' @param integration_control Optional named list of backend-specific
#'   parameters. With multiple methods, provide a list keyed by method, for
#'   example \code{list(harmony = list(theta = 3), coralysis = list(...))};
#'   an optional \code{.default} entry is merged into every method. For
#'   a complete executable template of every accepted field and its default,
#'   call \code{sn_get_integration_control_template()} or
#'   \code{sn_get_integration_control_template("scvi")}. For
#'   \code{"coralysis"}, use \code{icp_args} for
#'   \code{RunParallelDivisiveICP()} arguments, \code{pca_args} for
#'   \code{RunPCA()} arguments, and \code{store_sce = FALSE} only when the
#'   trained Coralysis SingleCellExperiment should not be kept under
#'   \code{object@misc$coralysis}. The default is \code{store_sce = TRUE} so
#'   native Coralysis references can be used directly by
#'   \code{sn_transfer_labels(method = "coralysis")}.
#'   For \code{"seurat_cca"} and \code{"seurat_rpca"}, values are forwarded
#'   to \code{Seurat::IntegrateLayers()}. For \code{"scvi"} and
#'   \code{"scanvi"}, common fields include \code{runtime_dir},
#'   \code{pixi_project}, \code{pixi_home}, \code{run_dir}, \code{pixi},
#'   \code{manifest_path}, \code{install_pixi}, \code{accelerator},
#'   \code{cuda_version}, \code{mirror}, \code{n_latent}, \code{max_epochs},
#'   \code{model_args}, \code{train_args}, and \code{write_h5ad};
#'   \code{"scanvi"} additionally requires \code{label_by} and accepts
#'   \code{unlabeled_category}. \code{"scpoli"} accepts optional
#'   \code{label_by}, \code{n_epochs}, \code{pretraining_epochs},
#'   \code{embedding_dims}, \code{latent_batch_size}, \code{model_args}, and
#'   \code{train_args}.
#'   \code{"bbknn"} accepts \code{bbknn_args} plus an optional
#'   \code{graph_name}; its imported graph is used instead of running
#'   \code{Seurat::FindNeighbors()}. \code{"totalvi"} additionally accepts
#'   \code{totalvi_model_args}, \code{totalvi_train_args}, and
#'   \code{protein_obsm_key}. \code{"mmochi"} additionally accepts
#'   \code{protein_layer}, \code{single_peaks}, \code{marker_bandwidths},
#'   \code{peak_overrides}, \code{inclusion_mask}, \code{landmark_args},
#'   \code{corrected_layer}, \code{store_corrected_layer},
#'   \code{single_sample_batch_key}, and \code{keep_single_sample_batch};
#'   Shennong runs MMoCHi's ADT landmark registration and imports the corrected
#'   protein matrix as a protein-derived reduction. When \code{batch = NULL},
#'   Shennong uses a constant internal backend batch key for single-sample
#'   registration. When Seurat accepts arbitrary assay layers, the corrected
#'   matrix is stored as \code{corrected_layer}; otherwise it is kept under
#'   \code{object@misc$mmochi$corrected_protein}.
#'   Use \code{sn_get_pixi_paths()} to inspect the
#'   generated directory layout, \code{sn_get_pixi_config_path()} to inspect the
#'   bundled \code{inst/pixi/} config, \code{sn_ensure_pixi()} to preinstall
#'   pixi, and \code{sn_configure_pixi_mirror()} to set Shennong-level mirrors.
#' @param nfeatures Number of variable features to select. Multiple values
#'   create separate preprocessing/embedding branches in a parameter-grid run.
#' @param hvg_features Optional character vector of user-supplied features to
#'   force into the selected backend feature set. For PCA-based workflows this
#'   is also the feature set used for scaling/PCA. These features are merged
#'   with internally selected HVGs and any rare-aware features after validating
#'   that they are present in \code{object}.
#' @param vars_to_regress Covariates to regress out in \code{ScaleData}.
#' @param resolution Resolution parameter for \code{FindClusters}. Multiple
#'   values create separate cluster columns while reusing the same graph and
#'   dimensional reductions.
#' @param cluster_algorithm Community-detection algorithm passed to
#'   \code{Seurat::FindClusters()}. Supported names are \code{"louvain"}
#'   (Seurat algorithm 1), \code{"louvain_multilevel"} (algorithm 2),
#'   \code{"slm"} (algorithm 3), and \code{"leiden"} (algorithm 4). Numeric
#'   values 1 through 4 are also accepted.
#'   Multiple explicitly supplied values form a parameter-grid axis.
#' @param cluster_name Optional metadata column name for the cluster labels.
#'   Defaults to Seurat's \code{"seurat_clusters"} behavior.
#' @param cluster_n_start,cluster_n_iter Number of starts and iterations passed
#'   to \code{Seurat::FindClusters()}.
#' @param cluster_random_seed Random seed passed to
#'   \code{Seurat::FindClusters()}.
#' @param seed Top-level reproducibility seed overriding
#'   \code{cluster_random_seed} (and any existing \code{integration_control$seed})
#'   when supplied. Precedence: \code{seed} > \code{cluster_random_seed} > default.
#' @param verbose Top-level progress logging switch forwarded to the clustering
#'   implementation; a \code{verbose} tail argument keeps precedence over it.
#' @param cluster_group_singletons Whether \code{Seurat::FindClusters()}
#'   should group singletons into the nearest cluster.
#' @param leiden_method Leiden implementation passed to
#'   \code{Seurat::FindClusters()} when \code{cluster_algorithm = "leiden"}.
#' @param leiden_objective_function Leiden objective function passed to
#'   \code{Seurat::FindClusters()}.
#' @param ... Additional clustering controls. Supported names include
#'   \code{cluster_control}, an optional named list of additional
#'   \code{Seurat::FindClusters()} arguments. Values here override Shennong's
#'   generated defaults.
#'   \code{reuse}: logical; when \code{TRUE}, reuse previously recorded
#'   \code{sn_run_cluster()} stages if their stored input signatures still match
#'   the current call. This lets resolution-only changes start at clustering,
#'   integration-method changes start at integration, and HVG changes start at
#'   feature selection instead of rerunning all earlier steps.
#'   \code{rerun_from}: optional stage name forcing recomputation from that stage
#'   onward while still allowing earlier matching stages to be reused. Supported
#'   values are \code{"normalize"}, \code{"cell_cycle"}, \code{"hvg"},
#'   \code{"pca"}, \code{"adt"}, \code{"integration"}, \code{"neighbors"},
#'   \code{"clusters"}, \code{"umap"}, and \code{"tsne"}.
#'   \code{auto_install}: logical; when \code{TRUE}, install missing optional
#'   clustering dependencies such as \pkg{leidenbase} before the relevant stage.
#'   \code{install_repos}: CRAN-like repositories used when \code{auto_install}
#'   installs CRAN packages.
#'   \code{install_ask}: passed to \code{BiocManager::install()} when
#'   \code{auto_install} installs Bioconductor packages through
#'   \code{sn_install_dependencies()}.
#'   \code{hvg_group_by}: optional metadata column used to compute highly variable
#'   genes within groups before merging and ranking them. When \code{NULL} and
#'   \code{batch} is supplied, Shennong reuses \code{batch} by default. Use
#'   \code{NULL} with \code{batch = NULL} to compute HVGs on the full object.
#'   \code{rare_feature_method}: optional rare-cell-aware feature methods appended
#'   to the base HVG set before PCA/clustering. Supported values are
#'   \code{"none"}, \code{"gini"}, and \code{"local_markers"}.
#'   \code{rare_feature_group_by}: optional metadata column used to define groups
#'   for \code{"local_markers"}. When \code{NULL}, Shennong builds a temporary
#'   coarse clustering from the base HVGs.
#'   \code{rare_feature_n}: number of rare-aware features to add per selected
#'   method. For example, \code{c("gini", "local_markers")} with
#'   \code{rare_feature_n = 50} can contribute up to 100 rare-aware features
#'   before de-duplication.
#'   Multiple values form a parameter-grid axis.
#'   \code{rare_feature_control}: named list of advanced rare-feature thresholds.
#'   Supported fields are \code{group_max_fraction}, \code{group_max_cells},
#'   \code{gene_max_fraction}, and \code{min_cells}.
#'   \code{block_genes}: character vector of bundled signature queries and/or
#'   custom gene symbols to exclude from internally selected HVGs. Signature
#'   queries can use leaf names such as \code{"ribo"} and
#'   \code{"cellCycle.G2M"} or full paths such as
#'   \code{"Programs/cellCycle.G1S"}; \code{"g1s"} and \code{"g2m"} are kept as
#'   short aliases for the cell-cycle signatures. Applies to both
#'   log-normalization and SCTransform workflows; explicit \code{hvg_features}
#'   are preserved even when they overlap a blocked signature.
#'   \code{theta}: the \code{theta} parameter for \code{harmony::RunHarmony}, controlling batch
#'   diversity preservation vs. correction. Used only when
#'   \code{integration_method = "harmony"}. Multiple values expand only the
#'   Harmony branch and do not duplicate other integration methods.
#'   \code{group_by_vars}: optional column name or character vector passed to
#'   \code{harmony::RunHarmony(group.by.vars = ...)}. Defaults to \code{batch}
#'   and is used only when \code{integration_method = "harmony"}.
#'   \code{npcs}: number of PCs to compute in \code{RunPCA}. Multiple values
#'   form a parameter-grid axis.
#'   \code{dims}: a numeric vector of PCs (dimensions) to use for neighbor search,
#'   clustering, and UMAP.
#'   \code{species}: optional species label. Used when block genes must be resolved
#'   from built-in signatures.
#'   \code{assay}: assay used for clustering. Defaults to \code{"RNA"}.
#'   \code{layer}: layer used as the input matrix. Defaults to \code{"counts"}.
#'   scVI, scANVI, scPoli, and the PCA upstream of BBKNN all honor this value,
#'   so a layer such as \code{"decontaminated_counts"} is used consistently.
#'   \code{modality}: workflow modality. \code{"rna"} runs the standard RNA-only
#'   workflow. \code{"cite_seq"} enables paired RNA+ADT workflows selected by
#'   \code{multimodal_method}.
#'   \code{multimodal_method}: CITE-seq backend used when
#'   \code{modality = "cite_seq"}. \code{"wnn"} combines RNA PCA with ADT PCA
#'   using Seurat's weighted nearest-neighbor workflow and clusters on
#'   \code{"wsnn"}. \code{"coralysis"} runs native Coralysis on the ADT assay
#'   as a log-normalized protein matrix. \code{"totalvi"} runs scvi-tools
#'   totalVI on RNA counts plus ADT counts and clusters on the imported totalVI
#'   latent representation.
#'   \code{"mmochi"} runs MMoCHi ADT landmark registration across batches, or
#'   in single-sample mode when \code{batch = NULL}, stores the corrected
#'   protein matrix when supported, computes a protein PCA reduction, and
#'   clusters on that reduction. When
#'   \code{NULL}, Shennong keeps the historical CITE-seq default
#'   \code{"wnn"} unless \code{integration_method} was explicitly set to one of
#'   the supported multimodal backends.
#'   \code{adt_assay}: assay containing antibody-derived tag counts for
#'   \code{modality = "cite_seq"}.
#'   \code{adt_layer}: layer in \code{adt_assay} used as ADT counts.
#'   \code{adt_features}: optional ADT/protein features used by CITE-seq backends.
#'   Defaults to all features in \code{adt_assay}.
#'   \code{adt_npcs}: number of ADT PCs to compute for \code{modality = "cite_seq"}.
#'   \code{adt_dims}: numeric vector of ADT PCs used in weighted nearest-neighbor
#'   graph construction. Defaults to \code{seq_len(min(18, adt_npcs))}.
#'   \code{wnn_control}: optional named list of additional
#'   \code{Seurat::FindMultiModalNeighbors()} arguments used only when
#'   \code{modality = "cite_seq"}. Values here override Shennong's generated
#'   defaults.
#'   \code{umap_control}: optional named list of additional
#'   \code{Seurat::RunUMAP()} arguments. Values here override Shennong's
#'   generated defaults, for example \code{n.neighbors}, \code{min.dist},
#'   \code{spread}, \code{metric}, \code{seed.use}, or \code{reduction.name}.
#'   \code{run_tsne}: logical; run t-SNE in addition to UMAP. It defaults to
#'   \code{FALSE}, including for multi-method and parameter-grid runs.
#'   \code{tsne_control}: optional named list of additional
#'   \code{Seurat::RunTSNE()} arguments. In multi-method mode the reduction name
#'   and key are generated per method and cannot overwrite another result.
#'   \code{checkpoint_dir}: optional directory for persistent parameter-grid
#'   checkpoints. After every completed run Shennong writes the current object,
#'   comparison manifest, performance records, and completed run IDs to a
#'   temporary file and atomically publishes it. Only the latest complete
#'   checkpoint for the call signature is retained. The signature includes a
#'   blockwise digest of the selected layer plus cell/feature identity and the
#'   metadata values used by batching, grouped HVGs, regression, or supervised
#'   integration.
#'   \code{resume}: logical; when \code{TRUE} (default), resume a matching
#'   checkpoint in \code{checkpoint_dir}. The same content/metadata-aware
#'   signature also covers package version, grid, and analysis arguments;
#'   incomplete \code{.partial} files are ignored.
#'   \code{checkpoint_compress}: logical; compress RDS checkpoints. It defaults
#'   to \code{FALSE} for faster writes at the cost of more disk space.
#'   \code{return_cluster}: if \code{TRUE}, return only the cluster assignments.
#'   Multi-method calls return a data frame with one cluster column per method.
#'   \code{verbose}: whether to print/log progress messages.
#'
#' @return A \code{Seurat} object with clustering results and embeddings, or a
#'   cluster_by vector if \code{return_cluster = TRUE}. Parameter-grid objects
#'   store per-run timing and memory fields in
#'   \code{object@misc$integration_comparison$performance}; native R backends
#'   report peak R heap usage, and Linux pixi backends additionally report the
#'   maximum child-process-tree RSS measured by GNU \code{time}. These values do
#'   not include GPU device memory.
#'
#' @examples
#' \dontrun{
#' seurat_obj <- sn_run_cluster(
#'   object = seurat_obj,
#'   normalization_method = "seurat",
#'   resolution = 0.8,
#'   cluster_algorithm = "leiden"
#' )
#'
#' seurat_obj <- sn_run_cluster(
#'   object = seurat_obj,
#'   batch = "sample_id",
#'   integration_method = "harmony",
#'   normalization_method = "seurat",
#'   hvg_group_by = "sample_id",
#'   nfeatures = 3000,
#'   resolution = 0.5,
#'   block_genes = c("ribo", "mito") # or a custom vector of gene symbols
#' )
#' }
#' @export
sn_run_cluster <- function(object,
                           batch = NULL,
                           normalization_method = c("seurat", "scran", "sctransform"),
                           integration_method = c("harmony", "unintegrated", "coralysis", "seurat_cca", "seurat_rpca", "scvi", "scanvi", "scpoli", "bbknn", "totalvi", "mmochi"),
                           integration_control = list(),
                           nfeatures = 3000,
                           hvg_features = NULL,
                           vars_to_regress = NULL,
                           resolution = 0.8,
                           cluster_algorithm = c("louvain", "louvain_multilevel", "slm", "leiden"),
                           cluster_name = NULL,
                           cluster_n_start = 10,
                           cluster_n_iter = 10,
                           cluster_random_seed = 717,
                           cluster_group_singletons = TRUE,
                           leiden_method = c("leidenbase", "igraph"),
                           leiden_objective_function = c("modularity", "CPM"),
                           seed = NULL,
                           verbose = TRUE,
                           ...) {
  integration_method_supplied <- !missing(integration_method)
  normalization_method_supplied <- !missing(normalization_method)
  nfeatures_supplied <- !missing(nfeatures)
  resolution_supplied <- !missing(resolution)
  cluster_algorithm_supplied <- !missing(cluster_algorithm)
  tail <- .sn_resolve_cluster_tail_args(list(...))
  block_genes_supplied <- "block_genes" %in% tail$supplied
  if (!is.null(seed)) {
    cluster_random_seed <- seed
    if (!is.null(integration_control$seed)) {
      integration_control$seed <- seed
    }
  }
  if (!"verbose" %in% tail$supplied) {
    tail$values$verbose <- verbose
  }

  cluster_args <- c(
    list(
      object = object,
      batch = batch,
      normalization_method = normalization_method,
      integration_method = integration_method,
      integration_control = integration_control,
      nfeatures = nfeatures,
      hvg_features = hvg_features,
      vars_to_regress = vars_to_regress,
      resolution = resolution,
      cluster_algorithm = cluster_algorithm,
      cluster_name = cluster_name,
      cluster_n_start = cluster_n_start,
      cluster_n_iter = cluster_n_iter,
      cluster_random_seed = cluster_random_seed,
      cluster_group_singletons = cluster_group_singletons,
      leiden_method = leiden_method,
      leiden_objective_function = leiden_objective_function
    ),
    tail$values,
    list(
      .integration_method_supplied = integration_method_supplied,
      .normalization_method_supplied = normalization_method_supplied,
      .nfeatures_supplied = nfeatures_supplied,
      .resolution_supplied = resolution_supplied,
      .cluster_algorithm_supplied = cluster_algorithm_supplied,
      .block_genes_supplied = block_genes_supplied,
      .run_tsne_supplied = "run_tsne" %in% tail$supplied,
      .tail_supplied = tail$supplied
    )
  )

  grid_requested <-
    (integration_method_supplied && length(integration_method) > 1L) ||
    (normalization_method_supplied && length(normalization_method) > 1L) ||
    (nfeatures_supplied && length(nfeatures) > 1L) ||
    (resolution_supplied && length(resolution) > 1L) ||
    (cluster_algorithm_supplied && length(cluster_algorithm) > 1L) ||
    any(vapply(
      intersect(c("rare_feature_n", "theta", "npcs"), tail$supplied),
      function(name) length(tail$values[[name]]) > 1L,
      logical(1)
    ))
  if (grid_requested) {
    return(.sn_run_cluster_multi(cluster_args))
  }
  .sn_run_cluster_impl(cluster_args)
}

.sn_multi_method_control <- function(integration_control, methods, method) {
  if (!is.list(integration_control)) {
    stop("`integration_control` must be a named list.", call. = FALSE)
  }
  control_names <- names(integration_control) %||% character(0)
  mapped <- any(control_names %in% c(.sn_supported_integration_methods(), ".default"))
  if (!mapped) {
    return(integration_control)
  }
  unknown <- setdiff(control_names, c(methods, ".default"))
  if (length(unknown) > 0L) {
    stop(
      "Unknown multi-method `integration_control` name(s): ",
      paste(unknown, collapse = ", "), ".",
      call. = FALSE
    )
  }
  default <- integration_control[[".default"]] %||% list()
  specific <- integration_control[[method]] %||% list()
  if (!is.list(default) || !is.list(specific)) {
    stop("Each per-method `integration_control` entry must be a named list.", call. = FALSE)
  }
  utils::modifyList(default, specific, keep.null = TRUE)
}

.sn_cluster_grid_id_value <- function(value) {
  value <- paste(value, collapse = "-")
  value <- gsub("-", "m", value, fixed = TRUE)
  value <- gsub("\\.", "p", value)
  value <- gsub("[^A-Za-z0-9]+", "-", value)
  gsub("(^-+|-+$)", "", value)
}

.sn_cluster_grid_axes <- function(args) {
  supplied <- args$.tail_supplied %||% character(0)
  choose <- function(value, was_supplied) {
    value <- unique(value)
    if (isTRUE(was_supplied)) value else value[[1L]]
  }
  list(
    normalization_method = choose(
      args$normalization_method,
      args$.normalization_method_supplied
    ),
    nfeatures = choose(args$nfeatures, args$.nfeatures_supplied),
    npcs = choose(args$npcs, "npcs" %in% supplied),
    resolution = choose(args$resolution, args$.resolution_supplied),
    cluster_algorithm = choose(
      args$cluster_algorithm,
      args$.cluster_algorithm_supplied
    ),
    rare_feature_n = choose(args$rare_feature_n, "rare_feature_n" %in% supplied)
  )
}

.sn_expand_cluster_grid <- function(args) {
  id_values <- function(values) {
    vapply(values, .sn_cluster_grid_id_value, character(1))
  }
  methods <- if (isTRUE(args$.integration_method_supplied)) {
    .sn_normalize_integration_methods(args$integration_method)
  } else {
    .sn_normalize_integration_methods(args$integration_method[[1L]])
  }
  axes <- .sn_cluster_grid_axes(args)
  base_grid <- do.call(
    expand.grid,
    c(axes, list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE))
  )
  theta_supplied <- "theta" %in% (args$.tail_supplied %||% character(0))
  rows <- list()
  index <- 0L
  for (method in methods) {
    method_control <- .sn_multi_method_control(
      integration_control = args$integration_control,
      methods = methods,
      method = method
    )
    theta_values <- if (identical(method, "harmony")) {
      unique(method_control$theta %||% if (theta_supplied) args$theta else args$theta[[1L]])
    } else {
      args$theta[[1L]]
    }
    for (theta in theta_values) {
      part <- base_grid
      part$method <- method
      part$theta <- theta
      part$.method_control <- I(rep(list(method_control), nrow(part)))
      index <- index + 1L
      rows[[index]] <- part
    }
  }
  grid <- do.call(rbind, rows)
  rownames(grid) <- NULL
  grid$preprocess_id <- paste0(
    "prep.", id_values(grid$normalization_method),
    ".hvg", id_values(grid$nfeatures),
    ".pc", id_values(grid$npcs),
    ".rare", id_values(grid$rare_feature_n)
  )
  grid$embedding_id <- paste0(
    grid$method, ".", grid$preprocess_id,
    ifelse(
      grid$method == "harmony",
      paste0(".theta", id_values(grid$theta)),
      ""
    )
  )
  axis_lengths <- vapply(axes, length, integer(1))
  expanded <- any(axis_lengths > 1L) ||
    (identical("harmony", methods) && length(unique(grid$theta)) > 1L) ||
    ("harmony" %in% methods && length(unique(grid$theta[grid$method == "harmony"])) > 1L)
  simple_multi <- !expanded && !anyDuplicated(grid$method)
  if (simple_multi) {
    grid$run_id <- grid$method
    grid$embedding_id <- grid$method
    grid$preprocess_id <- "shared"
  } else {
    grid$run_id <- paste0(
      grid$embedding_id,
      ".res", id_values(grid$resolution),
      ".", id_values(grid$cluster_algorithm),
      ifelse(
        length(unique(grid$rare_feature_n)) > 1L,
        paste0(".rare", id_values(grid$rare_feature_n)),
        ""
      )
    )
    grid$run_id <- make.unique(grid$run_id, sep = ".v")
  }
  grid <- grid[
    order(grid$preprocess_id, match(grid$method, methods), grid$embedding_id, grid$resolution, grid$run_id),
    , drop = FALSE
  ]
  rownames(grid) <- NULL
  attr(grid, "simple_multi") <- simple_multi
  grid
}

.sn_alias_cluster_reduction <- function(object, reduction, alias, key_prefix = "integrated") {
  if (is.null(reduction) || !reduction %in% names(object@reductions) || identical(reduction, alias)) {
    return(object)
  }
  aliased <- object[[reduction]]
  SeuratObject::Key(aliased) <- .sn_reduction_key(key_prefix, alias)
  object[[alias]] <- aliased
  object
}

.sn_cluster_existing_grid_graph <- function(object, args, graph_name, cluster_name) {
  algorithm <- .sn_resolve_find_clusters_algorithm(args$cluster_algorithm)
  leiden_method <- match.arg(args$leiden_method, c("leidenbase", "igraph"))
  leiden_objective_function <- match.arg(args$leiden_objective_function, c("modularity", "CPM"))
  .sn_ensure_cluster_algorithm_dependencies(
    cluster_algorithm_value = algorithm,
    leiden_method = leiden_method,
    auto_install = args$auto_install,
    repos = args$install_repos,
    ask = args$install_ask
  )
  cluster_args <- .sn_merge_control_args(
    defaults = list(
      object = object,
      graph.name = graph_name,
      resolution = args$resolution,
      algorithm = algorithm,
      n.start = args$cluster_n_start,
      n.iter = args$cluster_n_iter,
      random.seed = args$cluster_random_seed,
      group.singletons = args$cluster_group_singletons,
      leiden_method = leiden_method,
      leiden_objective_function = leiden_objective_function,
      cluster.name = cluster_name,
      verbose = args$verbose
    ),
    control = args$cluster_control
  )
  .sn_with_default_seurat_acceleration(
    .sn_call_with_symbolic_object(
      fun_call = quote(Seurat::FindClusters),
      object = object,
      args = cluster_args
    ),
    object = object,
    assay = args$assay
  )
}

.sn_cluster_checkpoint_signature <- function(args, grid) {
  signature_args <- args
  signature_args$object <- NULL
  signature_args$checkpoint_dir <- NULL
  signature_args$resume <- NULL
  signature_args$checkpoint_compress <- NULL
  signature_args$verbose <- NULL
  object <- args$object
  object_signature <- list(
    cells = colnames(object),
    features = rownames(object),
    assays = names(object@assays),
    assay = args$assay,
    layer = args$layer,
    analysis_input = .sn_cluster_analysis_input_signature(
      object = object,
      assay = args$assay,
      layer = args$layer
    ),
    layer_class = if (args$assay %in% names(object@assays)) {
      class(.sn_get_seurat_layer_data(object, assay = args$assay, layer = args$layer))
    } else {
      NULL
    },
    batch = if (!is.null(args$batch) && args$batch %in% colnames(object[[]])) {
      as.character(object[[args$batch, drop = TRUE]])
    } else {
      NULL
    },
    metadata = .sn_cluster_metadata_signature(
      object = object,
      columns = c(
        args$batch,
        args$hvg_group_by,
        args$rare_feature_group_by,
        args$vars_to_regress,
        args$group_by_vars,
        .sn_cluster_control_metadata_columns(args$integration_control)
      )
    )
  )
  digest::digest(
    list(
      schema_version = 1L,
      package_version = as.character(utils::packageVersion("Shennong")),
      object = object_signature,
      arguments = signature_args,
      grid = grid[, setdiff(colnames(grid), ".method_control"), drop = FALSE]
    ),
    algo = "sha256",
    serialize = TRUE
  )
}

.sn_cluster_checkpoint_files <- function(checkpoint_dir, signature) {
  pattern <- paste0("^sn-run-cluster-", substr(signature, 1L, 16L), "-[0-9]{6}\\.rds$")
  list.files(checkpoint_dir, pattern = pattern, full.names = TRUE)
}

.sn_read_cluster_checkpoint <- function(checkpoint_dir, signature) {
  files <- .sn_cluster_checkpoint_files(checkpoint_dir, signature)
  if (length(files) == 0L) return(NULL)
  indices <- suppressWarnings(as.integer(sub("^.*-([0-9]{6})\\.rds$", "\\1", files)))
  files <- files[order(indices, decreasing = TRUE)]
  for (path in files) {
    state <- tryCatch(readRDS(path), error = function(error) NULL)
    if (is.list(state) && identical(state$signature, signature)) {
      state$checkpoint_path <- normalizePath(path, winslash = "/", mustWork = TRUE)
      return(state)
    }
  }
  NULL
}

.sn_write_cluster_checkpoint <- function(state, checkpoint_dir, compress = FALSE) {
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  prefix <- substr(state$signature, 1L, 16L)
  target <- file.path(
    checkpoint_dir,
    sprintf("sn-run-cluster-%s-%06d.rds", prefix, as.integer(state$completed_index))
  )
  partial <- paste0(target, ".", Sys.getpid(), ".partial")
  on.exit(if (file.exists(partial)) unlink(partial), add = TRUE)
  saveRDS(state, file = partial, compress = compress)
  if (!file.rename(partial, target)) {
    stop("Could not atomically publish clustering checkpoint: ", target, call. = FALSE)
  }
  previous <- setdiff(.sn_cluster_checkpoint_files(checkpoint_dir, state$signature), target)
  if (length(previous) > 0L) unlink(previous)
  normalizePath(target, winslash = "/", mustWork = TRUE)
}

.sn_cluster_performance_table <- function(results) {
  rows <- lapply(results, function(result) {
    if (is.null(result)) return(NULL)
    workflow <- result$performance$workflow %||% list()
    integration <- result$performance$integration %||% list()
    data.frame(
      run_id = result$run_id,
      embedding_id = result$embedding_id,
      preprocess_id = result$preprocess_id,
      method = result$method,
      elapsed_seconds = workflow$elapsed_seconds %||% NA_real_,
      peak_memory_mb = workflow$peak_memory_mb %||% NA_real_,
      integration_elapsed_seconds = integration$elapsed_seconds %||% NA_real_,
      integration_peak_memory_mb = integration$peak_memory_mb %||% NA_real_,
      integration_peak_r_memory_mb = integration$peak_r_memory_mb %||% NA_real_,
      integration_backend_peak_rss_mb = integration$backend_peak_rss_mb %||% NA_real_,
      reused_embedding = isTRUE(result$performance$reused_embedding),
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) return(data.frame())
  do.call(rbind, rows)
}

.sn_run_cluster_multi <- function(args) {
  grid <- .sn_expand_cluster_grid(args)
  methods <- unique(grid$method)
  if (nrow(grid) > 100L) {
    warning(
      "The requested clustering parameter grid contains ", nrow(grid),
      " runs and may require substantial compute and storage.",
      call. = FALSE
    )
  }
  if (is.null(args$batch)) {
    stop("Parameter-grid and multi-method runs currently require `batch`.", call. = FALSE)
  }
  modality <- match.arg(args$modality, c("rna", "cite_seq"))
  if (!identical(modality, "rna")) {
    stop("Multiple `integration_method` values currently support `modality = \"rna\"` only.", call. = FALSE)
  }
  if (any(methods %in% c("totalvi", "mmochi"))) {
    stop("`totalvi` and `mmochi` cannot be used in an RNA multi-method comparison.", call. = FALSE)
  }

  return_cluster <- isTRUE(args$return_cluster)
  base_cluster_name <- args$cluster_name
  checkpoint_dir <- args$checkpoint_dir %||% NULL
  resume <- isTRUE(args$resume)
  checkpoint_compress <- args$checkpoint_compress
  if (!is.null(checkpoint_dir) && (!is.character(checkpoint_dir) || length(checkpoint_dir) != 1L || !nzchar(checkpoint_dir))) {
    stop("`checkpoint_dir` must be `NULL` or one non-empty directory path.", call. = FALSE)
  }
  if (!is.logical(args$resume) || length(args$resume) != 1L || is.na(args$resume)) {
    stop("`resume` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  if (!is.logical(checkpoint_compress) || length(checkpoint_compress) != 1L || is.na(checkpoint_compress)) {
    stop("`checkpoint_compress` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  checkpoint_signature <- .sn_cluster_checkpoint_signature(args, grid)
  checkpoint_dir <- if (!is.null(checkpoint_dir)) {
    dir.create(path.expand(checkpoint_dir), recursive = TRUE, showWarnings = FALSE)
    normalizePath(path.expand(checkpoint_dir), winslash = "/", mustWork = TRUE)
  } else {
    NULL
  }
  comparison <- list(
    schema_version = "2.1.0",
    run_ids = grid$run_id,
    methods = methods,
    batch_by = args$batch,
    assay = args$assay,
    layer = args$layer,
    grid = grid[, setdiff(colnames(grid), ".method_control"), drop = FALSE],
    results = stats::setNames(vector("list", nrow(grid)), grid$run_id),
    performance = data.frame(),
    checkpoint = list(
      enabled = !is.null(checkpoint_dir),
      directory = checkpoint_dir,
      signature = checkpoint_signature,
      resumed = FALSE,
      completed_run_ids = character(0),
      latest_path = NULL
    ),
    created_at = as.character(Sys.time())
  )
  current <- args$object
  simple_multi <- isTRUE(attr(grid, "simple_multi"))
  previous_preprocess_id <- NULL
  normalization_layers <- list()
  start_index <- 1L

  checkpoint_state <- if (!is.null(checkpoint_dir) && resume) {
    .sn_read_cluster_checkpoint(checkpoint_dir, checkpoint_signature)
  } else {
    NULL
  }
  if (!is.null(checkpoint_state)) {
    current <- checkpoint_state$object
    comparison <- checkpoint_state$comparison %||%
      current@misc$integration_comparison
    if (is.null(comparison)) {
      stop("The clustering checkpoint is missing its comparison manifest.", call. = FALSE)
    }
    comparison$checkpoint$resumed <- TRUE
    comparison$checkpoint$latest_path <- checkpoint_state$checkpoint_path
    previous_preprocess_id <- checkpoint_state$previous_preprocess_id %||% NULL
    normalization_layers <- checkpoint_state$normalization_layers %||% list()
    start_index <- as.integer(checkpoint_state$completed_index) + 1L
    if (isTRUE(args$verbose)) {
      .sn_log_info(
        "[sn_run_cluster] Resuming after {length(comparison$checkpoint$completed_run_ids)} completed run(s) from {checkpoint_state$checkpoint_path}."
      )
    }
  }

  if (isTRUE(args$verbose)) {
    .sn_log_info(
      "[sn_run_cluster] Expanding {nrow(grid)} run(s) across ",
      "{length(unique(grid$embedding_id))} unique embedding(s)."
    )
  }

  row_indices <- if (start_index <= nrow(grid)) seq.int(start_index, nrow(grid)) else integer(0)
  for (row_index in row_indices) {
    row <- grid[row_index, , drop = FALSE]
    method <- row$method[[1L]]
    run_id <- row$run_id[[1L]]
    embedding_id <- row$embedding_id[[1L]]
    preprocess_id <- row$preprocess_id[[1L]]
    if (!is.null(previous_preprocess_id) && !identical(previous_preprocess_id, preprocess_id)) {
      assay_layers <- SeuratObject::Layers(current[[args$assay]])
      if ("scale.data" %in% assay_layers) {
        current[[args$assay]]["scale.data"] <- NULL
      }
    }
    method_args <- args
    method_control <- row$.method_control[[1L]]
    method_args$object <- current
    method_args$normalization_method <- row$normalization_method[[1L]]
    method_args$integration_method <- method
    method_args$integration_control <- method_control
    method_args$nfeatures <- row$nfeatures[[1L]]
    method_args$npcs <- row$npcs[[1L]]
    method_args$resolution <- row$resolution[[1L]]
    method_args$cluster_algorithm <- row$cluster_algorithm[[1L]]
    method_args$rare_feature_n <- row$rare_feature_n[[1L]]
    method_args$theta <- row$theta[[1L]]
    method_args$group_by_vars <- method_control$group_by_vars %||% args$group_by_vars
    method_args$integration_control$theta <- NULL
    method_args$integration_control$group_by_vars <- NULL
    if (identical(method, "bbknn")) {
      method_args$integration_control$graph_name <- paste0(embedding_id, "_snn")
      method_args$integration_control$umap_reduction <- paste0("umap.", embedding_id)
    }
    method_args$cluster_name <- if (is.null(base_cluster_name)) {
      paste0(run_id, "_clusters")
    } else {
      paste0(base_cluster_name, ".", run_id)
    }
    method_args$umap_control <- utils::modifyList(
      args$umap_control,
      list(
        reduction.name = paste0("umap.", embedding_id),
        reduction.key = .sn_reduction_key("umap", embedding_id)
      ),
      keep.null = TRUE
    )
    method_args$run_tsne <- isTRUE(args$run_tsne)
    method_args$tsne_control <- utils::modifyList(
      args$tsne_control,
      list(
        reduction.name = paste0("tsne.", embedding_id),
        reduction.key = .sn_reduction_key("tsne", embedding_id)
      ),
      keep.null = TRUE
    )
    method_args$return_cluster <- return_cluster
    method_args$.return_object_for_multi <- return_cluster
    method_args$.result_namespace <- embedding_id

    existing_run <- names(comparison$results)[vapply(
      comparison$results,
      function(info) !is.null(info) && identical(info$embedding_id, embedding_id),
      logical(1)
    )]
    if (length(existing_run) > 0L) {
      existing <- comparison$results[[existing_run[[1L]]]]
      snn_candidates <- existing$graph_names[grepl("_snn$", existing$graph_names)]
      graph_name <- if (length(snn_candidates) > 0L) {
        snn_candidates[[1L]]
      } else {
        existing$graph_names[[1L]]
      }
      profiled_run <- .sn_cluster_performance(.sn_cluster_existing_grid_graph(
        object = current,
        args = method_args,
        graph_name = graph_name,
        cluster_name = method_args$cluster_name
      ))
      current <- profiled_run$value
      comparison$results[[run_id]] <- utils::modifyList(
        existing,
        list(
          run_id = run_id,
          cluster_column = method_args$cluster_name,
          parameters = as.list(row[, setdiff(colnames(row), c(".method_control", "run_id", "embedding_id", "preprocess_id", "method")), drop = FALSE]),
          performance = list(
            workflow = profiled_run$performance,
            integration = existing$performance$integration %||% list(),
            reused_embedding = TRUE,
            source_run_id = existing$run_id
          )
        ),
        keep.null = TRUE
      )
      previous_preprocess_id <- preprocess_id
      comparison$performance <- .sn_cluster_performance_table(comparison$results)
      comparison$checkpoint$completed_run_ids <- grid$run_id[seq_len(row_index)]
      current@misc$integration_comparison <- comparison
      if (!is.null(checkpoint_dir)) {
        checkpoint_path <- .sn_write_cluster_checkpoint(
          list(
            schema_version = 1L,
            signature = checkpoint_signature,
            completed_index = row_index,
            object = current,
            previous_preprocess_id = previous_preprocess_id,
            normalization_layers = normalization_layers,
            saved_at = as.character(Sys.time())
          ),
          checkpoint_dir = checkpoint_dir,
          compress = checkpoint_compress
        )
        comparison$checkpoint$latest_path <- checkpoint_path
        current@misc$integration_comparison <- comparison
      }
      next
    }

    if (isTRUE(args$verbose)) {
      .sn_log_info("[sn_run_cluster] Running grid entry '{run_id}'.")
    }
    profiled_run <- .sn_cluster_performance(.sn_run_cluster_impl(method_args))
    current <- profiled_run$value
    stages <- current@misc$sn_run_cluster$stages %||% list()
    normalization_method <- row$normalization_method[[1L]]
    normalized_layer <- "data"
    if (length(unique(grid$normalization_method)) > 1L) {
      normalized_layer <- normalization_layers[[normalization_method]] %||%
        paste0("data.sn.", .sn_cluster_grid_id_value(normalization_method))
      if (!normalized_layer %in% SeuratObject::Layers(current[[args$assay]])) {
        SeuratObject::LayerData(
          object = current,
          assay = args$assay,
          layer = normalized_layer
        ) <- .sn_get_seurat_layer_data(current, assay = args$assay, layer = "data")
      }
      normalization_layers[[normalization_method]] <- normalized_layer
    }
    native_reduction <- stages$integration$reduction %||%
      if (identical(method, "unintegrated")) "pca" else method
    stored_reduction <- if (simple_multi) native_reduction else paste0("integrated.", embedding_id)
    current <- .sn_alias_cluster_reduction(
      current,
      reduction = native_reduction,
      alias = stored_reduction,
      key_prefix = "integrated"
    )
    selected_features <- stages$integration$signature$features %||%
      SeuratObject::VariableFeatures(current[[args$assay]])
    comparison$results[[run_id]] <- list(
      run_id = run_id,
      embedding_id = embedding_id,
      preprocess_id = preprocess_id,
      method = method,
      parameters = as.list(row[, setdiff(colnames(row), c(".method_control", "run_id", "embedding_id", "preprocess_id", "method")), drop = FALSE]),
      integration_reduction = stored_reduction,
      cluster_column = method_args$cluster_name,
      graph_names = stages$neighbors$graph_names %||% character(0),
      umap_reduction = if (!return_cluster) stages$umap$reduction %||% NULL else NULL,
      tsne_reduction = if (!return_cluster && isTRUE(method_args$run_tsne)) {
        stages$tsne$reduction %||% NULL
      } else {
        NULL
      },
      input_features = selected_features,
      normalized_layer = normalized_layer,
      integration_control = method_control,
      integration = current@misc$integration %||% NULL,
      performance = list(
        workflow = profiled_run$performance,
        integration = stages$integration$performance %||%
          current@misc$integration$performance %||% list(),
        reused_embedding = FALSE,
        source_run_id = run_id
      )
    )
    previous_preprocess_id <- preprocess_id
    comparison$performance <- .sn_cluster_performance_table(comparison$results)
    comparison$checkpoint$completed_run_ids <- grid$run_id[seq_len(row_index)]
    current@misc$integration_comparison <- comparison
    if (!is.null(checkpoint_dir)) {
      checkpoint_path <- .sn_write_cluster_checkpoint(
        list(
          schema_version = 1L,
          signature = checkpoint_signature,
          completed_index = row_index,
          object = current,
          previous_preprocess_id = previous_preprocess_id,
          normalization_layers = normalization_layers,
          saved_at = as.character(Sys.time())
        ),
        checkpoint_dir = checkpoint_dir,
        compress = checkpoint_compress
      )
      comparison$checkpoint$latest_path <- checkpoint_path
      current@misc$integration_comparison <- comparison
    }
  }

  current@misc$integration_comparison <- comparison
  if (return_cluster) {
    columns <- vapply(comparison$results, `[[`, character(1), "cluster_column")
    return(current[[]][, columns, drop = FALSE])
  }
  .sn_log_seurat_command(object = current, assay = args$assay, name = "sn_run_cluster_multi")
}

.sn_run_cluster_impl <- function(args) {
  object <- args$object
  batch <- args$batch
  normalization_method <- args$normalization_method
  integration_method <- args$integration_method
  integration_control <- args$integration_control
  nfeatures <- args$nfeatures
  hvg_features <- args$hvg_features
  vars_to_regress <- args$vars_to_regress
  resolution <- args$resolution
  cluster_algorithm <- args$cluster_algorithm
  cluster_name <- args$cluster_name
  cluster_n_start <- args$cluster_n_start
  cluster_n_iter <- args$cluster_n_iter
  cluster_random_seed <- args$cluster_random_seed
  cluster_group_singletons <- args$cluster_group_singletons
  leiden_method <- args$leiden_method
  leiden_objective_function <- args$leiden_objective_function
  cluster_control <- args$cluster_control
  reuse <- args$reuse
  rerun_from <- args$rerun_from
  auto_install <- args$auto_install
  install_repos <- args$install_repos
  install_ask <- args$install_ask
  hvg_group_by <- args$hvg_group_by
  rare_feature_method <- args$rare_feature_method
  rare_feature_group_by <- args$rare_feature_group_by
  rare_feature_n <- args$rare_feature_n
  rare_feature_control <- args$rare_feature_control
  block_genes <- args$block_genes
  theta <- args$theta
  group_by_vars <- args$group_by_vars
  npcs <- args$npcs
  dims <- args$dims
  species <- args$species
  assay <- args$assay
  layer <- args$layer
  modality <- args$modality
  multimodal_method <- args$multimodal_method
  adt_assay <- args$adt_assay
  adt_layer <- args$adt_layer
  adt_features <- args$adt_features
  adt_npcs <- args$adt_npcs
  adt_dims <- args$adt_dims
  wnn_control <- args$wnn_control
  umap_control <- args$umap_control
  run_tsne <- args$run_tsne
  tsne_control <- args$tsne_control
  return_cluster <- args$return_cluster
  verbose <- args$verbose
  .integration_method_supplied <- args$.integration_method_supplied
  .block_genes_supplied <- args$.block_genes_supplied
  result_namespace <- args$.result_namespace %||% NULL
  return_object_for_multi <- isTRUE(args$.return_object_for_multi)

  check_installed("Seurat")

  if (!inherits(object, "Seurat")) {
    stop("Input must be a Seurat object.")
  }
  if (inherits(batch, "Seurat")) {
    stop(
      "`batch` received a Seurat object. ",
      "This usually means the input object was supplied twice, for example ",
      "`object %>% sn_run_cluster(object)`. Use `object %>% sn_run_cluster()` ",
      "or `sn_run_cluster(object)` instead.",
      call. = FALSE
    )
  }
  if (!is.null(batch) && (!is.character(batch) || length(batch) != 1L)) {
    stop("`batch` must be a single metadata column name or `NULL`.", call. = FALSE)
  }

  integration_method_supplied <- isTRUE(.integration_method_supplied)
  block_genes_supplied <- isTRUE(.block_genes_supplied)
  normalization_method <- match.arg(
    normalization_method,
    c("seurat", "scran", "sctransform")
  )
  if (!integration_method_supplied) {
    integration_method <- integration_method[[1L]]
  }
  integration_method <- .sn_normalize_integration_methods(integration_method)
  if (length(integration_method) != 1L) {
    stop("Internal clustering calls require one `integration_method`.", call. = FALSE)
  }
  modality <- match.arg(modality, c("rna", "cite_seq"))
  multimodal_method <- if (identical(modality, "cite_seq")) {
    .sn_resolve_multimodal_method(
      multimodal_method = multimodal_method,
      integration_method = integration_method,
      integration_method_supplied = integration_method_supplied
    )
  } else {
    if (!is.null(multimodal_method)) {
      stop("`multimodal_method` is used only when `modality = \"cite_seq\"`.", call. = FALSE)
    }
    NULL
  }
  cluster_algorithm_value <- .sn_resolve_find_clusters_algorithm(cluster_algorithm)
  leiden_method <- match.arg(leiden_method, c("leidenbase", "igraph"))
  leiden_objective_function <- match.arg(
    leiden_objective_function,
    c("modularity", "CPM")
  )
  rerun_from <- .sn_resolve_cluster_rerun_from(rerun_from)
  if (!is.list(integration_control)) {
    stop("`integration_control` must be a named list.", call. = FALSE)
  }
  if (identical(integration_method, "harmony")) {
    theta <- integration_control$theta %||% theta
    group_by_vars <- integration_control$group_by_vars %||% group_by_vars
  }
  if (!is.list(cluster_control)) {
    stop("`cluster_control` must be a named list.", call. = FALSE)
  }
  if (!is.list(wnn_control)) {
    stop("`wnn_control` must be a named list.", call. = FALSE)
  }
  if (!is.list(umap_control)) {
    stop("`umap_control` must be a named list.", call. = FALSE)
  }
  if (!is.logical(run_tsne) || length(run_tsne) != 1L || is.na(run_tsne)) {
    stop("`run_tsne` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  if (!is.list(tsne_control)) {
    stop("`tsne_control` must be a named list.", call. = FALSE)
  }
  rare_feature_method <- unique(match.arg(
    rare_feature_method,
    c("none", "gini", "local_markers"),
    several.ok = TRUE
  ))
  rare_feature_control <- .sn_resolve_rare_feature_control(control = rare_feature_control)
  if (!is_null(hvg_group_by) && !hvg_group_by %in% colnames(object[[]])) {
    stop(glue("`hvg_group_by` must be NULL or a metadata column name. '{hvg_group_by}' was not found."))
  }
  dims <- .sn_resolve_cluster_dims(dims = dims, npcs = npcs)
  adt_feature_info <- NULL
  adt_feature_set <- character(0)
  if (identical(modality, "cite_seq")) {
    if (identical(multimodal_method, "wnn") && !is.null(batch)) {
      stop(
        "`multimodal_method = \"wnn\"` currently supports single-object WNN clustering only. ",
        "For CITE-seq batch-aware workflows, use `multimodal_method = \"totalvi\"`, ",
        "`\"coralysis\"`, or `\"mmochi\"`.",
        call. = FALSE
      )
    }
    if (identical(multimodal_method, "coralysis") && is.null(batch)) {
      stop("`multimodal_method = \"coralysis\"` requires `batch`.", call. = FALSE)
    }
    if (identical(normalization_method, "sctransform") && identical(multimodal_method, "coralysis")) {
      stop("CITE-seq Coralysis protein workflows currently require `normalization_method = \"seurat\"` or `\"scran\"`.", call. = FALSE)
    }
    if (!is.character(adt_assay) || length(adt_assay) != 1L) {
      stop("`adt_assay` must be a single assay name.", call. = FALSE)
    }
    if (!is.character(adt_layer) || length(adt_layer) != 1L) {
      stop("`adt_layer` must be a single layer name.", call. = FALSE)
    }
    adt_npcs <- as.integer(adt_npcs)
    if (length(adt_npcs) != 1L || is.na(adt_npcs) || adt_npcs < 1L) {
      stop("`adt_npcs` must be a positive integer.", call. = FALSE)
    }
    adt_dims <- .sn_resolve_cluster_dims(dims = adt_dims %||% seq_len(min(18L, adt_npcs)), npcs = adt_npcs)
    .sn_validate_seurat_assay_layer(object = object, assay = adt_assay, layer = adt_layer)
    adt_feature_info <- .sn_resolve_assay_features(
      object = object,
      assay = adt_assay,
      features = adt_features,
      arg_name = "adt_features",
      verbose = verbose
    )
    adt_feature_set <- adt_feature_info$features
  } else if (integration_method %in% c("totalvi", "mmochi")) {
    stop("`integration_method = \"totalvi\"` or `\"mmochi\"` requires `modality = \"cite_seq\"`.", call. = FALSE)
  }

  needs_rna_workflow <- .sn_cluster_requires_rna_workflow(
    modality = modality,
    multimodal_method = multimodal_method
  )
  needs_rna_pca <- .sn_cluster_requires_rna_pca(
    modality = modality,
    multimodal_method = multimodal_method,
    batch = batch,
    integration_method = integration_method
  )
  needs_adt_data <- .sn_cluster_requires_adt_data(
    modality = modality,
    multimodal_method = multimodal_method,
    integration_control = integration_control
  )
  needs_adt_pca <- .sn_cluster_requires_adt_pca(
    modality = modality,
    multimodal_method = multimodal_method
  )

  analysis_input_signature <- .sn_cluster_analysis_input_signature(
    object = object,
    assay = assay,
    layer = layer
  )
  prepared <- .sn_prepare_seurat_analysis_input(
    object = object,
    assay = assay,
    layer = layer
  )
  object <- prepared$object
  adt_context <- NULL
  restore_analysis_inputs <- function(current_object) {
    if (!is.null(adt_context)) {
      current_object <- .sn_restore_seurat_analysis_input(object = current_object, context = adt_context)
    }
    .sn_restore_seurat_analysis_input(object = current_object, context = prepared$context)
  }

  if (!is_null(x = batch)) {
    if (!(batch %in% colnames(object@meta.data))) {
      stop(glue("Batch variable '{batch}' not found in metadata."))
    }
    if (normalization_method == "sctransform" && !integration_method %in% c("harmony", "unintegrated")) {
      stop(
        "SCTransform integration is currently supported only with `integration_method = \"harmony\"`; ",
        "use `\"unintegrated\"` only for the uncorrected PCA baseline.",
        call. = FALSE
      )
    }
    if (verbose) {
      .sn_log_info("[sn_run_cluster] Starting {integration_method} integration for batch = '{batch}'.")
    }
  }

  if (isTRUE(needs_rna_workflow) && is_null(hvg_group_by) && !is_null(batch)) {
    hvg_group_by <- batch
  }

  if (verbose) {
    .sn_log_info(
      "[sn_run_cluster] Normalization method = {normalization_method}; ",
      "batch = {batch %||% 'none'}; integration_method = {if (is.null(batch) || identical(modality, 'cite_seq')) 'none' else integration_method}; ",
      "modality = {modality}; multimodal_method = {multimodal_method %||% 'none'}."
    )
  }

  user_hvg_info <- .sn_resolve_user_hvg_features(
    features = hvg_features,
    object = object,
    verbose = verbose
  )
  user_hvg <- user_hvg_info$features

  if (!isTRUE(needs_rna_workflow)) {
    block_genes <- NULL
    rare_feature_method <- "none"
  }

  if (normalization_method == "sctransform" && !identical(rare_feature_method, "none")) {
    .sn_log_warn("`rare_feature_method` is only applied in log-normalization workflows; ignoring it for SCTransform.")
    rare_feature_method <- "none"
  }

  if (isTRUE(needs_rna_workflow) && !is_null(block_genes)) {
    species <- tryCatch(
      sn_get_species(object = object, species = species),
      error = function(error) {
        if (identical(normalization_method, "sctransform") && !isTRUE(block_genes_supplied)) {
          if (verbose) {
            .sn_log_warn(
              "Skipping default `block_genes` in the SCTransform workflow because species could not be inferred. ",
              "Provide `species` to enable blocked HVG filtering."
            )
          }
          return(NULL)
        }
        stop(error)
      }
    )
    if (is_null(species)) {
      block_genes <- NULL
    } else {
      if (verbose) .sn_log_info("Processing blocked genes.")
      block_genes <- .sn_resolve_block_genes(
        block_genes = block_genes,
        species = species,
        verbose = verbose
      )
    }
  }

  hvg_candidate_nfeatures <- nfeatures
  if (!is_null(block_genes)) {
    hvg_candidate_nfeatures <- min(
      nrow(object),
      nfeatures + length(intersect(block_genes, rownames(object)))
    )
    if (identical(normalization_method, "sctransform")) {
      hvg_candidate_nfeatures <- min(
        hvg_candidate_nfeatures,
        nrow(object),
        max(nfeatures, nfeatures * 2L)
      )
    }
  }

  hvg <- character(0)
  hvg_signature <- NULL
  pca_signature <- NULL

  normalization_signature <- list(
    method = normalization_method,
    assay = assay,
    layer = layer,
    layer_source_version = 2L,
    analysis_input = analysis_input_signature,
    nfeatures = if (identical(normalization_method, "sctransform")) nfeatures else NULL,
    hvg_candidate_nfeatures = if (identical(normalization_method, "sctransform")) hvg_candidate_nfeatures else NULL,
    vars_to_regress = if (identical(normalization_method, "sctransform")) vars_to_regress else NULL,
    vars_to_regress_metadata = if (identical(normalization_method, "sctransform")) {
      .sn_cluster_metadata_signature(object = object, columns = vars_to_regress)
    } else {
      NULL
    },
    user_hvg = if (identical(normalization_method, "sctransform")) user_hvg else NULL
  )
  can_reuse_normalization <- .sn_can_reuse_cluster_stage(
    object = object,
    stage = "normalize",
    signature = normalization_signature,
    reuse = reuse,
    rerun_from = rerun_from,
    required = function(current_object, stage_info) {
      if (identical(normalization_method, "sctransform")) {
        return("SCT" %in% names(current_object@assays))
      }
      .sn_has_seurat_layer(current_object, assay = assay, layer = "data")
    }
  )

  if (!isTRUE(needs_rna_workflow)) {
    if (verbose) {
      .sn_log_info(
        "[1/6] Skipping RNA normalization, feature selection, scaling, and PCA; ",
        "{multimodal_method} uses the ADT/protein assay directly."
      )
    }
  } else if (normalization_method == "sctransform") {
    if (can_reuse_normalization) {
      if (verbose) .sn_log_info("[1/5] Reusing SCTransform results.")
      SeuratObject::DefaultAssay(object = object) <- "SCT"
    } else {
      check_installed("glmGamPoi", reason = "for the SCTransform workflow.")

      if (verbose) .sn_log_info("[1/5] Running SCTransform.")
      sct_args <- list(
        object = object,
        variable.features.n = hvg_candidate_nfeatures,
        vars.to.regress = vars_to_regress,
        verbose = verbose,
        seed.use = 717
      )
      if (length(user_hvg) > 0L) {
        sct_args$return.only.var.genes <- FALSE
      }
      object <- .sn_with_default_seurat_acceleration(
        .sn_with_auto_future_globals(
          .sn_call_with_symbolic_object(
            fun_call = quote(Seurat::SCTransform),
            object = object,
            args = sct_args
          ),
          object = object,
          context = "SCTransform",
          verbose = verbose
        ),
        object = object,
        assay = assay
      )
      object <- .sn_record_cluster_stage(object, "normalize", normalization_signature)
    }

    hvg_signature <- list(
      method = "sctransform",
      nfeatures = nfeatures,
      hvg_candidate_nfeatures = hvg_candidate_nfeatures,
      block_genes = block_genes,
      user_hvg = user_hvg,
      missing_user_hvg = user_hvg_info$missing,
      normalization = normalization_signature
    )
    blocked_hvg <- character(0)
    if (.sn_can_reuse_cluster_stage(
      object = object,
      stage = "hvg",
      signature = hvg_signature,
      reuse = reuse,
      rerun_from = rerun_from,
      required = function(current_object, stage_info) length(stage_info$selected_features %||% character(0)) > 0L
    )) {
      if (verbose) .sn_log_info("[2/5] Reusing SCTransform feature set.")
      hvg <- object@misc$sn_run_cluster$stages$hvg$selected_features
    } else {
      base_hvg <- Seurat::VariableFeatures(object = object)
      blocked_hvg <- character(0)
      if (!is_null(block_genes)) {
        n_before <- length(base_hvg)
        blocked_hvg <- intersect(base_hvg, block_genes)
        base_hvg <- setdiff(base_hvg, block_genes)
        base_hvg <- utils::head(base_hvg, nfeatures)

        if (verbose) {
          .sn_log_info(
            "  Removed {length(blocked_hvg)} genes ({round(length(blocked_hvg)/n_before*100,1)}%) via block_genes"
          )
        }

        if (length(base_hvg) < nfeatures) {
          .sn_log_warn(
            "Only {length(base_hvg)} SCTransform HVGs left after block_genes filtering ",
            "(< requested {nfeatures}). Consider adjusting 'nfeatures' or 'block_genes'."
          )
        }
      } else {
        base_hvg <- utils::head(base_hvg, nfeatures)
      }
      hvg <- unique(c(base_hvg, user_hvg))
      object <- .sn_record_cluster_stage(
        object,
        "hvg",
        hvg_signature,
        selected_features = hvg,
        blocked_features = blocked_hvg
      )
    }
    Seurat::VariableFeatures(object = object) <- hvg
    object@misc$hvg_selection <- list(
      method = "sctransform",
      base_hvg_n = nfeatures,
      hvg_candidate_nfeatures = hvg_candidate_nfeatures,
      selected_features = hvg,
      blocked_features = object@misc$sn_run_cluster$stages$hvg$blocked_features %||% character(0),
      user_features = user_hvg,
      missing_user_features = user_hvg_info$missing
    )

    if (isTRUE(needs_rna_pca)) {
      pca_signature <- list(
        method = "sctransform",
        features = hvg,
        npcs = npcs,
        vars_to_regress = vars_to_regress,
        normalization = normalization_signature,
        hvg = hvg_signature
      )
      if (.sn_can_reuse_cluster_stage(
        object = object,
        stage = "pca",
        signature = pca_signature,
        reuse = reuse,
        rerun_from = rerun_from,
        required = function(current_object, stage_info) "pca" %in% names(current_object@reductions)
      )) {
        if (verbose) .sn_log_info("[3/5] Reusing PCA reduction.")
      } else {
        if (verbose) .sn_log_info("[3/5] Running PCA.")
        object <- .sn_with_default_seurat_acceleration(
          Seurat::RunPCA(
            object,
            npcs = npcs,
            features = hvg,
            verbose = verbose,
            seed.use = 717
          ),
          object = object,
          assay = "SCT"
        )
        object <- .sn_record_cluster_stage(object, "pca", pca_signature, reduction = "pca")
      }
    } else if (verbose) {
      .sn_log_info("[3/5] Skipping RNA PCA; the selected backend imports or computes its own latent representation.")
    }
  } else {
    if (can_reuse_normalization) {
      if (verbose) .sn_log_info("[1/6] Reusing normalized data.")
    } else if (normalization_method == "scran") {
      object <- sn_normalize_data(
        object = object,
        method = "scran",
        assay = assay,
        layer = layer
      )
      object <- .sn_record_cluster_stage(object, "normalize", normalization_signature)
    } else {
      standard_counts <- identical(layer, "counts") &&
        !isTRUE(prepared$context$needs_temp_counts)
      object <- if (standard_counts) {
        .sn_with_default_seurat_acceleration(
          Seurat::NormalizeData(
            object = object,
            assay = assay,
            verbose = verbose
          ),
          object = object,
          assay = assay
        )
      } else {
        .sn_with_acceleration_disabled(
          Seurat::NormalizeData(
            object = object,
            assay = assay,
            layer = layer,
            verbose = verbose
          )
        )
      }
      object <- .sn_record_cluster_stage(object, "normalize", normalization_signature)
    }

    if (!is_null(species)) {
      cell_cycle_signature <- list(species = species, normalization = normalization_signature)
      if (.sn_can_reuse_cluster_stage(
        object = object,
        stage = "cell_cycle",
        signature = cell_cycle_signature,
        reuse = reuse,
        rerun_from = rerun_from,
        required = function(current_object, stage_info) {
          all(c("S.Score", "G2M.Score", "Phase", "CC.Difference") %in% colnames(current_object[[]]))
        }
      )) {
        if (verbose) .sn_log_info("[2/6] Reusing cell-cycle scores.")
      } else {
        if (verbose) .sn_log_info("[2/6] Scoring cell cycle.")
        object <- sn_score_cell_cycle(object = object, species = species)
        object <- .sn_record_cluster_stage(object, "cell_cycle", cell_cycle_signature)
      }
    }

    hvg_signature <- list(
      method = normalization_method,
      nfeatures = nfeatures,
      hvg_candidate_nfeatures = hvg_candidate_nfeatures,
      hvg_group_by = hvg_group_by,
      block_genes = block_genes,
      rare_feature_method = rare_feature_method,
      rare_feature_group_by = rare_feature_group_by,
      rare_feature_n = rare_feature_n,
      rare_feature_control = rare_feature_control,
      rare_feature_resolution = if (identical(rare_feature_method, "none")) {
        NULL
      } else {
        min(resolution, 0.4)
      },
      grouping_metadata = .sn_cluster_metadata_signature(
        object = object,
        columns = c(hvg_group_by, rare_feature_group_by)
      ),
      user_hvg = user_hvg,
      missing_user_hvg = user_hvg_info$missing,
      normalization = normalization_signature
    )
    blocked_hvg <- character(0)
    if (.sn_can_reuse_cluster_stage(
      object = object,
      stage = "hvg",
      signature = hvg_signature,
      reuse = reuse,
      rerun_from = rerun_from,
      required = function(current_object, stage_info) length(stage_info$selected_features %||% character(0)) > 0L
    )) {
      if (verbose) .sn_log_info("[3/6] Reusing selected feature set.")
      stage_info <- object@misc$sn_run_cluster$stages$hvg
      hvg <- stage_info$selected_features
      blocked_hvg <- stage_info$blocked_features %||% character(0)
      selected_rare_features <- stage_info$rare_features %||% character(0)
      object@misc$rare_feature_selection <- stage_info$rare_feature_selection %||% object@misc$rare_feature_selection
    } else {
      if (verbose) {
        .sn_log_info("[3/6] Selecting highly variable features with hvg_group_by = {hvg_group_by %||% 'global'}.")
      }
      hvg_info <- .sn_select_variable_features(
        object = object,
        nfeatures = hvg_candidate_nfeatures,
        split_by = hvg_group_by,
        assay = assay,
        layer = "data",
        verbose = verbose
      )
      object <- hvg_info$object
      hvg <- hvg_info$features

      if (!is_null(block_genes)) {
        n_before <- length(hvg)
        blocked_hvg <- intersect(hvg, block_genes)
        hvg <- setdiff(hvg, block_genes)
        n_after_filter <- length(hvg)
        hvg <- utils::head(hvg, nfeatures)
        n_removed <- n_before - n_after_filter

        if (verbose) {
          .sn_log_info(
            "  Removed {n_removed} genes ({round(n_removed/n_before*100,1)}%) via block_genes"
          )
        }

        if (length(hvg) < nfeatures) {
          .sn_log_warn(
            "Only {length(hvg)} HVGs left (< requested {nfeatures}).\n",
            "Consider adjusting 'nfeatures' or 'block_genes'."
          )
        }
        hvg <- utils::head(hvg, nfeatures)
      }

      rare_feature_info <- .sn_select_rare_features(
        object = object,
        base_features = hvg,
        method = rare_feature_method,
        assay = assay,
        layer = "data",
        nfeatures = rare_feature_n,
        group_by = rare_feature_group_by,
        control = rare_feature_control,
        min_cells = rare_feature_control$min_cells,
        npcs = min(npcs, 20),
        dims = seq_len(min(10, npcs)),
        resolution = min(resolution, 0.4),
        verbose = verbose
      )
      selected_rare_features <- rare_feature_info$features
      hvg <- unique(c(hvg, selected_rare_features, user_hvg))
      rare_feature_store <- list(
        method = rare_feature_method,
        base_hvg_n = nfeatures,
        rare_feature_n = rare_feature_n,
        control = rare_feature_control,
        selected_rare_feature_n = length(selected_rare_features),
        selected_features = hvg,
        rare_features = selected_rare_features,
        rare_feature_table = rare_feature_info$metadata,
        rare_groups = if (!is.null(rare_feature_info$groups)) rare_feature_info$groups$rare_groups else character(0)
      )
      object@misc$rare_feature_selection <- rare_feature_store
      object <- .sn_record_cluster_stage(
        object,
        "hvg",
        hvg_signature,
        selected_features = hvg,
        blocked_features = blocked_hvg,
        rare_features = selected_rare_features,
        rare_feature_selection = rare_feature_store
      )
    }
    Seurat::VariableFeatures(object = object) <- hvg
    object@misc$hvg_selection <- list(
      method = normalization_method,
      hvg_group_by = hvg_group_by,
      base_hvg_n = nfeatures,
      selected_features = hvg,
      blocked_features = blocked_hvg,
      user_features = user_hvg,
      missing_user_features = user_hvg_info$missing
    )

    if (isTRUE(needs_rna_pca)) {
      pca_signature <- list(
        method = normalization_method,
        features = hvg,
        vars_to_regress = vars_to_regress,
        regression_metadata = .sn_cluster_metadata_signature(
          object = object,
          columns = vars_to_regress
        ),
        npcs = npcs,
        normalization = normalization_signature,
        hvg = hvg_signature
      )
      if (.sn_can_reuse_cluster_stage(
        object = object,
        stage = "pca",
        signature = pca_signature,
        reuse = reuse,
        rerun_from = rerun_from,
        required = function(current_object, stage_info) "pca" %in% names(current_object@reductions)
      )) {
        if (verbose) .sn_log_info("[4/6] Reusing scaled data and PCA reduction.")
      } else {
        if (verbose) .sn_log_info("[4/6] Scaling data.")
        object <- .sn_with_default_seurat_acceleration(
          Seurat::ScaleData(
            object = object,
            vars.to.regress = vars_to_regress,
            features = hvg,
            verbose = verbose
          ),
          object = object,
          assay = assay
        )

        if (verbose) .sn_log_info("[5/6] Running PCA.")
        object <- .sn_with_default_seurat_acceleration(
          Seurat::RunPCA(
            object,
            npcs = npcs,
            features = hvg,
            verbose = verbose,
            seed.use = 717
          ),
          object = object,
          assay = assay
        )
        object <- .sn_record_cluster_stage(object, "pca", pca_signature, reduction = "pca")
      }
    } else if (verbose) {
      .sn_log_info("[4/6] Skipping RNA scaling and PCA; the selected backend does not use a Seurat PCA reduction.")
    }
  }

  adt_signature <- NULL
  if (isTRUE(needs_adt_data)) {
    adt_signature <- list(
      modality = modality,
      multimodal_method = multimodal_method,
      assay = adt_assay,
      layer = adt_layer,
      layer_source_version = 2L,
      normalization_method = "CLR",
      margin = 2L,
      npcs = adt_npcs,
      run_pca = isTRUE(needs_adt_pca),
      features = adt_feature_set,
      missing_features = adt_feature_info$missing
    )
    can_reuse_adt <- .sn_can_reuse_cluster_stage(
      object = object,
      stage = "adt",
      signature = adt_signature,
      reuse = reuse,
      rerun_from = rerun_from,
      required = function(current_object, stage_info) {
        has_data <- .sn_has_seurat_layer(current_object, assay = adt_assay, layer = "data")
        if (!isTRUE(needs_adt_pca)) {
          return(has_data)
        }
        has_data &&
          "apca" %in% names(current_object@reductions) &&
          .sn_has_seurat_layer(current_object, assay = adt_assay, layer = "scale.data")
      }
    )

    adt_prepared <- .sn_prepare_seurat_analysis_input(
      object = object,
      assay = adt_assay,
      layer = adt_layer
    )
    object <- adt_prepared$object
    adt_context <- adt_prepared$context

    if (can_reuse_adt) {
      if (verbose) {
        if (isTRUE(needs_adt_pca)) {
          .sn_log_info("[5/6] Reusing ADT normalization and PCA.")
        } else {
          .sn_log_info("[5/6] Reusing ADT normalization.")
        }
      }
    } else {
      if (verbose) {
        if (isTRUE(needs_adt_pca)) {
          .sn_log_info("[5/6] Running ADT CLR normalization and PCA.")
        } else {
          .sn_log_info("[5/6] Running ADT CLR normalization.")
        }
      }
      object <- .sn_with_acceleration_disabled(
        Seurat::NormalizeData(
          object = object,
          assay = adt_assay,
          normalization.method = "CLR",
          margin = 2,
          verbose = verbose
        )
      )
      if (isTRUE(needs_adt_pca)) {
        object <- .sn_with_default_seurat_acceleration(
          Seurat::ScaleData(
            object = object,
            assay = adt_assay,
            features = adt_feature_set,
            verbose = verbose
          ),
          object = object,
          assay = adt_assay
        )
        object <- .sn_with_default_seurat_acceleration(
          Seurat::RunPCA(
            object = object,
            assay = adt_assay,
            features = adt_feature_set,
            npcs = adt_npcs,
            reduction.name = "apca",
            reduction.key = "apca_",
            verbose = verbose,
            seed.use = 717
          ),
          object = object,
          assay = adt_assay
        )
        object <- .sn_record_cluster_stage(object, "adt", adt_signature, reduction = "apca")
      } else {
        object <- .sn_record_cluster_stage(object, "adt", adt_signature)
      }
    }
    SeuratObject::DefaultAssay(object = object) <- if (isTRUE(needs_rna_workflow) && identical(normalization_method, "sctransform")) "SCT" else assay
  }

  integration_graph <- NULL
  integration_umap <- NULL
  integration_performance <- NULL
  backend_integration_control <- integration_control
  integration_metadata_signature <- .sn_cluster_metadata_signature(
    object = object,
    columns = c(
      batch,
      group_by_vars,
      .sn_cluster_control_metadata_columns(integration_control)
    )
  )
  if (identical(integration_method, "bbknn")) {
    backend_integration_control$umap_args <- umap_control
  }
  if (identical(modality, "cite_seq") && identical(multimodal_method, "wnn")) {
    reduction <- "wnn"
    integration_signature <- NULL
    reduction_signature <- list(
      method = "cite_seq_wnn",
      rna_reduction = "pca",
      adt_reduction = "apca",
      rna_dims = dims,
      adt_dims = adt_dims,
      pca = pca_signature,
      adt = adt_signature
    )
  } else if (identical(modality, "cite_seq") && identical(multimodal_method, "coralysis")) {
    integration_signature <- list(
      method = paste0("cite_seq_", multimodal_method),
      batch = batch,
      protein_assay = adt_assay,
      protein_layer = adt_layer,
      protein_features = adt_feature_set,
      dims = adt_dims,
      npcs = adt_npcs,
      integration_control = integration_control,
      store_sce = .sn_coralysis_store_sce(integration_control),
      adt = adt_signature
    )
    if (.sn_can_reuse_cluster_stage(
      object = object,
      stage = "integration",
      signature = integration_signature,
      reuse = reuse,
      rerun_from = rerun_from,
      required = function(current_object, stage_info) {
        !is.null(stage_info$reduction) && stage_info$reduction %in% names(current_object@reductions)
      }
    )) {
      reduction <- object@misc$sn_run_cluster$stages$integration$reduction
      if (verbose) .sn_log_info("[5/6] Reusing {reduction} protein integration reduction.")
    } else {
      if (verbose) .sn_log_info("[5/6] Running {multimodal_method} protein integration.")
      profiled_integration <- .sn_cluster_performance(.sn_run_batch_integration(
        object = object,
        method = multimodal_method,
        batch = batch,
        reduction = NULL,
        features = adt_feature_set,
        assay = adt_assay,
        layer = adt_layer,
        dims = adt_dims,
        npcs = adt_npcs,
        theta = theta,
        group_by_vars = group_by_vars,
        integration_control = integration_control,
        verbose = verbose
      ))
      integration <- profiled_integration$value
      object <- integration$object
      integration_performance <- .sn_complete_integration_performance(
        profiled_integration$performance,
        object
      )
      object@misc$integration$performance <- integration_performance
      reduction <- integration$reduction
      object@misc$integration$modality <- "cite_seq"
      object@misc$integration$protein_assay <- adt_assay
      object@misc$integration$protein_layer <- adt_layer
      object <- .sn_record_cluster_stage(
        object, "integration", integration_signature,
        reduction = reduction, performance = integration_performance
      )
    }
    reduction_signature <- integration_signature
  } else if (identical(modality, "cite_seq") && identical(multimodal_method, "mmochi")) {
    mmochi_control <- utils::modifyList(
      integration_control,
      list(
        protein_layer = integration_control$protein_layer %||% integration_control$adt_layer %||% "data"
      ),
      keep.null = TRUE
    )
    integration_signature <- list(
      method = "cite_seq_mmochi",
      batch = batch,
      protein_assay = adt_assay,
      protein_layer = mmochi_control$protein_layer,
      protein_features = adt_feature_set,
      dims = adt_dims,
      npcs = adt_npcs,
      integration_control = mmochi_control,
      adt = adt_signature
    )
    if (.sn_can_reuse_cluster_stage(
      object = object,
      stage = "integration",
      signature = integration_signature,
      reuse = reuse,
      rerun_from = rerun_from,
      required = function(current_object, stage_info) {
        !is.null(stage_info$reduction) && stage_info$reduction %in% names(current_object@reductions)
      }
    )) {
      reduction <- object@misc$sn_run_cluster$stages$integration$reduction
      if (verbose) .sn_log_info("[5/6] Reusing MMoCHi protein integration reduction.")
    } else {
      if (verbose) .sn_log_info("[5/6] Running MMoCHi ADT landmark registration.")
      profiled_integration <- .sn_cluster_performance(.sn_run_mmochi_integration(
        object = object,
        batch = batch,
        assay = adt_assay,
        protein_assay = adt_assay,
        protein_features = adt_feature_set,
        dims = adt_dims,
        npcs = adt_npcs,
        integration_control = mmochi_control,
        verbose = verbose
      ))
      integration <- profiled_integration$value
      object <- integration$object
      integration_performance <- .sn_complete_integration_performance(
        profiled_integration$performance,
        object
      )
      object@misc$integration$performance <- integration_performance
      reduction <- integration$reduction
      object@misc$integration$modality <- "cite_seq"
      object@misc$integration$protein_assay <- adt_assay
      object@misc$integration$protein_layer <- mmochi_control$protein_layer
      object <- .sn_record_cluster_stage(
        object, "integration", integration_signature,
        reduction = reduction, performance = integration_performance
      )
    }
    reduction_signature <- integration_signature
  } else if (identical(modality, "cite_seq") && identical(multimodal_method, "totalvi")) {
    totalvi_control <- utils::modifyList(
      integration_control,
      list(
        adt_assay = adt_assay,
        adt_layer = adt_layer,
        adt_features = adt_feature_set
      ),
      keep.null = TRUE
    )
    integration_signature <- list(
      method = "cite_seq_totalvi",
      batch = batch,
      features = hvg,
      protein_assay = adt_assay,
      protein_layer = adt_layer,
      protein_features = adt_feature_set,
      integration_control = totalvi_control,
      pca = pca_signature
    )
    if (.sn_can_reuse_cluster_stage(
      object = object,
      stage = "integration",
      signature = integration_signature,
      reuse = reuse,
      rerun_from = rerun_from,
      required = function(current_object, stage_info) {
        !is.null(stage_info$reduction) && stage_info$reduction %in% names(current_object@reductions)
      }
    )) {
      reduction <- object@misc$sn_run_cluster$stages$integration$reduction
      if (verbose) .sn_log_info("[5/6] Reusing totalVI integration reduction.")
    } else {
      if (verbose) .sn_log_info("[5/6] Running totalVI RNA+ADT integration.")
      profiled_integration <- .sn_cluster_performance(.sn_run_scvi_integration(
        object = object,
        method = "totalvi",
        batch = batch,
        features = hvg,
        assay = assay,
        layer = layer,
        integration_control = totalvi_control,
        verbose = verbose
      ))
      integration <- profiled_integration$value
      object <- integration$object
      integration_performance <- .sn_complete_integration_performance(
        profiled_integration$performance,
        object
      )
      object@misc$integration$performance <- integration_performance
      reduction <- integration$reduction
      object@misc$integration$modality <- "cite_seq"
      object@misc$integration$protein_assay <- adt_assay
      object@misc$integration$protein_layer <- adt_layer
      object <- .sn_record_cluster_stage(
        object, "integration", integration_signature,
        reduction = reduction, performance = integration_performance
      )
    }
    reduction_signature <- integration_signature
  } else {
    integration_signature <- list(
      method = integration_method,
      batch = batch,
      assay = assay,
      layer = layer,
      features = hvg,
      dims = dims,
      npcs = npcs,
      theta = theta,
      group_by_vars = group_by_vars,
      integration_metadata = integration_metadata_signature,
      integration_control = backend_integration_control,
      store_sce = if (identical(integration_method, "coralysis")) {
        .sn_coralysis_store_sce(integration_control)
      } else {
        NULL
      },
      pca = pca_signature
    )
    if (is_null(x = batch)) {
      reduction <- "pca"
    } else if (identical(integration_method, "unintegrated")) {
      reduction <- "pca"
      integration_performance <- list(
        elapsed_seconds = 0,
        peak_memory_mb = NA_real_,
        peak_r_memory_mb = NA_real_,
        backend_peak_rss_mb = NA_real_,
        memory_scope = "not_applicable",
        reused = TRUE
      )
      object@misc$integration <- list(
        method = "unintegrated",
        batch_by = batch,
        reduction = reduction,
        input_reduction = reduction,
        performance = integration_performance
      )
      object <- .sn_record_cluster_stage(
        object,
        "integration",
        integration_signature,
        reduction = reduction,
        performance = integration_performance
      )
    } else if (.sn_can_reuse_cluster_stage(
      object = object,
      stage = "integration",
      signature = integration_signature,
      reuse = reuse,
      rerun_from = rerun_from,
      required = function(current_object, stage_info) {
        has_reduction <- !is.null(stage_info$reduction) && stage_info$reduction %in% names(current_object@reductions)
        if (!identical(integration_method, "bbknn")) {
          return(has_reduction)
        }
        has_reduction &&
          !is.null(stage_info$graph_name) && stage_info$graph_name %in% names(current_object@graphs) &&
          !is.null(stage_info$umap_reduction) && stage_info$umap_reduction %in% names(current_object@reductions)
      }
    )) {
      reduction <- object@misc$sn_run_cluster$stages$integration$reduction
      integration_graph <- object@misc$sn_run_cluster$stages$integration$graph_name %||% NULL
      integration_umap <- object@misc$sn_run_cluster$stages$integration$umap_reduction %||% NULL
      if (verbose) .sn_log_info("[5/6] Reusing {reduction} integration reduction.")
    } else {
      if (verbose) .sn_log_info("[5/6] Running {integration_method} integration.")
      profiled_integration <- .sn_cluster_performance(.sn_run_batch_integration(
        object = object,
        method = integration_method,
        batch = batch,
        reduction = if (isTRUE(needs_rna_pca)) "pca" else NULL,
        features = hvg,
        assay = assay,
        layer = layer,
        dims = dims,
        npcs = npcs,
        theta = theta,
        group_by_vars = group_by_vars,
        integration_control = backend_integration_control,
        verbose = verbose
      ))
      integration <- profiled_integration$value
      object <- integration$object
      integration_performance <- .sn_complete_integration_performance(
        profiled_integration$performance,
        object
      )
      object@misc$integration$performance <- integration_performance
      reduction <- integration$reduction
      integration_graph <- integration$graph %||% NULL
      integration_umap <- integration$umap %||% NULL
      object <- .sn_record_cluster_stage(
        object,
        "integration",
        integration_signature,
        reduction = reduction,
        graph_name = integration_graph,
        umap_reduction = integration_umap,
        performance = integration_performance
      )
    }
    reduction_signature <- if (is_null(x = batch)) pca_signature else integration_signature
  }

  use_bbknn_graph <- identical(modality, "rna") && identical(integration_method, "bbknn") && !is.null(batch)
  if (verbose) .sn_log_info("[6/6] Clustering with integrated embeddings.")
  if (identical(modality, "cite_seq") && identical(multimodal_method, "wnn")) {
    dims <- .sn_valid_reduction_dims(object = object, reduction = "pca", dims = dims)
    adt_dims <- .sn_valid_reduction_dims(object = object, reduction = "apca", dims = adt_dims)
    neighbors_signature <- list(
      method = "weighted_nearest_neighbor",
      reduction_list = c("pca", "apca"),
      dims_list = list(rna = dims, adt = adt_dims),
      modality_weight_name = c("RNA.weight", "ADT.weight"),
      knn_range = max(1L, min(200L, ncol(object) - 1L)),
      wnn_control = wnn_control,
      upstream = reduction_signature
    )
  } else if (use_bbknn_graph) {
    if (is.null(integration_graph) || !integration_graph %in% names(object@graphs)) {
      stop("BBKNN integration did not provide a connectivity graph.", call. = FALSE)
    }
    dims <- .sn_valid_reduction_dims(object = object, reduction = reduction, dims = dims)
    neighbors_signature <- list(
      method = "bbknn",
      reduction = reduction,
      dims = dims,
      graph_name = integration_graph,
      upstream = reduction_signature
    )
  } else {
    reduction_dims <- if (identical(modality, "cite_seq") && multimodal_method %in% c("coralysis", "mmochi")) {
      adt_dims
    } else {
      dims
    }
    dims <- .sn_valid_reduction_dims(object = object, reduction = reduction, dims = reduction_dims)
    neighbors_signature <- list(
      reduction = reduction,
      dims = dims,
      upstream = reduction_signature
    )
  }
  if (.sn_can_reuse_cluster_stage(
    object = object,
    stage = "neighbors",
    signature = neighbors_signature,
    reuse = reuse,
    rerun_from = rerun_from,
    required = function(current_object, stage_info) {
      length(stage_info$graph_names %||% character(0)) > 0L &&
        all(stage_info$graph_names %in% names(current_object@graphs))
    }
  )) {
    if (verbose) .sn_log_info("[6/6] Reusing nearest-neighbor graph.")
  } else {
    if (use_bbknn_graph) {
      graph_names <- integration_graph
      snn_graph <- integration_graph
    } else {
      graph_names_before <- names(object@graphs)
      if (identical(modality, "cite_seq") && identical(multimodal_method, "wnn")) {
        wnn_args <- .sn_merge_control_args(
          defaults = list(
            object = object,
            reduction.list = list("pca", "apca"),
            dims.list = list(dims, adt_dims),
            modality.weight.name = c("RNA.weight", "ADT.weight"),
            knn.range = max(1L, min(200L, ncol(object) - 1L)),
            verbose = verbose
          ),
          control = wnn_control
        )
        object <- .sn_call_with_symbolic_object(
          fun_call = quote(Seurat::FindMultiModalNeighbors),
          object = object,
          args = wnn_args
        )
      } else {
        if (is.null(result_namespace)) {
          object <- .sn_with_default_seurat_acceleration(
            Seurat::FindNeighbors(
              object,
              reduction = reduction,
              dims = dims,
              verbose = verbose
            ),
            object = object,
            assay = assay
          )
        } else {
          object <- .sn_with_default_seurat_acceleration(
            Seurat::FindNeighbors(
              object,
              reduction = reduction,
              dims = dims,
              graph.name = c(
                paste0(result_namespace, "_nn"),
                paste0(result_namespace, "_snn")
              ),
              verbose = verbose
            ),
            object = object,
            assay = assay
          )
        }
      }
      graph_names_after <- names(object@graphs)
      if (identical(modality, "cite_seq") && identical(multimodal_method, "wnn") && "wsnn" %in% graph_names_after) {
        graph_names <- setdiff(graph_names_after, graph_names_before)
        if (length(graph_names) == 0L) {
          graph_names <- intersect(c("wknn", "wsnn"), graph_names_after)
        }
        snn_graph <- "wsnn"
      } else {
        reduction_assay <- tryCatch(
          SeuratObject::DefaultAssay(object[[reduction]]),
          error = function(error) assay
        )
        expected_graph_names <- if (is.null(result_namespace)) {
          paste0(reduction_assay, c("_nn", "_snn"))
        } else {
          paste0(result_namespace, c("_nn", "_snn"))
        }
        graph_names <- intersect(expected_graph_names, graph_names_after)
        expected_snn <- expected_graph_names[[2L]]
        if (expected_snn %in% graph_names_after) {
          snn_graph <- expected_snn
        } else {
          added_graphs <- setdiff(graph_names_after, graph_names_before)
          snn_candidates <- grep("_snn$", added_graphs, value = TRUE)
          if (length(snn_candidates) != 1L) {
            stop(
              "Could not identify the SNN graph created by `Seurat::FindNeighbors()`.",
              call. = FALSE
            )
          }
          snn_graph <- snn_candidates[[1L]]
          graph_names <- unique(c(graph_names, added_graphs))
        }
      }
    }
    object <- .sn_record_cluster_stage(
      object,
      "neighbors",
      neighbors_signature,
      graph_names = graph_names,
      snn_graph = snn_graph,
      nn_name = if (identical(modality, "cite_seq") && identical(multimodal_method, "wnn")) "weighted.nn" else NULL
    )
  }
  find_clusters_args <- .sn_merge_control_args(
    defaults = list(
      object = object,
      resolution = resolution,
      algorithm = cluster_algorithm_value,
      n.start = cluster_n_start,
      n.iter = cluster_n_iter,
      random.seed = cluster_random_seed,
      group.singletons = cluster_group_singletons,
      leiden_method = leiden_method,
      leiden_objective_function = leiden_objective_function,
      verbose = verbose
    ),
    control = cluster_control
  )
  if (!is.null(cluster_name)) {
    find_clusters_args$cluster.name <- cluster_name
  }
  cluster_column <- cluster_name %||% "seurat_clusters"
  neighbors_stage <- object@misc$sn_run_cluster$stages$neighbors %||% NULL
  if (is.null(find_clusters_args$graph.name) && !is.null(neighbors_stage$snn_graph)) {
    find_clusters_args$graph.name <- neighbors_stage$snn_graph
  }
  cluster_signature <- find_clusters_args
  cluster_signature$object <- NULL
  cluster_signature$verbose <- NULL
  cluster_signature$neighbors <- neighbors_signature
  if (.sn_can_reuse_cluster_stage(
    object = object,
    stage = "clusters",
    signature = cluster_signature,
    reuse = reuse,
    rerun_from = rerun_from,
    required = function(current_object, stage_info) cluster_column %in% colnames(current_object[[]])
  )) {
    if (verbose) .sn_log_info("[7/7] Reusing cluster assignments.")
  } else {
    .sn_ensure_cluster_algorithm_dependencies(
      cluster_algorithm_value = find_clusters_args$algorithm,
      leiden_method = find_clusters_args$leiden_method %||% leiden_method,
      auto_install = auto_install,
      repos = install_repos,
      ask = install_ask
    )
    object <- .sn_with_default_seurat_acceleration(
      .sn_call_with_symbolic_object(
        fun_call = quote(Seurat::FindClusters),
        object = object,
        args = find_clusters_args
      ),
      object = object,
      assay = assay
    )
    object <- .sn_record_cluster_stage(object, "clusters", cluster_signature, cluster_column = cluster_column)
  }

  if (return_cluster) {
    object <- restore_analysis_inputs(object)
    if (verbose) .sn_log_info("Integration completed successfully.")
    if (return_object_for_multi) {
      return(.sn_log_seurat_command(object = object, assay = assay, name = "sn_run_cluster"))
    }
    return(object@meta.data[, cluster_column])
  } else {
    if (identical(modality, "cite_seq") && identical(multimodal_method, "wnn")) {
      umap_args <- .sn_merge_control_args(
        defaults = list(
          nn.name = "weighted.nn",
          reduction.name = "wnn.umap",
          reduction.key = "wnnUMAP_",
          umap.method = "uwot",
          metric = "cosine",
          verbose = verbose,
          seed.use = 717
        ),
        control = umap_control
      )
      umap_signature <- umap_args
      umap_signature$verbose <- NULL
      umap_signature$method <- "weighted_nearest_neighbor"
      umap_signature$upstream <- neighbors_signature
    } else if (use_bbknn_graph) {
      if (is.null(integration_umap) || !integration_umap %in% names(object@reductions)) {
        stop("BBKNN integration did not provide a graph-derived UMAP reduction.", call. = FALSE)
      }
      umap_args <- list(reduction.name = integration_umap)
      umap_signature <- list(
        method = "bbknn_graph",
        graph = integration_graph,
        reduction = integration_umap,
        umap_control = umap_control
      )
      umap_signature$upstream <- neighbors_signature
    } else {
      umap_args <- .sn_merge_control_args(
        defaults = list(
          reduction = reduction,
          dims = dims,
          umap.method = "uwot",
          metric = "cosine",
          verbose = verbose,
          seed.use = 717
        ),
        control = umap_control
      )
      umap_signature <- umap_args
      umap_signature$verbose <- NULL
      umap_signature$upstream <- reduction_signature
    }
    if (.sn_can_reuse_cluster_stage(
      object = object,
      stage = "umap",
      signature = umap_signature,
      reuse = reuse,
      rerun_from = rerun_from,
      required = function(current_object, stage_info) {
        !is.null(stage_info$reduction) && stage_info$reduction %in% names(current_object@reductions)
      }
    )) {
      if (verbose) .sn_log_info("[7/7] Reusing UMAP reduction.")
    } else {
      if (verbose) .sn_log_info("[7/7] Running UMAP.")
      umap_reduction <- umap_args$reduction.name %||% "umap"
      if (use_bbknn_graph) {
        umap_reduction <- integration_umap
        object <- .sn_record_cluster_stage(object, "umap", umap_signature, reduction = umap_reduction)
      } else if (identical(modality, "cite_seq") && identical(multimodal_method, "wnn")) {
        object <- suppressWarnings(.sn_call_with_symbolic_object(
          fun_call = quote(Seurat::RunUMAP),
          object = object,
          args = umap_args
        ))
        object <- .sn_record_cluster_stage(object, "umap", umap_signature, reduction = umap_reduction)
      } else {
        object <- suppressWarnings(.sn_call_with_symbolic_object(
          fun_call = quote(Seurat::RunUMAP),
          object = object,
          args = umap_args
        ))
        object <- .sn_record_cluster_stage(object, "umap", umap_signature, reduction = umap_reduction)
      }
    }
    if (isTRUE(run_tsne)) {
      if (identical(modality, "cite_seq") && identical(multimodal_method, "wnn")) {
        stop("`run_tsne = TRUE` is not yet supported for WNN graph clustering.", call. = FALSE)
      }
      tsne_args <- .sn_merge_control_args(
        defaults = list(
          reduction = reduction,
          dims = dims,
          reduction.name = "tsne",
          reduction.key = "tSNE_",
          perplexity = max(1, min(30, floor((ncol(object) - 1L) / 3L))),
          check_duplicates = FALSE,
          seed.use = 717
        ),
        control = tsne_control
      )
      tsne_signature <- tsne_args
      tsne_signature$upstream <- reduction_signature
      if (.sn_can_reuse_cluster_stage(
        object = object,
        stage = "tsne",
        signature = tsne_signature,
        reuse = reuse,
        rerun_from = rerun_from,
        required = function(current_object, stage_info) {
          !is.null(stage_info$reduction) && stage_info$reduction %in% names(current_object@reductions)
        }
      )) {
        if (verbose) .sn_log_info("[8/8] Reusing t-SNE reduction.")
      } else {
        if (verbose) .sn_log_info("[8/8] Running t-SNE.")
        object <- suppressWarnings(.sn_call_with_symbolic_object(
          fun_call = quote(Seurat::RunTSNE),
          object = object,
          args = tsne_args
        ))
        object <- .sn_record_cluster_stage(
          object,
          "tsne",
          tsne_signature,
          reduction = tsne_args$reduction.name %||% "tsne"
        )
      }
    }
    object <- restore_analysis_inputs(object)
    if (verbose) .sn_log_info("Integration completed successfully.")
    return(.sn_log_seurat_command(object = object, assay = assay, name = "sn_run_cluster"))
  }
}

#' Run a Python analysis command through a managed Shennong pixi environment
#'
#' These are analysis-oriented wrappers around \code{sn_call_pixi_environment()}.
#' They prepare the corresponding package-bundled environment and run the
#' requested command. When \code{object} is supplied, method wrappers use a
#' Seurat object-level contract: export the object, run the packaged pixi
#' runner script, and import method outputs back into the object when the
#' backend produces cell-level metadata or embeddings.
#'
#' @param object Seurat object. Shennong writes the object to a Python
#'   interchange directory, runs the corresponding pixi script, and imports
#'   supported results.
#' @param reference_object Optional reference Seurat object for tools that map
#'   a query/spatial object against a single-cell reference, such as Tangram.
#' @param reference_assay,reference_layer Assay and layer used when exporting
#'   \code{reference_object}.
#' @param reference_signatures Optional file path or data frame of reference
#'   cell-state signatures for cell2location.
#' @param group_by Metadata column used by CellPhoneDB cell groups.
#' @param batch_by,label_by Metadata columns used by scArches/scPoli-style
#'   object workflows.
#' @param spatial_cols Two metadata columns containing spatial coordinates for
#'   spatial tools.
#' @param cell_type_by Reference metadata column containing cell-type labels
#'   for Tangram projection.
#' @param cluster_by Metadata column used by Squidpy neighborhood enrichment.
#' @param method_control Optional named list of backend-specific settings passed
#'   to the Python runner config.
#' @param assay Assay used for object-level infercnvpy input.
#' @param layer Assay layer used for object-level Python input. scPoli defaults
#'   to \code{"counts"}; infercnvpy and the other generic object wrappers
#'   default to \code{"data"} when present and otherwise \code{"counts"}.
#' @param species Species used to match bundled gene positions when
#'   \code{gene_order} and \code{gtf_file} are not supplied.
#' @param reference_by Metadata column containing normal/tumor annotations.
#' @param reference_cat One or more values in \code{reference_by} denoting
#'   normal reference cells.
#' @param gene_order Optional data frame with gene positions. It must contain a
#'   gene identifier column such as \code{feature}, \code{gene},
#'   \code{gene_name}, or \code{gene_id}, plus chromosome/start/end columns.
#' @param gtf_file Optional GTF file used by infercnvpy to annotate genomic
#'   positions instead of Shennong's bundled GENCODE table.
#' @param gtf_gene_id GTF attribute used by infercnvpy for matching.
#' @param adata_gene_id Optional AnnData var column used for matching a GTF.
#' @param output_dir Optional run directory. Defaults to
#'   \code{~/.shennong/runs/infercnvpy_*}.
#' @param runtime_dir Optional Shennong runtime directory.
#' @param key_added infercnvpy key used for the CNV representation.
#' @param window_size,step,dynamic_threshold,exclude_chromosomes,chunksize,n_jobs,calculate_gene_values,lfc_clip
#'   Parameters forwarded to \code{infercnvpy.tl.infercnv()}.
#' @param run_pca,run_neighbors,run_leiden,run_umap,score Logical flags for
#'   downstream infercnvpy analysis steps.
#' @param leiden_resolution Resolution passed to infercnvpy Leiden clustering.
#' @param cnv_score_group_by Optional grouping column for infercnvpy CNV scores.
#' @param metadata_prefix Prefix added to imported infercnvpy metadata columns.
#' @param result_name Name used under \code{object@misc$infercnvpy}.
#' @param return_object Whether to return the updated object. If \code{FALSE},
#'   return a run manifest list.
#' @param ... Additional arguments passed to \code{sn_call_pixi_environment()}.
#'
#' @return A Seurat object or a run manifest.
#'
#' @examples
#' \dontrun{
#' object <- sn_run_infercnvpy(
#'   object = object,
#'   reference_by = "cell_type",
#'   reference_cat = c("T cell", "Myeloid")
#' )
#' spatial <- sn_run_tangram(
#'   object = spatial,
#'   reference_object = reference,
#'   cell_type_by = "cell_type",
#'   spatial_cols = c("x", "y")
#' )
#' }
#'
#' @export
sn_run_scarches <- function(object,
                            assay = NULL,
                            layer = NULL,
                            batch_by = NULL,
                            label_by = NULL,
                            output_dir = NULL,
                            runtime_dir = NULL,
                            metadata_prefix = "scarches_",
                            result_name = "scarches",
                            return_object = TRUE,
                            method_control = list(),
                            ...) {
  .sn_run_python_object_method(
    object = object,
    environment = "scarches",
    script_name = "scarches_run.py",
    method = "scarches",
    assay = assay,
    layer = layer,
    output_dir = output_dir,
    runtime_dir = runtime_dir,
    metadata_prefix = metadata_prefix,
    result_name = result_name,
    return_object = return_object,
    config = c(list(batch_key = batch_by, labels_key = label_by), method_control),
    ...
  )
}

#' @rdname sn_run_scarches
#' @export
sn_run_scpoli <- function(object,
                          assay = NULL,
                          layer = NULL,
                          batch_by = NULL,
                          label_by = NULL,
                          output_dir = NULL,
                          runtime_dir = NULL,
                          metadata_prefix = "scpoli_",
                          result_name = "scpoli",
                          return_object = TRUE,
                          method_control = list(),
                          ...) {
  .sn_run_python_object_method(
    object = object,
    environment = "scarches",
    script_name = "scarches_run.py",
    method = "scpoli",
    assay = assay,
    layer = layer %||% "counts",
    output_dir = output_dir,
    runtime_dir = runtime_dir,
    metadata_prefix = metadata_prefix,
    result_name = result_name,
    return_object = return_object,
    config = c(list(batch_key = batch_by, labels_key = label_by), method_control),
    ...
  )
}

#' Deprecated alias of `sn_get_integration_control_template()`
#'
#' `sn_integration_control_template()` is a deprecated compatibility alias. Use [sn_get_integration_control_template()] directly;
#' the alias will be removed in a future release.
#'
#' @param ... Named arguments passed on to [sn_get_integration_control_template()].
#'
#' @return Result of `sn_get_integration_control_template(...)`.
#'
#' @export
sn_integration_control_template <- function(...) {
  .Deprecated("sn_get_integration_control_template", package = "Shennong")
  sn_get_integration_control_template(...)
}
