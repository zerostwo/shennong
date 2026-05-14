# Internal helper to normalize the clustering dimension arguments.
.sn_resolve_cluster_dims <- function(dims = NULL, npcs = 50) {
  dims <- dims %||% seq_len(min(20, npcs))
  dims <- as.integer(dims)
  dims[dims > 0]
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
  is.null(batch) || integration_method %in% c("harmony", "seurat_cca", "seurat_rpca")
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
  umap = 9L
)

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

.sn_merge_control_args <- function(defaults, control) {
  control <- control %||% list()
  if (!is.list(control)) {
    stop("`integration_control` must be a named list.", call. = FALSE)
  }
  utils::modifyList(defaults, control, keep.null = TRUE)
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

.sn_run_batch_integration <- function(object,
                                      method,
                                      batch,
                                      reduction,
                                      features,
                                      assay,
                                      dims,
                                      npcs,
                                      theta,
                                      group_by_vars,
                                      integration_control = list(),
                                      verbose = TRUE) {
  switch(
    method,
    harmony = .sn_run_harmony_integration(
      object = object,
      batch = batch,
      reduction = reduction,
      theta = theta,
      group_by_vars = group_by_vars,
      verbose = verbose
    ),
    coralysis = .sn_run_coralysis_integration(
      object = object,
      batch = batch,
      assay = assay,
      features = features,
      dims = dims,
      npcs = npcs,
      integration_control = integration_control,
      verbose = verbose
    ),
    seurat_cca = .sn_run_seurat_layer_integration(
      object = object,
      method = method,
      batch = batch,
      reduction = reduction,
      features = features,
      assay = assay,
      dims = dims,
      integration_control = integration_control,
      verbose = verbose
    ),
    seurat_rpca = .sn_run_seurat_layer_integration(
      object = object,
      method = method,
      batch = batch,
      reduction = reduction,
      features = features,
      assay = assay,
      dims = dims,
      integration_control = integration_control,
      verbose = verbose
    ),
    scvi = .sn_run_scvi_integration(
      object = object,
      method = method,
      batch = batch,
      features = features,
      assay = assay,
      integration_control = integration_control,
      verbose = verbose
    ),
    totalvi = .sn_run_scvi_integration(
      object = object,
      method = method,
      batch = batch,
      features = features,
      assay = assay,
      integration_control = integration_control,
      verbose = verbose
    ),
    scanvi = .sn_run_scvi_integration(
      object = object,
      method = method,
      batch = batch,
      features = features,
      assay = assay,
      integration_control = integration_control,
      verbose = verbose
    )
  )
}

.sn_run_harmony_integration <- function(object,
                                        batch,
                                        reduction,
                                        theta = 2,
                                        group_by_vars = NULL,
                                        verbose = TRUE) {
  check_installed("harmony")
  group_by_vars <- group_by_vars %||% batch
  object <- harmony::RunHarmony(
    object = object,
    group.by.vars = group_by_vars,
    theta = theta,
    reduction.use = reduction,
    verbose = verbose
  )
  object@misc$integration <- list(
    method = "harmony",
    batch_by = batch,
    group_by_vars = group_by_vars,
    reduction = "harmony",
    input_reduction = reduction,
    theta = theta
  )
  list(object = object, reduction = "harmony")
}

.sn_resolve_integration_assay_features <- function(object,
                                                   assay,
                                                   features,
                                                   backend = "Integration") {
  feature_set <- intersect(features, rownames(object[[assay]]))
  if (length(feature_set) < 2L) {
    stop(glue("{backend} integration requires at least two selected features present in the object."), call. = FALSE)
  }
  feature_set
}

.sn_mmochi_script_path <- function(script = NULL) {
  if (!is.null(script) && nzchar(script)) {
    script <- path.expand(script)
    if (!file.exists(script)) {
      stop("`integration_control$script` does not exist: ", script, call. = FALSE)
    }
    return(normalizePath(script, winslash = "/", mustWork = TRUE))
  }

  installed <- system.file("pixi", "mmochi", "scripts", "mmochi_run.py", package = "Shennong")
  if (nzchar(installed) && file.exists(installed)) {
    return(normalizePath(installed, winslash = "/", mustWork = TRUE))
  }

  source_path <- file.path(getwd(), "inst", "pixi", "mmochi", "scripts", "mmochi_run.py")
  if (file.exists(source_path)) {
    return(normalizePath(source_path, winslash = "/", mustWork = TRUE))
  }

  stop("Could not locate Shennong's MMoCHi Python runner.", call. = FALSE)
}

.sn_unique_mmochi_batch_key <- function(object, base = ".sn_mmochi_single_sample") {
  if (!is.character(base) || length(base) != 1L || !nzchar(base)) {
    stop("`integration_control$single_sample_batch_key` must be a non-empty string.", call. = FALSE)
  }
  metadata_columns <- colnames(object[[]])
  key <- base
  counter <- 1L
  while (key %in% metadata_columns) {
    counter <- counter + 1L
    key <- paste0(base, "_", counter)
  }
  key
}

.sn_write_mmochi_input <- function(object,
                                   input_dir,
                                   batch,
                                   protein_assay,
                                   protein_layer = "data",
                                   protein_features = NULL) {
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  protein <- .sn_get_seurat_layer_data(object = object, assay = protein_assay, layer = protein_layer)
  feature_set <- protein_features %||% rownames(protein)
  feature_set <- intersect(feature_set, rownames(protein))
  if (length(feature_set) < 2L) {
    stop("MMoCHi integration requires at least two ADT/protein features present in the object.", call. = FALSE)
  }
  protein <- protein[feature_set, colnames(object), drop = FALSE]
  protein_df <- as.data.frame(t(as.matrix(protein)), check.names = FALSE)
  protein_df <- cbind(cell_id = rownames(protein_df), protein_df)

  metadata <- object[[]]
  metadata <- metadata[colnames(object), , drop = FALSE]
  metadata <- data.frame(cell_id = rownames(metadata), metadata, check.names = FALSE)

  protein_path <- file.path(input_dir, "protein.csv")
  obs_path <- file.path(input_dir, "obs.csv")
  utils::write.csv(protein_df, protein_path, row.names = FALSE, quote = TRUE)
  utils::write.csv(metadata, obs_path, row.names = FALSE, quote = TRUE)

  list(
    input_dir = normalizePath(input_dir, winslash = "/", mustWork = TRUE),
    protein_path = normalizePath(protein_path, winslash = "/", mustWork = TRUE),
    obs_path = normalizePath(obs_path, winslash = "/", mustWork = TRUE),
    protein_features = feature_set,
    cells = colnames(object),
    batch = batch
  )
}

.sn_import_mmochi_results <- function(object,
                                      output_dir,
                                      assay,
                                      batch,
                                      backend_batch = batch,
                                      single_sample_batch = FALSE,
                                      protein_assay,
                                      protein_layer,
                                      protein_features,
                                      dims,
                                      npcs,
                                      reduction = "mmochi",
                                      corrected_layer = "mmochi.data",
                                      store_corrected_layer = TRUE) {
  corrected_path <- file.path(output_dir, "landmark_protein.csv")
  manifest_path <- file.path(output_dir, "manifest.json")
  if (!file.exists(corrected_path)) {
    stop("MMoCHi output is missing `landmark_protein.csv`: ", corrected_path, call. = FALSE)
  }

  corrected <- utils::read.csv(corrected_path, row.names = 1, check.names = FALSE)
  missing_cells <- setdiff(colnames(object), rownames(corrected))
  if (length(missing_cells) > 0L) {
    stop("MMoCHi landmark output is missing cells from the input object.", call. = FALSE)
  }
  missing_features <- setdiff(protein_features, colnames(corrected))
  if (length(missing_features) > 0L) {
    stop("MMoCHi landmark output is missing requested ADT/protein features.", call. = FALSE)
  }
  corrected <- as.matrix(corrected[colnames(object), protein_features, drop = FALSE])
  storage.mode(corrected) <- "numeric"
  corrected_features_by_cells <- t(corrected)

  stored_corrected_layer <- NULL
  stored_corrected_misc <- NULL
  if (isTRUE(store_corrected_layer)) {
    layer_error <- NULL
    object <- tryCatch(
      {
        SeuratObject::LayerData(object = object, assay = protein_assay, layer = corrected_layer) <- corrected_features_by_cells
        stored_corrected_layer <- corrected_layer
        object
      },
      error = function(e) {
        layer_error <<- conditionMessage(e)
        object
      }
    )
    if (!is.null(layer_error)) {
      object@misc$mmochi <- object@misc$mmochi %||% list()
      object@misc$mmochi$corrected_protein <- corrected_features_by_cells
      object@misc$mmochi$corrected_protein_layer_error <- layer_error
      stored_corrected_misc <- "object@misc$mmochi$corrected_protein"
    }
  }

  scaled <- t(scale(t(corrected_features_by_cells)))
  scaled[is.na(scaled)] <- 0
  available_pcs <- max(1L, min(as.integer(npcs), as.integer(max(dims)), nrow(scaled), ncol(scaled) - 1L))
  pca <- stats::prcomp(t(scaled), rank. = available_pcs, center = FALSE, scale. = FALSE)
  embeddings <- pca$x[, seq_len(available_pcs), drop = FALSE]
  embeddings <- embeddings[colnames(object), , drop = FALSE]
  reduction_key <- paste0(toupper(reduction), "_")
  colnames(embeddings) <- paste0(reduction_key, seq_len(ncol(embeddings)))
  object[[reduction]] <- Seurat::CreateDimReducObject(
    embeddings = embeddings,
    key = reduction_key,
    assay = assay
  )

  backend_manifest <- if (file.exists(manifest_path)) {
    jsonlite::read_json(manifest_path, simplifyVector = TRUE)
  } else {
    list()
  }
  object@misc$integration <- list(
    method = "mmochi",
    batch_by = batch,
    backend_batch_key = backend_batch,
    single_sample_batch = isTRUE(single_sample_batch),
    reduction = reduction,
    input_features = protein_features,
    protein_assay = protein_assay,
    protein_layer = protein_layer,
    corrected_layer = stored_corrected_layer,
    corrected_storage = if (!is.null(stored_corrected_layer)) "layer" else if (!is.null(stored_corrected_misc)) "misc" else NULL,
    corrected_misc = stored_corrected_misc,
    run_dir = normalizePath(output_dir, winslash = "/", mustWork = TRUE),
    mmochi_version = backend_manifest$mmochi_version %||% NULL
  )

  list(object = object, reduction = reduction)
}

.sn_run_mmochi_integration <- function(object,
                                       batch,
                                       assay,
                                       protein_assay,
                                       protein_features,
                                       dims,
                                       npcs,
                                       integration_control = list(),
                                       verbose = TRUE) {
  user_batch <- batch
  backend_batch <- batch
  single_sample_batch <- is.null(backend_batch) || !nzchar(backend_batch)
  if (single_sample_batch) {
    backend_batch <- .sn_unique_mmochi_batch_key(
      object = object,
      base = integration_control$single_sample_batch_key %||% ".sn_mmochi_single_sample"
    )
    object[[backend_batch]] <- rep("single_sample", ncol(object))
    if (verbose) {
      .sn_log_info("MMoCHi single-sample mode uses internal backend batch key '{backend_batch}'.")
    }
  } else if (!backend_batch %in% colnames(object[[]])) {
    stop("MMoCHi landmark registration requires `batch` to name a metadata column.", call. = FALSE)
  }

  runtime_dir <- .sn_shennong_runtime_dir(integration_control$runtime_dir %||% NULL)
  pixi_paths <- sn_pixi_paths(environment = "mmochi", runtime_dir = runtime_dir)
  pixi_home <- integration_control$pixi_home %||% pixi_paths$pixi_home
  mirror <- match.arg(integration_control$mirror %||% "default", c("default", "auto", "china", "tuna", "ustc", "bfsu"))
  resolved_mirror <- .sn_resolve_pixi_mirror(mirror)
  if (!identical(resolved_mirror, "default")) {
    sn_configure_pixi_mirror(
      mirror = mirror,
      pixi_home = pixi_home,
      runtime_dir = runtime_dir,
      append_original = integration_control$mirror_append_original %||% TRUE
    )
  }

  pixi_project <- integration_control$pixi_project %||%
    integration_control$pixi_project_dir %||%
    pixi_paths$project_dir
  run_dir <- integration_control$run_dir %||% .sn_default_scvi_run_dir(method = "mmochi", runtime_dir = runtime_dir)
  input_dir <- file.path(run_dir, "input")
  output_dir <- file.path(run_dir, "output")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  manifest_path <- .sn_prepare_scvi_pixi_project(
    project_dir = pixi_project,
    environment = "mmochi",
    manifest_path = integration_control$manifest_path %||% NULL,
    manifest_lines = integration_control$manifest_lines %||% NULL,
    overwrite = isTRUE(integration_control$overwrite_manifest),
    platforms = integration_control$platforms %||% NULL
  )
  protein_layer <- integration_control$protein_layer %||% integration_control$adt_layer %||% "data"
  input <- .sn_write_mmochi_input(
    object = object,
    input_dir = input_dir,
    batch = backend_batch,
    protein_assay = protein_assay,
    protein_layer = protein_layer,
    protein_features = protein_features
  )

  config <- list(
    method = "mmochi",
    batch_key = backend_batch,
    data_key = integration_control$data_key %||% "protein",
    key_added = integration_control$key_added %||% "landmark_protein",
    single_peaks = integration_control$single_peaks %||% list(),
    marker_bandwidths = integration_control$marker_bandwidths %||% list(),
    peak_overrides = integration_control$peak_overrides %||% list(),
    inclusion_mask = integration_control$inclusion_mask %||% NULL,
    landmark_args = integration_control$landmark_args %||% list(),
    show = integration_control$show %||% FALSE,
    protein_assay = protein_assay,
    protein_layer = protein_layer,
    protein_features = input$protein_features
  )
  config_path <- .sn_write_json_file(config, file.path(run_dir, "config.json"))
  script <- .sn_mmochi_script_path(integration_control$script %||% NULL)

  .sn_execute_scvi_pixi(
    pixi = integration_control$pixi %||% NULL,
    manifest_path = manifest_path,
    script = script,
    input_dir = input$input_dir,
    output_dir = output_dir,
    config_path = config_path,
    environment = integration_control$environment %||% "default",
    pixi_home = pixi_home,
    install_pixi = integration_control$install_pixi %||% TRUE,
    pixi_version = integration_control$pixi_version %||% "latest",
    pixi_download_url = integration_control$pixi_download_url %||% NULL,
    verbose = verbose,
    backend_label = "MMoCHi"
  )

  result <- .sn_import_mmochi_results(
    object = object,
    output_dir = output_dir,
    assay = assay,
    batch = user_batch,
    backend_batch = backend_batch,
    single_sample_batch = single_sample_batch,
    protein_assay = protein_assay,
    protein_layer = protein_layer,
    protein_features = input$protein_features,
    dims = dims,
    npcs = npcs,
    reduction = integration_control$reduction %||% "mmochi",
    corrected_layer = integration_control$corrected_layer %||% "mmochi.data",
    store_corrected_layer = integration_control$store_corrected_layer %||% TRUE
  )

  if (single_sample_batch && !isTRUE(integration_control$keep_single_sample_batch)) {
    result$object@meta.data[[backend_batch]] <- NULL
  }
  result
}

.sn_run_coralysis_integration <- function(object,
                                          batch,
                                          assay,
                                          features,
                                          dims,
                                          npcs,
                                          integration_control = list(),
                                          verbose = TRUE) {
  backend <- "Coralysis"
  method <- "coralysis"
  check_installed(c("Coralysis", "SingleCellExperiment", "SummarizedExperiment"))

  feature_set <- .sn_resolve_integration_assay_features(
    object = object,
    assay = assay,
    features = features,
    backend = backend
  )
  coralysis_ns <- asNamespace(backend)
  prepare_data <- get("PrepareData", envir = coralysis_ns, inherits = FALSE)
  run_icp <- get("RunParallelDivisiveICP", envir = coralysis_ns, inherits = FALSE)
  run_pca <- get("RunPCA", envir = coralysis_ns, inherits = FALSE)
  reduction_name <- "coralysis"
  reduction_key <- "CORALYSIS_"
  default_dimred_name <- "Coralysis"

  expr <- SeuratObject::LayerData(object = object, assay = assay, layer = "data")
  expr <- expr[feature_set, colnames(object), drop = FALSE]
  metadata <- object[[]]

  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = expr),
    colData = metadata
  )

  if (verbose) .sn_log_info("[sn_run_cluster] Preparing {backend} input.")
  sce <- prepare_data(object = sce)

  icp_args <- .sn_merge_control_args(
    defaults = {
      defaults <- list(
        object = sce,
        batch.label = batch,
        threads = 0,
        verbose = verbose,
        RNGseed = 717
      )
      batch_sizes <- if (!is.null(batch) && batch %in% colnames(metadata)) {
        table(metadata[[batch]])
      } else {
        ncol(object)
      }
      defaults$build.train.params <- list(
        nhvg = min(length(feature_set), 2000L),
        p = min(30L, max(1L, length(feature_set) - 1L), max(1L, min(batch_sizes) - 5L))
      )
      defaults
    },
    control = integration_control$icp_args
  )
  if (verbose) .sn_log_info("[sn_run_cluster] Running {backend} multi-level integration.")
  sce <- do.call(run_icp, icp_args)

  coralysis_dims <- max(dims)
  coralysis_p <- max(1L, min(as.integer(npcs), as.integer(coralysis_dims), ncol(object) - 1L))
  pca_args <- .sn_merge_control_args(
    defaults = list(
      object = sce,
      assay.name = "joint.probability",
      p = coralysis_p,
      dimred.name = default_dimred_name,
      return.model = TRUE
    ),
    control = integration_control$pca_args
  )
  if (verbose) .sn_log_info("[sn_run_cluster] Running PCA on {backend} joint probabilities.")
  sce <- do.call(run_pca, pca_args)

  dimred_name <- pca_args$dimred.name %||% default_dimred_name
  if (!dimred_name %in% SingleCellExperiment::reducedDimNames(sce)) {
    dimred_name <- utils::tail(SingleCellExperiment::reducedDimNames(sce), n = 1L)
  }
  embeddings <- SingleCellExperiment::reducedDim(sce, dimred_name)
  embeddings <- embeddings[colnames(object), , drop = FALSE]
  colnames(embeddings) <- paste0(reduction_key, seq_len(ncol(embeddings)))

  object[[reduction_name]] <- Seurat::CreateDimReducObject(
    embeddings = embeddings,
    key = reduction_key,
    assay = assay
  )
  store_sce <- .sn_coralysis_store_sce(integration_control)
  integration_metadata <- list(
    method = method,
    batch_by = batch,
    reduction = reduction_name,
    input_features = feature_set,
    coralysis_dimred = dimred_name,
    stored_sce = isTRUE(store_sce)
  )
  object@misc$integration <- integration_metadata
  if (isTRUE(store_sce)) {
    object@misc[[reduction_name]] <- sce
  } else {
    object@misc[[reduction_name]] <- NULL
  }

  list(object = object, reduction = reduction_name)
}

.sn_run_seurat_layer_integration <- function(object,
                                             method,
                                             batch,
                                             reduction,
                                             features,
                                             assay,
                                             dims,
                                             integration_control = list(),
                                             verbose = TRUE) {
  if (!exists("IntegrateLayers", envir = asNamespace("Seurat"), inherits = FALSE)) {
    stop("Seurat layer integration requires Seurat >= 5 with `IntegrateLayers()`.", call. = FALSE)
  }
  integration_fun <- switch(
    method,
    seurat_cca = Seurat::CCAIntegration,
    seurat_rpca = Seurat::RPCAIntegration
  )
  new_reduction <- switch(
    method,
    seurat_cca = "integrated.cca",
    seurat_rpca = "integrated.rpca"
  )

  old_default_assay <- SeuratObject::DefaultAssay(object = object)
  on.exit(SeuratObject::DefaultAssay(object = object) <- old_default_assay, add = TRUE)
  SeuratObject::DefaultAssay(object = object) <- assay
  object[[assay]] <- split(object[[assay]], f = object[[batch, drop = TRUE]])

  args <- .sn_merge_control_args(
    defaults = list(
      object = object,
      method = integration_fun,
      orig.reduction = reduction,
      assay = assay,
      features = features,
      dims = dims,
      new.reduction = new_reduction,
      verbose = verbose
    ),
    control = integration_control
  )
  if (verbose) .sn_log_info("[sn_run_cluster] Running Seurat layer integration with method = {method}.")
  object <- .sn_call_with_symbolic_object(
    fun_call = quote(Seurat::IntegrateLayers),
    object = object,
    args = args
  )
  object[[assay]] <- SeuratObject::JoinLayers(object[[assay]])
  SeuratObject::DefaultAssay(object = object) <- old_default_assay
  object@misc$integration <- list(
    method = method,
    batch_by = batch,
    reduction = new_reduction,
    input_reduction = reduction
  )

  list(object = object, reduction = new_reduction)
}

.sn_shennong_runtime_dir <- function(path = NULL) {
  candidates <- c(
    path,
    getOption("shennong.runtime_dir", NULL),
    Sys.getenv("SHENNONG_RUNTIME_DIR", unset = ""),
    Sys.getenv("SHENNONG_HOME", unset = ""),
    "~/.shennong"
  )
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  path <- candidates[[1]]
  path <- path.expand(path)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.sn_default_scvi_run_dir <- function(method, runtime_dir = NULL) {
  runtime_dir <- .sn_shennong_runtime_dir(runtime_dir)
  run_id <- paste0(
    method,
    "_",
    format(Sys.time(), "%Y%m%d_%H%M%S"),
    "_",
    Sys.getpid()
  )
  file.path(runtime_dir, "runs", run_id)
}

.sn_current_pixi_platform <- function() {
  sysname <- tolower(Sys.info()[["sysname"]] %||% .Platform$OS.type)
  machine <- tolower(Sys.info()[["machine"]] %||% "")
  if (grepl("darwin", sysname, fixed = TRUE)) {
    if (machine %in% c("arm64", "aarch64")) "osx-arm64" else "osx-64"
  } else if (.Platform$OS.type == "windows" || grepl("windows", sysname, fixed = TRUE)) {
    "win-64"
  } else if (machine %in% c("aarch64", "arm64")) {
    "linux-aarch64"
  } else {
    "linux-64"
  }
}

.sn_scvi_pixi_manifest_lines <- function(cuda_version = "12.6",
                                         platforms = NULL) {
  .sn_render_pixi_config(environment = "scvi", platforms = platforms, cuda_version = cuda_version)
}

.sn_prepare_scvi_pixi_project <- function(project_dir,
                                          environment = "scvi",
                                          manifest_path = NULL,
                                          manifest_lines = NULL,
                                          overwrite = FALSE,
                                          cuda_version = "12.6",
                                          platforms = NULL) {
  project_dir <- path.expand(project_dir)
  dir.create(project_dir, recursive = TRUE, showWarnings = FALSE)
  manifest_path <- manifest_path %||% file.path(project_dir, "pixi.toml")
  manifest_path <- path.expand(manifest_path)
  dir.create(dirname(manifest_path), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(manifest_path) || isTRUE(overwrite)) {
    writeLines(
      manifest_lines %||% .sn_render_pixi_config(
        environment = environment,
        cuda_version = cuda_version,
        platforms = platforms
      ),
      con = manifest_path,
      useBytes = TRUE
    )
  }
  normalizePath(manifest_path, winslash = "/", mustWork = TRUE)
}

.sn_scvi_script_path <- function(script = NULL) {
  if (!is.null(script) && nzchar(script)) {
    script <- path.expand(script)
    if (!file.exists(script)) {
      stop("`integration_control$script` does not exist: ", script, call. = FALSE)
    }
    return(normalizePath(script, winslash = "/", mustWork = TRUE))
  }

  installed <- system.file("pixi", "scvi", "scripts", "scvi_integration.py", package = "Shennong")
  if (nzchar(installed) && file.exists(installed)) {
    return(normalizePath(installed, winslash = "/", mustWork = TRUE))
  }

  source_path <- file.path(getwd(), "inst", "pixi", "scvi", "scripts", "scvi_integration.py")
  if (file.exists(source_path)) {
    return(normalizePath(source_path, winslash = "/", mustWork = TRUE))
  }

  stop("Could not locate Shennong's scVI Python runner.", call. = FALSE)
}

.sn_write_scvi_input <- function(object,
                                 input_dir,
                                 features,
                                 assay,
                                 batch,
                                 labels_key = NULL,
                                 protein_assay = NULL,
                                 protein_layer = "counts",
                                 protein_features = NULL) {
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

  feature_set <- intersect(features, rownames(object))
  if (length(feature_set) < 2L) {
    stop("scVI integration requires at least two selected features present in the object.", call. = FALSE)
  }

  counts <- .sn_get_seurat_layer_data(object = object, assay = assay, layer = "counts")
  counts <- counts[feature_set, colnames(object), drop = FALSE]
  counts <- .sn_as_sparse_matrix(counts)

  counts_path <- file.path(input_dir, "counts.mtx")
  features_path <- file.path(input_dir, "features.csv")
  cells_path <- file.path(input_dir, "cells.csv")
  obs_path <- file.path(input_dir, "obs.csv")

  Matrix::writeMM(obj = counts, file = counts_path)
  utils::write.csv(
    data.frame(feature_id = rownames(counts), stringsAsFactors = FALSE),
    file = features_path,
    row.names = FALSE,
    quote = TRUE
  )
  utils::write.csv(
    data.frame(cell_id = colnames(counts), stringsAsFactors = FALSE),
    file = cells_path,
    row.names = FALSE,
    quote = TRUE
  )

  obs_cols <- unique(c(batch, labels_key))
  obs_cols <- obs_cols[!is.na(obs_cols) & nzchar(obs_cols)]
  obs <- object[[]][colnames(counts), obs_cols, drop = FALSE]
  obs <- data.frame(cell_id = rownames(obs), obs, check.names = FALSE)
  utils::write.csv(obs, file = obs_path, row.names = FALSE, quote = TRUE)

  protein_counts_path <- NULL
  proteins_path <- NULL
  protein_feature_set <- character(0)
  if (!is.null(protein_assay) && nzchar(protein_assay)) {
    .sn_validate_seurat_assay_layer(object = object, assay = protein_assay, layer = protein_layer)
    protein_counts <- .sn_get_seurat_layer_data(object = object, assay = protein_assay, layer = protein_layer)
    protein_feature_set <- protein_features %||% rownames(protein_counts)
    protein_feature_set <- unique(protein_feature_set[!is.na(protein_feature_set) & nzchar(protein_feature_set)])
    protein_feature_set <- intersect(protein_feature_set, rownames(protein_counts))
    if (length(protein_feature_set) < 1L) {
      stop("totalVI integration requires at least one protein feature present in `adt_assay`.", call. = FALSE)
    }
    protein_counts <- protein_counts[protein_feature_set, colnames(counts), drop = FALSE]
    protein_counts <- .sn_as_sparse_matrix(protein_counts)
    protein_counts_path <- file.path(input_dir, "protein_counts.mtx")
    proteins_path <- file.path(input_dir, "proteins.csv")
    Matrix::writeMM(obj = protein_counts, file = protein_counts_path)
    utils::write.csv(
      data.frame(protein_id = rownames(protein_counts), stringsAsFactors = FALSE),
      file = proteins_path,
      row.names = FALSE,
      quote = TRUE
    )
  }

  list(
    input_dir = normalizePath(input_dir, winslash = "/", mustWork = TRUE),
    counts_path = normalizePath(counts_path, winslash = "/", mustWork = TRUE),
    features_path = normalizePath(features_path, winslash = "/", mustWork = TRUE),
    cells_path = normalizePath(cells_path, winslash = "/", mustWork = TRUE),
    obs_path = normalizePath(obs_path, winslash = "/", mustWork = TRUE),
    protein_counts_path = if (!is.null(protein_counts_path)) normalizePath(protein_counts_path, winslash = "/", mustWork = TRUE) else NULL,
    proteins_path = if (!is.null(proteins_path)) normalizePath(proteins_path, winslash = "/", mustWork = TRUE) else NULL,
    features = feature_set,
    protein_features = protein_feature_set
  )
}

.sn_write_json_file <- function(x, path) {
  jsonlite::write_json(x = x, path = path, auto_unbox = TRUE, pretty = TRUE, null = "null")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

.sn_execute_scvi_pixi <- function(pixi,
                                  manifest_path,
                                  script,
                                  input_dir,
                                  output_dir,
                                  config_path,
                                  environment = NULL,
                                  pixi_home = NULL,
                                  install_pixi = TRUE,
                                  pixi_version = "latest",
                                  pixi_download_url = NULL,
                                  verbose = TRUE,
                                  backend_label = "scVI") {
  pixi_info <- sn_ensure_pixi(
    pixi = pixi,
    install = install_pixi,
    version = pixi_version,
    pixi_home = path.expand("~/.pixi"),
    no_path_update = TRUE,
    download_url = pixi_download_url,
    quiet = !isTRUE(verbose)
  )
  pixi <- pixi_info$path

  args <- c(
    "run",
    "--manifest-path", shQuote(manifest_path),
    if (!is.null(environment) && nzchar(environment)) c("--environment", shQuote(environment)) else character(0),
    "python",
    shQuote(script),
    "--input-dir", shQuote(input_dir),
    "--output-dir", shQuote(output_dir),
    "--config", shQuote(config_path)
  )
  if (verbose) {
    .sn_log_info("[sn_run_cluster] Executing {backend_label} backend with pixi manifest: {manifest_path}.")
  }
  env <- if (!is.null(pixi_home) && nzchar(pixi_home)) {
    dir.create(pixi_home, recursive = TRUE, showWarnings = FALSE)
    paste0("PIXI_HOME=", normalizePath(pixi_home, winslash = "/", mustWork = TRUE))
  } else {
    character(0)
  }

  status <- tryCatch(
    system2(command = pixi, args = args, env = env, stdout = TRUE, stderr = TRUE),
    error = function(e) {
      stop(backend_label, " pixi execution failed. ", conditionMessage(e), call. = FALSE)
    }
  )
  exit_code <- attr(status, "status") %||% 0L
  if (!identical(exit_code, 0L)) {
    stop(backend_label, " pixi execution failed.\n", paste(status, collapse = "\n"), call. = FALSE)
  }
  invisible(status)
}

.sn_normalize_cuda_requirement <- function(cuda_version) {
  cuda_version <- as.character(cuda_version %||% "12.0")
  major_minor <- regmatches(cuda_version, regexpr("^[0-9]+(\\.[0-9]+)?", cuda_version))
  if (length(major_minor) == 0L || !nzchar(major_minor)) {
    return("12.0")
  }
  parts <- strsplit(major_minor, ".", fixed = TRUE)[[1]]
  if (length(parts) == 1L) {
    paste0(parts[[1]], ".0")
  } else {
    paste(parts[1:2], collapse = ".")
  }
}

.sn_default_scvi_cuda_version <- function(detected_cuda_version) {
  detected_cuda_version <- .sn_normalize_cuda_requirement(detected_cuda_version %||% "12.6")
  detected_major <- suppressWarnings(as.integer(strsplit(detected_cuda_version, ".", fixed = TRUE)[[1]][[1]]))
  if (!is.na(detected_major) && detected_major >= 12L) {
    return("12.6")
  }
  detected_cuda_version
}

.sn_resolve_scvi_accelerator <- function(accelerator = c("auto", "cpu", "gpu", "cuda")) {
  accelerator <- match.arg(accelerator)
  if (identical(accelerator, "gpu")) {
    accelerator <- "cuda"
  }
  detected <- sn_detect_accelerator(quiet = TRUE)
  environment <- if (identical(accelerator, "auto")) {
    if (identical(detected$backend, "cuda")) "gpu" else "cpu"
  } else if (identical(accelerator, "cuda")) {
    "gpu"
  } else {
    "cpu"
  }
  cuda_version <- .sn_default_scvi_cuda_version(detected$cuda_version)
  list(
    requested = accelerator,
    environment = environment,
    detected = detected,
    cuda_version = cuda_version
  )
}

.sn_import_scvi_results <- function(object,
                                    output_dir,
                                    method,
                                    assay,
                                    batch,
                                    features,
                                    reduction = NULL) {
  reduction <- reduction %||% method
  latent_path <- file.path(output_dir, "latent.csv")
  metadata_path <- file.path(output_dir, "obs.csv")
  manifest_path <- file.path(output_dir, "manifest.json")

  if (!file.exists(latent_path)) {
    stop("scVI output is missing `latent.csv`: ", latent_path, call. = FALSE)
  }
  latent <- utils::read.csv(latent_path, row.names = 1, check.names = FALSE)
  missing_cells <- setdiff(colnames(object), rownames(latent))
  if (length(missing_cells) > 0L) {
    stop("scVI latent output is missing cells from the input object.", call. = FALSE)
  }
  latent <- as.matrix(latent[colnames(object), , drop = FALSE])
  storage.mode(latent) <- "numeric"
  colnames(latent) <- paste0(toupper(reduction), "_", seq_len(ncol(latent)))
  object[[reduction]] <- Seurat::CreateDimReducObject(
    embeddings = latent,
    key = paste0(toupper(reduction), "_"),
    assay = assay
  )

  if (file.exists(metadata_path)) {
    metadata <- utils::read.csv(metadata_path, row.names = 1, check.names = FALSE)
    shared_cells <- intersect(colnames(object), rownames(metadata))
    if (length(shared_cells) > 0L && ncol(metadata) > 0L) {
      object <- Seurat::AddMetaData(object = object, metadata = metadata[shared_cells, , drop = FALSE])
    }
  }

  backend_manifest <- if (file.exists(manifest_path)) {
    jsonlite::read_json(manifest_path, simplifyVector = TRUE)
  } else {
    list()
  }
  object@misc$integration <- list(
    method = method,
    batch_by = batch,
    reduction = reduction,
    input_features = features,
    run_dir = normalizePath(output_dir, winslash = "/", mustWork = TRUE),
    output_h5ad = backend_manifest$output_h5ad %||% file.path(output_dir, "integrated.h5ad")
  )

  list(object = object, reduction = reduction)
}

.sn_run_scvi_integration <- function(object,
                                     method,
                                     batch,
                                     features,
                                     assay,
                                     integration_control = list(),
                                     verbose = TRUE) {
  protein_assay <- NULL
  protein_layer <- NULL
  protein_features <- NULL
  if (identical(method, "totalvi")) {
    protein_assay <- integration_control$protein_assay %||% integration_control$adt_assay %||% NULL
    if (is.null(protein_assay) || !nzchar(protein_assay)) {
      stop("totalVI integration requires `integration_control$adt_assay` or `protein_assay`.", call. = FALSE)
    }
    protein_layer <- integration_control$protein_layer %||% integration_control$adt_layer %||% "counts"
    protein_features <- integration_control$protein_features %||% integration_control$adt_features %||% NULL
  }

  if (identical(method, "scanvi")) {
    labels_key <- integration_control$label_by %||% NULL
    if (is.null(labels_key) || !nzchar(labels_key) || !labels_key %in% colnames(object[[]])) {
      stop(
        "`integration_control$label_by` must name a metadata column when ",
        "`integration_method = \"scanvi\"`.",
        call. = FALSE
      )
    }
  } else {
    labels_key <- integration_control$label_by %||% NULL
  }

  runtime_dir <- .sn_shennong_runtime_dir(integration_control$runtime_dir %||% NULL)
  pixi_paths <- sn_pixi_paths(environment = method, runtime_dir = runtime_dir)
  pixi_home <- integration_control$pixi_home %||% pixi_paths$pixi_home
  accelerator <- .sn_resolve_scvi_accelerator(integration_control$accelerator %||% "auto")
  pixi_environment <- integration_control$environment %||% accelerator$environment
  mirror <- match.arg(integration_control$mirror %||% "default", c("default", "auto", "china", "tuna", "ustc", "bfsu"))
  resolved_mirror <- .sn_resolve_pixi_mirror(mirror)
  if (!identical(resolved_mirror, "default")) {
    sn_configure_pixi_mirror(
      mirror = mirror,
      pixi_home = pixi_home,
      runtime_dir = runtime_dir,
      append_original = integration_control$mirror_append_original %||% TRUE
    )
  }
  pixi_project <- integration_control$pixi_project %||%
    integration_control$pixi_project_dir %||%
    pixi_paths$project_dir
  run_dir <- integration_control$run_dir %||% .sn_default_scvi_run_dir(method = method, runtime_dir = runtime_dir)
  input_dir <- file.path(run_dir, "input")
  output_dir <- file.path(run_dir, "output")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  manifest_path <- .sn_prepare_scvi_pixi_project(
    project_dir = pixi_project,
    environment = method,
    manifest_path = integration_control$manifest_path %||% NULL,
    manifest_lines = integration_control$manifest_lines %||% NULL,
    overwrite = isTRUE(integration_control$overwrite_manifest),
    cuda_version = integration_control$cuda_version %||% accelerator$cuda_version,
    platforms = integration_control$platforms %||% NULL
  )
  input <- .sn_write_scvi_input(
    object = object,
    input_dir = input_dir,
    features = features,
    assay = assay,
    batch = batch,
    labels_key = labels_key,
    protein_assay = protein_assay,
    protein_layer = protein_layer %||% "counts",
    protein_features = protein_features
  )

  config <- list(
    method = method,
    batch_key = batch,
    labels_key = labels_key,
    unlabeled_category = integration_control$unlabeled_category %||% "Unknown",
    reduction = integration_control$reduction %||% method,
    n_latent = integration_control$n_latent %||% 30L,
    seed = integration_control$seed %||% 717L,
    max_epochs = integration_control$max_epochs %||% NULL,
    scanvi_max_epochs = integration_control$scanvi_max_epochs %||% integration_control$max_epochs %||% NULL,
    model_args = integration_control$model_args %||% list(),
    train_args = integration_control$train_args %||% list(),
    scanvi_model_args = integration_control$scanvi_model_args %||% list(),
    scanvi_train_args = integration_control$scanvi_train_args %||% list(),
    totalvi_model_args = integration_control$totalvi_model_args %||% integration_control$model_args %||% list(),
    totalvi_train_args = integration_control$totalvi_train_args %||% integration_control$train_args %||% list(),
    protein_obsm_key = integration_control$protein_obsm_key %||% "protein_expression",
    protein_assay = protein_assay,
    protein_layer = protein_layer,
    protein_features = input$protein_features,
    write_h5ad = integration_control$write_h5ad %||% TRUE,
    accelerator = accelerator$requested,
    pixi_environment = pixi_environment
  )
  config_path <- .sn_write_json_file(config, file.path(run_dir, "config.json"))
  script <- .sn_scvi_script_path(integration_control$script %||% NULL)

  .sn_execute_scvi_pixi(
    pixi = integration_control$pixi %||% NULL,
    manifest_path = manifest_path,
    script = script,
    input_dir = input$input_dir,
    output_dir = output_dir,
    config_path = config_path,
    environment = pixi_environment,
    pixi_home = pixi_home,
    install_pixi = integration_control$install_pixi %||% TRUE,
    pixi_version = integration_control$pixi_version %||% "latest",
    pixi_download_url = integration_control$pixi_download_url %||% NULL,
    verbose = verbose
  )

  .sn_import_scvi_results(
    object = object,
    output_dir = output_dir,
    method = method,
    assay = assay,
    batch = batch,
    features = input$features,
    reduction = config$reduction
  )
}

#' Run scVI or scANVI integration through Shennong
#'
#' Convenience wrappers around \code{\link{sn_run_cluster}} for users who want
#' an explicit Python-method entry point. The underlying pixi environment is
#' prepared from the bundled scVI-family config under \code{inst/pixi/scvi/}
#' and materialized under \code{~/.shennong/pixi/scvi/}.
#'
#' @inheritParams sn_run_cluster
#' @param ... Additional arguments passed to \code{sn_run_cluster()}.
#'
#' @return A Seurat object returned by \code{sn_run_cluster()}.
#'
#' @examples
#' \dontrun{
#' obj <- sn_run_scvi(obj, batch = "sample_id")
#' obj <- sn_run_scanvi(
#'   obj,
#'   batch = "sample_id",
#'   integration_control = list(label_by = "cell_type")
#' )
#' }
#'
#' @export
sn_run_scvi <- function(object,
                        batch = NULL,
                        integration_control = list(),
                        ...) {
  sn_run_cluster(
    object = object,
    batch = batch,
    integration_method = "scvi",
    integration_control = integration_control,
    ...
  )
}

#' @rdname sn_run_scvi
#' @export
sn_run_scanvi <- function(object,
                          batch = NULL,
                          integration_control = list(),
                          ...) {
  sn_run_cluster(
    object = object,
    batch = batch,
    integration_method = "scanvi",
    integration_control = integration_control,
    ...
  )
}

.sn_gini_coefficient <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x) & x >= 0]
  if (length(x) == 0 || sum(x) <= 0) {
    return(0)
  }
  x <- sort(x, decreasing = FALSE)
  n <- length(x)
  (2 * sum(seq_len(n) * x) / (n * sum(x))) - ((n + 1) / n)
}

.sn_detect_rare_features_gini <- function(expr,
                                          nfeatures = 200,
                                          min_cells = 3,
                                          max_fraction = 0.1) {
  stopifnot(nrow(expr) > 0, ncol(expr) > 1)

  prevalence <- Matrix::rowSums(expr > 0)
  prevalence_fraction <- prevalence / ncol(expr)
  keep <- prevalence >= min_cells & prevalence_fraction <= max_fraction

  if (!any(keep)) {
    return(data.frame(
      feature = character(0),
      score = numeric(0),
      prevalence = integer(0),
      prevalence_fraction = numeric(0),
      mean_expression = numeric(0),
      stringsAsFactors = FALSE
    ))
  }

  kept_features <- rownames(expr)[keep]
  gini_score <- vapply(kept_features, function(feature) {
    .sn_gini_coefficient(expr[feature, ])
  }, numeric(1))
  mean_expression <- Matrix::rowMeans(expr[kept_features, , drop = FALSE])

  out <- data.frame(
    feature = kept_features,
    score = gini_score,
    prevalence = as.integer(prevalence[keep]),
    prevalence_fraction = as.numeric(prevalence_fraction[keep]),
    mean_expression = as.numeric(mean_expression),
    stringsAsFactors = FALSE
  )
  out <- out[order(-out$score, out$prevalence_fraction, -out$mean_expression, out$feature), , drop = FALSE]
  utils::head(out, nfeatures)
}

.sn_make_temporary_grouping <- function(object,
                                        features,
                                        npcs = 20,
                                        dims = 1:10,
                                        resolution = 0.2) {
  temp_object <- suppressWarnings(
    Seurat::ScaleData(
      object = object,
      features = features,
      verbose = FALSE
    )
  )
  temp_object <- suppressWarnings(
    Seurat::RunPCA(
      object = temp_object,
      features = features,
      npcs = npcs,
      verbose = FALSE,
      seed.use = 717
    )
  )
  temp_dims <- dims[dims <= npcs]
  temp_object <- Seurat::FindNeighbors(
    object = temp_object,
    reduction = "pca",
    dims = temp_dims,
    verbose = FALSE
  )
  temp_object <- Seurat::FindClusters(
    object = temp_object,
    resolution = resolution,
    random.seed = 717,
    verbose = FALSE
  )
  as.character(temp_object$seurat_clusters)
}

.sn_resolve_rare_groups <- function(object,
                                    group_by = NULL,
                                    features = NULL,
                                    npcs = 20,
                                    dims = 1:10,
                                    resolution = 0.2,
                                    max_fraction = 0.05,
                                    max_cells = 100) {
  groups <- if (!is.null(group_by)) {
    if (!group_by %in% colnames(object[[]])) {
      stop(glue("`rare_feature_group_by` column '{group_by}' was not found."), call. = FALSE)
    }
    as.character(object[[group_by, drop = TRUE]])
  } else {
    .sn_make_temporary_grouping(
      object = object,
      features = features,
      npcs = npcs,
      dims = dims,
      resolution = resolution
    )
  }

  group_sizes <- table(groups)
  rare_groups <- names(group_sizes)[
    group_sizes <= max_cells | (group_sizes / length(groups)) <= max_fraction
  ]

  list(
    groups = groups,
    rare_groups = rare_groups,
    group_sizes = group_sizes
  )
}

.sn_detect_rare_features_local_markers <- function(object,
                                                   groups,
                                                   rare_groups,
                                                   assay = "RNA",
                                                   layer = "data",
                                                   nfeatures = 200,
                                                   min_cells = 3) {
  if (length(rare_groups) == 0) {
    return(character(0))
  }

  expr <- .sn_get_seurat_layer_data(object = object, assay = assay, layer = layer)
  local_markers <- lapply(rare_groups, function(current_group) {
    in_group <- groups == current_group
    if (sum(in_group) < min_cells || sum(!in_group) < min_cells) {
      return(character(0))
    }
    avg_in <- Matrix::rowMeans(expr[, in_group, drop = FALSE])
    avg_out <- Matrix::rowMeans(expr[, !in_group, drop = FALSE])
    pct_in <- Matrix::rowSums(expr[, in_group, drop = FALSE] > 0) / sum(in_group)
    pct_out <- Matrix::rowSums(expr[, !in_group, drop = FALSE] > 0) / sum(!in_group)
    score <- avg_in - avg_out
    keep <- is.finite(score) & pct_in > pct_out & pct_in > 0
    ranked <- names(sort(score[keep], decreasing = TRUE))
    utils::head(ranked, nfeatures)
  })

  utils::head(unique(unlist(local_markers, use.names = FALSE)), nfeatures)
}

.sn_resolve_rare_feature_control <- function(control = list()) {
  if (is.null(control)) {
    control <- list()
  }
  if (!is.list(control)) {
    stop("`rare_feature_control` must be a named list.", call. = FALSE)
  }

  defaults <- list(
    group_max_fraction = 0.05,
    group_max_cells = 100,
    gene_max_fraction = 0.1,
    min_cells = 3
  )
  unknown <- setdiff(names(control), names(defaults))
  if (length(unknown) > 0L) {
    stop(
      glue("Unknown `rare_feature_control` field(s): {paste(unknown, collapse = ', ')}."),
      call. = FALSE
    )
  }
  resolved <- utils::modifyList(defaults, control)

  resolved$group_max_fraction <- as.numeric(resolved$group_max_fraction)
  resolved$group_max_cells <- as.integer(resolved$group_max_cells)
  resolved$gene_max_fraction <- as.numeric(resolved$gene_max_fraction)
  resolved$min_cells <- as.integer(resolved$min_cells)

  if (!is.finite(resolved$group_max_fraction) || resolved$group_max_fraction <= 0 || resolved$group_max_fraction > 1) {
    stop("`rare_feature_control$group_max_fraction` must be in (0, 1].", call. = FALSE)
  }
  if (!is.finite(resolved$gene_max_fraction) || resolved$gene_max_fraction <= 0 || resolved$gene_max_fraction > 1) {
    stop("`rare_feature_control$gene_max_fraction` must be in (0, 1].", call. = FALSE)
  }
  if (is.na(resolved$group_max_cells) || resolved$group_max_cells < 1L) {
    stop("`rare_feature_control$group_max_cells` must be a positive integer.", call. = FALSE)
  }
  if (is.na(resolved$min_cells) || resolved$min_cells < 1L) {
    stop("`rare_feature_control$min_cells` must be a positive integer.", call. = FALSE)
  }

  resolved
}

.sn_select_rare_features <- function(object,
                                     base_features,
                                     method = "none",
                                     assay = "RNA",
                                     layer = "data",
                                     nfeatures = 200,
                                     group_by = NULL,
                                     control = list(),
                                     min_cells = 3,
                                     npcs = 20,
                                     dims = 1:10,
                                     resolution = 0.2,
                                     verbose = TRUE) {
  control <- .sn_resolve_rare_feature_control(control = control)
  method <- unique(match.arg(method, c("none", "gini", "local_markers"), several.ok = TRUE))
  method <- setdiff(method, "none")
  if (length(method) == 0) {
    return(list(
      features = character(0),
      metadata = data.frame(),
      groups = NULL
    ))
  }

  expr <- .sn_get_seurat_layer_data(object = object, assay = assay, layer = layer)
  rare_features <- list()
  rare_meta <- list()
  rare_group_info <- NULL

  if ("local_markers" %in% method) {
    rare_group_info <- .sn_resolve_rare_groups(
      object = object,
      group_by = group_by,
      features = base_features,
      npcs = npcs,
      dims = dims,
      resolution = resolution,
      max_fraction = control$group_max_fraction,
      max_cells = control$group_max_cells
    )
  }

  if ("gini" %in% method) {
    gini_tbl <- .sn_detect_rare_features_gini(
      expr = expr,
      nfeatures = nfeatures,
      min_cells = control$min_cells %||% min_cells,
      max_fraction = control$gene_max_fraction
    )
    rare_features$gini <- gini_tbl$feature
    if (nrow(gini_tbl) > 0) {
      rare_meta[[length(rare_meta) + 1]] <- transform(gini_tbl, method = "gini")
    }
  }

  if ("local_markers" %in% method) {
    local_markers <- .sn_detect_rare_features_local_markers(
      object = object,
      groups = rare_group_info$groups,
      rare_groups = rare_group_info$rare_groups,
      assay = assay,
      layer = layer,
      nfeatures = nfeatures,
      min_cells = control$min_cells %||% min_cells
    )
    rare_features$local_markers <- local_markers
    if (length(local_markers) > 0) {
      rare_meta[[length(rare_meta) + 1]] <- data.frame(
        feature = local_markers,
        score = NA_real_,
        prevalence = NA_integer_,
        prevalence_fraction = NA_real_,
        mean_expression = NA_real_,
        method = "local_markers",
        stringsAsFactors = FALSE
      )
    }
  }

  combined <- unique(unlist(rare_features, use.names = FALSE))
  if (verbose) {
    .sn_log_info(
      "[rare_features] Added {length(combined)} rare-aware feature(s) from method(s): {paste(method, collapse = ', ')}"
    )
  }

  list(
    features = combined,
    metadata = .sn_bind_rows(rare_meta),
    groups = rare_group_info
  )
}

.sn_run_sccad <- function(expr,
                          cell_ids,
                          gene_ids,
                          python = NULL,
                          script = NULL,
                          normalization = FALSE,
                          seed = 2023,
                          rare_h = 0.01,
                          merge_h = 0.3,
                          overlap_h = 0.7,
                          save_full = FALSE) {
  python <- python %||% getOption("shennong.sccad_python", Sys.which("python"))
  if (!nzchar(python)) {
    stop("Could not find a Python executable for the scCAD backend.", call. = FALSE)
  }

  script <- script %||% getOption("shennong.sccad_script", Sys.getenv("SHENNONG_SCCAD_SCRIPT", unset = ""))
  if (!nzchar(script) || !file.exists(path.expand(script))) {
    stop(
      "Could not locate `scCAD.py`. Supply `sccad_script` or set `options(shennong.sccad_script = ...)`.",
      call. = FALSE
    )
  }
  script <- normalizePath(path.expand(script), winslash = "/", mustWork = TRUE)

  workdir <- tempfile("sccad_")
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(workdir, recursive = TRUE, force = TRUE), add = TRUE)

  expr_path <- file.path(workdir, "expression.csv")
  output_path <- file.path(workdir, "sccad_result.json")
  runner_path <- file.path(workdir, "run_sccad.py")

  expr_df <- as.data.frame(t(as.matrix(expr)))
  colnames(expr_df) <- gene_ids
  rownames(expr_df) <- cell_ids
  utils::write.csv(expr_df, file = expr_path, quote = FALSE)

  runner_lines <- c(
    "import json",
    "import os",
    "import sys",
    "import numpy as np",
    "import pandas as pd",
    sprintf("sys.path.insert(0, %s)", shQuote(dirname(script))),
    sprintf("import %s as sccad_module", tools::file_path_sans_ext(basename(script))),
    sprintf("expr = pd.read_csv(%s, index_col=0)", shQuote(expr_path)),
    "data = expr.to_numpy(dtype=float)",
    "cell_names = expr.index.to_numpy()",
    "gene_names = expr.columns.to_numpy()",
    sprintf(
      paste0(
        "result, score, sub_clusters, degs_list = sccad_module.scCAD(",
        "data=data, dataName='Shennong', cellNames=cell_names, geneNames=gene_names, ",
        "normalization=%s, seed=%s, rare_h=%s, merge_h=%s, overlap_h=%s, save_full=%s, save_path=%s)"
      ),
      if (isTRUE(normalization)) "True" else "False",
      as.integer(seed),
      as.numeric(rare_h),
      as.numeric(merge_h),
      as.numeric(overlap_h),
      if (isTRUE(save_full)) "True" else "False",
      shQuote(workdir)
    ),
    "def _to_str_list(x):",
    "    out = []",
    "    for item in x:",
    "        if isinstance(item, bytes):",
    "            out.append(item.decode('utf-8'))",
    "        else:",
    "            out.append(str(item))",
    "    return out",
    "rare_sets = [_to_str_list(cluster) for cluster_by in result]",
    "sub_clusters = [str(x) for x in sub_clusters]",
    "score = [float(x) for x in score]",
    "payload = {",
    "    'rare_sets': rare_sets,",
    "    'scores': score,",
    "    'sub_clusters': sub_clusters,",
    "    'degs_list': [[str(g) for g in genes] for genes in degs_list]",
    "}",
    sprintf("with open(%s, 'w') as handle:", shQuote(output_path)),
    "    json.dump(payload, handle)"
  )
  writeLines(runner_lines, con = runner_path, useBytes = TRUE)

  status <- tryCatch(
    system2(
      command = python,
      args = c(shQuote(runner_path)),
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(e) {
      stop(
        "SCA execution failed. Ensure the Python executable and the `shannonca` package are available. ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  exit_code <- attr(status, "status") %||% 0L
  if (!identical(exit_code, 0L) || !file.exists(output_path)) {
    stop(
      "scCAD execution failed. ",
      paste(status, collapse = "\n"),
      call. = FALSE
    )
  }

  jsonlite::read_json(output_path, simplifyVector = TRUE)
}

.sn_get_embedding_knn <- function(embeddings, k = 20, n_trees = 50) {
  k <- min(as.integer(k), nrow(embeddings) - 1L)
  if (k < 1L) {
    stop("At least two cells are required to score embedding rarity.", call. = FALSE)
  }

  if (rlang::is_installed("Seurat") && nrow(embeddings) > 1000L) {
    return(.sn_find_annoy_knn(
      embeddings = embeddings,
      k = k,
      n_trees = n_trees,
      include_distance = TRUE
    ))
  }

  .sn_exact_knn(
    embeddings = embeddings,
    k = k,
    include_distance = TRUE
  )
}

.sn_score_embedding_rarity <- function(embeddings, k = 20, n_trees = 50) {
  knn <- .sn_get_embedding_knn(embeddings = embeddings, k = k, n_trees = n_trees)
  rare_score <- rowMeans(knn$dist, na.rm = TRUE)
  rare_score[!is.finite(rare_score)] <- 0
  stats::setNames(as.numeric(rare_score), rownames(embeddings))
}

.sn_run_gapclust <- function(expr, k = 200) {
  check_installed_github("GapClust", "fabotao/GapClust", reason = "to detect rare cells with the GapClust backend.")

  result <- GapClust::GapClust(data = as.matrix(expr), k = as.integer(k))
  cell_ids <- colnames(expr)
  if (length(result) == 1 && is.na(result)) {
    return(data.frame(
      cell_id = cell_ids,
      rare_score = 0,
      rare_cell = FALSE,
      stringsAsFactors = FALSE
    ))
  }

  rare_membership <- sort(unique(unlist(result$rare_cell_indices, use.names = FALSE)))
  rare_score <- apply(result$rare_score, 1, function(x) {
    current <- x[is.finite(x)]
    if (length(current) == 0) {
      return(0)
    }
    max(current)
  })

  data.frame(
    cell_id = cell_ids,
    rare_score = as.numeric(rare_score),
    rare_cell = seq_along(cell_ids) %in% rare_membership,
    stringsAsFactors = FALSE
  )
}

.sn_run_sca <- function(expr,
                        python = NULL,
                        n_comps = 20,
                        iters = 3,
                        nbhd_size = 15,
                        model = "wilcoxon",
                        seed = 717) {
  python <- python %||% getOption("shennong.sca_python", Sys.which("python"))
  if (!nzchar(python)) {
    stop("Could not find a Python executable for the SCA backend.", call. = FALSE)
  }

  workdir <- tempfile("sca_")
  dir.create(workdir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(workdir, recursive = TRUE, force = TRUE), add = TRUE)

  expr_path <- file.path(workdir, "expression.csv")
  output_path <- file.path(workdir, "sca_result.json")
  runner_path <- file.path(workdir, "run_sca.py")

  expr_df <- as.data.frame(t(as.matrix(expr)))
  colnames(expr_df) <- rownames(expr)
  rownames(expr_df) <- colnames(expr)
  utils::write.csv(expr_df, file = expr_path, quote = FALSE)

  runner_lines <- c(
    "import json",
    "import pandas as pd",
    "from shannonca.dimred import reduce",
    sprintf("expr = pd.read_csv(%s, index_col=0)", shQuote(expr_path)),
    "reduction = reduce(",
    "    expr.to_numpy(dtype=float),",
    sprintf("    n_comps=%s,", min(as.integer(n_comps), ncol(expr_df) - 1L)),
    sprintf("    iters=%s,", as.integer(iters)),
    sprintf("    nbhd_size=%s,", as.integer(nbhd_size)),
    sprintf("    model=%s,", shQuote(model)),
    sprintf("    seed=%s", as.integer(seed)),
    ")",
    "payload = {'reduction': reduction.tolist()}",
    sprintf("with open(%s, 'w') as handle:", shQuote(output_path)),
    "    json.dump(payload, handle)"
  )
  writeLines(runner_lines, con = runner_path, useBytes = TRUE)

  status <- tryCatch(
    system2(
      command = python,
      args = c(shQuote(runner_path)),
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(e) {
      stop(
        "SCA execution failed. Ensure the Python executable and the `shannonca` package are available. ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  exit_code <- attr(status, "status") %||% 0L
  if (!identical(exit_code, 0L) || !file.exists(output_path)) {
    stop(
      "SCA execution failed. Ensure the Python package `shannonca` is installed. ",
      paste(status, collapse = "\n"),
      call. = FALSE
    )
  }

  result <- jsonlite::read_json(output_path, simplifyVector = TRUE)
  embedding <- as.matrix(result$reduction)
  rownames(embedding) <- colnames(expr)
  embedding
}

#' Detect rare cells with native or optional rare-cell backends
#'
#' @param object A \code{Seurat} object.
#' @param method Rare-cell method. Supported values are \code{"gini"},
#'   \code{"sccad"}, \code{"sca"}, \code{"gapclust"}, and
#'   \code{"challenging_groups"}.
#' @param group_by Optional metadata column used with
#'   \code{method = "challenging_groups"}.
#' @param reduction Reduction used by graph-based methods. Defaults to
#'   \code{"harmony"} when present, otherwise \code{"pca"}.
#' @param dims Optional embedding dimensions to use.
#' @param assay Assay used to extract expression values.
#' @param layer Layer used to extract expression values for score-based methods.
#' @param nfeatures Number of rare-aware genes to use for score construction.
#' @param min_cells Minimum number of cells a gene must be detected in before it
#'   is considered by score-based methods.
#' @param max_fraction Maximum expressing-cell fraction for score-based rare
#'   genes.
#' @param threshold Optional explicit threshold on the rare-cell score. When
#'   \code{NULL}, the function uses the upper-IQR rule.
#' @param k Number of neighbors for graph-based methods.
#' @param seed Random seed used by stochastic backends.
#' @param sccad_python Optional Python executable used by the scCAD backend.
#' @param sccad_script Optional path to the upstream \code{scCAD.py} script.
#' @param sccad_normalization Whether scCAD should normalize the provided
#'   matrix internally.
#' @param sccad_rare_h Rare threshold passed to scCAD.
#' @param sccad_merge_h Merge threshold passed to scCAD.
#' @param sccad_overlap_h Overlap threshold passed to scCAD.
#' @param gapclust_k Upper limit of the minor-cluster size used by GapClust.
#' @param sca_python Optional Python executable used by the SCA backend.
#' @param sca_n_comps Number of SCA components used before rarity scoring.
#' @param sca_iters Number of SCA iterations.
#' @param sca_nbhd_size Neighborhood size passed to SCA.
#' @param sca_model Scoring model passed to SCA.
#'
#' @return A data frame with one row per cell, including a \code{rare_score}
#'   column and a logical \code{rare_cell} flag.
#'
#' @examples
#' \dontrun{
#' data("pbmc_small", package = "Shennong")
#' rare_tbl <- sn_detect_rare_cells(pbmc_small, method = "gini")
#' head(rare_tbl)
#' }
#'
#' @export
sn_detect_rare_cells <- function(object,
                                 method = c("gini", "sccad", "sca", "gapclust", "challenging_groups"),
                                 group_by = NULL,
                                 reduction = .sn_default_metric_reduction(object),
                                 dims = NULL,
                                 assay = "RNA",
                                 layer = "data",
                                 nfeatures = 200,
                                 min_cells = 3,
                                 max_fraction = 0.1,
                                 threshold = NULL,
                                 k = 20,
                                 seed = 717,
                                 sccad_python = NULL,
                                 sccad_script = NULL,
                                 sccad_normalization = FALSE,
                                 sccad_rare_h = 0.01,
                                 sccad_merge_h = 0.3,
                                 sccad_overlap_h = 0.7,
                                 gapclust_k = 200,
                                 sca_python = NULL,
                                 sca_n_comps = 20,
                                 sca_iters = 3,
                                 sca_nbhd_size = 15,
                                 sca_model = "wilcoxon") {
  if (!inherits(object, "Seurat")) {
    stop("Input must be a Seurat object.", call. = FALSE)
  }

  method <- rlang::arg_match(method)
  expr <- .sn_get_seurat_layer_data(object = object, assay = assay, layer = layer)
  cell_ids <- colnames(object)

  if (identical(method, "gini")) {
    gini_tbl <- .sn_detect_rare_features_gini(
      expr = expr,
      nfeatures = nfeatures,
      min_cells = min_cells,
      max_fraction = max_fraction
    )
    if (nrow(gini_tbl) == 0) {
      return(data.frame(
        cell_id = cell_ids,
        method = method,
        rare_score = 0,
        rare_cell = FALSE,
        stringsAsFactors = FALSE
      ))
    }

    feature_mat <- expr[gini_tbl$feature, , drop = FALSE]
    feature_scale <- pmax(Matrix::rowMeans(feature_mat), 1e-8)
    normalized <- Matrix::Diagonal(x = 1 / feature_scale) %*% feature_mat
    rare_score <- Matrix::colMeans(normalized)
  } else if (identical(method, "sccad")) {
    sccad_result <- .sn_run_sccad(
      expr = expr,
      cell_ids = cell_ids,
      gene_ids = rownames(expr),
      python = sccad_python,
      script = sccad_script,
      normalization = sccad_normalization,
      seed = seed,
      rare_h = sccad_rare_h,
      merge_h = sccad_merge_h,
      overlap_h = sccad_overlap_h
    )
    rare_membership <- unique(unlist(sccad_result$rare_sets, use.names = FALSE))
    subcluster_score <- stats::setNames(
      as.numeric(sccad_result$scores),
      unique(sccad_result$sub_clusters)
    )
    rare_score <- unname(subcluster_score[as.character(sccad_result$sub_clusters)])
    rare_flag <- cell_ids %in% rare_membership
    return(data.frame(
      cell_id = cell_ids,
      method = method,
      rare_score = as.numeric(rare_score),
      rare_cell = rare_flag,
      subcluster = as.character(sccad_result$sub_clusters),
      stringsAsFactors = FALSE
    ))
  } else if (identical(method, "gapclust")) {
    return(transform(
      .sn_run_gapclust(expr = expr, k = gapclust_k),
      method = method
    ))
  } else if (identical(method, "sca")) {
    sca_embedding <- .sn_run_sca(
      expr = expr,
      python = sca_python,
      n_comps = sca_n_comps,
      iters = sca_iters,
      nbhd_size = sca_nbhd_size,
      model = sca_model,
      seed = seed
    )
    rare_score <- .sn_score_embedding_rarity(sca_embedding, k = k)
  } else {
    if (is.null(group_by)) {
      stop("`group_by` must be supplied when `method = \"challenging_groups\"`.", call. = FALSE)
    }
    group_tbl <- sn_identify_challenging_groups(
      x = object,
      group_by = group_by,
      reduction = reduction,
      dims = dims,
      k = k,
      neighbor_method = "auto",
      seed = seed
    )
    group_scores <- stats::setNames(group_tbl$challenge_score, group_tbl[[group_by]])
    rare_score <- unname(group_scores[as.character(object[[group_by, drop = TRUE]])])
  }

  score_threshold <- threshold %||% as.numeric(stats::quantile(rare_score, 0.75, na.rm = TRUE) + 1.5 * stats::IQR(rare_score, na.rm = TRUE))
  if (!is.finite(score_threshold)) {
    score_threshold <- Inf
  }

  data.frame(
    cell_id = cell_ids,
    method = method,
    rare_score = as.numeric(rare_score),
    rare_cell = as.numeric(rare_score) >= score_threshold,
    stringsAsFactors = FALSE
  )
}

.sn_select_variable_features <- function(object,
                                         nfeatures = 3000,
                                         split_by = NULL,
                                         assay = NULL,
                                         layer = NULL,
                                         verbose = TRUE) {
  if (is_null(split_by) || identical(split_by, "global")) {
    object <- Seurat::FindVariableFeatures(object, nfeatures = nfeatures, verbose = verbose)
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
    current_object <- Seurat::FindVariableFeatures(
      current_object,
      nfeatures = nfeatures,
      verbose = FALSE
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
#' @param integration_method Batch-integration backend used when \code{batch}
#'   is supplied. Supported values are \code{"harmony"},
#'   \code{"coralysis"}, \code{"seurat_cca"}, \code{"seurat_rpca"},
#'   \code{"scvi"}, \code{"scanvi"}, and \code{"totalvi"}. \code{"mmochi"} is
#'   accepted as a CITE-seq convenience alias and requires
#'   \code{modality = "cite_seq"}.
#'   \code{"harmony"} preserves the historical Shennong behavior.
#'   \code{"coralysis"} runs native Coralysis on the selected log-normalized
#'   feature set and stores the integrated embedding as the \code{"coralysis"}
#'   reduction. \code{"scvi"} and \code{"scanvi"} export the selected count
#'   matrix to a pixi-managed scverse environment under \code{~/.shennong/pixi/},
#'   run the Python backend, and import the latent representation as a Seurat
#'   reduction. \code{"totalvi"} is used for RNA+ADT CITE-seq workflows and is
#'   usually selected through \code{modality = "cite_seq"} and
#'   \code{multimodal_method = "totalvi"}.
#' @param integration_control Optional named list of backend-specific
#'   parameters. For \code{"coralysis"}, use \code{icp_args} for
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
#'   \code{unlabeled_category}. \code{"totalvi"} additionally accepts
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
#'   Use \code{sn_pixi_paths()} to inspect the
#'   generated directory layout, \code{sn_pixi_config_path()} to inspect the
#'   bundled \code{inst/pixi/} config, \code{sn_ensure_pixi()} to preinstall
#'   pixi, and \code{sn_configure_pixi_mirror()} to set Shennong-level mirrors.
#' @param nfeatures Number of variable features to select.
#' @param hvg_features Optional character vector of user-supplied features to
#'   force into the selected backend feature set. For PCA-based workflows this
#'   is also the feature set used for scaling/PCA. These features are merged
#'   with internally selected HVGs and any rare-aware features after validating
#'   that they are present in \code{object}.
#' @param vars_to_regress Covariates to regress out in \code{ScaleData}.
#' @param resolution Resolution parameter for \code{FindClusters}.
#' @param cluster_algorithm Community-detection algorithm passed to
#'   \code{Seurat::FindClusters()}. Supported names are \code{"louvain"}
#'   (Seurat algorithm 1), \code{"louvain_multilevel"} (algorithm 2),
#'   \code{"slm"} (algorithm 3), and \code{"leiden"} (algorithm 4). Numeric
#'   values 1 through 4 are also accepted.
#' @param cluster_name Optional metadata column name for the cluster labels.
#'   Defaults to Seurat's \code{"seurat_clusters"} behavior.
#' @param cluster_n_start,cluster_n_iter Number of starts and iterations passed
#'   to \code{Seurat::FindClusters()}.
#' @param cluster_random_seed Random seed passed to
#'   \code{Seurat::FindClusters()}.
#' @param cluster_group_singletons Whether \code{Seurat::FindClusters()}
#'   should group singletons into the nearest cluster.
#' @param leiden_method Leiden implementation passed to
#'   \code{Seurat::FindClusters()} when \code{cluster_algorithm = "leiden"}.
#' @param leiden_objective_function Leiden objective function passed to
#'   \code{Seurat::FindClusters()}.
#' @param cluster_control Optional named list of additional
#'   \code{Seurat::FindClusters()} arguments. Values here override Shennong's
#'   generated defaults.
#' @param reuse Logical; when \code{TRUE}, reuse previously recorded
#'   \code{sn_run_cluster()} stages if their stored input signatures still match
#'   the current call. This lets resolution-only changes start at clustering,
#'   integration-method changes start at integration, and HVG changes start at
#'   feature selection instead of rerunning all earlier steps.
#' @param rerun_from Optional stage name forcing recomputation from that stage
#'   onward while still allowing earlier matching stages to be reused. Supported
#'   values are \code{"normalize"}, \code{"cell_cycle"}, \code{"hvg"},
#'   \code{"pca"}, \code{"adt"}, \code{"integration"}, \code{"neighbors"},
#'   \code{"clusters"}, and \code{"umap"}.
#' @param auto_install Logical; when \code{TRUE}, install missing optional
#'   clustering dependencies such as \pkg{leidenbase} before the relevant stage.
#' @param install_repos CRAN-like repositories used when \code{auto_install}
#'   installs CRAN packages.
#' @param install_ask Passed to \code{BiocManager::install()} when
#'   \code{auto_install} installs Bioconductor packages through
#'   \code{sn_install_dependencies()}.
#' @param hvg_group_by Optional metadata column used to compute highly variable
#'   genes within groups before merging and ranking them. When \code{NULL} and
#'   \code{batch} is supplied, Shennong reuses \code{batch} by default. Use
#'   \code{NULL} with \code{batch = NULL} to compute HVGs on the full object.
#' @param rare_feature_method Optional rare-cell-aware feature methods appended
#'   to the base HVG set before PCA/clustering. Supported values are
#'   \code{"none"}, \code{"gini"}, and \code{"local_markers"}.
#' @param rare_feature_group_by Optional metadata column used to define groups
#'   for \code{"local_markers"}. When \code{NULL}, Shennong builds a temporary
#'   coarse clustering from the base HVGs.
#' @param rare_feature_n Number of rare-aware features to add per selected
#'   method. For example, \code{c("gini", "local_markers")} with
#'   \code{rare_feature_n = 50} can contribute up to 100 rare-aware features
#'   before de-duplication.
#' @param rare_feature_control Named list of advanced rare-feature thresholds.
#'   Supported fields are \code{group_max_fraction}, \code{group_max_cells},
#'   \code{gene_max_fraction}, and \code{min_cells}.
#' @param block_genes Either a character vector of predefined bundled signature
#'   categories (for example \code{c("ribo","mito")}) or a custom vector of
#'   gene symbols to exclude from HVGs.
#' @param theta The \code{theta} parameter for \code{harmony::RunHarmony}, controlling batch
#'   diversity preservation vs. correction. Used only when
#'   \code{integration_method = "harmony"}.
#' @param group_by_vars Optional column name or character vector passed to
#'   \code{harmony::RunHarmony(group.by.vars = ...)}. Defaults to \code{batch}
#'   and is used only when \code{integration_method = "harmony"}.
#' @param npcs Number of PCs to compute in \code{RunPCA}.
#' @param dims A numeric vector of PCs (dimensions) to use for neighbor search,
#'   clustering, and UMAP.
#' @param species Optional species label. Used when block genes must be resolved
#'   from built-in signatures.
#' @param assay Assay used for clustering. Defaults to \code{"RNA"}.
#' @param layer Layer used as the input count matrix. Defaults to \code{"counts"}.
#' @param modality Workflow modality. \code{"rna"} runs the standard RNA-only
#'   workflow. \code{"cite_seq"} enables paired RNA+ADT workflows selected by
#'   \code{multimodal_method}.
#' @param multimodal_method CITE-seq backend used when
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
#' @param adt_assay Assay containing antibody-derived tag counts for
#'   \code{modality = "cite_seq"}.
#' @param adt_layer Layer in \code{adt_assay} used as ADT counts.
#' @param adt_features Optional ADT/protein features used by CITE-seq backends.
#'   Defaults to all features in \code{adt_assay}.
#' @param adt_npcs Number of ADT PCs to compute for \code{modality = "cite_seq"}.
#' @param adt_dims Numeric vector of ADT PCs used in weighted nearest-neighbor
#'   graph construction. Defaults to \code{seq_len(min(18, adt_npcs))}.
#' @param wnn_control Optional named list of additional
#'   \code{Seurat::FindMultiModalNeighbors()} arguments used only when
#'   \code{modality = "cite_seq"}. Values here override Shennong's generated
#'   defaults.
#' @param umap_control Optional named list of additional
#'   \code{Seurat::RunUMAP()} arguments. Values here override Shennong's
#'   generated defaults, for example \code{n.neighbors}, \code{min.dist},
#'   \code{spread}, \code{metric}, \code{seed.use}, or \code{reduction.name}.
#' @param return_cluster If \code{TRUE}, return only the cluster_by assignments.
#' @param verbose Whether to print/log progress messages.
#'
#' @return A \code{Seurat} object with clustering results and embeddings, or a
#'   cluster_by vector if \code{return_cluster = TRUE}.
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
                           integration_method = c("harmony", "coralysis", "seurat_cca", "seurat_rpca", "scvi", "scanvi", "totalvi", "mmochi"),
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
                           return_cluster = FALSE,
                           verbose = TRUE) {
  check_installed("Seurat")
  check_installed("HGNChelper")

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

  integration_method_supplied <- !missing(integration_method)
  normalization_method <- match.arg(normalization_method)
  integration_method <- match.arg(integration_method)
  modality <- match.arg(modality)
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
  leiden_method <- match.arg(leiden_method)
  leiden_objective_function <- match.arg(leiden_objective_function)
  rerun_from <- .sn_resolve_cluster_rerun_from(rerun_from)
  if (!is.list(integration_control)) {
    stop("`integration_control` must be a named list.", call. = FALSE)
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
    if (normalization_method == "sctransform" && integration_method != "harmony") {
      stop("SCTransform integration is currently supported only with `integration_method = \"harmony\"`.", call. = FALSE)
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

  predefined_genesets <- c(
    "tcr", "immunoglobulins", "ribo", "mito",
    "heatshock", "noncoding", "pseudogenes", "g1s", "g2m"
  )

  if (!isTRUE(needs_rna_workflow)) {
    block_genes <- NULL
    rare_feature_method <- "none"
  }

  if (normalization_method == "sctransform" && !is_null(block_genes)) {
    .sn_log_warn("`block_genes` is only applied in log-normalization workflows; ignoring it for SCTransform.")
    block_genes <- NULL
  }
  if (normalization_method == "sctransform" && !identical(rare_feature_method, "none")) {
    .sn_log_warn("`rare_feature_method` is only applied in log-normalization workflows; ignoring it for SCTransform.")
    rare_feature_method <- "none"
  }

  if (isTRUE(needs_rna_workflow) && !is_null(block_genes)) {
    species <- sn_get_species(object = object, species = species)
    if (verbose) .sn_log_info("Processing blocked genes.")

    if (is_character(block_genes) && all(block_genes %in% predefined_genesets)) {
      block_genes <- sn_get_signatures(species = species, category = block_genes)
      if (verbose) {
        .sn_log_info("Loaded {length(block_genes)} blocked genes from built-in sets.")
      }
    } else {
      checked <- HGNChelper::checkGeneSymbols(block_genes, species = species)
      valid_genes <- checked$Suggested.Symbol[!is_na(checked$Suggested.Symbol)]
      invalid <- setdiff(block_genes, valid_genes)

      if (length(invalid) > 0) {
        .sn_log_warn(
          "Removed {length(invalid)} invalid genes from block list (e.g. {paste(utils::head(invalid, 3), collapse=', ')})."
        )
      }
      block_genes <- unique(valid_genes)
      if (verbose) .sn_log_info("Using {length(block_genes)} custom blocked genes.")
    }
  }

  hvg_candidate_nfeatures <- nfeatures
  if (!is_null(block_genes)) {
    hvg_candidate_nfeatures <- min(
      nrow(object),
      nfeatures + length(intersect(block_genes, rownames(object)))
    )
  }

  hvg <- character(0)
  hvg_signature <- NULL
  pca_signature <- NULL

  normalization_signature <- list(
    method = normalization_method,
    assay = assay,
    layer = layer,
    layer_source_version = 2L,
    nfeatures = if (identical(normalization_method, "sctransform")) nfeatures else NULL,
    vars_to_regress = if (identical(normalization_method, "sctransform")) vars_to_regress else NULL,
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
        variable.features.n = nfeatures,
        vars.to.regress = vars_to_regress,
        verbose = verbose,
        seed.use = 717
      )
      if (length(user_hvg) > 0L) {
        sct_args$return.only.var.genes <- FALSE
      }
      object <- .sn_with_auto_future_globals(
        .sn_call_with_symbolic_object(
          fun_call = quote(Seurat::SCTransform),
          object = object,
          args = sct_args
        ),
        object = object,
        context = "SCTransform",
        verbose = verbose
      )
      object <- .sn_record_cluster_stage(object, "normalize", normalization_signature)
    }

    hvg_signature <- list(
      method = "sctransform",
      nfeatures = nfeatures,
      user_hvg = user_hvg,
      missing_user_hvg = user_hvg_info$missing,
      normalization = normalization_signature
    )
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
      hvg <- unique(c(Seurat::VariableFeatures(object = object), user_hvg))
      object <- .sn_record_cluster_stage(object, "hvg", hvg_signature, selected_features = hvg)
    }
    Seurat::VariableFeatures(object = object) <- hvg
    object@misc$hvg_selection <- list(
      method = "sctransform",
      base_hvg_n = nfeatures,
      selected_features = hvg,
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
        object <- Seurat::RunPCA(
          object,
          npcs = npcs,
          features = hvg,
          verbose = verbose,
          seed.use = 717
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
      object <- Seurat::NormalizeData(object = object, assay = assay, layer = layer, verbose = verbose)
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
      rare_feature_resolution = min(resolution, 0.4),
      user_hvg = user_hvg,
      missing_user_hvg = user_hvg_info$missing,
      normalization = normalization_signature
    )
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
      user_features = user_hvg,
      missing_user_features = user_hvg_info$missing
    )

    if (isTRUE(needs_rna_pca)) {
      pca_signature <- list(
        method = normalization_method,
        features = hvg,
        vars_to_regress = vars_to_regress,
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
        object <- Seurat::ScaleData(
          object = object,
          vars.to.regress = vars_to_regress,
          features = hvg,
          verbose = verbose
        )

        if (verbose) .sn_log_info("[5/6] Running PCA.")
        object <- Seurat::RunPCA(
          object,
          npcs = npcs,
          features = hvg,
          verbose = verbose,
          seed.use = 717
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
      object <- Seurat::NormalizeData(
        object = object,
        assay = adt_assay,
        normalization.method = "CLR",
        margin = 2,
        verbose = verbose
      )
      if (isTRUE(needs_adt_pca)) {
        object <- Seurat::ScaleData(
          object = object,
          assay = adt_assay,
          features = adt_feature_set,
          verbose = verbose
        )
        object <- Seurat::RunPCA(
          object = object,
          assay = adt_assay,
          features = adt_feature_set,
          npcs = adt_npcs,
          reduction.name = "apca",
          reduction.key = "apca_",
          verbose = verbose,
          seed.use = 717
        )
        object <- .sn_record_cluster_stage(object, "adt", adt_signature, reduction = "apca")
      } else {
        object <- .sn_record_cluster_stage(object, "adt", adt_signature)
      }
    }
    SeuratObject::DefaultAssay(object = object) <- if (isTRUE(needs_rna_workflow) && identical(normalization_method, "sctransform")) "SCT" else assay
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
      integration <- .sn_run_batch_integration(
        object = object,
        method = multimodal_method,
        batch = batch,
        reduction = NULL,
        features = adt_feature_set,
        assay = adt_assay,
        dims = adt_dims,
        npcs = adt_npcs,
        theta = theta,
        group_by_vars = group_by_vars,
        integration_control = integration_control,
        verbose = verbose
      )
      object <- integration$object
      reduction <- integration$reduction
      object@misc$integration$modality <- "cite_seq"
      object@misc$integration$protein_assay <- adt_assay
      object@misc$integration$protein_layer <- adt_layer
      object <- .sn_record_cluster_stage(object, "integration", integration_signature, reduction = reduction)
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
      integration <- .sn_run_mmochi_integration(
        object = object,
        batch = batch,
        assay = adt_assay,
        protein_assay = adt_assay,
        protein_features = adt_feature_set,
        dims = adt_dims,
        npcs = adt_npcs,
        integration_control = mmochi_control,
        verbose = verbose
      )
      object <- integration$object
      reduction <- integration$reduction
      object@misc$integration$modality <- "cite_seq"
      object@misc$integration$protein_assay <- adt_assay
      object@misc$integration$protein_layer <- mmochi_control$protein_layer
      object <- .sn_record_cluster_stage(object, "integration", integration_signature, reduction = reduction)
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
      integration <- .sn_run_scvi_integration(
        object = object,
        method = "totalvi",
        batch = batch,
        features = hvg,
        assay = assay,
        integration_control = totalvi_control,
        verbose = verbose
      )
      object <- integration$object
      reduction <- integration$reduction
      object@misc$integration$modality <- "cite_seq"
      object@misc$integration$protein_assay <- adt_assay
      object@misc$integration$protein_layer <- adt_layer
      object <- .sn_record_cluster_stage(object, "integration", integration_signature, reduction = reduction)
    }
    reduction_signature <- integration_signature
  } else {
    integration_signature <- list(
      method = integration_method,
      batch = batch,
      features = hvg,
      dims = dims,
      npcs = npcs,
      theta = theta,
      group_by_vars = group_by_vars,
      integration_control = integration_control,
      store_sce = if (identical(integration_method, "coralysis")) {
        .sn_coralysis_store_sce(integration_control)
      } else {
        NULL
      },
      pca = pca_signature
    )
    if (is_null(x = batch)) {
      reduction <- "pca"
    } else if (.sn_can_reuse_cluster_stage(
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
      if (verbose) .sn_log_info("[5/6] Reusing {reduction} integration reduction.")
    } else {
      if (verbose) .sn_log_info("[5/6] Running {integration_method} integration.")
      integration <- .sn_run_batch_integration(
        object = object,
        method = integration_method,
        batch = batch,
        reduction = if (isTRUE(needs_rna_pca)) "pca" else NULL,
        features = hvg,
        assay = assay,
        dims = dims,
        npcs = npcs,
        theta = theta,
        group_by_vars = group_by_vars,
        integration_control = integration_control,
        verbose = verbose
      )
      object <- integration$object
      reduction <- integration$reduction
      object <- .sn_record_cluster_stage(object, "integration", integration_signature, reduction = reduction)
    }
    reduction_signature <- if (is_null(x = batch)) pca_signature else integration_signature
  }

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
      object <- Seurat::FindNeighbors(object, reduction = reduction, dims = dims, verbose = verbose)
    }
    graph_names_after <- names(object@graphs)
    graph_names <- setdiff(graph_names_after, graph_names_before)
    if (length(graph_names) == 0L) {
      graph_names <- graph_names_after
    }
    snn_graph <- if (identical(modality, "cite_seq") && identical(multimodal_method, "wnn") && "wsnn" %in% graph_names_after) {
      "wsnn"
    } else {
      snn_candidates <- grep("_snn$", graph_names, value = TRUE)
      if (length(snn_candidates) > 0L) {
        snn_candidates[[1]]
      } else {
        utils::tail(graph_names, n = 1L)
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
    object <- .sn_call_with_symbolic_object(
      fun_call = quote(Seurat::FindClusters),
      object = object,
      args = find_clusters_args
    )
    object <- .sn_record_cluster_stage(object, "clusters", cluster_signature, cluster_column = cluster_column)
  }

  if (return_cluster) {
    object <- restore_analysis_inputs(object)
    if (verbose) .sn_log_info("Integration completed successfully.")
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
      if (identical(modality, "cite_seq") && identical(multimodal_method, "wnn")) {
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
    object <- restore_analysis_inputs(object)
    if (verbose) .sn_log_info("Integration completed successfully.")
    return(.sn_log_seurat_command(object = object, assay = assay, name = "sn_run_cluster"))
  }
}

.sn_find_transfer_anchors_backend <- function(...) {
  Seurat::FindTransferAnchors(...)
}

.sn_transfer_data_backend <- function(...) {
  Seurat::TransferData(...)
}

.sn_prepare_label_transfer_name <- function(prefix, cells) {
  paste0(prefix, "_", .sn_metadata_suffix(cells))
}

.sn_transfer_labels_scanvi <- function(object,
                                       reference,
                                       label_by,
                                       prediction_prefix,
                                       assay = NULL,
                                       batch_by = NULL,
                                       features = NULL,
                                       transfer_control = list(),
                                       return_anchors = FALSE,
                                       verbose = TRUE,
                                       method_name = "scanvi") {
  if (!inherits(reference, "Seurat")) {
    stop("`reference` must be a Seurat object for scANVI/scArches label_by transfer.", call. = FALSE)
  }
  if (!label_by %in% colnames(reference[[]])) {
    stop(glue("`label_by` column '{label_by}' was not found in `reference` metadata."), call. = FALSE)
  }
  if (!is.list(transfer_control)) {
    stop("`transfer_control` must be a named list.", call. = FALSE)
  }

  assay <- assay %||% SeuratObject::DefaultAssay(object = object)
  reference_assay <- transfer_control$reference_assay %||% assay
  query_assay <- transfer_control$query_assay %||% assay
  if (!reference_assay %in% names(reference@assays)) {
    stop(glue("Reference assay '{reference_assay}' was not found."), call. = FALSE)
  }
  if (!query_assay %in% names(object@assays)) {
    stop(glue("Query assay '{query_assay}' was not found."), call. = FALSE)
  }

  reference_prefix <- transfer_control$reference_prefix %||% "reference"
  query_prefix <- transfer_control$query_prefix %||% "query"
  reference_cells <- .sn_prepare_label_transfer_name(reference_prefix, colnames(reference))
  query_cells <- .sn_prepare_label_transfer_name(query_prefix, colnames(object))

  reference <- Seurat::RenameCells(reference, new.names = reference_cells)
  object_for_transfer <- Seurat::RenameCells(object, new.names = query_cells)
  reference$.sn_transfer_role <- "reference"
  object_for_transfer$.sn_transfer_role <- "query"
  transfer_label_by <- transfer_control$label_by %||% transfer_control$label_by %||% ".sn_transfer_label"
  unlabeled_category <- transfer_control$unlabeled_category %||% "Unknown"
  reference[[transfer_label_by]] <- as.character(reference[[label_by, drop = TRUE]])
  object_for_transfer[[transfer_label_by]] <- unlabeled_category

  batch_by <- batch_by %||% transfer_control$batch_by %||% ".sn_transfer_batch"
  if (!batch_by %in% colnames(reference[[]])) {
    reference[[batch_by]] <- "reference"
  }
  if (!batch_by %in% colnames(object_for_transfer[[]])) {
    object_for_transfer[[batch_by]] <- "query"
  }

  common_features <- intersect(rownames(reference[[reference_assay]]), rownames(object_for_transfer[[query_assay]]))
  feature_set <- features %||% common_features
  feature_set <- intersect(feature_set, common_features)
  if (length(feature_set) < 2L) {
    stop("scANVI/scArches label_by transfer requires at least two shared features.", call. = FALSE)
  }

  old_reference_assay <- SeuratObject::DefaultAssay(reference)
  old_query_assay <- SeuratObject::DefaultAssay(object_for_transfer)
  on.exit(SeuratObject::DefaultAssay(reference) <- old_reference_assay, add = TRUE)
  on.exit(SeuratObject::DefaultAssay(object_for_transfer) <- old_query_assay, add = TRUE)
  SeuratObject::DefaultAssay(reference) <- reference_assay
  SeuratObject::DefaultAssay(object_for_transfer) <- query_assay
  combined <- merge(reference, y = object_for_transfer, merge.data = FALSE)
  combined_assay <- SeuratObject::DefaultAssay(combined)

  integration_control <- transfer_control
  integration_control$label_by <- transfer_label_by
  integration_control$label_by <- transfer_label_by
  integration_control$unlabeled_category <- unlabeled_category
  integration_control$reduction <- integration_control$reduction %||% method_name

  fit <- .sn_run_scvi_integration(
    object = combined,
    method = "scanvi",
    batch = batch_by,
    features = feature_set,
    assay = combined_assay,
    integration_control = integration_control,
    verbose = verbose
  )
  combined <- fit$object
  prediction_col <- transfer_control$prediction_col %||% "scanvi_prediction"
  if (!prediction_col %in% colnames(combined[[]])) {
    stop("scANVI/scArches backend did not return a prediction column.", call. = FALSE)
  }

  query_predictions <- combined[[prediction_col, drop = TRUE]][query_cells]
  metadata <- data.frame(row.names = colnames(object))
  metadata[[paste0(prediction_prefix, "_label")]] <- as.character(query_predictions)
  object <- Seurat::AddMetaData(object = object, metadata = metadata)
  object@misc$label_transfer[[prediction_prefix]] <- list(
    method = method_name,
    label_by = label_by,
    label_col = label_by,
    batch_by = batch_by,
    transfer_label_by = transfer_label_by,
    unlabeled_category = unlabeled_category,
    prediction_columns = colnames(metadata),
    run_dir = combined@misc$integration$run_dir %||% NULL,
    output_h5ad = combined@misc$integration$output_h5ad %||% NULL,
    transfer_control = transfer_control
  )

  if (isTRUE(return_anchors)) {
    return(list(query = object, combined = combined, prediction_col = prediction_col))
  }
  object
}

.sn_coralysis_reference_mapping_backend <- function(...) {
  Coralysis::ReferenceMapping(...)
}

.sn_metadata_suffix <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", as.character(x))
  x <- gsub("^_+|_+$", "", x)
  x[!nzchar(x)] <- "value"
  x
}

.sn_get_seurat_logcounts_sce <- function(object,
                                         assay = NULL,
                                         layer = "data",
                                         verbose = TRUE) {
  check_installed(c("SingleCellExperiment", "SummarizedExperiment"))

  assay <- assay %||% SeuratObject::DefaultAssay(object = object)
  if (identical(layer, "data") && !"data" %in% SeuratObject::Layers(object = object[[assay]])) {
    if (verbose) .sn_log_info("Normalizing query object before Coralysis reference mapping.")
    object <- Seurat::NormalizeData(object = object, assay = assay, verbose = verbose)
  }

  expr <- SeuratObject::LayerData(object = object, assay = assay, layer = layer)
  SingleCellExperiment::SingleCellExperiment(
    assays = list(logcounts = expr),
    colData = object[[]]
  )
}

.sn_sync_coralysis_reference_label <- function(ref_sce,
                                              reference = NULL,
                                              label_col = NULL) {
  if (
    inherits(reference, "Seurat") &&
      is.character(label_col) &&
      length(label_col) == 1L &&
      label_col %in% colnames(reference[[]])
  ) {
    labels <- reference[[label_col, drop = TRUE]]
    names(labels) <- colnames(reference)
    SummarizedExperiment::colData(ref_sce)[[label_col]] <- labels[colnames(ref_sce)]
  }
  ref_sce
}

.sn_get_coralysis_reference_sce <- function(reference,
                                           reference_assay = NULL,
                                           reference_layer = "data",
                                           verbose = TRUE,
                                           label_col = NULL) {
  if (inherits(reference, "SingleCellExperiment")) {
    return(reference)
  }
  if (!inherits(reference, "Seurat")) {
    stop("`reference` must be a Seurat or SingleCellExperiment object for Coralysis label_by transfer.", call. = FALSE)
  }
  if (inherits(reference@misc$coralysis, "SingleCellExperiment")) {
    return(.sn_sync_coralysis_reference_label(
      ref_sce = reference@misc$coralysis,
      reference = reference,
      label_col = label_col
    ))
  }
  stop(
    "Coralysis label_by transfer requires a native Coralysis-trained reference stored under `reference@misc$coralysis`.\n",
    "Run `sn_run_cluster(reference, batch = ..., integration_method = \"coralysis\")` first, ",
    "or avoid `integration_control = list(store_sce = FALSE)` if the object should be used as a reference.",
    call. = FALSE
  )
}

.sn_minimize_coralysis_pca_model <- function(pca_model) {
  if (!is.list(pca_model)) {
    return(pca_model)
  }
  keep <- intersect(c("x", "center", "scale", "rotation"), names(pca_model))
  pca_model[keep]
}

.sn_prepare_coralysis_label_transfer_reference <- function(object,
                                                          label_by,
                                                          metadata_columns = NULL,
                                                          keep_umap_model = FALSE,
                                                          verbose = TRUE) {
  check_installed(c("SingleCellExperiment", "SummarizedExperiment", "S4Vectors"))

  ref_sce <- .sn_get_coralysis_reference_sce(
    reference = object,
    verbose = verbose,
    label_col = label_by
  )
  coldata <- SummarizedExperiment::colData(ref_sce)
  if (!label_by %in% colnames(coldata)) {
    stop(glue("`label_by` column '{label_by}' was not found in the Coralysis reference colData."), call. = FALSE)
  }

  coralysis <- S4Vectors::metadata(ref_sce)$coralysis
  required <- c("models", "pca.model", "pca.params")
  missing_required <- required[vapply(coralysis[required], is.null, logical(1))]
  if (length(missing_required) > 0L) {
    stop(
      "The Coralysis reference is missing required field(s): ",
      paste(missing_required, collapse = ", "),
      call. = FALSE
    )
  }
  if (is.null(coralysis$pca.params$select.icp.tables)) {
    stop("The Coralysis reference is missing `pca.params$select.icp.tables`.", call. = FALSE)
  }

  keep_cols <- unique(c(label_by, metadata_columns))
  keep_cols <- intersect(keep_cols, colnames(coldata))
  coldata <- coldata[, keep_cols, drop = FALSE]

  reference <- SingleCellExperiment::SingleCellExperiment(
    assays = list(),
    rowData = S4Vectors::DataFrame(row.names = rownames(ref_sce)),
    colData = coldata
  )
  reference_metadata <- S4Vectors::metadata(reference)
  reference_metadata$coralysis <- list(
    models = coralysis$models,
    pca.model = .sn_minimize_coralysis_pca_model(coralysis$pca.model),
    pca.params = list(select.icp.tables = coralysis$pca.params$select.icp.tables)
  )
  if (isTRUE(keep_umap_model) && !is.null(coralysis$umap.model)) {
    reference_metadata$coralysis$umap.model <- coralysis$umap.model
  }
  reference_metadata$shennong_reference <- list(
    method = "coralysis",
    label_by = label_by,
    n_features = nrow(reference),
    n_cells = ncol(reference),
    assays_dropped = TRUE,
    reduced_dims_dropped = TRUE,
    joint_probability_dropped = TRUE
  )
  S4Vectors::metadata(reference) <- reference_metadata
  reference
}

.sn_prepare_seurat_label_transfer_reference <- function(object,
                                                        label_by,
                                                        method,
                                                        assay = NULL,
                                                        layers = NULL,
                                                        features = NULL,
                                                        reduction = NULL,
                                                        metadata_columns = NULL) {
  check_installed("Seurat")
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object for this label-transfer method.", call. = FALSE)
  }
  if (!label_by %in% colnames(object[[]])) {
    stop(glue("`label_by` column '{label_by}' was not found in `object` metadata."), call. = FALSE)
  }

  assay <- assay %||% SeuratObject::DefaultAssay(object = object)
  if (!assay %in% names(object@assays)) {
    stop(glue("Assay '{assay}' was not found."), call. = FALSE)
  }

  available_layers <- SeuratObject::Layers(object = object[[assay]])
  if (is.null(layers)) {
    layers <- switch(
      method,
      seurat = intersect(c("counts", "data"), available_layers),
      scanvi = intersect("counts", available_layers),
      scarches = intersect("counts", available_layers)
    )
    if (length(layers) == 0L) {
      layers <- available_layers[1]
    }
  } else if (identical(layers, "all")) {
    layers <- available_layers
  }
  missing_layers <- setdiff(layers, available_layers)
  if (length(missing_layers) > 0L) {
    stop(
      "Layer(s) not found in assay '", assay, "': ",
      paste(missing_layers, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(reduction) && identical(method, "seurat") && "pca" %in% names(object@reductions)) {
    reduction <- "pca"
  }
  dimreducs <- reduction %||% character(0)

  old_assay <- SeuratObject::DefaultAssay(object = object)
  on.exit(SeuratObject::DefaultAssay(object = object) <- old_assay, add = TRUE)
  SeuratObject::DefaultAssay(object = object) <- assay
  reference <- Seurat::DietSeurat(
    object = object,
    assays = assay,
    layers = layers,
    features = features,
    dimreducs = dimreducs,
    graphs = character(0),
    misc = FALSE
  )

  keep_cols <- unique(c(label_by, metadata_columns))
  keep_cols <- intersect(keep_cols, colnames(reference[[]]))
  reference@meta.data <- reference@meta.data[, keep_cols, drop = FALSE]
  reference@misc$label_transfer_reference <- list(
    method = method,
    label_by = label_by,
    assay = assay,
    layers = layers,
    features = features,
    reduction = reduction
  )
  reference
}

.sn_write_label_transfer_reference <- function(reference,
                                               path = NULL,
                                               overwrite = FALSE,
                                               ...) {
  if (is.null(path)) {
    return(reference)
  }
  if (file.exists(path) && !isTRUE(overwrite)) {
    stop("`path` already exists. Set `overwrite = TRUE` to replace it: ", path, call. = FALSE)
  }
  sn_write(reference, path = path, ...)
  reference
}

#' Prepare a compact label-transfer reference
#'
#' \code{sn_prepare_label_transfer_reference()} converts a full analysis object
#' into a smaller reference object for \code{\link{sn_transfer_labels}}. For
#' native Coralysis, it returns a minimal \code{SingleCellExperiment} containing
#' the trained Coralysis models, PCA model, feature names, and selected
#' reference labels, while dropping the large reference assays, reduced
#' dimensions, and stored joint probabilities. For Seurat, scANVI, and scArches
#' workflows, it returns a slim Seurat object with only the selected assay
#' layers, labels, optional features, and optional reduction.
#'
#' @param object A Seurat object. For \code{method = "coralysis"}, an existing
#'   Coralysis-trained \code{SingleCellExperiment} is also accepted.
#' @param label_by Metadata column containing reference labels.
#' @param method Label-transfer backend the reference should support.
#' @param assay Assay to keep for Seurat/scANVI/scArches references.
#' @param layers Layers to keep for Seurat/scANVI/scArches references. Defaults
#'   to \code{counts} and \code{data} for Seurat transfer, and \code{counts} for
#'   scANVI/scArches. Use \code{"all"} to retain all layers in \code{assay}.
#' @param features Optional features to keep for Seurat/scANVI/scArches
#'   references.
#' @param reduction Optional dimensional reduction to keep. Defaults to
#'   \code{"pca"} for Seurat references when available.
#' @param metadata_columns Additional metadata columns to keep alongside
#'   \code{label_by}.
#' @param keep_umap_model For Coralysis references, keep the UMAP projection
#'   model if present. This is only needed when using
#'   \code{transfer_control = list(project.umap = TRUE)}.
#' @param path Optional output path. Use \code{.qs2} for serialized reference
#'   objects.
#' @param overwrite Logical; overwrite an existing \code{path}.
#' @param verbose Whether to print progress messages.
#' @param ... Additional arguments passed to \code{\link{sn_write}} when
#'   \code{path} is supplied.
#'
#' @return A compact Seurat or SingleCellExperiment reference object.
#'
#' @examples
#' \dontrun{
#' coral_ref <- sn_prepare_label_transfer_reference(
#'   reference,
#'   label_by = "cell_type",
#'   method = "coralysis",
#'   path = "data/processed/pbmc_coralysis_reference.qs2",
#'   overwrite = TRUE
#' )
#'
#' query <- sn_transfer_labels(
#'   query,
#'   reference = coral_ref,
#'   label_by = "cell_type",
#'   method = "coralysis"
#' )
#' }
#' @export
sn_prepare_label_transfer_reference <- function(object,
                                                label_by,
                                                method = c("coralysis", "seurat", "scanvi", "scarches"),
                                                assay = NULL,
                                                layers = NULL,
                                                features = NULL,
                                                reduction = NULL,
                                                metadata_columns = NULL,
                                                keep_umap_model = FALSE,
                                                path = NULL,
                                                overwrite = FALSE,
                                                verbose = TRUE,
                                                ...) {
  method <- match.arg(method)
  if (!is.character(label_by) || length(label_by) != 1L || !nzchar(label_by)) {
    stop("`label_by` must be a non-empty metadata column name.", call. = FALSE)
  }

  reference <- if (identical(method, "coralysis")) {
    .sn_prepare_coralysis_label_transfer_reference(
      object = object,
      label_by = label_by,
      metadata_columns = metadata_columns,
      keep_umap_model = keep_umap_model,
      verbose = verbose
    )
  } else {
    .sn_prepare_seurat_label_transfer_reference(
      object = object,
      label_by = label_by,
      method = method,
      assay = assay,
      layers = layers,
      features = features,
      reduction = reduction,
      metadata_columns = metadata_columns
    )
  }

  .sn_write_label_transfer_reference(
    reference = reference,
    path = path,
    overwrite = overwrite,
    ...
  )
}

.sn_transfer_labels_coralysis <- function(object,
                                          reference,
                                          label_col,
                                          prediction_prefix,
                                          reference_assay = NULL,
                                          query_assay = NULL,
                                          reference_layer = "data",
                                          query_layer = "data",
                                          transfer_control = list(),
                                          return_anchors = FALSE,
                                          verbose = TRUE) {
  check_installed("Coralysis")
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat query object for Coralysis label_by transfer.", call. = FALSE)
  }
  if (!is.list(transfer_control)) {
    stop("`transfer_control` must be a named list.", call. = FALSE)
  }

  ref_sce <- .sn_get_coralysis_reference_sce(
    reference = reference,
    reference_assay = reference_assay,
    reference_layer = reference_layer,
    verbose = verbose,
    label_col = label_col
  )
  if (!label_col %in% colnames(SummarizedExperiment::colData(ref_sce))) {
    stop(glue("`label_col` column '{label_col}' was not found in Coralysis reference colData."), call. = FALSE)
  }
  query_sce <- .sn_get_seurat_logcounts_sce(
    object = object,
    assay = query_assay,
    layer = query_layer,
    verbose = verbose
  )
  query_sce <- Coralysis::PrepareData(object = query_sce)

  args <- .sn_merge_control_args(
    defaults = list(
      ref = ref_sce,
      query = query_sce,
      ref.label = label_col
    ),
    control = transfer_control
  )
  mapped <- do.call(.sn_coralysis_reference_mapping_backend, args)
  mapped_coldata <- as.data.frame(SummarizedExperiment::colData(mapped))
  mapped_coldata <- mapped_coldata[colnames(object), , drop = FALSE]

  label_source <- if ("pruned_coral_labels" %in% colnames(mapped_coldata)) {
    "pruned_coral_labels"
  } else {
    "coral_labels"
  }
  if (!label_source %in% colnames(mapped_coldata)) {
    stop("Coralysis::ReferenceMapping() did not return `coral_labels`.", call. = FALSE)
  }

  metadata <- data.frame(row.names = colnames(object))
  metadata[[paste0(prediction_prefix, "_label")]] <- mapped_coldata[[label_source]]
  if ("coral_probability" %in% colnames(mapped_coldata)) {
    metadata[[paste0(prediction_prefix, "_score")]] <- mapped_coldata$coral_probability
  }
  if ("coral_labels" %in% colnames(mapped_coldata)) {
    metadata[[paste0(prediction_prefix, "_raw_label")]] <- mapped_coldata$coral_labels
  }

  object <- Seurat::AddMetaData(object, metadata = metadata)
  object@misc$label_transfer[[prediction_prefix]] <- list(
    method = "coralysis",
    label_col = label_col,
    prediction_columns = colnames(metadata),
    transfer_control = transfer_control
  )

  if (isTRUE(return_anchors)) {
    return(list(query = object, mapping = mapped))
  }
  object
}

#' Transfer labels from a Seurat reference to a query object
#'
#' \code{sn_transfer_labels()} is a Shennong wrapper for reference mapping. It
#' keeps the common path compact: transfer one metadata label, add the predicted
#' label_by and confidence score back to the query, and store a small provenance
#' record in \code{query@misc$label_transfer}. The default \code{method =
#' "seurat"} wraps Seurat's \code{FindTransferAnchors()} and
#' \code{TransferData()} workflow. \code{method = "coralysis"} projects the
#' query onto a native Coralysis-trained reference with
#' \code{Coralysis::ReferenceMapping()}. \code{method = "scanvi"} and
#' \code{method = "scarches"} use the managed scVI-family pixi backend to train
#' a semi-supervised scANVI model with reference labels and query cells marked
#' as unlabeled, then import the predicted query labels.
#'
#' @param object A Seurat query object to annotate. This argument comes first
#'   so the function can be used in pipes.
#' @param reference A labeled Seurat reference object.
#' @param label_by Metadata column in \code{reference} to transfer.
#' @param method Label-transfer backend. \code{"seurat"} uses Seurat anchors;
#'   \code{"coralysis"} uses native \code{Coralysis::ReferenceMapping()} and
#'   requires a trained Coralysis reference stored under
#'   \code{reference@misc$coralysis}; \code{"scanvi"} and \code{"scarches"} use
#'   the scVI-family pixi backend.
#' @param prediction_prefix Prefix for metadata columns added to
#'   \code{query}. Defaults to \code{paste0(label_by, "_transfer")}.
#' @param normalization_method Normalization method passed to
#'   \code{Seurat::FindTransferAnchors()}.
#' @param reference_assay,query_assay Assays passed to
#'   \code{Seurat::FindTransferAnchors()} and \code{Seurat::TransferData()}.
#'   For \code{method = "coralysis"}, \code{query_assay} controls the assay
#'   converted to query \code{logcounts}; the reference assay is ignored when a
#'   stored Coralysis reference is available.
#' @param reference_layer,query_layer Layers used as log-normalized expression
#'   for \code{method = "coralysis"}. The query defaults to the Seurat
#'   \code{"data"} layer and is normalized first if that layer is absent.
#' @param reduction Dimensional reduction strategy passed to
#'   \code{Seurat::FindTransferAnchors()}.
#' @param reference_reduction Optional reference reduction passed to
#'   \code{Seurat::FindTransferAnchors()}.
#' @param features Optional features used to find transfer anchors.
#' @param dims Dimensions used for anchor scoring and label_by transfer.
#' @param npcs Number of PCs used by \code{Seurat::FindTransferAnchors()}.
#' @param k_anchor,k_filter,k_score,k_weight Seurat anchor/weighting
#'   parameters.
#' @param store_prediction_scores If \code{TRUE}, also store per-label
#'   prediction scores as query metadata columns.
#' @param return_anchors If \code{TRUE}, return a list containing the annotated
#'   query and backend artifacts. For Coralysis, the artifact is the mapped
#'   SingleCellExperiment.
#' @param transfer_control Optional backend-specific list. For
#'   \code{method = "coralysis"}, values are forwarded to
#'   \code{Coralysis::ReferenceMapping()}. For \code{method = "scanvi"} or
#'   \code{"scarches"}, common values include \code{batch_by},
#'   \code{runtime_dir}, \code{pixi_project}, \code{max_epochs},
#'   \code{scanvi_max_epochs}, \code{accelerator}, \code{mirror}, and
#'   \code{install_pixi}.
#' @param verbose Whether to print Seurat progress messages.
#' @param ... Additional arguments passed to \code{Seurat::FindTransferAnchors()}.
#'
#' @return A Seurat query object with transferred labels, or a list when
#'   \code{return_anchors = TRUE}.
#'
#' @examples
#' \dontrun{
#' query <- sn_transfer_labels(
#'   object = query,
#'   reference = reference,
#'   label_by = "cell_type",
#'   dims = 1:30
#' )
#' }
#' @export
sn_transfer_labels <- function(object = NULL,
                               reference,
                               label_by = NULL,
                               method = c("seurat", "coralysis", "scanvi", "scarches"),
                               prediction_prefix = NULL,
                               normalization_method = "LogNormalize",
                               reference_assay = NULL,
                               query_assay = NULL,
                               reference_layer = "data",
                               query_layer = "data",
                               reduction = "pcaproject",
                               reference_reduction = NULL,
                               features = NULL,
                               dims = 1:30,
                               npcs = 30,
                               k_anchor = 5,
                               k_filter = NA,
                               k_score = 30,
                               k_weight = 50,
                               store_prediction_scores = FALSE,
                               return_anchors = FALSE,
                               transfer_control = list(),
                               verbose = TRUE,
                               ...) {
  check_installed("Seurat")
  method <- match.arg(method)
  if (is.null(object)) {
    stop("`object` must be supplied as the query Seurat object.", call. = FALSE)
  }

  if (!inherits(reference, "Seurat") && !(identical(method, "coralysis") && inherits(reference, "SingleCellExperiment"))) {
    stop("`reference` must be a Seurat object.", call. = FALSE)
  }
  if (!inherits(object, "Seurat")) {
    stop("`object` must be a Seurat query object.", call. = FALSE)
  }
  if (!is.character(label_by) || length(label_by) != 1L || !nzchar(label_by)) {
    stop("`label_by` must be a non-empty metadata column name.", call. = FALSE)
  }
  if (method == "coralysis") {
    prediction_prefix <- prediction_prefix %||% paste0(label_by, "_coralysis")
    return(.sn_transfer_labels_coralysis(
      object = object,
      reference = reference,
      label_col = label_by,
      prediction_prefix = prediction_prefix,
      reference_assay = reference_assay,
      query_assay = query_assay,
      reference_layer = reference_layer,
      query_layer = query_layer,
      transfer_control = transfer_control,
      return_anchors = return_anchors,
      verbose = verbose
    ))
  }
  if (method %in% c("scanvi", "scarches")) {
    prediction_prefix <- prediction_prefix %||% paste0(label_by, "_", method)
    return(.sn_transfer_labels_scanvi(
      object = object,
      reference = reference,
      label_by = label_by,
      prediction_prefix = prediction_prefix,
      assay = query_assay %||% reference_assay,
      batch_by = transfer_control$batch_by %||% NULL,
      features = features,
      transfer_control = transfer_control,
      return_anchors = return_anchors,
      verbose = verbose,
      method_name = method
    ))
  }

  if (!label_by %in% colnames(reference[[]])) {
    stop(glue("`label_by` column '{label_by}' was not found in `reference` metadata."), call. = FALSE)
  }

  ref_labels <- reference[[label_by, drop = TRUE]]
  names(ref_labels) <- colnames(reference)
  prediction_prefix <- prediction_prefix %||% paste0(label_by, "_transfer")

  anchors <- .sn_find_transfer_anchors_backend(
    reference = reference,
    query = object,
    normalization.method = normalization_method,
    reference.assay = reference_assay,
    query.assay = query_assay,
    reduction = reduction,
    reference.reduction = reference_reduction,
    features = features,
    npcs = npcs,
    dims = dims,
    k.anchor = k_anchor,
    k.filter = k_filter,
    k.score = k_score,
    verbose = verbose,
    ...
  )

  predictions <- .sn_transfer_data_backend(
    anchorset = anchors,
    refdata = ref_labels,
    reference = reference,
    query = object,
    query.assay = query_assay,
    weight.reduction = reduction,
    dims = dims,
    k.weight = k_weight,
    verbose = verbose
  )
  predictions <- as.data.frame(predictions)
  if (!all(colnames(object) %in% rownames(predictions)) && nrow(predictions) == ncol(object)) {
    rownames(predictions) <- colnames(object)
  }
  predictions <- predictions[colnames(object), , drop = FALSE]
  if (!"predicted.id" %in% colnames(predictions)) {
    stop("Seurat::TransferData() did not return a `predicted.id` column.", call. = FALSE)
  }

  metadata <- data.frame(row.names = colnames(object))
  metadata[[paste0(prediction_prefix, "_label")]] <- predictions$predicted.id
  if ("prediction.score.max" %in% colnames(predictions)) {
    metadata[[paste0(prediction_prefix, "_score")]] <- predictions$prediction.score.max
  }
  if (isTRUE(store_prediction_scores)) {
    score_cols <- grep("^prediction\\.score\\.", colnames(predictions), value = TRUE)
    score_cols <- setdiff(score_cols, "prediction.score.max")
    for (score_col in score_cols) {
      suffix <- sub("^prediction\\.score\\.", "", score_col)
      metadata[[paste0(prediction_prefix, "_score_", .sn_metadata_suffix(suffix))]] <- predictions[[score_col]]
    }
  }

  object <- Seurat::AddMetaData(object, metadata = metadata)
  object@misc$label_transfer[[prediction_prefix]] <- list(
    label_by = label_by,
    label_col = label_by,
    normalization_method = normalization_method,
    reduction = reduction,
    dims = dims,
    features = features,
    prediction_columns = colnames(metadata)
  )

  if (isTRUE(return_anchors)) {
    return(list(query = object, anchors = anchors, predictions = predictions))
  }
  object
}

.sn_run_scdesign3_backend <- function(...) {
  scDesign3::scdesign3(...)
}

.sn_seurat_to_sce_for_scdesign3 <- function(object,
                                            assay = "RNA",
                                            layer = "counts") {
  counts <- .sn_get_seurat_layer_data(object = object, assay = assay, layer = layer)
  SingleCellExperiment::SingleCellExperiment(
    list(counts = .sn_as_sparse_matrix(counts)),
    colData = object[[]]
  )
}

.sn_validate_scdesign3_columns <- function(sce,
                                           celltype = NULL,
                                           pseudotime = NULL,
                                           spatial = NULL,
                                           other_covariates = NULL) {
  available <- colnames(SummarizedExperiment::colData(sce))
  requested <- unique(c(celltype, pseudotime, spatial, other_covariates))
  requested <- requested[!is.na(requested) & nzchar(requested)]
  missing <- setdiff(requested, available)
  if (length(missing) > 0L) {
    stop(glue("Missing scDesign3 covariate column(s): {paste(missing, collapse = ', ')}."), call. = FALSE)
  }
  invisible(TRUE)
}

.sn_default_scdesign3_formula <- function(celltype = NULL,
                                          pseudotime = NULL,
                                          spatial = NULL,
                                          other_covariates = NULL,
                                          k = 10) {
  if (!is.null(pseudotime) && length(pseudotime) > 0L) {
    return(glue("s({pseudotime[[1]]}, bs = 'cr', k = {k})"))
  }
  if (!is.null(spatial) && length(spatial) == 2L) {
    return(glue("s({spatial[[1]]}, {spatial[[2]]}, bs = 'gp', k = {k})"))
  }
  terms <- unique(c(celltype, other_covariates))
  terms <- terms[!is.na(terms) & nzchar(terms)]
  if (length(terms) == 0L) {
    return("1")
  }
  paste(terms, collapse = " + ")
}

.sn_default_scdesign3_corr_formula <- function(celltype = NULL,
                                               pseudotime = NULL,
                                               spatial = NULL,
                                               other_covariates = NULL) {
  terms <- unique(c(pseudotime, spatial, celltype, other_covariates))
  terms <- terms[!is.na(terms) & nzchar(terms)]
  if (length(terms) == 0L) {
    return("1")
  }
  paste(terms, collapse = " + ")
}

.sn_extract_scdesign3_counts <- function(result, expected_genes = NULL) {
  if (!is.list(result) || !"new_count" %in% names(result)) {
    stop("scDesign3 did not return a `new_count` result.", call. = FALSE)
  }
  counts <- result$new_count
  if (is.list(counts)) {
    counts <- counts[[1]]
  }
  counts <- .sn_as_sparse_matrix(counts)
  if (is.null(rownames(counts)) && !is.null(expected_genes) && length(expected_genes) == nrow(counts)) {
    rownames(counts) <- expected_genes
  }
  if (is.null(colnames(counts))) {
    colnames(counts) <- paste0("sim_cell_", seq_len(ncol(counts)))
  }
  counts
}

#' Simulate single-cell counts with scDesign3
#'
#' \code{sn_simulate()} provides a method-based simulation entry point.
#' Currently \code{method = "scdesign3"} delegates to
#' \code{sn_simulate_scdesign3()}.
#'
#' @param object A Seurat or SingleCellExperiment object.
#' @param method Simulation backend. Currently supports \code{"scdesign3"}.
#' @param ... Additional arguments passed to the selected backend.
#'
#' @return Simulated data in the backend's requested format.
#'
#' @examples
#' \dontrun{
#' sim <- sn_simulate(seurat_obj, method = "scdesign3", celltype = "cell_type")
#' }
#' @export
sn_simulate <- function(object,
                        method = c("scdesign3"),
                        ...) {
  method <- match.arg(method)
  switch(
    method,
    scdesign3 = sn_simulate_scdesign3(object = object, ...)
  )
}

#' Simulate single-cell counts with scDesign3
#'
#' \code{sn_simulate_scdesign3()} prepares a Seurat or SingleCellExperiment
#' object for \code{scDesign3::scdesign3()}, chooses simple default formulas
#' from the supplied covariates, and returns simulated counts as a Seurat
#' object, SingleCellExperiment, sparse count matrix, or raw scDesign3 result.
#' Prefer \code{sn_simulate(method = "scdesign3")} for new code.
#'
#' @param object A Seurat or SingleCellExperiment object.
#' @param celltype Column in \code{colData(object)} or Seurat metadata used as
#'   the cell-type covariate. If \code{NULL}, \code{"seurat_clusters"} is used
#'   when available.
#' @param pseudotime Optional pseudotime covariate column name.
#' @param spatial Optional length-two vector naming spatial coordinate columns.
#' @param other_covariates Optional additional covariate columns.
#' @param ncell Number of cells to simulate. Defaults to the input cell count.
#' @param mu_formula,sigma_formula,corr_formula scDesign3 model formulas. When
#'   \code{mu_formula} or \code{corr_formula} is \code{NULL}, Shennong builds a
#'   simple default from \code{pseudotime}, \code{spatial}, \code{celltype}, and
#'   \code{other_covariates}.
#' @param family_use Marginal distribution passed to scDesign3.
#' @param n_cores Number of cores passed to scDesign3.
#' @param assay,layer Assay/layer used when \code{object} is a Seurat object.
#' @param assay_use Assay name used when \code{object} is already a
#'   SingleCellExperiment. Defaults to \code{"counts"}.
#' @param return One of \code{"seurat"}, \code{"sce"}, \code{"counts"}, or
#'   \code{"result"}.
#' @param project Project name for returned Seurat objects.
#' @param combine_original If \code{TRUE}, return a merged Seurat object
#'   containing original and simulated cells with a \code{simulation_source}
#'   metadata column. Only applies when \code{return = "seurat"} and
#'   \code{object} is a Seurat object.
#' @param seed Optional random seed.
#' @param ... Additional arguments passed to \code{scDesign3::scdesign3()}.
#'
#' @return Simulated data in the requested format.
#'
#' @examples
#' \dontrun{
#' sim <- sn_simulate_scdesign3(
#'   object = seurat_obj,
#'   celltype = "cell_type",
#'   ncell = 1000,
#'   n_cores = 4
#' )
#' }
#' @export
sn_simulate_scdesign3 <- function(object,
                                  celltype = NULL,
                                  pseudotime = NULL,
                                  spatial = NULL,
                                  other_covariates = NULL,
                                  ncell = NULL,
                                  mu_formula = NULL,
                                  sigma_formula = "1",
                                  corr_formula = NULL,
                                  family_use = "nb",
                                  n_cores = 2,
                                  assay = "RNA",
                                  layer = "counts",
                                  assay_use = "counts",
                                  return = c("seurat", "sce", "counts", "result"),
                                  project = "scdesign3",
                                  combine_original = FALSE,
                                  seed = 717,
                                  ...) {
  check_installed("scDesign3", reason = "to simulate single-cell data with scDesign3.")
  check_installed(c("SingleCellExperiment", "SummarizedExperiment"))

  return <- match.arg(return)
  if (inherits(object, "Seurat")) {
    sce <- .sn_seurat_to_sce_for_scdesign3(object = object, assay = assay, layer = layer)
  } else if (inherits(object, "SingleCellExperiment")) {
    sce <- object
  } else {
    stop("`object` must be a Seurat or SingleCellExperiment object.", call. = FALSE)
  }

  coldata <- SummarizedExperiment::colData(sce)
  if (is.null(celltype) && "seurat_clusters" %in% colnames(coldata)) {
    celltype <- "seurat_clusters"
  }
  .sn_validate_scdesign3_columns(
    sce = sce,
    celltype = celltype,
    pseudotime = pseudotime,
    spatial = spatial,
    other_covariates = other_covariates
  )

  ncell <- ncell %||% ncol(sce)
  formula_k <- max(3L, min(10L, floor(ncell / 5)))
  mu_formula <- mu_formula %||% .sn_default_scdesign3_formula(
    celltype = celltype,
    pseudotime = pseudotime,
    spatial = spatial,
    other_covariates = other_covariates,
    k = formula_k
  )
  corr_formula <- corr_formula %||% .sn_default_scdesign3_corr_formula(
    celltype = celltype,
    pseudotime = pseudotime,
    spatial = spatial,
    other_covariates = other_covariates
  )

  if (!is.null(seed)) {
    set.seed(seed)
  }
  result <- .sn_run_scdesign3_backend(
    sce = sce,
    assay_use = assay_use,
    celltype = celltype,
    pseudotime = pseudotime,
    spatial = spatial,
    other_covariates = other_covariates,
    ncell = ncell,
    mu_formula = mu_formula,
    sigma_formula = sigma_formula,
    family_use = family_use,
    n_cores = n_cores,
    corr_formula = corr_formula,
    ...
  )
  result$shennong <- list(
    celltype = celltype,
    pseudotime = pseudotime,
    spatial = spatial,
    other_covariates = other_covariates,
    mu_formula = mu_formula,
    sigma_formula = sigma_formula,
    corr_formula = corr_formula,
    family_use = family_use,
    ncell = ncell
  )

  if (return == "result") {
    return(result)
  }

  sim_counts <- .sn_extract_scdesign3_counts(
    result = result,
    expected_genes = rownames(sce)
  )
  sim_meta <- result$new_covariate
  if (is.null(sim_meta)) {
    sim_meta <- data.frame(row.names = colnames(sim_counts))
  } else {
    sim_meta <- as.data.frame(sim_meta)
    if (nrow(sim_meta) == ncol(sim_counts)) {
      rownames(sim_meta) <- colnames(sim_counts)
    }
  }

  if (return == "counts") {
    return(sim_counts)
  }
  if (return == "sce") {
    return(SingleCellExperiment::SingleCellExperiment(
      list(counts = sim_counts),
      colData = sim_meta
    ))
  }

  check_installed("Seurat")
  sim_object <- Seurat::CreateSeuratObject(
    counts = sim_counts,
    meta.data = sim_meta,
    project = project
  )
  sim_object@misc$scdesign3 <- result$shennong
  if (isTRUE(combine_original) && inherits(object, "Seurat")) {
    original <- object
    original$simulation_source <- "original"
    sim_object$simulation_source <- "simulated"
    combined <- merge(
      x = original,
      y = sim_object,
      add.cell.ids = c("original", "simulated"),
      project = project
    )
    combined@misc$scdesign3 <- result$shennong
    return(combined)
  }
  sim_object
}


#' Run CellTypist for automated cell type annotation
#'
#' @param x A Seurat object or a path to a count matrix / AnnData file that CellTypist can consume.
#' @param celltypist Path to the `celltypist` binary. Defaults to "/opt/mambaforge/envs/scverse/bin/celltypist".
#' @param model Model used for predictions. Defaults to "Immune_All_Low.pkl".
#' @param outdir Directory to store the output files. If NULL, use a temporary directory.
#' @param prefix Prefix for the output files. By default, use the model name plus a dot.
#' @param mode Choose the cell type with the largest score/probability (`"best_match"`) or enable multi-label classification (`"prob_match"`).
#' @param p_thres Probability threshold for the multi-label classification. Ignored if `mode = "best_match"`.
#' @param majority_voting Logical. Whether to refine labels using majority voting after over-clustering.
#' @param over_clustering Input file or a string key specifying an existing metadata column in the AnnData object, or "auto".
#' @param min_prop For the dominant cell type within a subcluster, the minimum proportion of cells required to name the subcluster by this cell type.
#' @param transpose_input Logical. If `TRUE`, add the `--transpose-input` argument when calling `celltypist`.
#' @param gene_file If the provided input is in the `mtx` format, path to the file storing gene information. Otherwise ignored.
#' @param cell_file If the provided input is in the `mtx` format, path to the file storing cell information. Otherwise ignored.
#' @param assay Assay used when exporting Seurat counts to CellTypist. Defaults to \code{"RNA"}.
#' @param layer Layer used as the input count matrix. Defaults to \code{"counts"}.
#' @param xlsx Logical. If `TRUE`, merge output tables into a single Excel (.xlsx). Defaults to `FALSE`.
#' @param plot_results Logical. If `TRUE`, plot the prediction results. Defaults to `FALSE`.
#' @param quiet Logical. If `TRUE`, hide the banner and config info from `celltypist`. Defaults to `FALSE`.
#'
#' @return When \code{x} is a Seurat object, a Seurat object with prediction
#'   columns added to metadata. When \code{x} is a path, the CellTypist
#'   prediction table is returned.
#'
#' @examples
#' \dontrun{
#' data("pbmc_small", package = "Shennong")
#' pbmc <- sn_run_cluster(pbmc, normalization_method = "seurat", verbose = FALSE)
#' pbmc <- sn_run_celltypist(pbmc, model = "Immune_All_Low.pkl")
#' head(colnames(pbmc[[]]))
#' }
#' @export
sn_run_celltypist <- function(x,
                              celltypist = NULL,
                              model = "Immune_All_Low.pkl",
                              outdir = NULL,
                              prefix = NULL,
                              mode = c("best_match", "prob_match"),
                              p_thres = 0.5,
                              majority_voting = TRUE,
                              over_clustering = "auto",
                              min_prop = 0,
                              transpose_input = TRUE,
                              gene_file = NULL,
                              cell_file = NULL,
                              assay = "RNA",
                              layer = "counts",
                              xlsx = FALSE,
                              plot_results = FALSE,
                              quiet = FALSE) {
  check_installed(c("Seurat", "data.table", "logger", "glue"),
    reason = "to run CellTypist analysis."
  )

  mode <- match.arg(mode)
  celltypist <- celltypist %||% getOption("shennong.celltypist_path", Sys.which("celltypist"))
  sn_check_file(celltypist)

  .sn_log_info("Starting CellTypist analysis with model = {model}.")
  tictoc::tic("Total CellTypist runtime")

  if (is_null(outdir)) {
    outdir <- tempfile("celltypist_")
    dir.create(outdir, recursive = TRUE)
    .sn_log_info("Using a temporary output directory: {outdir}.")
  } else {
    outdir <- sn_set_path(outdir)
    .sn_log_info("Using the user-specified output directory: {outdir}.")
  }

  on.exit(
    {
      if (dir.exists(outdir) && grepl("^celltypist_", basename(outdir))) {
        unlink(outdir, recursive = TRUE)
        log_debug("Cleaned temporary directory: {outdir}")
      }
    },
    add = TRUE
  )

  if (inherits(x, "Seurat")) {
    x_is_seurat <- TRUE
    .sn_log_info("Converting the Seurat object to CellTypist input format.")
    input_data <- file.path(outdir, "counts.csv")

    counts <- tryCatch(
      .sn_get_seurat_layer_data(object = x, assay = assay, layer = layer),
      error = function(e) {
        .sn_log_error("Failed to extract the counts layer: {e$message}.")
        stop("Counts layer extraction failed")
      }
    )

    data.table::fwrite(
      x = as.data.frame(counts),
      file = input_data,
      quote = FALSE,
      row.names = TRUE,
      showProgress = FALSE
    )
    log_debug("Count matrix written to {input_data} ({file.size(input_data)} bytes)")
  } else {
    x_is_seurat <- FALSE
    input_data <- x
    .sn_log_info("Using the precomputed input matrix: {input_data}.")
  }

  over_clustering_path <- NULL

  if (over_clustering != "auto") {
    if (isTRUE(x_is_seurat)) {
      if (over_clustering %in% colnames(x@meta.data)) {
        .sn_log_info("Using the existing clustering column: {over_clustering}.")
        over_clustering_path <- file.path(outdir, "over_clustering.txt")
        writeLines(
          as.character(x[[over_clustering, drop = TRUE]]),
          over_clustering_path
        )
      } else if (file.exists(over_clustering)) {
        over_clustering_path <- over_clustering
      } else {
        stop(
          "`over_clustering` must be a metadata column or existing file when `x` is a Seurat object.",
          call. = FALSE
        )
      }
    } else {
      over_clustering_path <- over_clustering
    }
  } else if (over_clustering == "auto") {
    .sn_log_warn("Automatic over-clustering is not implemented yet.")
  }

  model_name <- tools::file_path_sans_ext(basename(model))
  prefix <- prefix %||% glue("{model_name}.")
  predicted_labels_path <- file.path(outdir, glue("{prefix}predicted_labels.csv"))

  cmd_args <- c(
    "--indata", shQuote(input_data),
    "--model", shQuote(model),
    "--mode", mode,
    "--outdir", shQuote(outdir),
    "--prefix", prefix
  )

  add_arg <- function(args, flag, value, condition = TRUE) {
    if (condition && !is_null(value)) c(args, flag, shQuote(value)) else args
  }

  cmd_args <- add_arg(cmd_args, "--gene-file", gene_file, !is_null(gene_file))
  cmd_args <- add_arg(cmd_args, "--cell-file", cell_file, !is_null(cell_file))
  cmd_args <- add_arg(cmd_args, "--p-thres", p_thres, mode == "prob_match")
  cmd_args <- add_arg(cmd_args, "--over-clustering", over_clustering_path, !is_null(over_clustering_path))
  cmd_args <- add_arg(cmd_args, "--min-prop", min_prop, majority_voting)

  if (transpose_input) cmd_args <- c(cmd_args, "--transpose-input")
  if (xlsx) cmd_args <- c(cmd_args, "--xlsx")
  if (plot_results) cmd_args <- c(cmd_args, "--plot-results")
  if (quiet) cmd_args <- c(cmd_args, "--quiet")
  if (majority_voting) cmd_args <- c(cmd_args, "--majority-voting")
  .sn_log_info("Executing CellTypist with command:\n{celltypist} {paste(cmd_args, collapse = ' ')}")

  exit_code <- system2(
    command = celltypist,
    args = cmd_args,
    stdout = if (quiet) FALSE else "",
    stderr = if (quiet) FALSE else ""
  )

  if (exit_code != 0) {
    .sn_log_error("CellTypist failed with exit code {exit_code}.")
    stop("CellTypist execution failed. Check logs for details.")
  }
  log_debug("CellTypist output written to {outdir}")

  if (!file.exists(predicted_labels_path)) {
    .sn_log_error("Prediction output is missing: {predicted_labels_path}.")
    stop("CellTypist did not generate expected output files")
  }

  predicted_labels <- utils::read.csv(predicted_labels_path, row.names = 1)
  if (majority_voting) {
    colnames(predicted_labels) <- c(
      paste0(gsub("\\.pkl$", "", model_name), "_predicted_labels"),
      paste0(gsub("\\.pkl$", "", model_name), "_over_clustering"),
      paste0(gsub("\\.pkl$", "", model_name), "_majority_voting")
    )
  } else {
    colnames(predicted_labels) <- c(
      paste0(gsub("\\.pkl$", "", model_name), "_predicted_labels")
    )
  }

  if (!isTRUE(x_is_seurat)) {
    tictoc::toc()
    .sn_log_info("CellTypist analysis completed successfully.")
    return(tibble::as_tibble(predicted_labels, rownames = "cell"))
  }

  .sn_log_info("Adding {ncol(predicted_labels)} metadata columns to the Seurat object.")
  x <- SeuratObject::AddMetaData(x, metadata = predicted_labels)

  tictoc::toc()
  .sn_log_info("CellTypist analysis completed successfully.")

  .sn_log_seurat_command(object = x, assay = assay, name = "sn_run_celltypist")
}
