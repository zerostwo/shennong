# Integration backends for `sn_run_cluster()`.
#
# Extracted from analysis_clustering.R: batch/Harmony/MMOCHI/Coralysis/
# Seurat-layer/scVI/scPoli/BBKNN adapters plus their pixi runtime helpers,
# the scArches/scPoli object-method wrappers, and deprecated shims.

.sn_run_batch_integration <- function(object,
                                      method,
                                      batch,
                                      reduction,
                                      features,
                                      assay,
                                      layer,
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
      layer = layer,
      integration_control = integration_control,
      verbose = verbose
    ),
    totalvi = .sn_run_scvi_integration(
      object = object,
      method = method,
      batch = batch,
      features = features,
      assay = assay,
      layer = layer,
      integration_control = integration_control,
      verbose = verbose
    ),
    scanvi = .sn_run_scvi_integration(
      object = object,
      method = method,
      batch = batch,
      features = features,
      assay = assay,
      layer = layer,
      integration_control = integration_control,
      verbose = verbose
    ),
    scpoli = .sn_run_scpoli_integration(
      object = object,
      batch = batch,
      features = features,
      assay = assay,
      layer = layer,
      integration_control = integration_control,
      verbose = verbose
    ),
    bbknn = .sn_run_bbknn_integration(
      object = object,
      batch = batch,
      reduction = reduction,
      assay = assay,
      layer = layer,
      dims = dims,
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
                                   protein_features = NULL,
                                   max_dense_gb = 2) {
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  protein <- .sn_get_seurat_layer_data(object = object, assay = protein_assay, layer = protein_layer)
  feature_set <- protein_features %||% rownames(protein)
  feature_set <- intersect(feature_set, rownames(protein))
  if (length(feature_set) < 2L) {
    stop("MMoCHi integration requires at least two ADT/protein features present in the object.", call. = FALSE)
  }
  protein <- protein[feature_set, colnames(object), drop = FALSE]
  .sn_assert_dense_materialization_budget(
    protein,
    max_dense_gb = max_dense_gb,
    name = "MMoCHi protein export",
    peak_copies = 4
  )
  protein_df <- as.data.frame(t(as.matrix(protein)), check.names = FALSE)
  protein_df <- cbind(cell_id = rownames(protein_df), protein_df)

  if (is.null(batch) || !nzchar(batch) || !batch %in% colnames(object[[]])) {
    stop("MMoCHi export requires `batch` to name a metadata column.", call. = FALSE)
  }
  metadata <- object[[]][colnames(object), batch, drop = FALSE]
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
                                      store_corrected_layer = TRUE,
                                      max_artifact_import_gb = 0.5) {
  corrected_path <- file.path(output_dir, "landmark_protein.csv")
  backend_manifest <- .sn_read_integration_backend_manifest(
    output_dir = output_dir,
    method = "mmochi",
    n_cells = ncol(object)
  )
  if (!file.exists(corrected_path)) {
    stop("MMoCHi output is missing `landmark_protein.csv`: ", corrected_path, call. = FALSE)
  }

  .sn_assert_python_artifact_budget(
    corrected_path,
    max_artifact_import_gb,
    "MMoCHi corrected protein output"
  )
  corrected <- utils::read.csv(corrected_path, row.names = 1, check.names = FALSE)
  if (anyDuplicated(rownames(corrected)) || !setequal(colnames(object), rownames(corrected))) {
    stop("MMoCHi landmark output cell identifiers do not exactly match the input object.", call. = FALSE)
  }
  if (anyDuplicated(colnames(corrected)) || !setequal(protein_features, colnames(corrected))) {
    stop("MMoCHi landmark output features do not exactly match the requested ADT/protein features.", call. = FALSE)
  }
  corrected <- as.matrix(corrected[colnames(object), protein_features, drop = FALSE])
  storage.mode(corrected) <- "numeric"
  if (length(corrected) == 0L || any(!is.finite(corrected))) {
    stop("MMoCHi landmark output must be a non-empty finite numeric matrix.", call. = FALSE)
  }
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
  .sn_with_integration_python_run(
    method = "mmochi",
    integration_control = integration_control,
    code = function(control) {
      .sn_run_mmochi_integration_impl(
        object = object,
        batch = batch,
        assay = assay,
        protein_assay = protein_assay,
        protein_features = protein_features,
        dims = dims,
        npcs = npcs,
        integration_control = control,
        verbose = verbose
      )
    }
  )
}

.sn_run_mmochi_integration_impl <- function(object,
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
  pixi_paths <- sn_get_pixi_paths(environment = "mmochi", runtime_dir = runtime_dir)
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
  run_dir <- integration_control$run_dir %||% .sn_default_python_run_dir(method = "mmochi", runtime_dir = runtime_dir)
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
    protein_features = protein_features,
    max_dense_gb = integration_control$max_dense_gb %||% 2
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

  backend_run <- .sn_execute_scvi_pixi(
    pixi = integration_control$pixi %||% NULL,
    manifest_path = manifest_path,
    script = script,
    input_dir = input$input_dir,
    output_dir = output_dir,
    config_path = config_path,
    environment = integration_control$environment %||% "default",
    pixi_home = pixi_home,
    install_pixi = integration_control$install_pixi %||% TRUE,
    pixi_version = integration_control$pixi_version %||% "0.69.0",
    pixi_download_url = integration_control$pixi_download_url %||% NULL,
    pixi_sha256 = integration_control$pixi_sha256 %||% NULL,
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
    store_corrected_layer = integration_control$store_corrected_layer %||% TRUE,
    max_artifact_import_gb = integration_control$max_artifact_import_gb %||% 0.5
  )

  if (single_sample_batch && !isTRUE(integration_control$keep_single_sample_batch)) {
    result$object@meta.data[[backend_batch]] <- NULL
  }
  result$object@misc$integration$backend_performance <- backend_run$performance
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
  .sn_with_explicit_acceleration_or_disabled(
    .sn_run_coralysis_integration_impl(
      object = object,
      batch = batch,
      assay = assay,
      features = features,
      dims = dims,
      npcs = npcs,
      integration_control = integration_control,
      verbose = verbose
    ),
    patches = "coralysis"
  )
}

.sn_run_coralysis_integration_impl <- function(object,
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
  expr <- .sn_as_sparse_matrix(expr)
  if (!inherits(expr, "dgCMatrix")) {
    expr <- methods::as(expr, "dgCMatrix")
  }
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
        threads = .sn_coralysis_default_threads(),
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

.sn_coralysis_default_threads <- function() 1L

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
  integration_fun <- .sn_with_default_seurat_acceleration(
    switch(
      method,
      seurat_cca = Seurat::CCAIntegration,
      seurat_rpca = Seurat::RPCAIntegration
    ),
    object = object,
    assay = assay
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
  resolved_new_reduction <- args$new.reduction %||% new_reduction
  if (!is.character(resolved_new_reduction) ||
      length(resolved_new_reduction) != 1L ||
      is.na(resolved_new_reduction) ||
      !nzchar(resolved_new_reduction)) {
    stop("`integration_control$new.reduction` must be one non-empty reduction name.", call. = FALSE)
  }
  args$new.reduction <- resolved_new_reduction
  if (verbose) .sn_log_info("[sn_run_cluster] Running Seurat layer integration with method = {method}.")
  object <- .sn_with_acceleration_disabled(
    .sn_call_with_symbolic_object(
      fun_call = quote(Seurat::IntegrateLayers),
      object = object,
      args = args
    )
  )
  object[[assay]] <- .sn_with_acceleration_disabled(
    SeuratObject::JoinLayers(object[[assay]])
  )
  SeuratObject::DefaultAssay(object = object) <- old_default_assay
  if (!resolved_new_reduction %in% names(object@reductions)) {
    stop(
      "Seurat layer integration did not create the requested reduction `",
      resolved_new_reduction,
      "`.",
      call. = FALSE
    )
  }
  object@misc$integration <- list(
    method = method,
    batch_by = batch,
    reduction = resolved_new_reduction,
    input_reduction = args$orig.reduction %||% reduction
  )

  list(object = object, reduction = resolved_new_reduction)
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
  manifest_path <- normalizePath(manifest_path, winslash = "/", mustWork = TRUE)
  .sn_prepare_pixi_lock(
    environment = environment,
    manifest_path = manifest_path,
    platforms = unique(as.character(platforms %||% .sn_current_pixi_platform())),
    cuda_version = cuda_version,
    overwrite = overwrite
  )
  manifest_path
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
                                 layer = "counts",
                                 batch,
                                 labels_key = NULL,
                                 protein_assay = NULL,
                                 protein_layer = "counts",
                                 protein_features = NULL) {
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)

  if (!.sn_name_declares_count_scale(layer)) {
    stop(
      "scVI-family integration requires a raw/count-like RNA layer; `", layer,
      "` is not declared count-scale.",
      call. = FALSE
    )
  }
  counts <- .sn_get_seurat_layer_data(object = object, assay = assay, layer = layer)
  counts <- counts[, colnames(object), drop = FALSE]
  available_features <- rownames(counts)
  if (is.null(available_features) || anyNA(available_features) ||
      any(!nzchar(available_features)) || anyDuplicated(available_features)) {
    stop("The target RNA assay layer must have unique, non-empty feature identifiers.", call. = FALSE)
  }
  requested_features <- unique(as.character(features))
  requested_features <- requested_features[!is.na(requested_features) & nzchar(requested_features)]
  feature_set <- intersect(requested_features, available_features)
  if (length(feature_set) < 2L) {
    stop("scVI integration requires at least two selected features present in the target assay.", call. = FALSE)
  }

  counts <- counts[feature_set, , drop = FALSE]
  counts <- .sn_as_sparse_matrix(counts)
  counts <- .sn_validate_scvi_raw_counts(counts, label = "RNA counts")

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
    if (!.sn_name_declares_count_scale(protein_layer)) {
      stop(
        "totalVI requires a raw/count-like ADT/protein layer; `", protein_layer,
        "` is not declared count-scale.",
        call. = FALSE
      )
    }
    protein_counts <- .sn_get_seurat_layer_data(object = object, assay = protein_assay, layer = protein_layer)
    protein_feature_set <- protein_features %||% rownames(protein_counts)
    protein_feature_set <- unique(protein_feature_set[!is.na(protein_feature_set) & nzchar(protein_feature_set)])
    protein_feature_set <- intersect(protein_feature_set, rownames(protein_counts))
    if (length(protein_feature_set) < 1L) {
      stop("totalVI integration requires at least one protein feature present in `adt_assay`.", call. = FALSE)
    }
    protein_counts <- protein_counts[protein_feature_set, colnames(counts), drop = FALSE]
    protein_counts <- .sn_as_sparse_matrix(protein_counts)
    protein_counts <- .sn_validate_scvi_raw_counts(protein_counts, label = "ADT/protein counts")
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
    layer = layer,
    protein_features = protein_feature_set
  )
}

.sn_validate_scvi_raw_counts <- function(counts, label) {
  .sn_validate_python_raw_counts(counts, label = label)
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
                                  pixi_version = "0.69.0",
                                  pixi_download_url = NULL,
                                  pixi_sha256 = NULL,
                                  verbose = TRUE,
                                  backend_label = "scVI") {
  pixi_info <- sn_ensure_pixi(
    pixi = pixi,
    install = install_pixi,
    version = pixi_version,
    pixi_home = pixi_home %||% path.expand("~/.pixi"),
    no_path_update = TRUE,
    download_url = pixi_download_url,
    sha256 = pixi_sha256,
    quiet = !isTRUE(verbose)
  )
  pixi <- pixi_info$path

  args <- c(
    "run",
    "--locked",
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

  resource_path <- file.path(output_dir, "resource-usage.txt")
  use_gnu_time <- identical(Sys.info()[["sysname"]], "Linux") &&
    file.exists("/usr/bin/time")
  command <- pixi
  command_args <- args
  if (use_gnu_time) {
    command <- "/usr/bin/time"
    command_args <- c("-v", "-o", shQuote(resource_path), shQuote(pixi), args)
  }
  started_at <- Sys.time()
  elapsed_start <- proc.time()[["elapsed"]]

  status <- tryCatch(
    system2(command = command, args = command_args, env = env, stdout = TRUE, stderr = TRUE),
    error = function(e) {
      stop(backend_label, " pixi execution failed. ", conditionMessage(e), call. = FALSE)
    }
  )
  exit_code <- attr(status, "status") %||% 0L
  if (!identical(exit_code, 0L)) {
    stop(backend_label, " pixi execution failed.\n", paste(status, collapse = "\n"), call. = FALSE)
  }
  peak_rss_mb <- NA_real_
  if (file.exists(resource_path)) {
    resource <- readLines(resource_path, warn = FALSE)
    rss_line <- grep("Maximum resident set size", resource, value = TRUE)
    if (length(rss_line) > 0L) {
      peak_rss_kb <- suppressWarnings(as.numeric(sub("^.*:[[:space:]]*", "", rss_line[[1L]])))
      if (is.finite(peak_rss_kb)) peak_rss_mb <- peak_rss_kb / 1024
    }
  }
  invisible(list(
    output = status,
    performance = list(
      elapsed_seconds = unname(proc.time()[["elapsed"]] - elapsed_start),
      backend_peak_rss_mb = peak_rss_mb,
      memory_scope = if (is.finite(peak_rss_mb)) "backend_process_tree_rss" else "unavailable",
      started_at = as.character(started_at),
      completed_at = as.character(Sys.time()),
      resource_path = if (file.exists(resource_path)) normalizePath(resource_path, winslash = "/", mustWork = TRUE) else NULL
    )
  ))
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
                                    layer,
                                    batch,
                                    features,
                                    protein_features = NULL,
                                    reduction = NULL,
                                    max_artifact_import_gb = 0.5) {
  reduction <- reduction %||% method
  latent_path <- file.path(output_dir, "latent.csv")
  metadata_path <- file.path(output_dir, "obs.csv")
  backend_manifest <- .sn_read_integration_backend_manifest(
    output_dir = output_dir,
    method = method,
    n_cells = ncol(object)
  )
  manifest_features <- suppressWarnings(as.integer(backend_manifest$n_features %||% NA_integer_))
  if (length(manifest_features) != 1L || is.na(manifest_features) ||
      manifest_features != length(features)) {
    stop("scVI manifest `n_features` does not match the exported target-assay features.", call. = FALSE)
  }
  if (identical(method, "totalvi")) {
    manifest_proteins <- suppressWarnings(as.integer(backend_manifest$n_proteins %||% NA_integer_))
    if (length(manifest_proteins) != 1L || is.na(manifest_proteins) ||
        manifest_proteins != length(protein_features)) {
      stop("totalVI manifest `n_proteins` does not match the exported ADT features.", call. = FALSE)
    }
  }

  if (!file.exists(latent_path)) {
    stop("scVI output is missing `latent.csv`: ", latent_path, call. = FALSE)
  }
  .sn_assert_python_artifact_budget(latent_path, max_artifact_import_gb, "scVI latent output")
  latent <- .sn_read_embedding_csv(latent_path, cells = colnames(object))
  colnames(latent) <- paste0(toupper(reduction), "_", seq_len(ncol(latent)))
  object[[reduction]] <- Seurat::CreateDimReducObject(
    embeddings = latent,
    key = paste0(toupper(reduction), "_"),
    assay = assay
  )

  if (file.exists(metadata_path)) {
    .sn_assert_python_artifact_budget(metadata_path, max_artifact_import_gb, "scVI metadata output")
    metadata <- utils::read.csv(metadata_path, row.names = 1, check.names = FALSE)
    .sn_validate_exact_python_ids(
      rownames(metadata), colnames(object),
      "scVI metadata output cell"
    )
    if (anyDuplicated(colnames(metadata))) {
      stop("scVI metadata output columns must be unique.", call. = FALSE)
    }
    if (ncol(metadata) > 0L) {
      object <- Seurat::AddMetaData(object = object, metadata = metadata[colnames(object), , drop = FALSE])
    }
  }

  output_h5ad <- .sn_integration_manifest_artifact(
    backend_manifest,
    field = "output_h5ad",
    output_dir = output_dir
  )
  object@misc$integration <- list(
    method = method,
    batch_by = batch,
    reduction = reduction,
    assay = assay,
    source_layer = layer,
    input_features = features,
    run_dir = normalizePath(output_dir, winslash = "/", mustWork = TRUE),
    output_h5ad = output_h5ad
  )

  list(object = object, reduction = reduction)
}

.sn_run_scvi_integration <- function(object,
                                     method,
                                     batch,
                                     features,
                                     assay,
                                     layer = "counts",
                                     integration_control = list(),
                                     verbose = TRUE) {
  .sn_with_integration_python_run(
    method = method,
    integration_control = integration_control,
    code = function(control) {
      .sn_run_scvi_integration_impl(
        object = object,
        method = method,
        batch = batch,
        features = features,
        assay = assay,
        layer = layer,
        integration_control = control,
        verbose = verbose
      )
    }
  )
}

.sn_run_scvi_integration_impl <- function(object,
                                     method,
                                     batch,
                                     features,
                                     assay,
                                     layer = "counts",
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
  pixi_paths <- sn_get_pixi_paths(environment = method, runtime_dir = runtime_dir)
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
  run_dir <- integration_control$run_dir %||% .sn_default_python_run_dir(method = method, runtime_dir = runtime_dir)
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
    layer = layer,
    batch = batch,
    labels_key = labels_key,
    protein_assay = protein_assay,
    protein_layer = protein_layer %||% "counts",
    protein_features = protein_features
  )

  config <- list(
    method = method,
    batch_key = batch,
    source_layer = layer,
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
    write_h5ad = integration_control$write_h5ad %||% FALSE,
    accelerator = accelerator$requested,
    pixi_environment = pixi_environment
  )
  config_path <- .sn_write_json_file(config, file.path(run_dir, "config.json"))
  script <- .sn_scvi_script_path(integration_control$script %||% NULL)

  backend_run <- .sn_execute_scvi_pixi(
    pixi = integration_control$pixi %||% NULL,
    manifest_path = manifest_path,
    script = script,
    input_dir = input$input_dir,
    output_dir = output_dir,
    config_path = config_path,
    environment = pixi_environment,
    pixi_home = pixi_home,
    install_pixi = integration_control$install_pixi %||% TRUE,
    pixi_version = integration_control$pixi_version %||% "0.69.0",
    pixi_download_url = integration_control$pixi_download_url %||% NULL,
    pixi_sha256 = integration_control$pixi_sha256 %||% NULL,
    verbose = verbose
  )

  result <- .sn_import_scvi_results(
    object = object,
    output_dir = output_dir,
    method = method,
    assay = assay,
    layer = layer,
    batch = batch,
    features = input$features,
    protein_features = input$protein_features,
    reduction = config$reduction,
    max_artifact_import_gb = integration_control$max_artifact_import_gb %||% 0.5
  )
  result$object@misc$integration$backend_performance <- backend_run$performance
  result
}

.sn_scpoli_script_path <- function(script = NULL) {
  if (!is.null(script) && nzchar(script)) {
    script <- path.expand(script)
    if (!file.exists(script)) {
      stop("`integration_control$script` does not exist: ", script, call. = FALSE)
    }
    return(normalizePath(script, winslash = "/", mustWork = TRUE))
  }

  installed <- system.file("pixi", "scarches", "scripts", "scpoli_integration.py", package = "Shennong")
  if (nzchar(installed) && file.exists(installed)) {
    return(normalizePath(installed, winslash = "/", mustWork = TRUE))
  }
  source_path <- file.path(getwd(), "inst", "pixi", "scarches", "scripts", "scpoli_integration.py")
  if (file.exists(source_path)) {
    return(normalizePath(source_path, winslash = "/", mustWork = TRUE))
  }
  stop("Could not locate Shennong's scPoli integration runner.", call. = FALSE)
}

.sn_run_scpoli_integration <- function(object,
                                       batch,
                                       features,
                                       assay,
                                       layer = "counts",
                                       integration_control = list(),
                                       verbose = TRUE) {
  .sn_with_integration_python_run(
    method = "scpoli",
    integration_control = integration_control,
    code = function(control) {
      .sn_run_scpoli_integration_impl(
        object = object,
        batch = batch,
        features = features,
        assay = assay,
        layer = layer,
        integration_control = control,
        verbose = verbose
      )
    }
  )
}

.sn_run_scpoli_integration_impl <- function(object,
                                       batch,
                                       features,
                                       assay,
                                       layer = "counts",
                                       integration_control = list(),
                                       verbose = TRUE) {
  labels_key <- integration_control$label_by %||% NULL
  if (!is.null(labels_key) && (!nzchar(labels_key) || !labels_key %in% colnames(object[[]]))) {
    stop("`integration_control$label_by` must name a metadata column when supplied for scPoli.", call. = FALSE)
  }

  runtime_dir <- .sn_shennong_runtime_dir(integration_control$runtime_dir %||% NULL)
  pixi_paths <- sn_get_pixi_paths(environment = "scpoli", runtime_dir = runtime_dir)
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
  run_dir <- integration_control$run_dir %||% .sn_default_python_run_dir(method = "scpoli", runtime_dir = runtime_dir)
  input_dir <- file.path(run_dir, "input")
  output_dir <- file.path(run_dir, "output")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  manifest_path <- .sn_prepare_scvi_pixi_project(
    project_dir = pixi_project,
    environment = "scpoli",
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
    layer = layer,
    batch = batch,
    labels_key = labels_key
  )
  n_epochs <- as.integer(integration_control$n_epochs %||% integration_control$max_epochs %||% 100L)
  config <- list(
    method = "scpoli",
    batch_key = batch,
    labels_key = labels_key,
    source_layer = layer,
    reduction = integration_control$reduction %||% "scpoli",
    n_latent = integration_control$n_latent %||% 10L,
    embedding_dims = integration_control$embedding_dims %||% 5L,
    latent_batch_size = integration_control$latent_batch_size %||% 2048L,
    seed = integration_control$seed %||% 717L,
    n_epochs = n_epochs,
    pretraining_epochs = integration_control$pretraining_epochs %||% floor(n_epochs * 0.9),
    model_args = integration_control$model_args %||% list(),
    train_args = integration_control$train_args %||% list(),
    write_h5ad = integration_control$write_h5ad %||% FALSE,
    save_model = integration_control$save_model %||% FALSE,
    accelerator = accelerator$requested,
    pixi_environment = pixi_environment
  )
  config_path <- .sn_write_json_file(config, file.path(run_dir, "config.json"))
  backend_run <- .sn_execute_scvi_pixi(
    pixi = integration_control$pixi %||% NULL,
    manifest_path = manifest_path,
    script = .sn_scpoli_script_path(integration_control$script %||% NULL),
    input_dir = input$input_dir,
    output_dir = output_dir,
    config_path = config_path,
    environment = pixi_environment,
    pixi_home = pixi_home,
    install_pixi = integration_control$install_pixi %||% TRUE,
    pixi_version = integration_control$pixi_version %||% "0.69.0",
    pixi_download_url = integration_control$pixi_download_url %||% NULL,
    pixi_sha256 = integration_control$pixi_sha256 %||% NULL,
    verbose = verbose,
    backend_label = "scPoli"
  )
  result <- .sn_import_scvi_results(
    object = object,
    output_dir = output_dir,
    method = "scpoli",
    assay = assay,
    layer = layer,
    batch = batch,
    features = input$features,
    reduction = config$reduction,
    max_artifact_import_gb = integration_control$max_artifact_import_gb %||% 0.5
  )
  result$object@misc$integration$backend_performance <- backend_run$performance
  result
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
