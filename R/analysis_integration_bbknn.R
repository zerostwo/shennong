# BBKNN integration backend.
#
# This module owns BBKNN script discovery, managed-runtime execution, graph
# import, and the object-level integration adapter.

.sn_bbknn_script_path <- function(script = NULL) {
  if (!is.null(script) && nzchar(script)) {
    script <- path.expand(script)
    if (!file.exists(script)) {
      stop("`integration_control$script` does not exist: ", script, call. = FALSE)
    }
    return(normalizePath(script, winslash = "/", mustWork = TRUE))
  }
  installed <- system.file("pixi", "bbknn", "scripts", "bbknn_integration.py", package = "Shennong")
  if (nzchar(installed) && file.exists(installed)) {
    return(normalizePath(installed, winslash = "/", mustWork = TRUE))
  }
  source_path <- file.path(getwd(), "inst", "pixi", "bbknn", "scripts", "bbknn_integration.py")
  if (file.exists(source_path)) {
    return(normalizePath(source_path, winslash = "/", mustWork = TRUE))
  }
  stop("Could not locate Shennong's BBKNN integration runner.", call. = FALSE)
}

.sn_run_bbknn_integration <- function(object,
                                      batch,
                                      reduction = "pca",
                                      assay,
                                      layer = "counts",
                                      dims,
                                      integration_control = list(),
                                      verbose = TRUE) {
  .sn_with_integration_python_run(
    method = "bbknn",
    integration_control = integration_control,
    code = function(control) {
      .sn_run_bbknn_integration_impl(
        object = object,
        batch = batch,
        reduction = reduction,
        assay = assay,
        layer = layer,
        dims = dims,
        integration_control = control,
        verbose = verbose
      )
    }
  )
}

.sn_run_bbknn_integration_impl <- function(object,
                                      batch,
                                      reduction = "pca",
                                      assay,
                                      layer = "counts",
                                      dims,
                                      integration_control = list(),
                                      verbose = TRUE) {
  if (is.null(batch) || !nzchar(batch) || !batch %in% colnames(object[[]])) {
    stop("BBKNN integration requires `batch` to name a metadata column.", call. = FALSE)
  }
  dims <- .sn_valid_reduction_dims(object = object, reduction = reduction, dims = dims)
  embedding <- Seurat::Embeddings(object = object, reduction = reduction)[, dims, drop = FALSE]
  runtime_dir <- .sn_shennong_runtime_dir(integration_control$runtime_dir %||% NULL)
  pixi_paths <- sn_get_pixi_paths(environment = "bbknn", runtime_dir = runtime_dir)
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
  run_dir <- integration_control$run_dir %||% .sn_default_python_run_dir(method = "bbknn", runtime_dir = runtime_dir)
  input_dir <- file.path(run_dir, "input")
  output_dir <- file.path(run_dir, "output")
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(embedding, file.path(input_dir, "pca.csv"), quote = TRUE)
  utils::write.csv(
    data.frame(cell_id = rownames(embedding), batch = as.character(object[[]][rownames(embedding), batch]), stringsAsFactors = FALSE),
    file.path(input_dir, "obs.csv"),
    row.names = FALSE,
    quote = TRUE
  )
  manifest_path <- .sn_prepare_scvi_pixi_project(
    project_dir = pixi_project,
    environment = "bbknn",
    manifest_path = integration_control$manifest_path %||% NULL,
    manifest_lines = integration_control$manifest_lines %||% NULL,
    overwrite = isTRUE(integration_control$overwrite_manifest),
    platforms = integration_control$platforms %||% NULL
  )
  graph_name <- integration_control$graph_name %||% "bbknn_snn"
  config <- list(
    method = "bbknn",
    batch_key = "batch",
    source_layer = layer,
    graph_name = graph_name,
    seed = integration_control$seed %||% 717L,
    bbknn_args = integration_control$bbknn_args %||% list(),
    umap_args = integration_control$umap_args %||% list()
  )
  config_path <- .sn_write_json_file(config, file.path(run_dir, "config.json"))
  backend_run <- .sn_execute_scvi_pixi(
    pixi = integration_control$pixi %||% NULL,
    manifest_path = manifest_path,
    script = .sn_bbknn_script_path(integration_control$script %||% NULL),
    input_dir = normalizePath(input_dir, winslash = "/", mustWork = TRUE),
    output_dir = output_dir,
    config_path = config_path,
    environment = integration_control$environment %||% "default",
    pixi_home = pixi_home,
    install_pixi = integration_control$install_pixi %||% TRUE,
    pixi_version = integration_control$pixi_version %||% "0.69.0",
    pixi_download_url = integration_control$pixi_download_url %||% NULL,
    pixi_sha256 = integration_control$pixi_sha256 %||% NULL,
    verbose = verbose,
    backend_label = "BBKNN"
  )
  graph_path <- file.path(output_dir, "connectivities.mtx")
  distances_path <- file.path(output_dir, "distances.mtx")
  cells_path <- file.path(output_dir, "cells.csv")
  required_paths <- c(graph_path, distances_path, cells_path)
  if (any(!file.exists(required_paths))) {
    stop(
      "BBKNN output is missing required file(s): ",
      paste(basename(required_paths[!file.exists(required_paths)]), collapse = ", "),
      call. = FALSE
    )
  }
  max_import_gb <- integration_control$max_artifact_import_gb %||% 0.5
  .sn_assert_python_artifact_budget(
    required_paths,
    max_import_gb,
    "BBKNN graph outputs"
  )
  graph_cells <- utils::read.csv(cells_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!identical(colnames(graph_cells), "cell_id")) {
    stop("BBKNN `cells.csv` must contain exactly one `cell_id` column.", call. = FALSE)
  }
  .sn_validate_exact_python_ids(
    graph_cells$cell_id,
    colnames(object),
    "BBKNN graph cell",
    ordered = TRUE
  )
  graph <- .sn_as_sparse_matrix(Matrix::readMM(graph_path))
  if (length(dim(graph)) != 2L || any(dim(graph) != c(ncol(object), ncol(object)))) {
    stop("BBKNN connectivity output must contain one row and column per input cell.", call. = FALSE)
  }
  distances <- .sn_as_sparse_matrix(Matrix::readMM(distances_path))
  if (length(dim(distances)) != 2L || any(dim(distances) != c(ncol(object), ncol(object)))) {
    stop("BBKNN distance output must contain one row and column per input cell.", call. = FALSE)
  }
  for (entry in list(connectivity = graph, distance = distances)) {
    values <- entry@x
    if (anyNA(values) || any(!is.finite(values)) || any(values < 0)) {
      stop("BBKNN graph outputs must contain finite, non-negative weights.", call. = FALSE)
    }
  }
  if (!isTRUE(Matrix::isSymmetric(graph, tol = 1e-7))) {
    stop("BBKNN connectivity output must be symmetric.", call. = FALSE)
  }
  if (any(Matrix::rowSums(graph != 0) == 0)) {
    stop("BBKNN connectivity output must connect every input cell.", call. = FALSE)
  }
  dimnames(graph) <- list(colnames(object), colnames(object))
  object[[graph_name]] <- SeuratObject::as.Graph(graph)
  umap_path <- file.path(output_dir, "umap.csv")
  if (!file.exists(umap_path)) {
    stop("BBKNN output is missing `umap.csv`: ", umap_path, call. = FALSE)
  }
  .sn_assert_python_artifact_budget(umap_path, max_import_gb, "BBKNN UMAP output")
  umap <- .sn_read_embedding_csv(umap_path, cells = colnames(object))
  umap_reduction <- integration_control$umap_reduction %||% "umap"
  colnames(umap) <- paste0("UMAP_", seq_len(ncol(umap)))
  object[[umap_reduction]] <- Seurat::CreateDimReducObject(
    embeddings = umap,
    key = "UMAP_",
    assay = assay
  )
  backend_manifest <- .sn_read_integration_backend_manifest(
    output_dir = output_dir,
    method = "bbknn",
    n_cells = ncol(object)
  )
  manifest_pcs <- suppressWarnings(as.integer(backend_manifest$n_pcs %||% NA_integer_))
  manifest_batches <- suppressWarnings(as.integer(backend_manifest$n_batches %||% NA_integer_))
  manifest_umap_dims <- suppressWarnings(as.integer(backend_manifest$umap_dimensions %||% NA_integer_))
  if (length(manifest_pcs) != 1L || is.na(manifest_pcs) || manifest_pcs != length(dims)) {
    stop("BBKNN manifest `n_pcs` does not match the exported PCA dimensions.", call. = FALSE)
  }
  expected_batches <- length(unique(as.character(object[[batch, drop = TRUE]])))
  if (length(manifest_batches) != 1L || is.na(manifest_batches) ||
      manifest_batches != expected_batches) {
    stop("BBKNN manifest `n_batches` does not match the exported batch metadata.", call. = FALSE)
  }
  if (length(manifest_umap_dims) != 1L || is.na(manifest_umap_dims) ||
      manifest_umap_dims != ncol(umap)) {
    stop("BBKNN manifest `umap_dimensions` does not match `umap.csv`.", call. = FALSE)
  }
  object@misc$integration <- list(
    method = "bbknn",
    batch_by = batch,
    reduction = reduction,
    assay = assay,
    source_layer = layer,
    input_dims = dims,
    graph = graph_name,
    umap_reduction = umap_reduction,
    run_dir = normalizePath(output_dir, winslash = "/", mustWork = TRUE),
    bbknn_version = backend_manifest$bbknn_version %||% NULL,
    backend_performance = backend_run$performance
  )
  list(object = object, reduction = reduction, graph = graph_name, umap = umap_reduction)
}
