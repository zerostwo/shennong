.sn_trajectory_embedding <- function(object, reduction = NULL, dims = NULL) {
  reductions <- SeuratObject::Reductions(object)
  if (length(reductions) == 0L) {
    stop("`object` must contain a dimensional reduction before trajectory inference.", call. = FALSE)
  }
  if (is_null(reduction)) {
    preferred <- intersect(c("pca", "umap", "tsne"), reductions)
    reduction <- if (length(preferred) > 0L) preferred[[1]] else reductions[[1]]
  }
  if (!reduction %in% reductions) {
    stop("Reduction '", reduction, "' was not found in `object`.", call. = FALSE)
  }
  embedding <- Seurat::Embeddings(object, reduction = reduction)
  if (is_null(dims)) {
    dims <- if (identical(reduction, "pca")) seq_len(min(10L, ncol(embedding))) else seq_len(min(2L, ncol(embedding)))
  }
  dims <- as.integer(dims)
  if (length(dims) == 0L || anyNA(dims) || any(dims < 1L | dims > ncol(embedding))) {
    stop("`dims` must select columns present in the requested reduction.", call. = FALSE)
  }
  if (length(dims) < 2L) {
    stop("Trajectory inference requires at least two reduction dimensions.", call. = FALSE)
  }
  embedding <- as.matrix(embedding[, dims, drop = FALSE])
  list(matrix = embedding, reduction = reduction, dims = dims)
}

.sn_trajectory_clusters <- function(object, cluster_by = NULL) {
  if (is_null(cluster_by)) {
    cluster_by <- if ("seurat_clusters" %in% colnames(object[[]])) "seurat_clusters" else "ident"
  }
  clusters <- if (identical(cluster_by, "ident")) {
    as.character(SeuratObject::Idents(object))
  } else {
    if (!cluster_by %in% colnames(object[[]])) {
      stop("`cluster_by` column '", cluster_by, "' was not found in object metadata.", call. = FALSE)
    }
    as.character(object[[cluster_by, drop = TRUE]])
  }
  names(clusters) <- colnames(object)
  if (anyNA(clusters) || any(!nzchar(clusters))) {
    stop("Trajectory cluster labels cannot contain missing or empty values.", call. = FALSE)
  }
  if (length(unique(clusters)) < 2L) {
    stop("Trajectory inference requires at least two clusters.", call. = FALSE)
  }
  list(values = clusters, cluster_by = cluster_by)
}

.sn_validate_requested_lineages <- function(lineages, clusters) {
  if (is_null(lineages)) return(NULL)
  if (!is.list(lineages) || length(lineages) == 0L) {
    stop("`lineages` must be NULL or a non-empty list of cluster paths.", call. = FALSE)
  }
  lineages <- lapply(lineages, as.character)
  invalid <- unique(setdiff(unlist(lineages, use.names = FALSE), unique(clusters)))
  if (length(invalid) > 0L) {
    stop("Requested lineage cluster(s) were not found: ", paste(invalid, collapse = ", "), ".", call. = FALSE)
  }
  if (any(lengths(lineages) < 2L)) {
    stop("Every requested lineage must contain at least two clusters.", call. = FALSE)
  }
  if (is.null(names(lineages)) || any(!nzchar(names(lineages)))) {
    names(lineages) <- paste0("Lineage", seq_along(lineages))
  }
  lineages
}

.sn_trajectory_curve_table <- function(curves) {
  dplyr::bind_rows(lapply(names(curves), function(lineage) {
    curve <- curves[[lineage]]
    coordinates <- as.matrix(curve$s)
    order <- curve$ord %||% seq_len(nrow(coordinates))
    coordinates <- coordinates[order, , drop = FALSE]
    table <- as.data.frame(coordinates)
    names(table) <- paste0("dimension_", seq_len(ncol(table)))
    table$lineage <- lineage
    table$point <- seq_len(nrow(table))
    tibble::as_tibble(table[, c("lineage", "point", paste0("dimension_", seq_len(ncol(coordinates)))), drop = FALSE])
  }))
}

.sn_trajectory_cell_table <- function(embedding, clusters, pseudotime, weights) {
  lineage_names <- colnames(pseudotime) %||% paste0("Lineage", seq_len(ncol(pseudotime)))
  colnames(pseudotime) <- lineage_names
  colnames(weights) <- lineage_names
  primary_index <- max.col(weights, ties.method = "first")
  no_assignment <- rowSums(weights, na.rm = TRUE) <= 0
  primary_lineage <- lineage_names[primary_index]
  primary_lineage[no_assignment] <- NA_character_
  primary_pseudotime <- pseudotime[cbind(seq_len(nrow(pseudotime)), primary_index)]
  primary_pseudotime[no_assignment] <- NA_real_

  table <- tibble::tibble(
    cell = rownames(embedding),
    cluster = unname(clusters[rownames(embedding)]),
    primary_lineage = primary_lineage,
    primary_pseudotime = as.numeric(primary_pseudotime)
  )
  for (index in seq_along(lineage_names)) {
    safe <- make.names(lineage_names[[index]])
    table[[paste0("pseudotime_", safe)]] <- as.numeric(pseudotime[, index])
    table[[paste0("weight_", safe)]] <- as.numeric(weights[, index])
  }
  table
}

.sn_trajectory_terminal_table <- function(lineages, cells) {
  dplyr::bind_rows(lapply(names(lineages), function(name) {
    terminal <- utils::tail(lineages[[name]], 1L)
    lineage_weight <- cells[[paste0("weight_", make.names(name))]]
    tibble::tibble(
      lineage = name,
      terminal_cluster = terminal,
      n_terminal_cells = sum(cells$cluster == terminal),
      effective_lineage_cells = sum(lineage_weight, na.rm = TRUE)
    )
  }))
}

.sn_trajectory_backend_matrix <- function(value, cells, value_name) {
  if (is.vector(value) && !is.list(value)) {
    names_value <- names(value)
    matrix_value <- matrix(as.numeric(value), ncol = 1L, dimnames = list(names_value %||% cells, "Lineage1"))
  } else if (is.data.frame(value) && all(c("cell", "lineage", value_name) %in% names(value))) {
    lineages <- unique(as.character(value$lineage))
    matrix_value <- matrix(NA_real_, nrow = length(cells), ncol = length(lineages), dimnames = list(cells, lineages))
    row <- match(as.character(value$cell), cells); column <- match(as.character(value$lineage), lineages)
    valid <- !is.na(row) & !is.na(column)
    matrix_value[cbind(row[valid], column[valid])] <- as.numeric(value[[value_name]][valid])
  } else if (is.data.frame(value) && "cell" %in% names(value)) {
    rownames(value) <- as.character(value$cell); value$cell <- NULL
    matrix_value <- as.matrix(value)
  } else {
    matrix_value <- as.matrix(value)
  }
  storage.mode(matrix_value) <- "double"
  if (is_null(rownames(matrix_value))) rownames(matrix_value) <- cells
  missing <- setdiff(cells, rownames(matrix_value))
  if (length(missing) > 0L) stop("Trajectory backend output is missing cell(s): ", paste(utils::head(missing, 5), collapse = ", "), ".", call. = FALSE)
  matrix_value <- matrix_value[cells, , drop = FALSE]
  if (is_null(colnames(matrix_value))) colnames(matrix_value) <- paste0("Lineage", seq_len(ncol(matrix_value)))
  matrix_value
}

.sn_standardize_trajectory_backend <- function(output, embedding, clusters, requested, start, end, method) {
  if (!is.list(output)) stop("Trajectory backend output must be a list.", call. = FALSE)
  pseudotime_value <- output$pseudotime %||% output$pseudotimes %||% output$cell_pseudotime
  if (is_null(pseudotime_value)) stop("Trajectory backend output requires `pseudotime`.", call. = FALSE)
  pseudotime <- .sn_trajectory_backend_matrix(pseudotime_value, rownames(embedding), "pseudotime")
  weights_value <- output$weights %||% output$lineage_weights %||% output$probabilities
  weights <- if (is_null(weights_value)) {
    matrix(as.numeric(is.finite(pseudotime)), nrow = nrow(pseudotime), dimnames = dimnames(pseudotime))
  } else .sn_trajectory_backend_matrix(weights_value, rownames(embedding), "weight")
  if (ncol(weights) != ncol(pseudotime)) stop("Trajectory weights and pseudotime must have the same number of lineages.", call. = FALSE)
  if (any(weights < 0, na.rm = TRUE)) stop("Trajectory lineage weights cannot be negative.", call. = FALSE)
  weights[!is.finite(weights)] <- 0
  colnames(weights) <- colnames(pseudotime)
  lineages <- output$lineages %||% requested
  if (is_null(lineages)) {
    lineages <- lapply(seq_len(ncol(pseudotime)), function(index) {
      scores <- tapply(pseudotime[, index], clusters[rownames(pseudotime)], stats::median, na.rm = TRUE)
      names(sort(scores[is.finite(scores)]))
    })
    names(lineages) <- colnames(pseudotime)
  }
  lineages <- lapply(lineages, as.character)
  if (length(lineages) != ncol(pseudotime) || any(lengths(lineages) == 0L)) {
    stop("Trajectory lineages must provide one non-empty cluster path per pseudotime column.", call. = FALSE)
  }
  if (is_null(names(lineages)) || any(!nzchar(names(lineages)))) names(lineages) <- colnames(pseudotime)
  if (!identical(names(lineages), colnames(pseudotime)) && length(lineages) == ncol(pseudotime)) names(lineages) <- colnames(pseudotime)
  cells <- .sn_trajectory_cell_table(embedding, clusters, pseudotime, weights)
  curves <- output$curves %||% tibble::tibble()
  if (!is.data.frame(curves)) stop("Trajectory backend `curves` must be a data frame when supplied.", call. = FALSE)
  terminals <- output$terminal_states %||% .sn_trajectory_terminal_table(lineages, cells)
  list(
    pseudotime = pseudotime, weights = weights, lineages = lineages,
    cells = cells, terminals = tibble::as_tibble(terminals), curves = tibble::as_tibble(curves),
    backend_cells = tibble::as_tibble(output$cell_metadata %||% tibble::tibble()),
    graphs = output$graphs %||% list(),
    model = output$model %||% NULL, warnings = as.character(output$warnings %||% character()),
    diagnostics = output$diagnostics %||% list(),
    backend = output$backend %||% method
  )
}

.sn_run_monocle3_trajectory <- function(object, embedding, clusters, start, assay, counts_layer, backend_control) {
  check_installed("monocle3", reason = "to run Monocle 3 trajectory inference.")
  check_installed("SingleCellExperiment", reason = "to store Monocle 3 reduced dimensions.")
  assay <- assay %||% SeuratObject::DefaultAssay(object)
  counts <- .sn_get_seurat_layer_data(object, assay = assay, layer = counts_layer)
  cell_metadata <- object[[]][colnames(counts), , drop = FALSE]
  gene_metadata <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))
  cds <- monocle3::new_cell_data_set(counts, cell_metadata = cell_metadata, gene_metadata = gene_metadata)
  umap <- embedding[, seq_len(2L), drop = FALSE]
  rownames(umap) <- rownames(embedding)
  SingleCellExperiment::reducedDim(cds, "UMAP") <- umap
  cluster_args <- utils::modifyList(list(cds = cds, reduction_method = "UMAP", verbose = FALSE), backend_control$monocle3$cluster %||% list(), keep.null = TRUE)
  cds <- do.call(monocle3::cluster_cells, cluster_args)
  graph_args <- utils::modifyList(list(cds = cds, use_partition = backend_control$monocle3$use_partition %||% TRUE, verbose = FALSE), backend_control$monocle3$learn_graph %||% list(), keep.null = TRUE)
  cds <- do.call(monocle3::learn_graph, graph_args)
  root_cells <- backend_control$monocle3$root_cells %||% if (length(start) > 0L) names(clusters)[clusters %in% start] else NULL
  if (is_null(root_cells) || length(root_cells) == 0L) {
    stop("Monocle 3 requires `start` cluster(s) or `backend_control$monocle3$root_cells`.", call. = FALSE)
  }
  cds <- monocle3::order_cells(cds, reduction_method = "UMAP", root_cells = root_cells)
  pseudotime <- monocle3::pseudotime(cds)
  partitions <- monocle3::partitions(cds, reduction_method = "UMAP")
  monocle_clusters <- monocle3::clusters(cds, reduction_method = "UMAP")
  cluster_medians <- tapply(pseudotime, clusters[names(pseudotime)], stats::median, na.rm = TRUE)
  cluster_order <- names(sort(cluster_medians[is.finite(cluster_medians)]))
  list(
    pseudotime = pseudotime,
    weights = stats::setNames(as.numeric(is.finite(pseudotime)), names(pseudotime)),
    lineages = list(Lineage1 = cluster_order),
    curves = tibble::tibble(),
    cell_metadata = tibble::tibble(
      cell = names(pseudotime),
      monocle_cluster = as.character(monocle_clusters[names(pseudotime)]),
      partition = as.character(partitions[names(pseudotime)])
    ),
    graphs = list(principal_graph = monocle3::principal_graph(cds)[["UMAP"]]),
    model = list(class = class(cds), principal_graph = TRUE),
    diagnostics = list(
      partitions = length(unique(partitions)),
      monocle_clusters = length(unique(monocle_clusters)),
      finite_pseudotime_cells = sum(is.finite(pseudotime))
    ),
    backend = "monocle3"
  )
}

.sn_trajectory_slingshot <- function(embedding, clusters, start, end, backend_control = list()) {
  check_installed("slingshot", reason = "to run trajectory inference.")
  defaults <- list(
    data = embedding,
    clusterLabels = unname(clusters[rownames(embedding)]),
    start.clus = start,
    end.clus = end
  )
  fit <- do.call(
    slingshot::slingshot,
    utils::modifyList(defaults, backend_control, keep.null = TRUE)
  )
  pseudotime <- slingshot::slingPseudotime(fit, na = FALSE)
  weights <- slingshot::slingCurveWeights(fit)
  rownames(pseudotime) <- rownames(weights) <- rownames(embedding)
  lineages <- slingshot::slingLineages(fit)
  list(
    pseudotime = pseudotime,
    weights = weights,
    lineages = lineages,
    curves = .sn_trajectory_curve_table(slingshot::slingCurves(fit)),
    model = list(class = class(fit), lineage_names = names(lineages)),
    backend = "slingshot"
  )
}

.sn_trajectory_test_table <- function(table, test) {
  table <- as.data.frame(table)
  table$feature <- rownames(table)
  rownames(table) <- NULL
  if ("pvalue" %in% names(table)) {
    table$adjusted_p_value <- stats::p.adjust(table$pvalue, method = "BH")
  }
  table$test <- test
  tibble::as_tibble(table[, c("feature", "test", setdiff(names(table), c("feature", "test"))), drop = FALSE])
}

.sn_fit_trajectory_dynamics <- function(object,
                                        pseudotime,
                                        weights,
                                        assay,
                                        counts_layer,
                                        features,
                                        max_features,
                                        nknots,
                                        trend_features,
                                        trend_points,
                                        backend_control) {
  check_installed("tradeSeq", reason = "to test dynamic genes along a trajectory.")
  if (!.sn_has_seurat_layer(object, assay = assay, layer = counts_layer)) {
    stop("Counts layer '", counts_layer, "' was not found in assay '", assay, "'.", call. = FALSE)
  }
  counts <- .sn_get_seurat_layer_data(object, assay = assay, layer = counts_layer)
  if (is_null(features)) {
    features <- SeuratObject::VariableFeatures(object[[assay]])
    if (length(features) == 0L) features <- rownames(counts)
  }
  features <- intersect(unique(as.character(features)), rownames(counts))
  max_features <- as.integer(max_features)
  if (!is.finite(max_features) || max_features < 1L) stop("`max_dynamic_features` must be positive.", call. = FALSE)
  features <- utils::head(features, max_features)
  if (length(features) < 3L) stop("Dynamic-gene testing requires at least three retained count features.", call. = FALSE)
  counts <- counts[features, rownames(pseudotime), drop = FALSE]

  defaults <- list(
    counts = counts,
    pseudotime = pseudotime,
    cellWeights = weights,
    nknots = as.integer(nknots),
    verbose = FALSE,
    parallel = FALSE
  )
  fit <- do.call(tradeSeq::fitGAM, utils::modifyList(defaults, backend_control, keep.null = TRUE))
  association <- .sn_trajectory_test_table(tradeSeq::associationTest(fit), "association")
  branch <- tibble::tibble()
  branch_warnings <- character()
  if (ncol(pseudotime) > 1L) {
    branch <- tryCatch(
      dplyr::bind_rows(
        .sn_trajectory_test_table(tradeSeq::patternTest(fit), "pattern"),
        .sn_trajectory_test_table(tradeSeq::diffEndTest(fit), "differential_end")
      ),
      error = function(error) {
        branch_warnings <<- paste0("Branch-specific tradeSeq tests failed: ", conditionMessage(error))
        tibble::tibble()
      }
    )
  } else {
    branch_warnings <- "Branch-specific tests require at least two inferred lineages."
  }

  ranked <- association$feature[order(association$adjusted_p_value, association$pvalue, na.last = TRUE)]
  selected_trends <- intersect(trend_features %||% utils::head(ranked, 20L), rownames(fit))
  trends <- if (length(selected_trends) > 0L) {
    tibble::as_tibble(tradeSeq::predictSmooth(
      fit,
      gene = selected_trends,
      nPoints = as.integer(trend_points),
      tidy = TRUE
    ))
  } else {
    tibble::tibble()
  }
  model_data <- as.data.frame(SummarizedExperiment::rowData(fit)$tradeSeq)
  convergence <- tibble::tibble(
    feature = rownames(fit),
    converged = as.logical(model_data$converged)
  )
  list(
    association = association,
    branch = branch,
    trends = trends,
    convergence = convergence,
    warnings = branch_warnings,
    diagnostics = list(
      tested_features = length(features),
      converged_features = sum(convergence$converged, na.rm = TRUE),
      nknots = as.integer(nknots),
      trend_features = selected_trends
    )
  )
}

#' Infer trajectories and test dynamic genes
#'
#' Runs a cluster-aware Slingshot trajectory, stores per-cell pseudotime and
#' lineage probabilities, and optionally fits tradeSeq negative-binomial GAMs.
#' The returned result follows the unified Shennong analysis-result contract.
#'
#' @param object A Seurat object with a dimensional reduction and cluster labels.
#' @param method Trajectory backend. Slingshot and Monocle 3 run directly;
#'   Palantir accepts a standardized external runner or result.
#' @param reduction Reduction used for inference. Defaults to PCA, then UMAP or
#'   the first available reduction.
#' @param cluster_by Metadata column containing cluster labels. Defaults to
#'   `seurat_clusters`, then active identities.
#' @param start,end Optional start and terminal cluster labels.
#' @param lineages Optional named list of expected cluster paths. Their endpoints
#'   constrain Slingshot and inferred paths are checked against them.
#' @param store_name Name used under the `trajectory` result type.
#' @param dims Reduction dimensions used for inference.
#' @param assay Assay used for dynamic-gene counts.
#' @param counts_layer Raw-count layer used by tradeSeq.
#' @param test_dynamic Fit tradeSeq dynamic and branch tests.
#' @param dynamic_features Optional features tested by tradeSeq.
#' @param max_dynamic_features Maximum number of dynamic features fitted by default.
#' @param nknots Number of tradeSeq spline knots.
#' @param trend_features Features for which fitted trends are retained.
#' @param trend_points Number of fitted points per lineage and feature.
#' @param backend_control Backend controls. Use named `slingshot`, `monocle3`,
#'   and `tradeSeq` argument lists for direct backends. Palantir accepts
#'   `runner` or `result`; the same explicit adapter boundary can override
#'   Monocle 3 for externally managed execution.
#' @param return_object Return the updated object instead of the result.
#'
#' @return A Seurat object or a unified trajectory result.
#'
#' @examples
#' \dontrun{
#' object <- sn_run_trajectory(
#'   object, reduction = "pca", cluster_by = "seurat_clusters",
#'   start = "0", store_name = "development"
#' )
#' result <- sn_get_result(object, "trajectory", "development")
#' result$tables$cells
#' }
#'
#' @export
sn_run_trajectory <- function(object,
                              method = c("slingshot", "monocle3", "palantir"),
                              reduction = NULL,
                              cluster_by = NULL,
                              start = NULL,
                              end = NULL,
                              lineages = NULL,
                              store_name = "trajectory",
                              dims = NULL,
                              assay = NULL,
                              counts_layer = "counts",
                              test_dynamic = TRUE,
                              dynamic_features = NULL,
                              max_dynamic_features = 2000L,
                              nknots = 6L,
                              trend_features = NULL,
                              trend_points = 100L,
                              backend_control = list(),
                              return_object = TRUE) {
  .sn_validate_result_object(object)
  method <- match.arg(method)
  embedding <- .sn_trajectory_embedding(object, reduction = reduction, dims = dims)
  clusters <- .sn_trajectory_clusters(object, cluster_by = cluster_by)
  requested <- .sn_validate_requested_lineages(lineages, clusters$values)
  if (!is_null(requested)) {
    start <- unique(c(start, vapply(requested, `[[`, character(1), 1L)))
    end <- unique(c(end, vapply(requested, function(path) utils::tail(path, 1L), character(1))))
  }
  invalid_endpoints <- setdiff(c(as.character(start), as.character(end)), unique(clusters$values))
  if (length(invalid_endpoints) > 0L) {
    stop("Trajectory endpoint cluster(s) were not found: ", paste(invalid_endpoints, collapse = ", "), ".", call. = FALSE)
  }

  output <- if (identical(method, "slingshot")) {
    .sn_trajectory_slingshot(
      embedding = embedding$matrix,
      clusters = clusters$values,
      start = start,
      end = end,
      backend_control = backend_control$slingshot %||% list()
    )
  } else if (is.function(backend_control$runner)) {
      backend_control$runner(
        object = object, method = method, embedding = embedding$matrix,
        clusters = clusters$values, start = start, end = end,
        lineages = requested, backend_control = backend_control
      )
  } else if (!is_null(backend_control$result)) {
    backend_control$result
  } else if (identical(method, "monocle3")) {
    .sn_run_monocle3_trajectory(object, embedding$matrix, clusters$values, start, assay, counts_layer, backend_control)
  } else {
    stop("Trajectory method '", method, "' requires `backend_control$runner` or `backend_control$result`.", call. = FALSE)
  }
  standardized <- .sn_standardize_trajectory_backend(output, embedding$matrix, clusters$values, requested, start, end, method)
  pseudotime <- standardized$pseudotime; weights <- standardized$weights
  inferred_lineages <- standardized$lineages; cells <- standardized$cells
  terminals <- standardized$terminals; curve_table <- standardized$curves
  backend_model <- standardized$model; backend_name <- standardized$backend
  backend_warnings <- standardized$warnings
  backend_cells <- standardized$backend_cells
  backend_graphs <- standardized$graphs
  backend_diagnostics <- standardized$diagnostics

  requested_warning <- character()
  if (!is_null(requested)) {
    encoded <- function(paths) vapply(paths, paste, collapse = ">", character(1))
    missing_paths <- setdiff(encoded(requested), encoded(inferred_lineages))
    if (length(missing_paths) > 0L) {
      requested_warning <- paste0("Expected lineage path(s) were not recovered exactly: ", paste(missing_paths, collapse = ", "))
    }
  }
  assay <- assay %||% Seurat::DefaultAssay(object)
  dynamics <- list(
    association = tibble::tibble(), branch = tibble::tibble(), trends = tibble::tibble(),
    convergence = tibble::tibble(), warnings = character(),
    diagnostics = list(tested_features = 0L, converged_features = 0L, nknots = as.integer(nknots), trend_features = character())
  )
  if (isTRUE(test_dynamic)) {
    dynamics <- .sn_fit_trajectory_dynamics(
      object = object,
      pseudotime = pseudotime,
      weights = weights,
      assay = assay,
      counts_layer = counts_layer,
      features = dynamic_features,
      max_features = max_dynamic_features,
      nknots = nknots,
      trend_features = trend_features,
      trend_points = trend_points,
      backend_control = backend_control$tradeSeq %||% list()
    )
  }

  object[[paste0(store_name, "_pseudotime")]] <- cells$primary_pseudotime
  object[[paste0(store_name, "_lineage")]] <- cells$primary_lineage
  result <- list(
    schema_version = "1.0",
    analysis_type = "trajectory",
    name = store_name,
    method = method,
    backend = if (isTRUE(test_dynamic)) paste0(backend_name, "+tradeSeq") else backend_name,
    input = list(
      cells = nrow(embedding$matrix),
      reduction = embedding$reduction,
      dimensions = embedding$dims,
      cluster_by = clusters$cluster_by,
      assay = assay,
      counts_layer = if (isTRUE(test_dynamic)) counts_layer else NULL
    ),
    parameters = list(
      start = start,
      end = end,
      requested_lineages = requested,
      test_dynamic = isTRUE(test_dynamic),
      max_dynamic_features = as.integer(max_dynamic_features),
      nknots = as.integer(nknots),
      trend_points = as.integer(trend_points)
    ),
    tables = list(
      cells = cells,
      terminal_states = terminals,
      curves = curve_table,
      dynamic_genes = dynamics$association,
      branch_genes = dynamics$branch,
      fitted_trends = dynamics$trends,
      convergence = dynamics$convergence,
      backend_cells = backend_cells
    ),
    embeddings = list(reduction = embedding$matrix),
    graphs = c(list(lineages = inferred_lineages), backend_graphs),
    models = list(
      trajectory = backend_model,
      slingshot = if (identical(method, "slingshot")) backend_model else NULL,
      tradeSeq = dynamics$diagnostics
    ),
    diagnostics = list(
      n_lineages = ncol(pseudotime),
      assigned_cells = sum(!is.na(cells$primary_lineage)),
      unassigned_cells = sum(is.na(cells$primary_lineage)),
      dynamic = dynamics$diagnostics,
      backend = backend_diagnostics
    ),
    warnings = c(requested_warning, backend_warnings, dynamics$warnings),
    provenance = .sn_analysis_provenance(random_seed = backend_control$seed %||% NA_integer_)
  )
  sn_validate_result(result)
  object <- sn_store_result(object, "trajectory", store_name, result)
  object <- .sn_log_seurat_command(object = object, assay = assay, name = "sn_run_trajectory")
  if (isTRUE(return_object)) object else sn_get_result(object, "trajectory", store_name)
}
