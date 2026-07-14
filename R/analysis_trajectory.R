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
#' @param method Trajectory backend. Slingshot is currently implemented;
#'   Monocle 3 and Palantir remain discoverable planned backends.
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
#' @param backend_control Named `slingshot` and `tradeSeq` argument lists.
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
  if (!identical(method, "slingshot")) {
    stop("Trajectory method '", method, "' is registered but not implemented yet; use `sn_method_status()` for setup details.", call. = FALSE)
  }
  check_installed("slingshot", reason = "to run trajectory inference.")
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

  defaults <- list(
    data = embedding$matrix,
    clusterLabels = unname(clusters$values[rownames(embedding$matrix)]),
    start.clus = start,
    end.clus = end
  )
  fit <- do.call(slingshot::slingshot, utils::modifyList(defaults, backend_control$slingshot %||% list(), keep.null = TRUE))
  pseudotime <- slingshot::slingPseudotime(fit, na = FALSE)
  weights <- slingshot::slingCurveWeights(fit)
  inferred_lineages <- slingshot::slingLineages(fit)
  curves <- slingshot::slingCurves(fit)
  rownames(pseudotime) <- rownames(weights) <- rownames(embedding$matrix)
  cells <- .sn_trajectory_cell_table(embedding$matrix, clusters$values, pseudotime, weights)
  terminals <- .sn_trajectory_terminal_table(inferred_lineages, cells)
  curve_table <- .sn_trajectory_curve_table(curves)

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
    backend = if (isTRUE(test_dynamic)) "slingshot+tradeSeq" else "slingshot",
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
      convergence = dynamics$convergence
    ),
    embeddings = list(reduction = embedding$matrix),
    graphs = list(lineages = inferred_lineages),
    models = list(
      slingshot = list(class = class(fit), lineage_names = names(inferred_lineages)),
      tradeSeq = dynamics$diagnostics
    ),
    diagnostics = list(
      n_lineages = ncol(pseudotime),
      assigned_cells = sum(!is.na(cells$primary_lineage)),
      unassigned_cells = sum(is.na(cells$primary_lineage)),
      dynamic = dynamics$diagnostics
    ),
    warnings = c(requested_warning, dynamics$warnings),
    provenance = .sn_analysis_provenance(random_seed = backend_control$seed %||% NA_integer_)
  )
  sn_validate_result(result)
  object <- sn_store_result(object, "trajectory", store_name, result)
  object <- .sn_log_seurat_command(object = object, assay = assay, name = "sn_run_trajectory")
  if (isTRUE(return_object)) object else sn_get_result(object, "trajectory", store_name)
}
