.sn_spatial_coordinates <- function(object, spatial_cols = NULL) {
  columns <- .sn_resolve_spatial_cols(object, spatial_cols = spatial_cols, required = TRUE)
  metadata <- object[[]]
  coordinates <- tibble::tibble(
    cell = rownames(metadata),
    spatial_x = suppressWarnings(as.numeric(metadata[[columns[[1]]]])),
    spatial_y = suppressWarnings(as.numeric(metadata[[columns[[2]]]]))
  )
  if (any(!is.finite(coordinates$spatial_x)) || any(!is.finite(coordinates$spatial_y))) {
    stop("Spatial coordinates must be finite numeric values.", call. = FALSE)
  }
  list(table = coordinates, columns = columns)
}

.sn_spatial_column <- function(candidates, names) {
  hit <- intersect(candidates, names)
  if (length(hit) == 0L) NULL else hit[[1]]
}

.sn_spatial_knn_graph <- function(coordinates, k = 6L) {
  coordinates <- tibble::as_tibble(coordinates)
  n <- nrow(coordinates)
  k <- min(as.integer(k), n - 1L)
  if (n < 3L || k < 1L) stop("Spatial graph construction requires at least three locations.", call. = FALSE)
  xy <- as.matrix(coordinates[, c("spatial_x", "spatial_y")])
  rows <- lapply(seq_len(n), function(index) {
    delta <- sweep(xy, 2, xy[index, ], "-")
    distance <- sqrt(rowSums(delta ^ 2))
    distance[[index]] <- Inf
    neighbors <- order(distance)[seq_len(k)]
    neighbor_distance <- distance[neighbors]
    tibble::tibble(
      source = coordinates$cell[[index]], target = coordinates$cell[neighbors],
      distance = neighbor_distance,
      weight = 1 / pmax(neighbor_distance, .Machine$double.eps)
    )
  })
  dplyr::bind_rows(rows)
}

.sn_spatial_expression <- function(object, assay, layer, features = NULL, max_features = 2000L) {
  expression <- .sn_annotation_expression(object, assay = assay, layer = layer)
  matrix <- expression$matrix
  if (is_null(features)) {
    means <- Matrix::rowMeans(matrix)
    variances <- Matrix::rowMeans(matrix ^ 2) - means ^ 2
    features <- names(utils::head(sort(variances, decreasing = TRUE), as.integer(max_features)))
  }
  features <- intersect(as.character(features), rownames(matrix))
  if (length(features) == 0L) stop("No requested spatial features were found.", call. = FALSE)
  expression$matrix <- matrix[features, , drop = FALSE]
  expression
}

.sn_morans_i <- function(values, graph, cells) {
  values <- as.numeric(values[cells])
  centered <- values - mean(values, na.rm = TRUE)
  denominator <- sum(centered ^ 2, na.rm = TRUE)
  if (!is.finite(denominator) || denominator <= .Machine$double.eps) return(NA_real_)
  source <- match(graph$source, cells)
  target <- match(graph$target, cells)
  weights <- graph$weight / mean(graph$weight)
  length(cells) / sum(weights) * sum(weights * centered[source] * centered[target], na.rm = TRUE) / denominator
}

.sn_spatial_morans_table <- function(expression, coordinates, graph, n_permutations, seed) {
  cells <- intersect(coordinates$cell, colnames(expression$matrix))
  graph <- graph[graph$source %in% cells & graph$target %in% cells, , drop = FALSE]
  set.seed(as.integer(seed))
  rows <- lapply(rownames(expression$matrix), function(feature) {
    values <- stats::setNames(as.numeric(expression$matrix[feature, cells]), cells)
    observed <- .sn_morans_i(values, graph, cells)
    null <- if (n_permutations > 0L && is.finite(observed)) {
      replicate(as.integer(n_permutations), .sn_morans_i(sample(values), graph, cells))
    } else {
      numeric()
    }
    p_value <- if (length(null) > 0L) (1 + sum(abs(null) >= abs(observed), na.rm = TRUE)) / (length(null) + 1) else NA_real_
    tibble::tibble(
      feature = feature, statistic = "morans_i", score = observed,
      p_value = p_value, null_mean = if (length(null)) mean(null, na.rm = TRUE) else NA_real_,
      null_sd = if (length(null) > 1L) stats::sd(null, na.rm = TRUE) else NA_real_
    )
  })
  table <- dplyr::bind_rows(rows)
  table$adjusted_p_value <- stats::p.adjust(table$p_value, method = "BH")
  table$rank <- rank(-table$score, ties.method = "first", na.last = "keep")
  table[order(table$adjusted_p_value, -table$score, na.last = TRUE), , drop = FALSE]
}

.sn_standardize_spatial_features <- function(output, method) {
  table <- tibble::as_tibble(output$table %||% output$features %||% output$svg)
  feature <- .sn_spatial_column(c("feature", "gene", "symbol"), names(table))
  score <- .sn_spatial_column(c("score", "morans_i", "I", "LR_stat", "statistic"), names(table))
  p_value <- .sn_spatial_column(c("p_value", "pval", "pval_norm"), names(table))
  adjusted <- .sn_spatial_column(c("adjusted_p_value", "padj", "FDR", "pval_norm_fdr_bh"), names(table))
  if (is_null(feature) || is_null(score)) stop("Spatial feature output requires feature/gene and score/statistic columns.", call. = FALSE)
  out <- tibble::tibble(
    feature = as.character(table[[feature]]), statistic = method,
    score = suppressWarnings(as.numeric(table[[score]])),
    p_value = if (is_null(p_value)) NA_real_ else suppressWarnings(as.numeric(table[[p_value]])),
    adjusted_p_value = if (is_null(adjusted)) NA_real_ else suppressWarnings(as.numeric(table[[adjusted]]))
  )
  if (all(is.na(out$adjusted_p_value)) && any(is.finite(out$p_value))) out$adjusted_p_value <- stats::p.adjust(out$p_value, "BH")
  out$rank <- rank(-out$score, ties.method = "first", na.last = "keep")
  out
}

#' Find spatially variable features
#'
#' @param object A Seurat object with coordinate metadata.
#' @param method Moran's I, nnSVG, or an explicit SPARK-X adapter.
#' @param spatial_cols Coordinate metadata columns.
#' @param assay,layer Expression assay and layer.
#' @param features Features to test.
#' @param store_name Stored result name.
#' @param backend_control Method controls or an explicit `runner`/`result`.
#' @param return_object Return the modified object or result.
#' @return A Seurat object or unified spatial-feature result.
#' @export
sn_find_spatial_features <- function(object,
                                     method = c("morans_i", "nnsvg", "sparkx"),
                                     spatial_cols = NULL,
                                     assay = NULL,
                                     layer = "data",
                                     features = NULL,
                                     store_name = "spatial_features",
                                     backend_control = list(),
                                     return_object = TRUE) {
  .sn_validate_result_object(object)
  method <- match.arg(method)
  coordinates <- .sn_spatial_coordinates(object, spatial_cols)
  expression <- .sn_spatial_expression(object, assay, layer, features, backend_control$max_features %||% 2000L)
  output <- if (is.function(backend_control$runner)) {
    backend_control$runner(object = object, method = method, coordinates = coordinates$table, expression = expression, backend_control = backend_control)
  } else if (!is_null(backend_control$result)) {
    backend_control$result
  } else if (identical(method, "morans_i")) {
    graph <- .sn_spatial_knn_graph(coordinates$table, backend_control$k %||% 6L)
    list(table = .sn_spatial_morans_table(
      expression, coordinates$table, graph,
      backend_control$n_permutations %||% 99L, backend_control$seed %||% 717L
    ), graph = graph)
  } else if (identical(method, "nnsvg")) {
    check_installed("nnSVG", reason = "to identify spatial features with nnSVG.")
    check_installed("SpatialExperiment", reason = "to construct nnSVG input.")
    counts <- expression$matrix[, coordinates$table$cell, drop = FALSE]
    spe <- SpatialExperiment::SpatialExperiment(
      assays = list(expression = counts),
      spatialCoords = as.matrix(coordinates$table[, c("spatial_x", "spatial_y")])
    )
    fit <- nnSVG::nnSVG(spe, assay_name = "expression", n_neighbors = backend_control$n_neighbors %||% 10L)
    list(table = data.frame(feature = rownames(fit), SummarizedExperiment::rowData(fit), check.names = FALSE), model = fit)
  } else {
    stop("SPARK-X requires `backend_control$runner` or `backend_control$result`.", call. = FALSE)
  }
  table <- if (identical(method, "morans_i") && all(c("feature", "score") %in% names(output$table))) tibble::as_tibble(output$table) else .sn_standardize_spatial_features(output, method)
  graph <- tibble::as_tibble(output$graph %||% tibble::tibble())
  result <- list(
    schema_version = "1.0", analysis_type = "spatial_features", name = store_name,
    method = method, backend = method,
    input = list(assay = expression$assay, layer = expression$layer, coordinate_columns = coordinates$columns, locations = nrow(coordinates$table)),
    parameters = list(k = backend_control$k %||% 6L, n_permutations = backend_control$n_permutations %||% 99L),
    tables = list(primary = table, features = table, coordinates = coordinates$table),
    embeddings = list(spatial = as.matrix(coordinates$table[, c("spatial_x", "spatial_y")])),
    graphs = list(spatial = graph), models = list(backend = output$model %||% NULL),
    diagnostics = list(tested_features = nrow(table), finite_scores = sum(is.finite(table$score))),
    warnings = as.character(output$warnings %||% character()),
    provenance = .sn_analysis_provenance(random_seed = backend_control$seed %||% 717L)
  )
  rownames(result$embeddings$spatial) <- coordinates$table$cell
  sn_validate_result(result)
  object <- sn_store_result(object, "spatial_features", store_name, result)
  object <- .sn_log_seurat_command(object, assay = expression$assay, name = "sn_find_spatial_features")
  if (isTRUE(return_object)) object else sn_get_result(object, "spatial_features", store_name)
}

.sn_standardize_spatial_domains <- function(output, object) {
  domains <- tibble::as_tibble(output$domains %||% output$table)
  cell <- .sn_spatial_column(c("cell", "cell_id", "entity"), names(domains))
  domain <- .sn_spatial_column(c("domain", "cluster", "label", "spatial_domain"), names(domains))
  if (is_null(cell) || is_null(domain)) stop("Spatial domain output requires cell and domain/cluster columns.", call. = FALSE)
  domains <- tibble::tibble(cell = as.character(domains[[cell]]), domain = as.character(domains[[domain]]))
  domains <- domains[domains$cell %in% colnames(object) & !is.na(domains$domain) & nzchar(domains$domain), , drop = FALSE]
  if (nrow(domains) == 0L) stop("No spatial domain assignments matched `object`.", call. = FALSE)
  domains
}

.sn_run_banksy_domains <- function(object, coordinates, assay, layer, control) {
  check_installed("Banksy", reason = "to identify spatial domains with BANKSY.")
  check_installed("SpatialExperiment", reason = "to construct BANKSY input.")
  expression <- .sn_annotation_expression(object, assay, layer)
  matrix <- expression$matrix[, coordinates$cell, drop = FALSE]
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(expression = matrix),
    spatialCoords = as.matrix(coordinates[, c("spatial_x", "spatial_y")])
  )
  original <- colnames(SummarizedExperiment::colData(spe))
  spe <- Banksy::computeBanksy(spe, assay_name = "expression", M = control$M %||% 1L, k_geom = control$k_geom %||% c(15L, 30L))
  spe <- Banksy::runBanksyPCA(spe, M = control$M %||% 1L, lambda = control$lambda %||% 0.8, npcs = control$npcs %||% 20L)
  spe <- Banksy::clusterBanksy(
    spe, M = control$M %||% 1L, lambda = control$lambda %||% 0.8,
    resolution = control$resolution %||% 1, algo = control$algorithm %||% "leiden",
    seed = control$seed %||% 717L
  )
  added <- setdiff(colnames(SummarizedExperiment::colData(spe)), original)
  cluster_column <- utils::tail(added[grepl("clust|cluster", added, ignore.case = TRUE)], 1L)
  if (length(cluster_column) == 0L) stop("BANKSY completed without a discoverable cluster column.", call. = FALSE)
  list(
    domains = tibble::tibble(cell = colnames(spe), domain = as.character(SummarizedExperiment::colData(spe)[[cluster_column]])),
    model = list(class = class(spe), cluster_column = cluster_column)
  )
}

#' Identify spatial domains
#'
#' @param object A Seurat object with coordinates.
#' @param method BANKSY or an explicit stLearn/BayesSpace/CellCharter adapter.
#' @param spatial_cols Coordinate metadata columns.
#' @param assay,layer Expression assay and layer.
#' @param store_name Stored result name.
#' @param backend_control Backend controls or an explicit `runner`/`result`.
#' @param return_object Return the modified object or result.
#' @return A Seurat object or spatial-domain result.
#' @export
sn_find_spatial_domains <- function(object,
                                    method = c("banksy", "stlearn", "bayesspace", "cellcharter"),
                                    spatial_cols = NULL,
                                    assay = NULL,
                                    layer = "counts",
                                    store_name = "spatial_domains",
                                    backend_control = list(),
                                    return_object = TRUE) {
  .sn_validate_result_object(object)
  method <- match.arg(method)
  coordinates <- .sn_spatial_coordinates(object, spatial_cols)
  output <- if (is.function(backend_control$runner)) {
    backend_control$runner(object = object, method = method, coordinates = coordinates$table, assay = assay, layer = layer, backend_control = backend_control)
  } else if (!is_null(backend_control$result)) {
    backend_control$result
  } else if (identical(method, "banksy")) {
    .sn_run_banksy_domains(object, coordinates$table, assay, layer, backend_control)
  } else {
    stop("The ", method, " domain backend requires `backend_control$runner` or `backend_control$result`.", call. = FALSE)
  }
  domains <- .sn_standardize_spatial_domains(output, object)
  values <- stats::setNames(domains$domain, domains$cell)
  object[[store_name]] <- values[colnames(object)]
  result <- list(
    schema_version = "1.0", analysis_type = "spatial_domains", name = store_name,
    method = method, backend = method,
    input = list(assay = assay %||% SeuratObject::DefaultAssay(object), layer = layer, coordinate_columns = coordinates$columns, locations = nrow(coordinates$table)),
    parameters = list(lambda = backend_control$lambda %||% NULL, resolution = backend_control$resolution %||% NULL),
    tables = list(primary = domains, domains = domains, coordinates = coordinates$table),
    embeddings = list(spatial = stats::setNames(as.matrix(coordinates$table[, c("spatial_x", "spatial_y")]), NULL)),
    graphs = list(), models = list(backend = output$model %||% NULL),
    diagnostics = list(domains = length(unique(domains$domain)), assigned_locations = nrow(domains)),
    warnings = as.character(output$warnings %||% character()),
    provenance = .sn_analysis_provenance(random_seed = backend_control$seed %||% 717L)
  )
  rownames(result$embeddings$spatial) <- coordinates$table$cell
  sn_validate_result(result)
  object <- sn_store_result(object, "spatial_domains", store_name, result)
  object <- .sn_log_seurat_command(object, assay = result$input$assay, name = "sn_find_spatial_domains")
  if (isTRUE(return_object)) object else sn_get_result(object, "spatial_domains", store_name)
}

.sn_spatial_enrichment <- function(graph, labels, n_permutations, seed) {
  levels <- sort(unique(labels[!is.na(labels) & nzchar(labels)]))
  pairs <- expand.grid(source_group = levels, target_group = levels, stringsAsFactors = FALSE)
  count_pairs <- function(current) {
    source <- current[graph$source]
    target <- current[graph$target]
    key <- paste(source, target, sep = "\r")
    counts <- table(key)
    values <- as.numeric(counts[paste(pairs$source_group, pairs$target_group, sep = "\r")])
    values[is.na(values)] <- 0
    values
  }
  observed <- count_pairs(labels)
  set.seed(as.integer(seed))
  null <- if (n_permutations > 0L) replicate(as.integer(n_permutations), count_pairs(stats::setNames(sample(labels), names(labels)))) else matrix(numeric(), nrow = nrow(pairs))
  null_mean <- if (length(null)) rowMeans(null) else rep(NA_real_, nrow(pairs))
  null_sd <- if (ncol(null) > 1L) apply(null, 1, stats::sd) else rep(NA_real_, nrow(pairs))
  pairs$observed <- observed
  pairs$expected <- null_mean
  pairs$z_score <- ifelse(is.finite(null_sd) & null_sd > 0, (observed - null_mean) / null_sd, 0)
  pairs$p_value <- if (length(null)) vapply(seq_len(nrow(pairs)), function(i) (1 + sum(abs(null[i, ] - null_mean[[i]]) >= abs(observed[[i]] - null_mean[[i]]))) / (ncol(null) + 1), numeric(1)) else NA_real_
  pairs$adjusted_p_value <- stats::p.adjust(pairs$p_value, "BH")
  tibble::as_tibble(pairs)
}

.sn_spatial_cooccurrence <- function(graph, labels, bins = 4L) {
  breaks <- unique(stats::quantile(graph$distance, probs = seq(0, 1, length.out = as.integer(bins) + 1L), na.rm = TRUE))
  if (length(breaks) < 2L) breaks <- range(graph$distance) + c(-1, 1) * .Machine$double.eps
  graph$distance_bin <- cut(graph$distance, breaks = breaks, include.lowest = TRUE, ordered_result = TRUE)
  graph$source_group <- unname(labels[graph$source])
  graph$target_group <- unname(labels[graph$target])
  table <- as.data.frame(table(graph$source_group, graph$target_group, graph$distance_bin), stringsAsFactors = FALSE)
  names(table) <- c("source_group", "target_group", "distance_bin", "count")
  totals <- tapply(table$count, table$distance_bin, sum)
  table$proportion <- table$count / pmax(1, unname(totals[as.character(table$distance_bin)]))
  tibble::as_tibble(table)
}

#' Analyze spatial neighborhoods
#'
#' @param object A Seurat object with coordinates and labels.
#' @param method Local k-nearest-neighbor analysis or an explicit Squidpy adapter.
#' @param group_by Metadata labels used for enrichment and co-occurrence.
#' @param spatial_cols Coordinate metadata columns.
#' @param store_name Stored result name.
#' @param backend_control Graph/permutation controls or `runner`/`result`.
#' @param return_object Return the modified object or result.
#' @return A Seurat object or spatial-neighborhood result.
#' @export
sn_run_spatial_neighborhood <- function(object,
                                        method = c("knn", "squidpy"),
                                        group_by,
                                        spatial_cols = NULL,
                                        store_name = "spatial_neighborhood",
                                        backend_control = list(),
                                        return_object = TRUE) {
  .sn_validate_result_object(object)
  method <- match.arg(method)
  if (missing(group_by) || !group_by %in% colnames(object[[]])) stop("`group_by` must name object metadata.", call. = FALSE)
  coordinates <- .sn_spatial_coordinates(object, spatial_cols)
  labels <- stats::setNames(as.character(object[[group_by, drop = TRUE]]), colnames(object))
  output <- if (is.function(backend_control$runner)) {
    backend_control$runner(object = object, method = method, coordinates = coordinates$table, labels = labels, backend_control = backend_control)
  } else if (!is_null(backend_control$result)) {
    backend_control$result
  } else if (identical(method, "knn")) {
    graph <- .sn_spatial_knn_graph(coordinates$table, backend_control$k %||% 6L)
    list(
      graph = graph,
      enrichment = .sn_spatial_enrichment(graph, labels, backend_control$n_permutations %||% 99L, backend_control$seed %||% 717L),
      cooccurrence = .sn_spatial_cooccurrence(graph, labels, backend_control$distance_bins %||% 4L)
    )
  } else {
    stop("Squidpy neighborhood analysis requires `backend_control$runner` or `backend_control$result`.", call. = FALSE)
  }
  graph <- tibble::as_tibble(output$graph %||% output$edges)
  enrichment <- tibble::as_tibble(output$enrichment %||% output$nhood_enrichment)
  cooccurrence <- tibble::as_tibble(output$cooccurrence %||% tibble::tibble())
  if (!all(c("source", "target", "distance") %in% names(graph))) stop("Spatial graph requires source, target, and distance columns.", call. = FALSE)
  result <- list(
    schema_version = "1.0", analysis_type = "spatial_neighborhood", name = store_name,
    method = method, backend = method,
    input = list(group_by = group_by, coordinate_columns = coordinates$columns, locations = nrow(coordinates$table)),
    parameters = list(k = backend_control$k %||% 6L, n_permutations = backend_control$n_permutations %||% 99L),
    tables = list(primary = enrichment, enrichment = enrichment, cooccurrence = cooccurrence, coordinates = coordinates$table),
    embeddings = list(), graphs = list(spatial = graph), models = list(backend = output$model %||% NULL),
    diagnostics = list(edges = nrow(graph), label_groups = length(unique(labels)), median_distance = stats::median(graph$distance)),
    warnings = as.character(output$warnings %||% character()),
    provenance = .sn_analysis_provenance(random_seed = backend_control$seed %||% 717L)
  )
  sn_validate_result(result)
  object <- sn_store_result(object, "spatial_neighborhood", store_name, result)
  object <- .sn_log_seurat_command(object, name = "sn_run_spatial_neighborhood")
  if (isTRUE(return_object)) object else sn_get_result(object, "spatial_neighborhood", store_name)
}

.sn_spatial_group_distance <- function(coordinates, labels) {
  xy <- as.matrix(coordinates[, c("spatial_x", "spatial_y")])
  rownames(xy) <- coordinates$cell
  groups <- split(names(labels), labels)
  dplyr::bind_rows(lapply(names(groups), function(source) {
    dplyr::bind_rows(lapply(names(groups), function(target) {
      source_xy <- xy[intersect(groups[[source]], rownames(xy)), , drop = FALSE]
      target_xy <- xy[intersect(groups[[target]], rownames(xy)), , drop = FALSE]
      nearest <- apply(source_xy, 1, function(point) min(sqrt(rowSums(sweep(target_xy, 2, point, "-") ^ 2))))
      tibble::tibble(source = source, target = target, spatial_distance = mean(nearest), source_locations = nrow(source_xy), target_locations = nrow(target_xy))
    }))
  }))
}

#' Add spatial distance evidence to a communication result
#'
#' @param object A Seurat object.
#' @param communication_name Stored communication result name.
#' @param communication Optional communication result supplied directly.
#' @param group_by Metadata column matching communication source/target labels.
#' @param spatial_cols Coordinate metadata columns.
#' @param max_distance Optional maximum mean nearest-group distance.
#' @param store_name Stored result name.
#' @param return_object Return the modified object or result.
#' @return A Seurat object or spatial-communication result.
#' @export
sn_run_spatial_communication <- function(object,
                                         communication_name = "communication",
                                         communication = NULL,
                                         group_by,
                                         spatial_cols = NULL,
                                         max_distance = NULL,
                                         store_name = "spatial_communication",
                                         return_object = TRUE) {
  .sn_validate_result_object(object)
  if (missing(group_by) || !group_by %in% colnames(object[[]])) stop("`group_by` must name object metadata.", call. = FALSE)
  communication <- communication %||% sn_get_result(object, "communication", communication_name)
  sn_validate_result(communication)
  interactions <- tibble::as_tibble(communication$tables$primary)
  if (!all(c("source", "target") %in% names(interactions))) stop("Communication result requires source and target columns.", call. = FALSE)
  coordinates <- .sn_spatial_coordinates(object, spatial_cols)
  labels <- stats::setNames(as.character(object[[group_by, drop = TRUE]]), colnames(object))
  distances <- .sn_spatial_group_distance(coordinates$table, labels)
  spatial <- dplyr::left_join(interactions, distances, by = c("source", "target"))
  spatial$within_distance <- if (is_null(max_distance)) is.finite(spatial$spatial_distance) else is.finite(spatial$spatial_distance) & spatial$spatial_distance <= max_distance
  filtered <- spatial[spatial$within_distance, , drop = FALSE]
  result <- list(
    schema_version = "1.0", analysis_type = "spatial_communication", name = store_name,
    method = paste0(communication$method, "+distance"), backend = communication$backend,
    input = list(communication_name = communication$name, group_by = group_by, coordinate_columns = coordinates$columns, locations = nrow(coordinates$table)),
    parameters = list(max_distance = max_distance),
    tables = list(primary = filtered, all_interactions = spatial, group_distances = distances, coordinates = coordinates$table),
    embeddings = list(), graphs = list(), models = list(source_result = list(name = communication$name, method = communication$method)),
    diagnostics = list(input_interactions = nrow(interactions), retained_interactions = nrow(filtered), distance_groups = nrow(distances)),
    warnings = character(), provenance = .sn_analysis_provenance()
  )
  sn_validate_result(result)
  object <- sn_store_result(object, "spatial_communication", store_name, result)
  object <- .sn_log_seurat_command(object, name = "sn_run_spatial_communication")
  if (isTRUE(return_object)) object else sn_get_result(object, "spatial_communication", store_name)
}

#' Run spatial deconvolution through cell2location
#' @inheritParams sn_run_cell2location
#' @export
sn_run_spatial_deconvolution <- function(...) sn_run_cell2location(...)

#' Map single cells to space through Tangram
#' @inheritParams sn_run_tangram
#' @export
sn_run_spatial_mapping <- function(...) sn_run_tangram(...)

#' Integrate spatial samples with an explicit backend adapter
#'
#' @param object A Seurat object.
#' @param method Integration backend label.
#' @param spatial_cols Coordinate columns.
#' @param store_name Stored result name.
#' @param backend_control Required `runner` or `result` returning an embedding
#'   table with a `cell` column.
#' @param return_object Return the modified object or result.
#' @return A Seurat object or spatial-integration result.
#' @export
sn_integrate_spatial <- function(object,
                                 method = c("staligner", "harmony", "custom"),
                                 spatial_cols = NULL,
                                 store_name = "spatial_integration",
                                 backend_control = list(),
                                 return_object = TRUE) {
  .sn_validate_result_object(object)
  method <- match.arg(method)
  coordinates <- .sn_spatial_coordinates(object, spatial_cols)
  output <- if (is.function(backend_control$runner)) backend_control$runner(object = object, method = method, coordinates = coordinates$table, backend_control = backend_control) else backend_control$result
  if (is_null(output) || !is.list(output)) stop("Spatial integration requires `backend_control$runner` or `backend_control$result`.", call. = FALSE)
  embedding <- tibble::as_tibble(output$embedding %||% output$latent)
  if (!"cell" %in% names(embedding) || ncol(embedding) < 3L) stop("Spatial integration embedding requires cell plus at least two dimensions.", call. = FALSE)
  matrix <- as.matrix(embedding[, setdiff(names(embedding), "cell"), drop = FALSE])
  storage.mode(matrix) <- "numeric"
  rownames(matrix) <- embedding$cell
  result <- list(
    schema_version = "1.0", analysis_type = "spatial_integration", name = store_name,
    method = method, backend = method,
    input = list(coordinate_columns = coordinates$columns, locations = nrow(coordinates$table)),
    parameters = list(), tables = list(primary = embedding, coordinates = coordinates$table),
    embeddings = list(integrated = matrix), graphs = list(), models = list(backend = output$model %||% NULL),
    diagnostics = list(dimensions = ncol(matrix), locations = nrow(matrix)),
    warnings = as.character(output$warnings %||% character()), provenance = .sn_analysis_provenance(random_seed = backend_control$seed %||% NA_integer_)
  )
  sn_validate_result(result)
  object <- sn_store_result(object, "spatial_integration", store_name, result)
  object <- .sn_log_seurat_command(object, name = "sn_integrate_spatial")
  if (isTRUE(return_object)) object else sn_get_result(object, "spatial_integration", store_name)
}

#' Unified spatial workflow dispatcher
#'
#' @param object A Seurat object.
#' @param task Spatial task.
#' @param method Backend method.
#' @param ... Arguments forwarded to the task-specific function.
#' @return The task-specific result.
#' @export
sn_run_spatial <- function(object,
                           task = c("qc", "svg", "domain", "neighborhood", "deconvolution", "mapping", "integration", "communication"),
                           method = "auto",
                           ...) {
  task <- match.arg(task)
  if (identical(task, "qc")) {
    coordinates <- .sn_spatial_coordinates(object, list(...)$spatial_cols %||% NULL)
    graph <- .sn_spatial_knn_graph(coordinates$table, list(...)$k %||% 6L)
    return(list(coordinates = coordinates$table, graph = graph, diagnostics = list(locations = nrow(coordinates$table), median_neighbor_distance = stats::median(graph$distance))))
  }
  dispatch <- switch(
    task,
    svg = sn_find_spatial_features,
    domain = sn_find_spatial_domains,
    neighborhood = sn_run_spatial_neighborhood,
    deconvolution = sn_run_spatial_deconvolution,
    mapping = sn_run_spatial_mapping,
    integration = sn_integrate_spatial,
    communication = sn_run_spatial_communication
  )
  args <- list(...)
  if (!identical(method, "auto") && task %in% c("svg", "domain", "neighborhood", "integration")) args$method <- method
  do.call(dispatch, c(list(object = object), args))
}
