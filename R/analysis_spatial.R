.sn_spatial_coordinates <- function(object, spatial_cols = NULL, sample_by = NULL) {
  columns <- .sn_resolve_spatial_cols(object, spatial_cols = spatial_cols, required = TRUE)
  metadata <- object[[]]
  if (!is_null(sample_by) && !sample_by %in% colnames(metadata)) {
    stop("`sample_by` column '", sample_by, "' was not found in object metadata.", call. = FALSE)
  }
  coordinates <- tibble::tibble(
    cell = rownames(metadata),
    spatial_x = suppressWarnings(as.numeric(metadata[[columns[[1]]]])),
    spatial_y = suppressWarnings(as.numeric(metadata[[columns[[2]]]])),
    spatial_sample = if (is_null(sample_by)) "__single_section__" else as.character(metadata[[sample_by]])
  )
  if (any(!is.finite(coordinates$spatial_x)) || any(!is.finite(coordinates$spatial_y))) {
    stop("Spatial coordinates must be finite numeric values.", call. = FALSE)
  }
  if (anyNA(coordinates$spatial_sample) || any(!nzchar(coordinates$spatial_sample))) {
    stop("Spatial sample/section labels cannot be missing or empty.", call. = FALSE)
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
  k <- as.integer(k)
  if (n < 3L || length(k) != 1L || is.na(k) || k < 1L) {
    stop("Spatial graph construction requires at least three locations and a positive `k`.", call. = FALSE)
  }
  section <- coordinates$spatial_sample %||% rep("__single_section__", n)
  xy <- as.matrix(coordinates[, c("spatial_x", "spatial_y")])
  rows <- lapply(split(seq_len(n), section), function(indices) {
    if (length(indices) < 2L) return(tibble::tibble())
    section_k <- min(k, length(indices) - 1L)
    dplyr::bind_rows(lapply(indices, function(index) {
      delta <- sweep(xy[indices, , drop = FALSE], 2, xy[index, ], "-")
      distance <- sqrt(rowSums(delta ^ 2))
      distance[indices == index] <- Inf
      neighbors <- indices[order(distance)[seq_len(section_k)]]
      neighbor_distance <- sqrt(rowSums(sweep(xy[neighbors, , drop = FALSE], 2, xy[index, ], "-") ^ 2))
      tibble::tibble(
        source = coordinates$cell[[index]], target = coordinates$cell[neighbors],
        distance = neighbor_distance,
        weight = 1 / pmax(neighbor_distance, .Machine$double.eps),
        spatial_sample = section[[index]]
      )
    }))
  })
  graph <- dplyr::bind_rows(rows)
  if (nrow(graph) == 0L) {
    stop("No within-sample spatial neighbors could be constructed.", call. = FALSE)
  }
  graph
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
  if (!is_null(names(values))) {
    missing <- setdiff(cells, names(values))
    if (length(missing) > 0L) stop("Moran values are missing requested cells.", call. = FALSE)
    values <- as.numeric(values[cells])
  } else {
    if (length(values) != length(cells)) stop("Unnamed Moran values must match the cell count.", call. = FALSE)
    values <- as.numeric(values)
  }
  centered <- values - mean(values, na.rm = TRUE)
  denominator <- sum(centered ^ 2, na.rm = TRUE)
  if (!is.finite(denominator) || denominator <= .Machine$double.eps) return(NA_real_)
  source <- match(graph$source, cells)
  target <- match(graph$target, cells)
  weights <- graph$weight / mean(graph$weight)
  length(cells) / sum(weights) * sum(weights * centered[source] * centered[target], na.rm = TRUE) / denominator
}

.sn_spatial_permute_values <- function(values, cells, samples = NULL) {
  if (!is_null(names(values))) {
    missing <- setdiff(cells, names(values))
    if (length(missing) > 0L) stop("Spatial permutation values are missing cells.", call. = FALSE)
    values <- values[cells]
  } else {
    if (length(values) != length(cells)) stop("Spatial permutation values must match the cell count.", call. = FALSE)
    names(values) <- cells
  }
  samples <- samples %||% stats::setNames(rep("__single_section__", length(cells)), cells)
  samples <- as.character(samples[cells])
  if (anyNA(samples)) stop("Spatial permutation strata are missing cells.", call. = FALSE)
  permuted <- values
  for (indices in split(seq_along(cells), samples)) {
    permuted[indices] <- unname(values[indices])[sample.int(length(indices))]
  }
  permuted
}

.sn_validate_spatial_permutation_count <- function(n_permutations) {
  if (!is.numeric(n_permutations) || length(n_permutations) != 1L ||
      is.na(n_permutations) || !is.finite(n_permutations) ||
      n_permutations < 0 || n_permutations > .Machine$integer.max ||
      n_permutations != floor(n_permutations)) {
    stop("`n_permutations` must be one non-negative integer.", call. = FALSE)
  }
  as.integer(n_permutations)
}

.sn_spatial_morans_table <- function(expression, coordinates, graph, n_permutations, seed) {
  n_permutations <- .sn_validate_spatial_permutation_count(n_permutations)
  cells <- intersect(coordinates$cell, colnames(expression$matrix))
  graph <- graph[graph$source %in% cells & graph$target %in% cells, , drop = FALSE]
  samples <- stats::setNames(coordinates$spatial_sample, coordinates$cell)
  rows <- .sn_with_seed(seed, lapply(rownames(expression$matrix), function(feature) {
    values <- stats::setNames(as.numeric(expression$matrix[feature, cells]), cells)
    observed <- .sn_morans_i(values, graph, cells)
    null <- if (n_permutations > 0L && is.finite(observed)) {
      replicate(
        as.integer(n_permutations),
        .sn_morans_i(.sn_spatial_permute_values(values, cells, samples), graph, cells)
      )
    } else {
      numeric()
    }
    finite_null <- null[is.finite(null)]
    null_mean <- if (length(finite_null) > 0L) mean(finite_null) else NA_real_
    p_value <- if (length(finite_null) > 0L) {
      (1 + sum(abs(finite_null - null_mean) >= abs(observed - null_mean))) /
        (length(finite_null) + 1)
    } else {
      NA_real_
    }
    tibble::tibble(
      feature = feature, statistic = "morans_i", score = observed,
      p_value = p_value, null_mean = null_mean,
      null_sd = if (length(finite_null) > 1L) stats::sd(finite_null) else NA_real_
    )
  }))
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
#' @param sample_by Optional metadata column defining independent samples or
#'   tissue sections. Spatial graphs and permutations are restricted within
#'   these boundaries.
#' @param assay,layer Expression assay and layer.
#' @param features Features to test.
#' @param result_id Stored result name.
#' @param backend_control Method controls or an explicit `runner`/`result`.
#' @param return_object Return the modified object or result.
#' @param seed Top-level reproducibility seed. Precedence: \code{seed} >
#'   \code{backend_control$seed} > task default; stamped into provenance.
#' @param verbose Top-level progress switch forwarded through
#'   \code{backend_control$verbose} when explicitly supplied.
#' @return A Seurat object or unified spatial-feature result.
#' @export
sn_find_spatial_features <- function(object,
                                     method = c("morans_i", "nnsvg", "sparkx"),
                                     spatial_cols = NULL,
                                     assay = NULL,
                                     layer = "data",
                                     features = NULL,
                                     result_id = "spatial_features",
                                     backend_control = list(),
                                     return_object = TRUE,
                                     seed = NULL,
                                     verbose = TRUE,
                                     sample_by = NULL) {
  result_id <- .sn_validate_result_id(result_id)
  .sn_validate_seurat_object(object)
  method <- match.arg(method)
  backend_control$seed <- seed %||% backend_control$seed
  if (!missing(verbose)) {
    backend_control$verbose <- isTRUE(verbose)
  }
  coordinates <- .sn_spatial_coordinates(object, spatial_cols, sample_by = sample_by)
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
    if (!is_null(sample_by) && length(unique(coordinates$table$spatial_sample)) > 1L) {
      stop(
        "The built-in nnSVG path accepts one tissue section at a time. ",
        "Subset by `sample_by` or supply an explicit sample-aware runner/result.",
        call. = FALSE
      )
    }
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
    schema_version = .sn_analysis_result_schema_version(), analysis_type = "spatial_features", result_id = result_id,
    method = method, backend = method,
    input = list(assay = expression$assay, layer = expression$layer, coordinate_columns = coordinates$columns, sample_by = sample_by, locations = nrow(coordinates$table)),
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
  object <- sn_store_result(object, "spatial_features", result_id, result)
  object <- .sn_log_seurat_command(object, assay = expression$assay, name = "sn_find_spatial_features")
  if (isTRUE(return_object)) object else sn_get_result(object, "spatial_features", result_id)
}

.sn_standardize_spatial_domains <- function(output, object) {
  domains <- tibble::as_tibble(output$domains %||% output$table)
  cell <- .sn_spatial_column(c("cell", "cell_id", "entity"), names(domains))
  domain <- .sn_spatial_column(c("domain", "cluster", "label", "spatial_domain"), names(domains))
  if (is_null(cell) || is_null(domain)) stop("Spatial domain output requires cell and domain/cluster columns.", call. = FALSE)
  domains <- tibble::tibble(cell = as.character(domains[[cell]]), domain = as.character(domains[[domain]]))
  if (anyNA(domains$cell) || any(!nzchar(domains$cell)) || anyDuplicated(domains$cell)) {
    stop("Spatial domain cell identifiers must be unique and non-empty.", call. = FALSE)
  }
  expected <- colnames(object)
  unknown <- setdiff(domains$cell, expected)
  missing <- setdiff(expected, domains$cell)
  if (length(unknown) > 0L || length(missing) > 0L) {
    stop(
      "Spatial domain assignments must match every object cell exactly; unknown = ",
      length(unknown), ", missing = ", length(missing), ".",
      call. = FALSE
    )
  }
  if (anyNA(domains$domain) || any(!nzchar(domains$domain))) {
    stop("Spatial domain labels must be non-missing and non-empty.", call. = FALSE)
  }
  domains[match(expected, domains$cell), , drop = FALSE]
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
#' @param sample_by Optional metadata column defining independent tissue
#'   sections. The built-in BANKSY path fails closed for pooled sections.
#' @param assay,layer Expression assay and layer.
#' @param result_id Stored result name.
#' @param backend_control Backend controls or an explicit `runner`/`result`.
#' @param return_object Return the modified object or result.
#' @param seed Top-level reproducibility seed. Precedence: \code{seed} >
#'   \code{backend_control$seed} > task default; stamped into provenance.
#' @param verbose Top-level progress switch forwarded through
#'   \code{backend_control$verbose} when explicitly supplied.
#' @return A Seurat object or spatial-domain result.
#' @export
sn_find_spatial_domains <- function(object,
                                    method = c("banksy", "stlearn", "bayesspace", "cellcharter"),
                                    spatial_cols = NULL,
                                    assay = NULL,
                                    layer = "counts",
                                    result_id = "spatial_domains",
                                    backend_control = list(),
                                    return_object = TRUE,
                                    seed = NULL,
                                    verbose = TRUE,
                                    sample_by = NULL) {
  result_id <- .sn_validate_result_id(result_id)
  .sn_validate_seurat_object(object)
  method <- match.arg(method)
  backend_control$seed <- seed %||% backend_control$seed
  if (!missing(verbose)) {
    backend_control$verbose <- isTRUE(verbose)
  }
  coordinates <- .sn_spatial_coordinates(object, spatial_cols, sample_by = sample_by)
  output <- if (is.function(backend_control$runner)) {
    backend_control$runner(object = object, method = method, coordinates = coordinates$table, assay = assay, layer = layer, backend_control = backend_control)
  } else if (!is_null(backend_control$result)) {
    backend_control$result
  } else if (identical(method, "banksy")) {
    if (!is_null(sample_by) && length(unique(coordinates$table$spatial_sample)) > 1L) {
      stop(
        "The built-in BANKSY path accepts one tissue section at a time. ",
        "Subset by `sample_by` or supply an explicit sample-aware runner/result.",
        call. = FALSE
      )
    }
    .sn_run_banksy_domains(object, coordinates$table, assay, layer, backend_control)
  } else {
    stop("The ", method, " domain backend requires `backend_control$runner` or `backend_control$result`.", call. = FALSE)
  }
  domains <- .sn_standardize_spatial_domains(output, object)
  values <- stats::setNames(domains$domain, domains$cell)
  object[[result_id]] <- values[colnames(object)]
  result <- list(
    schema_version = .sn_analysis_result_schema_version(), analysis_type = "spatial_domains", result_id = result_id,
    method = method, backend = method,
    input = list(assay = assay %||% SeuratObject::DefaultAssay(object), layer = layer, coordinate_columns = coordinates$columns, sample_by = sample_by, locations = nrow(coordinates$table)),
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
  object <- sn_store_result(object, "spatial_domains", result_id, result)
  object <- .sn_log_seurat_command(object, assay = result$input$assay, name = "sn_find_spatial_domains")
  if (isTRUE(return_object)) object else sn_get_result(object, "spatial_domains", result_id)
}

.sn_spatial_enrichment <- function(graph, labels, n_permutations, seed, samples = NULL) {
  n_permutations <- .sn_validate_spatial_permutation_count(n_permutations)
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
  null <- if (n_permutations > 0L) {
    values <- .sn_with_seed(seed, replicate(
      n_permutations,
      count_pairs(.sn_spatial_permute_values(labels, names(labels), samples))
    ))
    matrix(values, nrow = nrow(pairs), ncol = n_permutations)
  } else {
    matrix(numeric(), nrow = nrow(pairs))
  }
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
#' @param sample_by Optional metadata column defining independent samples or
#'   tissue sections. Neighbors and label permutations never cross boundaries.
#' @param result_id Stored result name.
#' @param backend_control Graph/permutation controls or `runner`/`result`.
#' @param return_object Return the modified object or result.
#' @return A Seurat object or spatial-neighborhood result.
#' @export
sn_run_spatial_neighborhood <- function(object,
                                        method = c("knn", "squidpy"),
                                        group_by,
                                        spatial_cols = NULL,
                                        result_id = "spatial_neighborhood",
                                        backend_control = list(),
                                        return_object = TRUE,
                                        sample_by = NULL) {
  result_id <- .sn_validate_result_id(result_id)
  .sn_validate_seurat_object(object)
  method <- match.arg(method)
  if (missing(group_by) || !group_by %in% colnames(object[[]])) stop("`group_by` must name object metadata.", call. = FALSE)
  coordinates <- .sn_spatial_coordinates(object, spatial_cols, sample_by = sample_by)
  labels <- stats::setNames(as.character(object[[group_by, drop = TRUE]]), colnames(object))
  samples <- stats::setNames(coordinates$table$spatial_sample, coordinates$table$cell)
  output <- if (is.function(backend_control$runner)) {
    backend_control$runner(object = object, method = method, coordinates = coordinates$table, labels = labels, backend_control = backend_control)
  } else if (!is_null(backend_control$result)) {
    backend_control$result
  } else if (identical(method, "knn")) {
    graph <- .sn_spatial_knn_graph(coordinates$table, backend_control$k %||% 6L)
    list(
      graph = graph,
      enrichment = .sn_spatial_enrichment(graph, labels, backend_control$n_permutations %||% 99L, backend_control$seed %||% 717L, samples = samples),
      cooccurrence = .sn_spatial_cooccurrence(graph, labels, backend_control$distance_bins %||% 4L)
    )
  } else {
    stop("Squidpy neighborhood analysis requires `backend_control$runner` or `backend_control$result`.", call. = FALSE)
  }
  graph <- tibble::as_tibble(output$graph %||% output$edges)
  enrichment <- tibble::as_tibble(output$enrichment %||% output$nhood_enrichment)
  cooccurrence <- tibble::as_tibble(output$cooccurrence %||% tibble::tibble())
  if (!all(c("source", "target", "distance") %in% names(graph))) stop("Spatial graph requires source, target, and distance columns.", call. = FALSE)
  graph$source <- as.character(graph$source)
  graph$target <- as.character(graph$target)
  if (anyNA(graph$source) || anyNA(graph$target) ||
      any(!nzchar(graph$source)) || any(!nzchar(graph$target))) {
    stop("Spatial graph endpoints must be non-missing cell identifiers.", call. = FALSE)
  }
  known_cells <- coordinates$table$cell
  if (!all(graph$source %in% known_cells) || !all(graph$target %in% known_cells)) {
    stop("Spatial graph endpoints must all belong to the analyzed object.", call. = FALSE)
  }
  edge_key <- paste(graph$source, graph$target, sep = "\r")
  if (anyDuplicated(edge_key)) {
    stop("Spatial graph must not contain duplicate directed edges.", call. = FALSE)
  }
  raw_distance <- graph$distance
  graph$distance <- suppressWarnings(as.numeric(as.character(raw_distance)))
  if (length(graph$distance) != length(raw_distance) ||
      any(!is.finite(graph$distance)) || any(graph$distance < 0)) {
    stop("Spatial graph distances must be finite non-negative numbers.", call. = FALSE)
  }
  if (!is_null(sample_by) && nrow(graph) > 0L) {
    cross_section <- unname(samples[graph$source]) != unname(samples[graph$target])
    if (anyNA(cross_section) || any(cross_section)) {
      stop("Spatial graph edges must remain within `sample_by` boundaries.", call. = FALSE)
    }
  }
  result <- list(
    schema_version = .sn_analysis_result_schema_version(), analysis_type = "spatial_neighborhood", result_id = result_id,
    method = method, backend = method,
    input = list(group_by = group_by, sample_by = sample_by, coordinate_columns = coordinates$columns, locations = nrow(coordinates$table)),
    parameters = list(k = backend_control$k %||% 6L, n_permutations = backend_control$n_permutations %||% 99L),
    tables = list(primary = enrichment, enrichment = enrichment, cooccurrence = cooccurrence, coordinates = coordinates$table),
    embeddings = list(), graphs = list(spatial = graph), models = list(backend = output$model %||% NULL),
    diagnostics = list(edges = nrow(graph), label_groups = length(unique(labels)), median_distance = stats::median(graph$distance)),
    warnings = as.character(output$warnings %||% character()),
    provenance = .sn_analysis_provenance(random_seed = backend_control$seed %||% 717L)
  )
  sn_validate_result(result)
  object <- sn_store_result(object, "spatial_neighborhood", result_id, result)
  object <- .sn_log_seurat_command(object, name = "sn_run_spatial_neighborhood")
  if (isTRUE(return_object)) object else sn_get_result(object, "spatial_neighborhood", result_id)
}

.sn_spatial_group_distance_by_sample <- function(coordinates, labels, samples = NULL) {
  xy <- as.matrix(coordinates[, c("spatial_x", "spatial_y")])
  rownames(xy) <- coordinates$cell
  samples <- samples %||% stats::setNames(rep("__single_section__", length(labels)), names(labels))
  samples <- as.character(samples[names(labels)])
  dplyr::bind_rows(lapply(unique(samples), function(current_sample) {
    sample_cells <- names(labels)[samples == current_sample]
    sample_labels <- labels[sample_cells]
    groups <- split(names(sample_labels), sample_labels)
    dplyr::bind_rows(lapply(names(groups), function(source) {
      dplyr::bind_rows(lapply(names(groups), function(target) {
        source_cells <- intersect(groups[[source]], rownames(xy))
        target_cells <- intersect(groups[[target]], rownames(xy))
        nearest <- vapply(source_cells, function(source_cell) {
          candidates <- if (identical(source, target)) setdiff(target_cells, source_cell) else target_cells
          if (length(candidates) == 0L) return(NA_real_)
          delta <- sweep(xy[candidates, , drop = FALSE], 2, xy[source_cell, ], "-")
          min(sqrt(rowSums(delta ^ 2)))
        }, numeric(1))
        tibble::tibble(
          source = source, target = target, spatial_sample = current_sample,
          spatial_distance = if (any(is.finite(nearest))) mean(nearest[is.finite(nearest)]) else NA_real_,
          source_locations = length(source_cells), target_locations = length(target_cells),
          contributing_source_locations = sum(is.finite(nearest))
        )
      }))
    }))
  }))
}

.sn_spatial_group_distance <- function(coordinates, labels, samples = NULL) {
  by_sample <- .sn_spatial_group_distance_by_sample(coordinates, labels, samples)
  groups <- split(by_sample, paste(by_sample$source, by_sample$target, sep = "\r"))
  dplyr::bind_rows(lapply(groups, function(current) {
    finite <- is.finite(current$spatial_distance) & current$contributing_source_locations > 0L
    tibble::tibble(
      source = current$source[[1]], target = current$target[[1]],
      spatial_distance = if (any(finite)) stats::weighted.mean(
        current$spatial_distance[finite], current$contributing_source_locations[finite]
      ) else NA_real_,
      source_locations = sum(current$source_locations),
      target_locations = sum(current$target_locations),
      contributing_source_locations = sum(current$contributing_source_locations),
      contributing_samples = sum(finite)
    )
  }))
}

#' Add spatial distance evidence to a communication result
#'
#' @param object A Seurat object.
#' @param source_result_id Stored communication result name.
#' @param communication Optional communication result supplied directly.
#' @param group_by Metadata column matching communication source/target labels.
#' @param spatial_cols Coordinate metadata columns.
#' @param sample_by Optional metadata column defining independent samples or
#'   tissue sections. Distances are calculated within sections and then
#'   aggregated, never between sections. When communication rows contain a
#'   non-missing `sample` column, distances are matched by source, target, and
#'   sample instead of using the cross-sample aggregate.
#' @param max_distance Optional finite non-negative maximum mean nearest-group
#'   distance.
#' @param result_id Stored result name.
#' @param return_object Return the modified object or result.
#' @return A Seurat object or spatial-communication result.
#' @export
sn_run_spatial_communication <- function(object,
                                         source_result_id = "communication",
                                         communication = NULL,
                                         group_by,
                                         spatial_cols = NULL,
                                         max_distance = NULL,
                                         result_id = "spatial_communication",
                                         return_object = TRUE,
                                         sample_by = NULL) {
  result_id <- .sn_validate_result_id(result_id)
  .sn_validate_seurat_object(object)
  if (missing(group_by) || !group_by %in% colnames(object[[]])) stop("`group_by` must name object metadata.", call. = FALSE)
  if (!is_null(max_distance) &&
      (!is.numeric(max_distance) || length(max_distance) != 1L ||
       is.na(max_distance) || !is.finite(max_distance) || max_distance < 0)) {
    stop("`max_distance` must be NULL or one finite non-negative number.", call. = FALSE)
  }
  communication <- communication %||% sn_get_result(object, "cell_communication", source_result_id)
  sn_validate_result(communication)
  interactions <- tibble::as_tibble(communication$tables$primary)
  if (!all(c("source", "target") %in% names(interactions))) stop("Communication result requires source and target columns.", call. = FALSE)
  coordinates <- .sn_spatial_coordinates(object, spatial_cols, sample_by = sample_by)
  labels <- stats::setNames(as.character(object[[group_by, drop = TRUE]]), colnames(object))
  samples <- stats::setNames(coordinates$table$spatial_sample, coordinates$table$cell)
  distances_by_sample <- .sn_spatial_group_distance_by_sample(coordinates$table, labels, samples)
  distances <- .sn_spatial_group_distance(coordinates$table, labels, samples)
  interaction_samples <- if ("sample" %in% names(interactions)) as.character(interactions$sample) else character()
  has_sample <- length(interaction_samples) > 0L & !is.na(interaction_samples) & nzchar(interaction_samples)
  if (any(has_sample) && !all(has_sample)) {
    stop("Communication interactions must either all provide `sample` or all omit it for spatial distance matching.", call. = FALSE)
  }
  sample_aware_join <- !is_null(sample_by) && length(has_sample) > 0L && all(has_sample)
  if (isTRUE(sample_aware_join)) {
    missing_samples <- setdiff(unique(interaction_samples), unique(as.character(distances_by_sample$spatial_sample)))
    if (length(missing_samples) > 0L) {
      stop(
        "Communication interaction sample(s) were not found in spatial metadata: ",
        paste(missing_samples, collapse = ", "), ".",
        call. = FALSE
      )
    }
    distances_for_join <- distances_by_sample
    names(distances_for_join)[names(distances_for_join) == "spatial_sample"] <- "sample"
    spatial <- dplyr::left_join(interactions, distances_for_join, by = c("source", "target", "sample"))
  } else {
    spatial <- dplyr::left_join(interactions, distances, by = c("source", "target"))
  }
  spatial$within_distance <- if (is_null(max_distance)) is.finite(spatial$spatial_distance) else is.finite(spatial$spatial_distance) & spatial$spatial_distance <= max_distance
  filtered <- spatial[spatial$within_distance, , drop = FALSE]
  result <- list(
    schema_version = .sn_analysis_result_schema_version(), analysis_type = "spatial_communication", result_id = result_id,
    method = paste0(communication$method, "+distance"), backend = communication$backend,
    input = list(source_result_id = communication$result_id, group_by = group_by, sample_by = sample_by, coordinate_columns = coordinates$columns, locations = nrow(coordinates$table)),
    parameters = list(max_distance = max_distance),
    tables = list(primary = filtered, all_interactions = spatial, group_distances = distances, group_distances_by_sample = distances_by_sample, coordinates = coordinates$table),
    embeddings = list(), graphs = list(), models = list(source_result = list(result_id = communication$result_id, method = communication$method)),
    diagnostics = list(
      input_interactions = nrow(interactions), retained_interactions = nrow(filtered),
      distance_groups = nrow(distances), distance_matching = if (isTRUE(sample_aware_join)) "source_target_sample" else "source_target_aggregate"
    ),
    warnings = character(), provenance = .sn_analysis_provenance()
  )
  sn_validate_result(result)
  object <- sn_store_result(object, "spatial_communication", result_id, result)
  object <- .sn_log_seurat_command(object, name = "sn_run_spatial_communication")
  if (isTRUE(return_object)) object else sn_get_result(object, "spatial_communication", result_id)
}

#' Run spatial deconvolution through cell2location
#'
#' `sn_run_spatial_deconvolution()` is a deprecated compatibility alias. Use
#' [sn_run_cell2location()] directly.
#' @inheritParams sn_run_cell2location
#' @export
sn_run_spatial_deconvolution <- function(...) {
  .Deprecated("sn_run_cell2location", package = "Shennong")
  sn_run_cell2location(...)
}

#' Map single cells to space through Tangram
#'
#' `sn_run_spatial_mapping()` is a deprecated compatibility alias. Use
#' [sn_run_tangram()] directly.
#' @inheritParams sn_run_tangram
#' @export
sn_run_spatial_mapping <- function(...) {
  .Deprecated("sn_run_tangram", package = "Shennong")
  sn_run_tangram(...)
}

#' Integrate spatial samples with an explicit backend adapter
#'
#' @param object A Seurat object.
#' @param method Integration backend label.
#' @param spatial_cols Coordinate columns.
#' @param result_id Stored result name.
#' @param backend_control Required `runner` or `result` returning an embedding
#'   table with a `cell` column.
#' @param return_object Return the modified object or result.
#' @return A Seurat object or spatial-integration result.
#' @export
sn_integrate_spatial <- function(object,
                                 method = c("staligner", "harmony", "custom"),
                                 spatial_cols = NULL,
                                 result_id = "spatial_integration",
                                 backend_control = list(),
                                 return_object = TRUE) {
  result_id <- .sn_validate_result_id(result_id)
  .sn_validate_seurat_object(object)
  method <- match.arg(method)
  coordinates <- .sn_spatial_coordinates(object, spatial_cols)
  output <- if (is.function(backend_control$runner)) backend_control$runner(object = object, method = method, coordinates = coordinates$table, backend_control = backend_control) else backend_control$result
  if (is_null(output) || !is.list(output)) stop("Spatial integration requires `backend_control$runner` or `backend_control$result`.", call. = FALSE)
  embedding <- tibble::as_tibble(output$embedding %||% output$latent)
  if (!"cell" %in% names(embedding) || ncol(embedding) < 3L) stop("Spatial integration embedding requires cell plus at least two dimensions.", call. = FALSE)
  embedding$cell <- as.character(embedding$cell)
  expected_cells <- coordinates$table$cell
  if (anyNA(embedding$cell) || any(!nzchar(embedding$cell)) || anyDuplicated(embedding$cell) ||
      !setequal(embedding$cell, expected_cells)) {
    stop("Spatial integration embedding cells must uniquely match every analyzed object cell.", call. = FALSE)
  }
  embedding <- embedding[match(expected_cells, embedding$cell), , drop = FALSE]
  dimensions <- setdiff(names(embedding), "cell")
  for (dimension in dimensions) {
    values <- suppressWarnings(as.numeric(as.character(embedding[[dimension]])))
    if (length(values) != nrow(embedding) || any(!is.finite(values))) {
      stop(
        "Spatial integration dimension '", dimension,
        "' must contain only finite numeric values.",
        call. = FALSE
      )
    }
    embedding[[dimension]] <- values
  }
  matrix <- as.matrix(embedding[, dimensions, drop = FALSE])
  rownames(matrix) <- embedding$cell
  result <- list(
    schema_version = .sn_analysis_result_schema_version(), analysis_type = "spatial_integration", result_id = result_id,
    method = method, backend = method,
    input = list(coordinate_columns = coordinates$columns, locations = nrow(coordinates$table)),
    parameters = list(), tables = list(primary = embedding, coordinates = coordinates$table),
    embeddings = list(integrated = matrix), graphs = list(), models = list(backend = output$model %||% NULL),
    diagnostics = list(dimensions = ncol(matrix), locations = nrow(matrix)),
    warnings = as.character(output$warnings %||% character()), provenance = .sn_analysis_provenance(random_seed = backend_control$seed %||% NA_integer_)
  )
  sn_validate_result(result)
  object <- sn_store_result(object, "spatial_integration", result_id, result)
  object <- .sn_log_seurat_command(object, name = "sn_integrate_spatial")
  if (isTRUE(return_object)) object else sn_get_result(object, "spatial_integration", result_id)
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
    controls <- list(...)
    coordinates <- .sn_spatial_coordinates(
      object,
      controls$spatial_cols %||% NULL,
      sample_by = controls$sample_by %||% NULL
    )
    graph <- .sn_spatial_knn_graph(coordinates$table, controls$k %||% 6L)
    return(list(
      coordinates = coordinates$table,
      graph = graph,
      diagnostics = list(
        locations = nrow(coordinates$table),
        samples = length(unique(coordinates$table$spatial_sample)),
        median_neighbor_distance = stats::median(graph$distance)
      )
    ))
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
#' @rdname sn_run_scarches
#' @export
sn_run_cell2location <- function(object,
                                 assay = NULL,
                                 layer = "counts",
                                 reference_signatures = NULL,
                                 spatial_cols = NULL,
                                 output_dir = NULL,
                                 runtime_dir = NULL,
                                 metadata_prefix = "cell2location_",
                                 artifact_id = "cell2location",
                                 return_object = TRUE,
                                 method_control = list(),
                                 keep_run_dir = NULL,
                                 max_artifact_import_gb = 0.5,
                                 ...) {
  .sn_run_python_object_method(
    object = object,
    environment = "cell2location",
    script_name = "cell2location_run.py",
    method = "cell2location",
    assay = assay,
    layer = layer,
    spatial_cols = spatial_cols,
    output_dir = output_dir,
    runtime_dir = runtime_dir,
    metadata_prefix = metadata_prefix,
    result_name = artifact_id,
    return_object = return_object,
    keep_run_dir = keep_run_dir,
    max_artifact_import_gb = max_artifact_import_gb,
    config = c(list(reference_signatures = reference_signatures), method_control),
    ...
  )
}

#' @rdname sn_run_scarches
#' @export
sn_run_tangram <- function(object,
                           reference_object = NULL,
                           assay = NULL,
                           layer = NULL,
                           reference_assay = NULL,
                           reference_layer = NULL,
                           spatial_cols = NULL,
                           cell_type_by = NULL,
                           output_dir = NULL,
                           runtime_dir = NULL,
                           metadata_prefix = "tangram_",
                           artifact_id = "tangram",
                           return_object = TRUE,
                           method_control = list(),
                           keep_run_dir = NULL,
                           max_artifact_import_gb = 0.5,
                           ...) {
  if (is.null(reference_object)) {
    stop("`reference_object` is required for `sn_run_tangram()`.", call. = FALSE)
  }
  .sn_run_python_object_method(
    object = object,
    reference_object = reference_object,
    environment = "tangram",
    script_name = "tangram_run.py",
    method = "tangram",
    assay = assay,
    layer = layer,
    reference_assay = reference_assay,
    reference_layer = reference_layer,
    spatial_cols = spatial_cols,
    output_dir = output_dir,
    runtime_dir = runtime_dir,
    metadata_prefix = metadata_prefix,
    result_name = artifact_id,
    return_object = return_object,
    keep_run_dir = keep_run_dir,
    max_artifact_import_gb = max_artifact_import_gb,
    config = c(list(cell_type_key = cell_type_by), method_control),
    ...
  )
}

#' @rdname sn_run_scarches
#' @export
sn_run_squidpy <- function(object,
                           assay = NULL,
                           layer = NULL,
                           spatial_cols = NULL,
                           cluster_by = NULL,
                           output_dir = NULL,
                           runtime_dir = NULL,
                           metadata_prefix = "squidpy_",
                           artifact_id = "squidpy",
                           return_object = TRUE,
                           method_control = list(),
                           keep_run_dir = NULL,
                           max_artifact_import_gb = 0.5,
                           ...) {
  .sn_run_python_object_method(
    object = object,
    environment = "squidpy",
    script_name = "squidpy_run.py",
    method = "squidpy",
    assay = assay,
    layer = layer,
    spatial_cols = spatial_cols,
    output_dir = output_dir,
    runtime_dir = runtime_dir,
    metadata_prefix = metadata_prefix,
    result_name = artifact_id,
    return_object = return_object,
    keep_run_dir = keep_run_dir,
    max_artifact_import_gb = max_artifact_import_gb,
    config = c(list(cluster_key = cluster_by), method_control),
    ...
  )
}

#' @rdname sn_run_scarches
#' @export
sn_run_spatialdata <- function(object,
                               assay = NULL,
                               layer = NULL,
                               spatial_cols = NULL,
                               output_dir = NULL,
                               runtime_dir = NULL,
                               metadata_prefix = "spatialdata_",
                               artifact_id = "spatialdata",
                               return_object = TRUE,
                               method_control = list(),
                               keep_run_dir = NULL,
                               max_artifact_import_gb = 0.5,
                               ...) {
  .sn_run_python_object_method(
    object = object,
    environment = "spatialdata",
    script_name = "spatialdata_run.py",
    method = "spatialdata",
    assay = assay,
    layer = layer,
    spatial_cols = spatial_cols,
    output_dir = output_dir,
    runtime_dir = runtime_dir,
    metadata_prefix = metadata_prefix,
    result_name = artifact_id,
    return_object = return_object,
    keep_run_dir = keep_run_dir,
    max_artifact_import_gb = max_artifact_import_gb,
    config = method_control,
    ...
  )
}

#' @rdname sn_run_scarches
#' @export
sn_run_stlearn <- function(object,
                           assay = NULL,
                           layer = NULL,
                           spatial_cols = NULL,
                           output_dir = NULL,
                           runtime_dir = NULL,
                           metadata_prefix = "stlearn_",
                           artifact_id = "stlearn",
                           return_object = TRUE,
                           method_control = list(),
                           keep_run_dir = NULL,
                           max_artifact_import_gb = 0.5,
                           ...) {
  .sn_run_python_object_method(
    object = object,
    environment = "stlearn",
    script_name = "stlearn_run.py",
    method = "stlearn",
    assay = assay,
    layer = layer,
    spatial_cols = spatial_cols,
    output_dir = output_dir,
    runtime_dir = runtime_dir,
    metadata_prefix = metadata_prefix,
    result_name = artifact_id,
    return_object = return_object,
    keep_run_dir = keep_run_dir,
    max_artifact_import_gb = max_artifact_import_gb,
    config = method_control,
    ...
  )
}
