.sn_nmf_rank <- function(n_programs, n_cells, n_features) {
  if (is.character(n_programs)) {
    if (!identical(n_programs, "auto")) stop("`n_programs` must be a positive integer or \"auto\".", call. = FALSE)
    return(max(2L, min(10L, as.integer(round(sqrt(min(n_cells, n_features)) / 2)))))
  }
  rank <- as.integer(n_programs)
  if (length(rank) != 1L || !is.finite(rank) || rank < 2L || rank >= min(n_cells, n_features)) {
    stop("`n_programs` must be at least 2 and smaller than both cells and features.", call. = FALSE)
  }
  rank
}

.sn_nmf_fit <- function(matrix, rank, seed, max_iter = 200L, tolerance = 1e-5) {
  set.seed(seed)
  matrix <- as.matrix(matrix)
  matrix[!is.finite(matrix)] <- 0
  if (any(matrix < 0)) stop("NMF requires non-negative expression values.", call. = FALSE)
  epsilon <- .Machine$double.eps
  weights <- matrix(stats::runif(nrow(matrix) * rank, min = 0.01, max = 1), nrow = nrow(matrix), ncol = rank)
  activity <- matrix(stats::runif(rank * ncol(matrix), min = 0.01, max = 1), nrow = rank, ncol = ncol(matrix))
  previous <- Inf
  converged <- FALSE
  iteration <- 0L
  for (iteration in seq_len(as.integer(max_iter))) {
    activity <- activity * (crossprod(weights, matrix) / (crossprod(weights, weights) %*% activity + epsilon))
    weights <- weights * ((matrix %*% t(activity)) / (weights %*% (activity %*% t(activity)) + epsilon))
    scale <- sqrt(colSums(weights ^ 2))
    scale[!is.finite(scale) | scale <= epsilon] <- 1
    weights <- sweep(weights, 2, scale, "/")
    activity <- sweep(activity, 1, scale, "*")
    if (iteration == 1L || iteration %% 10L == 0L || iteration == max_iter) {
      error <- sqrt(sum((matrix - weights %*% activity) ^ 2)) / max(sqrt(sum(matrix ^ 2)), epsilon)
      if (is.finite(previous) && abs(previous - error) <= tolerance * max(1, previous)) {
        converged <- TRUE
        break
      }
      previous <- error
    }
  }
  error <- sqrt(sum((matrix - weights %*% activity) ^ 2)) / max(sqrt(sum(matrix ^ 2)), epsilon)
  rownames(weights) <- rownames(matrix)
  colnames(activity) <- colnames(matrix)
  list(weights = weights, activity = activity, error = error, iterations = iteration, converged = converged, seed = seed)
}

.sn_nmf_stability <- function(best, fits) {
  if (length(fits) < 2L) return(rep(NA_real_, ncol(best$weights)))
  comparisons <- lapply(fits, function(fit) {
    correlation <- suppressWarnings(stats::cor(best$weights, fit$weights, method = "spearman", use = "pairwise.complete.obs"))
    apply(correlation, 1, function(values) {
      values <- values[is.finite(values)]
      if (length(values) == 0L) NA_real_ else max(values)
    })
  })
  rowMeans(do.call(cbind, comparisons), na.rm = TRUE)
}

.sn_program_discovery_tables <- function(fit, program_prefix = "program", top_genes = 50L) {
  program_names <- paste0(program_prefix, "_", seq_len(ncol(fit$weights)))
  colnames(fit$weights) <- program_names
  rownames(fit$activity) <- program_names
  weights <- dplyr::bind_rows(lapply(program_names, function(program) {
    values <- fit$weights[, program]
    order <- order(values, decreasing = TRUE)
    tibble::tibble(
      program = program, feature = rownames(fit$weights)[order],
      weight = as.numeric(values[order]), rank = seq_along(order),
      selected = seq_along(order) <= as.integer(top_genes)
    )
  }))
  activity <- dplyr::bind_rows(lapply(program_names, function(program) {
    values <- as.numeric(fit$activity[program, ])
    tibble::tibble(
      cell = colnames(fit$activity), program = program,
      score = values
    )
  }))
  list(weights = weights, activity = activity, weight_matrix = fit$weights, activity_matrix = fit$activity)
}

.sn_discover_programs_nmf <- function(object,
                                      n_programs,
                                      group_by,
                                      assay,
                                      layer,
                                      features,
                                      backend_control) {
  expression <- .sn_annotation_expression(object, assay = assay, layer = layer)
  matrix <- expression$matrix
  means <- Matrix::rowMeans(matrix)
  variances <- Matrix::rowMeans(matrix ^ 2) - means ^ 2
  features <- features %||% names(utils::head(sort(variances, decreasing = TRUE), backend_control$n_features %||% min(2000L, nrow(matrix))))
  features <- intersect(features, rownames(matrix))
  if (length(features) < 3L) stop("Program discovery requires at least three matched features.", call. = FALSE)
  groups <- if (is_null(group_by)) {
    stats::setNames(rep("all", ncol(object)), colnames(object))
  } else {
    if (!group_by %in% colnames(object[[]])) stop("`group_by` column was not found.", call. = FALSE)
    stats::setNames(as.character(object[[group_by, drop = TRUE]]), colnames(object))
  }
  nrun <- as.integer(backend_control$nrun %||% 5L)
  max_iter <- as.integer(backend_control$max_iter %||% 200L)
  tolerance <- backend_control$tolerance %||% 1e-5
  seed <- as.integer(backend_control$seed %||% 717L)
  max_dense <- backend_control$max_dense_values %||% 5e6
  outputs <- lapply(unique(groups[!is.na(groups) & nzchar(groups)]), function(group) {
    cells <- names(groups)[groups == group]
    current <- matrix[features, cells, drop = FALSE]
    if (length(current) > max_dense) {
      stop("NMF input exceeds `backend_control$max_dense_values`; reduce features/cells or use cNMF.", call. = FALSE)
    }
    rank <- .sn_nmf_rank(n_programs, n_cells = ncol(current), n_features = nrow(current))
    fits <- lapply(seq_len(nrun), function(run) {
      .sn_nmf_fit(current, rank, seed = seed + run - 1L, max_iter = max_iter, tolerance = tolerance)
    })
    errors <- vapply(fits, `[[`, numeric(1), "error")
    best_index <- which.min(errors)
    best <- fits[[best_index]]
    best$stability <- .sn_nmf_stability(best, fits[-best_index])
    prefix <- if (identical(group, "all")) "program" else paste0(gsub("[^[:alnum:]_]+", "_", group), "_program")
    tables <- .sn_program_discovery_tables(best, prefix, top_genes = backend_control$top_genes %||% 50L)
    tables$activity$group <- group
    tables$weights$group <- group
    tables$diagnostics <- tibble::tibble(
      group = group, run = seq_len(nrun), seed = seed + seq_len(nrun) - 1L,
      reconstruction_error = errors,
      converged = vapply(fits, `[[`, logical(1), "converged"),
      iterations = vapply(fits, `[[`, integer(1), "iterations"),
      selected = seq_len(nrun) == best_index
    )
    tables$stability <- tibble::tibble(
      group = group, program = colnames(tables$weight_matrix), stability = best$stability
    )
    tables
  })
  list(
    object = object,
    weights = dplyr::bind_rows(lapply(outputs, `[[`, "weights")),
    activity = dplyr::bind_rows(lapply(outputs, `[[`, "activity")),
    diagnostics = dplyr::bind_rows(lapply(outputs, `[[`, "diagnostics")),
    stability = dplyr::bind_rows(lapply(outputs, `[[`, "stability")),
    expression = expression
  )
}

.sn_standardize_program_discovery <- function(output, object, method) {
  if (!is.list(output)) stop("Program discovery adapter output must be a list.", call. = FALSE)
  weights <- tibble::as_tibble(output$weights %||% output$gene_weights)
  activity <- tibble::as_tibble(output$activity %||% output$scores)
  required_weights <- c("program", "feature", "weight")
  required_activity <- c("cell", "program", "score")
  if (!all(required_weights %in% names(weights))) stop("Program weights require program, feature, and weight columns.", call. = FALSE)
  if (!all(required_activity %in% names(activity))) stop("Program activity requires cell, program, and score columns.", call. = FALSE)
  activity <- activity[activity$cell %in% colnames(object), , drop = FALSE]
  if (!"method" %in% names(weights)) weights$method <- method
  if (!"method" %in% names(activity)) activity$method <- method
  if (!"group" %in% names(weights)) weights$group <- "all"
  if (!"group" %in% names(activity)) activity$group <- "all"
  list(
    object = output$object %||% object, weights = weights, activity = activity,
    diagnostics = tibble::as_tibble(output$diagnostics %||% tibble::tibble()),
    stability = tibble::as_tibble(output$stability %||% tibble::tibble()),
    expression = output$expression %||% list(assay = SeuratObject::DefaultAssay(object), layer = NA_character_)
  )
}

#' Discover latent gene programs
#'
#' @param object A Seurat object.
#' @param method Program-discovery backend. NMF is implemented locally; cNMF
#'   and Hotspot use explicit runner/result adapters.
#' @param n_programs Positive rank or `"auto"`.
#' @param group_by Optional metadata column for independent within-group discovery.
#' @param name Stored result name.
#' @param assay,layer Expression assay and non-negative layer.
#' @param features Optional features used for discovery.
#' @param backend_control NMF controls (`nrun`, `max_iter`, `seed`,
#'   `n_features`, `top_genes`) or an external `runner`/`result`.
#' @param return_object Return the modified object or unified result.
#' @return A Seurat object or unified program-discovery result.
#' @examples
#' \dontrun{
#' object <- sn_discover_programs(object, method = "nmf", n_programs = 6)
#' programs <- sn_get_result(object, "program_discovery", "programs_nmf")
#' programs$tables$gene_weights
#' }
#' @export
sn_discover_programs <- function(object,
                                 method = c("nmf", "cnmf", "hotspot"),
                                 n_programs = "auto",
                                 group_by = NULL,
                                 name = NULL,
                                 assay = NULL,
                                 layer = "data",
                                 features = NULL,
                                 backend_control = list(),
                                 return_object = TRUE) {
  .sn_validate_result_object(object)
  method <- match.arg(method)
  name <- name %||% paste0("programs_", method)
  if (identical(method, "nmf")) {
    output <- .sn_discover_programs_nmf(object, n_programs, group_by, assay, layer, features, backend_control)
  } else if (is.function(backend_control$runner)) {
    output <- backend_control$runner(
      object = object, method = method, n_programs = n_programs,
      group_by = group_by, assay = assay, layer = layer, features = features
    )
  } else if (!is_null(backend_control$result)) {
    output <- backend_control$result
  } else {
    stop("The ", method, " adapter requires `backend_control$runner` or `backend_control$result`.", call. = FALSE)
  }
  output <- .sn_standardize_program_discovery(output, object, method)
  object <- output$object
  activity <- output$activity
  prefix <- gsub("[^[:alnum:]_]+", "_", name)
  metadata <- data.frame(row.names = colnames(object))
  for (program in unique(activity$program)) {
    values <- stats::setNames(activity$score[activity$program == program], activity$cell[activity$program == program])
    metadata[[paste(prefix, gsub("[^[:alnum:]_]+", "_", program), sep = "_")]] <- as.numeric(values[colnames(object)])
  }
  object <- SeuratObject::AddMetaData(object, metadata = metadata)
  result <- list(
    schema_version = "1.0", analysis_type = "program_discovery", name = name,
    method = method, backend = method,
    input = list(
      assay = output$expression$assay %||% assay %||% SeuratObject::DefaultAssay(object),
      layer = output$expression$layer %||% layer, group_by = group_by,
      cells = ncol(object), features = features
    ),
    parameters = list(n_programs = n_programs),
    tables = list(
      primary = activity, activity = activity, gene_weights = output$weights,
      stability = output$stability, fit_diagnostics = output$diagnostics
    ),
    embeddings = list(), graphs = list(), models = list(),
    diagnostics = list(
      programs = length(unique(activity$program)),
      groups = length(unique(activity$group)),
      converged_runs = if ("converged" %in% names(output$diagnostics)) sum(output$diagnostics$converged) else NA_integer_
    ),
    warnings = character(), provenance = .sn_analysis_provenance(random_seed = backend_control$seed %||% NA_integer_)
  )
  sn_validate_result(result)
  object <- sn_store_result(object, "program_discovery", name, result)
  object <- .sn_log_seurat_command(object, assay = result$input$assay, name = "sn_discover_programs")
  if (isTRUE(return_object)) object else sn_get_result(object, "program_discovery", name)
}
