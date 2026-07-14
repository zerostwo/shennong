.sn_binary_auc <- function(labels, scores, positive) {
  labels <- as.character(labels)
  keep <- is.finite(scores) & !is.na(labels)
  labels <- labels[keep]
  scores <- scores[keep]
  positives <- labels == positive
  n_positive <- sum(positives)
  n_negative <- sum(!positives)
  if (n_positive == 0L || n_negative == 0L) return(NA_real_)
  ranks <- rank(scores, ties.method = "average")
  (sum(ranks[positives]) - n_positive * (n_positive + 1) / 2) / (n_positive * n_negative)
}

.sn_priority_folds <- function(labels, samples = NULL, folds = 3L, seed = 717L) {
  if (!is_null(samples)) {
    sample_ids <- unique(as.character(samples))
    return(lapply(sample_ids, function(sample) which(as.character(samples) == sample)))
  }
  set.seed(seed)
  assignments <- integer(length(labels))
  for (level in unique(labels)) {
    indices <- which(labels == level)
    assignments[indices] <- sample(rep(seq_len(folds), length.out = length(indices)))
  }
  split(seq_along(labels), assignments)
}

.sn_centroid_predictions <- function(expression, labels, folds, positive) {
  predictions <- rep(NA_real_, nrow(expression))
  negative <- setdiff(unique(labels), positive)
  if (length(negative) != 1L) return(predictions)
  for (test_indices in folds) {
    train_indices <- setdiff(seq_len(nrow(expression)), test_indices)
    train_labels <- labels[train_indices]
    if (length(unique(train_labels)) != 2L) next
    center <- colMeans(expression[train_indices, , drop = FALSE])
    spread <- apply(expression[train_indices, , drop = FALSE], 2, stats::sd)
    spread[!is.finite(spread) | spread == 0] <- 1
    train <- sweep(sweep(expression[train_indices, , drop = FALSE], 2, center, "-"), 2, spread, "/")
    test <- sweep(sweep(expression[test_indices, , drop = FALSE], 2, center, "-"), 2, spread, "/")
    direction <- colMeans(train[train_labels == positive, , drop = FALSE]) -
      colMeans(train[train_labels == negative, , drop = FALSE])
    raw_train <- as.numeric(train %*% direction)
    midpoint <- mean(c(mean(raw_train[train_labels == positive]), mean(raw_train[train_labels == negative])))
    scale_value <- stats::sd(raw_train)
    if (!is.finite(scale_value) || scale_value == 0) scale_value <- 1
    predictions[test_indices] <- stats::plogis((as.numeric(test %*% direction) - midpoint) / scale_value)
  }
  predictions
}

.sn_priority_augur <- function(object,
                               phenotype,
                               state_by,
                               sample_by,
                               contrast,
                               assay,
                               layer,
                               features,
                               max_features,
                               max_cells_per_state,
                               folds,
                               permutations,
                               seed) {
  metadata <- object[[]]
  required <- c(phenotype, state_by, sample_by)
  missing <- setdiff(required, colnames(metadata))
  if (length(missing) > 0L) stop("Missing state-priority metadata column(s): ", paste(missing, collapse = ", "), ".", call. = FALSE)
  labels_all <- as.character(metadata[[phenotype]])
  if (!is_null(sample_by)) {
    .sn_validate_constant_within_sample(metadata, sample_col = sample_by, group_col = phenotype)
  }
  levels <- unique(labels_all[!is.na(labels_all)])
  if (is_null(contrast)) {
    if (length(levels) != 2L) stop("Binary state prioritization requires `contrast = c(case, control)`.", call. = FALSE)
    contrast <- levels
  }
  if (length(contrast) != 2L) stop("`contrast` must be `c(case, control)`.", call. = FALSE)
  expression <- .sn_annotation_expression(object, assay = assay, layer = layer)
  if (is_null(features)) {
    features <- SeuratObject::VariableFeatures(object[[expression$assay]])
    if (length(features) == 0L) features <- rownames(expression$matrix)
  }
  features <- utils::head(intersect(features, rownames(expression$matrix)), as.integer(max_features))
  if (length(features) < 2L) stop("State prioritization requires at least two expression features.", call. = FALSE)
  states <- unique(as.character(metadata[[state_by]]))
  states <- states[!is.na(states)]
  rows <- list()
  cell_rows <- list()
  sample_rows <- list()
  null_rows <- list()
  set.seed(seed)
  for (state in states) {
    indices <- which(as.character(metadata[[state_by]]) == state & labels_all %in% contrast)
    if (length(indices) > max_cells_per_state) {
      by_label <- split(indices, labels_all[indices])
      indices <- unlist(lapply(by_label, function(current) sample(current, min(length(current), floor(max_cells_per_state / 2)))), use.names = FALSE)
    }
    labels <- labels_all[indices]
    if (any(table(factor(labels, levels = contrast)) < 10L)) next
    samples <- if (is_null(sample_by)) NULL else as.character(metadata[[sample_by]][indices])
    matrix <- t(as.matrix(expression$matrix[features, indices, drop = FALSE]))
    fold_list <- .sn_priority_folds(labels, samples = samples, folds = folds, seed = seed)
    prediction <- .sn_centroid_predictions(matrix, labels, fold_list, positive = contrast[[1]])
    auc <- .sn_binary_auc(labels, prediction, positive = contrast[[1]])
    permutations <- as.integer(permutations)
    if (length(permutations) != 1L || is.na(permutations) || permutations < 9L) {
      stop("`permutations` must be a single integer of at least 9.", call. = FALSE)
    }
    permutation_auc <- rep(NA_real_, permutations)
    for (iteration in seq_len(permutations)) {
      permuted <- if (is_null(samples)) {
        sample(labels)
      } else {
        sample_labels <- tapply(labels, samples, function(value) value[[1]])
        mapping <- stats::setNames(sample(unname(sample_labels)), names(sample_labels))
        unname(mapping[samples])
      }
      permutation_auc[[iteration]] <- .sn_binary_auc(
        permuted,
        .sn_centroid_predictions(matrix, permuted, fold_list, positive = contrast[[1]]),
        positive = contrast[[1]]
      )
    }
    empirical <- (1 + sum(permutation_auc >= auc, na.rm = TRUE)) / (1 + sum(is.finite(permutation_auc)))
    rows[[state]] <- tibble::tibble(
      state = state,
      priority_score = auc,
      auc = auc,
      phenotype_association = auc - 0.5,
      uncertainty = stats::sd(permutation_auc, na.rm = TRUE),
      null_mean = mean(permutation_auc, na.rm = TRUE),
      p_value = empirical,
      n_cells = length(indices),
      n_samples = if (is_null(samples)) NA_integer_ else length(unique(samples))
    )
    cell_rows[[state]] <- tibble::tibble(
      cell = colnames(object)[indices],
      state = state,
      phenotype = labels,
      predicted_probability = prediction
    )
    if (!is_null(samples)) {
      grouped <- split(seq_along(samples), samples)
      sample_rows[[state]] <- dplyr::bind_rows(lapply(names(grouped), function(sample) {
        current <- grouped[[sample]]
        tibble::tibble(
          sample = sample,
          state = state,
          phenotype = labels[current[[1]]],
          n_cells = length(current),
          contribution = mean(prediction[current], na.rm = TRUE)
        )
      }))
    }
    null_rows[[state]] <- tibble::tibble(state = state, permutation = seq_len(permutations), null_score = permutation_auc)
  }
  ranking <- dplyr::bind_rows(rows)
  if (nrow(ranking) == 0L) stop("No state retained enough cells from both phenotype groups.", call. = FALSE)
  ranking$adjusted_p_value <- stats::p.adjust(ranking$p_value, method = "BH")
  ranking <- ranking[order(ranking$priority_score, decreasing = TRUE), , drop = FALSE]
  list(
    ranking = ranking,
    cells = dplyr::bind_rows(cell_rows),
    samples = dplyr::bind_rows(sample_rows),
    null = dplyr::bind_rows(null_rows),
    assay = expression$assay,
    layer = expression$layer,
    contrast = contrast,
    features = features,
    warnings = character()
  )
}

.sn_priority_rareq <- function(object,
                               phenotype,
                               sample_by,
                               contrast,
                               reduction,
                               dims,
                               backend_control) {
  check_installed("RareQ", reason = "to run RareQ state discovery.")
  if (is_null(sample_by)) stop("RareQ phenotype prioritization requires `sample_by` for sample-level association.", call. = FALSE)
  working <- object
  available_dims <- ncol(Seurat::Embeddings(working, reduction = reduction))
  dims <- dims[dims >= 1L & dims <= available_dims]
  if (length(dims) < 2L) stop("RareQ requires at least two available reduction dimensions.", call. = FALSE)
  if (length(working@neighbors) == 0L) {
    working <- Seurat::FindNeighbors(
      working,
      reduction = reduction,
      dims = dims,
      k.param = backend_control$k.param %||% 20L,
      compute.SNN = FALSE,
      return.neighbor = TRUE,
      verbose = FALSE
    )
  }
  clusters <- do.call(RareQ::FindRare, utils::modifyList(
    list(sc_object = working, assay = Seurat::DefaultAssay(working)),
    backend_control$find_rare %||% list(), keep.null = TRUE
  ))
  working$sn_rareq_state <- as.character(clusters)
  association <- sn_test_abundance(
    working,
    method = "permutation",
    sample_by = sample_by,
    condition_by = phenotype,
    cell_type_by = "sn_rareq_state",
    contrast = contrast,
    permutations = backend_control$permutations %||% 999L,
    seed = backend_control$seed %||% 717L,
    return_object = FALSE
  )
  ranking <- association$tables$primary
  fractions <- prop.table(table(clusters))
  ranking$state <- ranking$feature
  ranking$state_fraction <- unname(fractions[ranking$state])
  ranking$phenotype_association <- ranking$estimate
  ranking$priority_score <- abs(ranking$estimate) * (1 - ranking$state_fraction)
  ranking$uncertainty <- ranking$null_sd
  list(
    ranking = ranking[order(ranking$priority_score, decreasing = TRUE), , drop = FALSE],
    cells = tibble::tibble(cell = colnames(object), state = as.character(clusters)),
    samples = association$tables$sample_contributions,
    null = association$tables$permutation_null,
    contrast = association$parameters$contrast,
    warnings = character()
  )
}

.sn_priority_scissor <- function(object,
                                 state_by,
                                 sample_by,
                                 bulk_expression,
                                 bulk_phenotype,
                                 family,
                                 backend_control) {
  check_installed("Scissor", reason = "to run phenotype-guided Scissor selection.")
  if (is_null(bulk_expression) || is_null(bulk_phenotype)) {
    stop("Scissor requires `bulk_expression` and `bulk_phenotype`; a cell metadata phenotype is not a substitute.", call. = FALSE)
  }
  bulk_expression <- as.matrix(bulk_expression)
  if (ncol(bulk_expression) != NROW(bulk_phenotype)) stop("`bulk_phenotype` must describe every bulk-expression column.", call. = FALSE)
  assay <- backend_control$assay %||% Seurat::DefaultAssay(object)
  expression <- .sn_annotation_expression(object, assay = assay, layer = backend_control$layer %||% "data")
  counts <- .sn_get_seurat_layer_data(object, assay = assay, layer = "counts")
  working <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = object[[]])
  suppressWarnings(working[["RNA"]] <- SeuratObject::CreateAssayObject(counts = counts))
  working <- Seurat::NormalizeData(working, verbose = FALSE)
  working <- Seurat::FindVariableFeatures(
    working,
    nfeatures = min(backend_control$nfeatures %||% 2000L, nrow(working)),
    verbose = FALSE
  )
  working <- Seurat::ScaleData(working, verbose = FALSE)
  npcs <- min(backend_control$npcs %||% 20L, length(SeuratObject::VariableFeatures(working)), ncol(working) - 1L)
  if (npcs < 2L) stop("Scissor requires at least two usable PCA dimensions.", call. = FALSE)
  working <- suppressWarnings(Seurat::RunPCA(working, npcs = npcs, verbose = FALSE))
  working <- Seurat::FindNeighbors(working, dims = seq_len(min(10L, npcs)), verbose = FALSE)
  save_file <- tempfile("shennong-scissor-", fileext = ".RData")
  on.exit(unlink(save_file), add = TRUE)
  defaults <- list(
    bulk_dataset = bulk_expression,
    sc_dataset = working,
    phenotype = bulk_phenotype,
    tag = backend_control$tag %||% names(table(bulk_phenotype)),
    alpha = backend_control$alpha %||% 0.05,
    cutoff = backend_control$cutoff %||% 0.2,
    family = family,
    Save_file = save_file
  )
  backend_warnings <- character()
  fit <- withCallingHandlers(
    do.call(Scissor::Scissor, defaults),
    warning = function(warning) {
      backend_warnings <<- c(backend_warnings, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )
  coefficients <- tibble::tibble(
    cell = colnames(expression$matrix),
    coefficient = as.numeric(fit$Coefs),
    selection = ifelse(fit$Coefs > 0, "Scissor+", ifelse(fit$Coefs < 0, "Scissor-", "Unselected"))
  )
  metadata <- object[[]]
  coefficients$state <- as.character(metadata[[state_by]][match(coefficients$cell, rownames(metadata))])
  if (!is_null(sample_by)) coefficients$sample <- as.character(metadata[[sample_by]][match(coefficients$cell, rownames(metadata))])
  ranking <- dplyr::bind_rows(lapply(split(coefficients, coefficients$state), function(current) {
    tibble::tibble(
      state = current$state[[1]],
      priority_score = mean(abs(current$coefficient)),
      phenotype_association = mean(current$coefficient),
      selected_fraction = mean(current$selection != "Unselected"),
      uncertainty = stats::sd(current$coefficient) / sqrt(nrow(current)),
      n_cells = nrow(current)
    )
  }))
  ranking <- ranking[order(ranking$priority_score, decreasing = TRUE), , drop = FALSE]
  samples <- if (is_null(sample_by)) tibble::tibble() else dplyr::bind_rows(lapply(
    split(coefficients, interaction(coefficients$state, coefficients$sample, drop = TRUE)),
    function(current) tibble::tibble(sample = current$sample[[1]], state = current$state[[1]], contribution = mean(current$coefficient), n_cells = nrow(current))
  ))
  list(
    ranking = ranking,
    cells = coefficients,
    samples = samples,
    null = tibble::tibble(),
    contrast = NULL,
    warnings = c(
      backend_warnings,
      "Scissor does not return a phenotype-permutation null by default; inspect coefficient uncertainty and run its reliability workflow for confirmatory analysis."
    ),
    model = fit$para
  )
}

#' Prioritize phenotype-responsive or rare cell states
#'
#' @param object A Seurat object.
#' @param method State-priority backend.
#' @param phenotype Cell metadata phenotype for Augur/RareQ, or a descriptive
#'   label when Scissor receives explicit bulk inputs.
#' @param sample_by Optional biological sample column. Strongly recommended and
#'   required for RareQ phenotype association.
#' @param state_by Existing state/cell-type column for Augur and Scissor.
#' @param contrast Binary phenotype labels ordered as `c(case, control)`.
#' @param assay,layer Expression source for Augur.
#' @param features,max_features Features used for separability scoring.
#' @param max_cells_per_state Maximum balanced cells per state.
#' @param folds Cross-validation folds when samples are unavailable.
#' @param permutations Number of label permutations.
#' @param reduction,dims Reduction used to construct RareQ neighbors.
#' @param bulk_expression Gene-by-bulk-sample matrix required by Scissor.
#' @param bulk_phenotype Bulk phenotype vector or survival matrix for Scissor.
#' @param family Scissor phenotype family.
#' @param store_name Stored result name.
#' @param seed Random seed.
#' @param backend_control Backend-specific options.
#' @param return_object Return the updated object instead of the result.
#' @return A Seurat object or unified state-priority result.
#' @examples
#' \dontrun{
#' object <- sn_prioritize_states(
#'   object, phenotype = "condition", sample_by = "sample",
#'   state_by = "cell_type", contrast = c("treated", "control")
#' )
#' sn_get_result(object, "state_priority", "priority")
#' }
#' @export
sn_prioritize_states <- function(object,
                                 method = c("augur", "scissor", "rareq"),
                                 phenotype,
                                 sample_by = NULL,
                                 state_by = NULL,
                                 contrast = NULL,
                                 assay = NULL,
                                 layer = "data",
                                 features = NULL,
                                 max_features = 500L,
                                 max_cells_per_state = 500L,
                                 folds = 3L,
                                 permutations = 100L,
                                 reduction = "pca",
                                 dims = 1:20,
                                 bulk_expression = NULL,
                                 bulk_phenotype = NULL,
                                 family = c("binomial", "gaussian", "cox"),
                                 store_name = "priority",
                                 seed = 717L,
                                 backend_control = list(),
                                 return_object = TRUE) {
  .sn_validate_result_object(object)
  method <- match.arg(method)
  family <- match.arg(family)
  if (!identical(method, "rareq") && (is_null(state_by) || !state_by %in% colnames(object[[]]))) {
    stop("`state_by` must name object metadata for Augur or Scissor prioritization.", call. = FALSE)
  }
  backend <- switch(
    method,
    augur = .sn_priority_augur(
      object, phenotype, state_by, sample_by, contrast, assay, layer, features,
      max_features, max_cells_per_state, folds, as.integer(permutations), seed
    ),
    rareq = .sn_priority_rareq(object, phenotype, sample_by, contrast, reduction, dims, backend_control$rareq %||% list()),
    scissor = .sn_priority_scissor(object, state_by, sample_by, bulk_expression, bulk_phenotype, family, backend_control$scissor %||% list())
  )
  result <- list(
    schema_version = "1.0",
    analysis_type = "state_priority",
    name = store_name,
    method = method,
    backend = switch(method, augur = "Shennong sample-aware Augur", rareq = "RareQ", scissor = "Scissor"),
    input = list(cells = ncol(object), phenotype = phenotype, sample_by = sample_by, state_by = state_by),
    parameters = list(contrast = backend$contrast, permutations = as.integer(permutations), family = if (method == "scissor") family else NULL),
    tables = list(primary = backend$ranking, cells = backend$cells, sample_contributions = backend$samples, permutation_null = backend$null),
    embeddings = list(),
    graphs = list(),
    models = list(summary = backend$model %||% list()),
    diagnostics = list(states = nrow(backend$ranking), sample_aware = !is_null(sample_by)),
    warnings = as.character(backend$warnings),
    provenance = .sn_analysis_provenance(random_seed = seed)
  )
  sn_validate_result(result)
  object <- sn_store_result(object, "state_priority", store_name, result)
  object <- .sn_log_seurat_command(object = object, assay = assay, name = "sn_prioritize_states")
  if (isTRUE(return_object)) object else sn_get_result(object, "state_priority", store_name)
}
