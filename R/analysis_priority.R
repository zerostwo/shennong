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
  .sn_priority_scissor_impl(
    object = object,
    state_by = state_by,
    sample_by = sample_by,
    bulk_expression = bulk_expression,
    bulk_phenotype = bulk_phenotype,
    family = family,
    backend_control = backend_control
  )
}

.sn_scissor_backend_control <- function(backend_control, assay = NULL, layer = NULL) {
  if (is_null(backend_control)) backend_control <- list()
  if (!is.list(backend_control)) {
    stop("`backend_control` must be a list.", call. = FALSE)
  }
  control <- if (is.list(backend_control$scissor)) backend_control$scissor else backend_control
  if (!is_null(assay) && is_null(control$assay)) control$assay <- assay
  if (!is_null(layer) && is_null(control$layer)) control$layer <- layer
  control
}

.sn_scissor_validate_control <- function(control) {
  alpha <- control$alpha
  if (!is_null(alpha) && (!is.numeric(alpha) || length(alpha) < 1L ||
      any(!is.finite(alpha)) || any(alpha <= 0 | alpha > 1))) {
    stop("Scissor `alpha` must contain finite values in (0, 1].", call. = FALSE)
  }
  cutoff <- control$cutoff %||% 0.2
  if (!is.numeric(cutoff) || length(cutoff) != 1L || !is.finite(cutoff) ||
      cutoff < 0 || cutoff > 1) {
    stop("Scissor `cutoff` must be one finite value in [0, 1].", call. = FALSE)
  }
  for (field in c("nfeatures", "npcs", "min_shared_features", "max_correlation_rows")) {
    value <- control[[field]]
    if (!is_null(value) && (!is.numeric(value) || length(value) != 1L ||
        !is.finite(value) || value < 1 || value != as.integer(value))) {
      stop(sprintf("Scissor `%s` must be one positive integer.", field), call. = FALSE)
    }
  }
  for (field in c("reliability", "retain_correlation_matrix")) {
    value <- control[[field]]
    if (!is_null(value) && (!is.logical(value) || length(value) != 1L || is.na(value))) {
      stop(sprintf("Scissor `%s` must be TRUE or FALSE.", field), call. = FALSE)
    }
  }
  for (field in c("runner", "reliability_runner")) {
    value <- control[[field]]
    if (!is_null(value) && !is.function(value)) {
      stop(sprintf("Scissor `%s` must be a function.", field), call. = FALSE)
    }
  }
  reliability_control <- control$reliability_control
  if (!is_null(reliability_control) && !is.list(reliability_control)) {
    stop("Scissor `reliability_control` must be a list.", call. = FALSE)
  }
  if (!is_null(reliability_control)) {
    fdr_cutoff <- reliability_control$fdr_cutoff
    if (!is_null(fdr_cutoff) && (!is.numeric(fdr_cutoff) || length(fdr_cutoff) != 1L ||
        !is.finite(fdr_cutoff) || fdr_cutoff < 0 || fdr_cutoff > 1)) {
      stop("Scissor reliability `fdr_cutoff` must be one finite value in [0, 1].", call. = FALSE)
    }
    bootstrap_n <- reliability_control$bootstrap_n
    if (!is_null(bootstrap_n) && (!is.numeric(bootstrap_n) || length(bootstrap_n) != 1L ||
        !is.finite(bootstrap_n) || bootstrap_n < 1 || bootstrap_n != as.integer(bootstrap_n))) {
      stop("Scissor reliability `bootstrap_n` must be one positive integer.", call. = FALSE)
    }
  }
  invisible(control)
}

.sn_scissor_validate_seed <- function(seed) {
  if (!is.numeric(seed) || length(seed) != 1L || !is.finite(seed) ||
      seed < 0 || seed > .Machine$integer.max || seed != as.integer(seed)) {
    stop("Scissor `seed` must be one non-negative integer.", call. = FALSE)
  }
  as.integer(seed)
}

.sn_scissor_validate_bulk <- function(bulk_expression) {
  if (is.data.frame(bulk_expression) &&
      !all(vapply(bulk_expression, is.numeric, logical(1)))) {
    stop("`bulk_expression` must contain only numeric values.", call. = FALSE)
  }
  bulk_expression <- as.matrix(bulk_expression)
  if (!is.numeric(bulk_expression) || length(dim(bulk_expression)) != 2L) {
    stop("`bulk_expression` must be a numeric gene-by-sample matrix.", call. = FALSE)
  }
  if (nrow(bulk_expression) < 2L || ncol(bulk_expression) < 2L) {
    stop("`bulk_expression` must contain at least two genes and two samples.", call. = FALSE)
  }
  if (any(!is.finite(bulk_expression))) {
    stop("`bulk_expression` must not contain missing or non-finite values.", call. = FALSE)
  }
  genes <- rownames(bulk_expression)
  samples <- colnames(bulk_expression)
  if (is_null(genes) || anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes)) {
    stop("`bulk_expression` must have unique, non-empty gene row names.", call. = FALSE)
  }
  if (is_null(samples) || anyNA(samples) || any(!nzchar(samples)) || anyDuplicated(samples)) {
    stop("`bulk_expression` must have unique, non-empty sample column names.", call. = FALSE)
  }
  storage.mode(bulk_expression) <- "double"
  bulk_expression
}

.sn_scissor_has_explicit_row_names <- function(x) {
  rows <- rownames(x)
  if (is_null(rows)) return(FALSE)
  if (is.data.frame(x)) return(.row_names_info(x, type = 1L) >= 0L)
  TRUE
}

.sn_scissor_align_vector <- function(values, bulk_samples, label) {
  if (length(values) != length(bulk_samples)) {
    stop(sprintf("`%s` must contain exactly one entry per bulk sample.", label), call. = FALSE)
  }
  value_names <- names(values)
  has_any_names <- !is_null(value_names) && any(!is.na(value_names) & nzchar(value_names))
  if (!has_any_names) {
    names(values) <- bulk_samples
    return(list(values = values, mode = "positional"))
  }
  if (anyNA(value_names) || any(!nzchar(value_names)) || anyDuplicated(value_names)) {
    stop(sprintf("Named `%s` entries must have unique, non-empty sample names.", label), call. = FALSE)
  }
  if (!setequal(value_names, bulk_samples)) {
    stop(sprintf("Named `%s` entries must match `bulk_expression` sample names exactly.", label), call. = FALSE)
  }
  list(values = values[match(bulk_samples, value_names)], mode = "named")
}

.sn_scissor_align_phenotype <- function(bulk_phenotype, bulk_samples, family, control) {
  if (family == "cox") {
    if (!is.matrix(bulk_phenotype) && !is.data.frame(bulk_phenotype)) {
      stop("For `family = \"cox\"`, `bulk_phenotype` must be a two-column time/event matrix.", call. = FALSE)
    }
    if (NCOL(bulk_phenotype) != 2L || NROW(bulk_phenotype) != length(bulk_samples)) {
      stop("For `family = \"cox\"`, `bulk_phenotype` must have one row per bulk sample and exactly two columns.", call. = FALSE)
    }
    alignment_mode <- "positional"
    if (.sn_scissor_has_explicit_row_names(bulk_phenotype)) {
      phenotype_rows <- rownames(bulk_phenotype)
      if (anyNA(phenotype_rows) || any(!nzchar(phenotype_rows)) || anyDuplicated(phenotype_rows)) {
        stop("Cox `bulk_phenotype` row names must be unique and non-empty.", call. = FALSE)
      }
      if (!setequal(phenotype_rows, bulk_samples)) {
        stop("Cox `bulk_phenotype` row names must match `bulk_expression` sample names exactly.", call. = FALSE)
      }
      bulk_phenotype <- bulk_phenotype[match(bulk_samples, phenotype_rows), , drop = FALSE]
      alignment_mode <- "named"
    }
    phenotype <- as.matrix(bulk_phenotype)
    if (!is.numeric(phenotype) || any(!is.finite(phenotype))) {
      stop("Cox `bulk_phenotype` time and event columns must be finite numeric values.", call. = FALSE)
    }
    if (any(phenotype[, 1L] <= 0)) {
      stop("Cox survival times must be greater than zero.", call. = FALSE)
    }
    events <- phenotype[, 2L]
    if (any(!events %in% c(0, 1)) || length(unique(events)) < 2L) {
      stop("Cox event status must contain both 0 (censored) and 1 (event).", call. = FALSE)
    }
    rownames(phenotype) <- bulk_samples
    if (is_null(colnames(phenotype))) colnames(phenotype) <- c("time", "event")
    return(list(
      values = phenotype,
      tag = control$tag %||% NULL,
      mode = alignment_mode,
      levels = c("censored", "event")
    ))
  }

  phenotype <- bulk_phenotype
  if (is.matrix(phenotype) || is.data.frame(phenotype)) {
    if (NCOL(phenotype) != 1L) {
      stop(sprintf("For `family = \"%s\"`, `bulk_phenotype` must be a vector or one-column table.", family), call. = FALSE)
    }
    row_names <- if (.sn_scissor_has_explicit_row_names(phenotype)) rownames(phenotype) else NULL
    phenotype <- phenotype[, 1L]
    if (!is_null(row_names)) names(phenotype) <- row_names
  }
  aligned <- .sn_scissor_align_vector(phenotype, bulk_samples, "bulk_phenotype")
  phenotype <- aligned$values
  if (anyNA(phenotype)) {
    stop("`bulk_phenotype` must not contain missing values.", call. = FALSE)
  }
  if (family == "binomial") {
    labels <- as.character(phenotype)
    observed <- if (is.factor(phenotype)) levels(droplevels(phenotype)) else unique(labels)
    observed <- observed[observed %in% labels]
    if (length(observed) != 2L) {
      stop("For `family = \"binomial\"`, `bulk_phenotype` must contain exactly two classes.", call. = FALSE)
    }
    values <- match(labels, observed) - 1L
    names(values) <- bulk_samples
    tag <- control$tag %||% observed
    if (length(tag) != 2L || anyNA(tag) || any(!nzchar(as.character(tag)))) {
      stop("Scissor `tag` must contain exactly two non-empty labels for a binomial phenotype.", call. = FALSE)
    }
    return(list(values = values, tag = as.character(tag), mode = aligned$mode, levels = observed))
  }
  values <- if (is.numeric(phenotype)) as.numeric(phenotype) else suppressWarnings(as.numeric(as.character(phenotype)))
  if (any(!is.finite(values))) {
    stop("For `family = \"gaussian\"`, `bulk_phenotype` must contain finite numeric values.", call. = FALSE)
  }
  if (length(unique(values)) < 2L) {
    stop("For `family = \"gaussian\"`, `bulk_phenotype` must vary across samples.", call. = FALSE)
  }
  names(values) <- bulk_samples
  observed <- names(table(values))
  tag <- control$tag %||% observed
  if (length(tag) != length(observed)) {
    stop("For a Gaussian phenotype, Scissor `tag` must have one label per observed phenotype value.", call. = FALSE)
  }
  list(values = values, tag = as.character(tag), mode = aligned$mode, levels = observed)
}

.sn_scissor_align_inputs <- function(object,
                                     state_by,
                                     sample_by,
                                     bulk_expression,
                                     bulk_phenotype,
                                     family,
                                     control) {
  .sn_scissor_validate_control(control)
  metadata <- object[[]]
  required_columns <- c(state_by, sample_by)
  required_columns <- unique(required_columns[!vapply(required_columns, is_null, logical(1))])
  missing_columns <- setdiff(required_columns, colnames(metadata))
  if (length(missing_columns) > 0L) {
    stop(sprintf("Scissor metadata column(s) not found: %s.", paste(missing_columns, collapse = ", ")), call. = FALSE)
  }
  for (column in required_columns) {
    values <- as.character(metadata[[column]])
    if (anyNA(values) || any(!nzchar(values))) {
      stop(sprintf("Scissor metadata column `%s` must not contain missing or empty values.", column), call. = FALSE)
    }
  }
  bulk <- .sn_scissor_validate_bulk(bulk_expression)
  phenotype <- .sn_scissor_align_phenotype(bulk_phenotype, colnames(bulk), family, control)
  assay <- control$assay %||% Seurat::DefaultAssay(object)
  layer <- control$layer %||% "data"
  expression_matrix <- .sn_get_seurat_layer_data(object, assay = assay, layer = layer)
  expression_cells <- colnames(expression_matrix)
  object_cells <- colnames(object)
  if (is_null(expression_cells) || !setequal(expression_cells, object_cells)) {
    stop(
      "The selected Scissor expression layer must contain every input cell exactly once.",
      call. = FALSE
    )
  }
  expression_matrix <- expression_matrix[, object_cells, drop = FALSE]
  expression_values <- if (inherits(expression_matrix, "sparseMatrix")) {
    expression_matrix@x
  } else {
    as.numeric(expression_matrix)
  }
  if (any(!is.finite(expression_values))) {
    stop("The selected Scissor expression layer contains non-finite values.", call. = FALSE)
  }
  expression <- list(matrix = expression_matrix, assay = assay, layer = layer)
  shared_features <- intersect(rownames(bulk), rownames(expression_matrix))
  min_shared_features <- as.integer(control$min_shared_features %||% 10L)
  if (length(shared_features) < min_shared_features) {
    stop(sprintf(
      "Scissor requires at least %d shared feature(s) between bulk and single-cell expression; found %d.",
      min_shared_features, length(shared_features)
    ), call. = FALSE)
  }
  list(
    bulk_expression = bulk,
    phenotype = phenotype$values,
    tag = phenotype$tag,
    phenotype_alignment = phenotype$mode,
    phenotype_levels = phenotype$levels,
    shared_features = shared_features,
    assay = assay,
    layer = layer,
    expression = expression
  )
}

.sn_scissor_standardize_reliability <- function(value, cells) {
  if (is_null(value) || NROW(value) == 0L) return(tibble::tibble())
  value <- as.data.frame(value, stringsAsFactors = FALSE)
  clean_names <- tolower(gsub("[^[:alnum:]]+", "_", trimws(colnames(value))))
  clean_names <- gsub("^_+|_+$", "", clean_names)
  colnames(value) <- make.unique(clean_names, sep = "_")
  if (!"cell" %in% colnames(value)) {
    row_ids <- rownames(value)
    if (!is_null(row_ids) && !identical(row_ids, as.character(seq_len(nrow(value))))) {
      value$cell <- row_ids
    }
  }
  if (!"cell" %in% colnames(value)) {
    stop("Scissor reliability output must identify cells by a `cell` column or row names.", call. = FALSE)
  }
  if ("cell" %in% colnames(value)) {
    value$cell <- as.character(value$cell)
    if (anyDuplicated(value$cell) || any(!value$cell %in% cells)) {
      stop("Scissor reliability rows must identify unique cells from the input object.", call. = FALSE)
    }
    value <- value[, c("cell", setdiff(colnames(value), "cell")), drop = FALSE]
  }
  tibble::as_tibble(value)
}

.sn_scissor_correlation_tables <- function(value,
                                           bulk_samples,
                                           cell_table,
                                           retain_matrix,
                                           max_rows) {
  empty <- list(summary = tibble::tibble(), matrix = tibble::tibble(), warnings = character())
  if (is_null(value)) return(empty)
  value <- as.matrix(value)
  if (!is.numeric(value) || nrow(value) != length(bulk_samples) || ncol(value) != nrow(cell_table)) {
    stop("The Scissor correlation matrix must be bulk-sample by single-cell.", call. = FALSE)
  }
  if (is_null(rownames(value))) rownames(value) <- bulk_samples
  if (is_null(colnames(value))) colnames(value) <- cell_table$cell
  if (!setequal(rownames(value), bulk_samples) || !setequal(colnames(value), cell_table$cell)) {
    stop("The Scissor correlation matrix names must align with bulk samples and single-cell barcodes.", call. = FALSE)
  }
  value <- value[bulk_samples, cell_table$cell, drop = FALSE]

  finite_summary <- function(x, fun) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else fun(x)
  }
  correlation_summary <- tibble::tibble(
    cell = colnames(value),
    mean_correlation = apply(value, 2L, finite_summary, fun = mean),
    median_correlation = apply(value, 2L, finite_summary, fun = stats::median),
    min_correlation = apply(value, 2L, finite_summary, fun = min),
    max_correlation = apply(value, 2L, finite_summary, fun = max),
    positive_fraction = apply(value, 2L, function(x) {
      x <- x[is.finite(x)]
      if (length(x) == 0L) NA_real_ else mean(x > 0)
    }),
    negative_fraction = apply(value, 2L, function(x) {
      x <- x[is.finite(x)]
      if (length(x) == 0L) NA_real_ else mean(x < 0)
    })
  )
  correlation_summary <- dplyr::left_join(correlation_summary, cell_table, by = "cell")

  long <- tibble::tibble()
  correlation_warnings <- character()
  if (isTRUE(retain_matrix)) {
    rows <- length(value)
    if (rows > max_rows) {
      correlation_warnings <- sprintf(
        "The full Scissor correlation matrix has %d rows and was not retained because `max_correlation_rows` is %d.",
        rows, max_rows
      )
    } else {
      long <- tibble::tibble(
        sample = rep(rownames(value), times = ncol(value)),
        cell = rep(colnames(value), each = nrow(value)),
        correlation = as.vector(value)
      )
    }
  }
  list(summary = correlation_summary, matrix = long, warnings = correlation_warnings)
}

.sn_scissor_run_reliability <- function(fit,
                                        save_file,
                                        correlation_matrix,
                                        aligned,
                                        object,
                                        cell_table,
                                        control) {
  warnings <- character()
  reliability <- control$reliability_result
  if (is_null(reliability) && is.function(control$reliability_runner)) {
    reliability <- tryCatch(
      control$reliability_runner(
        fit = fit,
        load_file = save_file,
        correlation_matrix = correlation_matrix,
        phenotype = aligned$phenotype,
        object = object
      ),
      error = function(error) {
        warnings <<- paste0("Scissor reliability runner failed: ", conditionMessage(error))
        NULL
      }
    )
  }
  if (is_null(reliability) && isTRUE(control$reliability)) {
    selected_cells <- cell_table$cell[cell_table$selection != "Unselected"]
    if (length(selected_cells) == 0L) {
      warnings <- "Scissor reliability was requested but no cells were selected."
    } else if (!file.exists(save_file)) {
      warnings <- "Scissor reliability was requested, but the backend did not create the required saved model data."
    } else {
      check_installed("Scissor", reason = "to evaluate Scissor cell-level reliability.")
      reliability_control <- control$reliability_control %||% list()
      reliability <- tryCatch(
        withCallingHandlers(
          Scissor::evaluate.cell(
            Load_file = save_file,
            Scissor_result = fit,
            FDR_cutoff = reliability_control$fdr_cutoff %||% 0.05,
            bootstrap_n = as.integer(reliability_control$bootstrap_n %||% 100L)
          ),
          warning = function(warning) {
            warnings <<- c(warnings, conditionMessage(warning))
            invokeRestart("muffleWarning")
          }
        ),
        error = function(error) {
          warnings <<- c(warnings, paste0("Scissor reliability evaluation failed: ", conditionMessage(error)))
          NULL
        }
      )
    }
  }
  reliability <- .sn_scissor_standardize_reliability(reliability, cells = cell_table$cell)
  if (nrow(reliability) > 0L && "cell" %in% colnames(reliability)) {
    context_columns <- c(
      "cell",
      setdiff(c("state", "selection", intersect("sample", colnames(cell_table))), colnames(reliability))
    )
    reliability <- dplyr::left_join(reliability, cell_table[, context_columns, drop = FALSE], by = "cell")
  }
  list(table = reliability, warnings = warnings)
}

.sn_priority_scissor_impl <- function(object,
                                      state_by,
                                      sample_by,
                                      bulk_expression,
                                      bulk_phenotype,
                                      family,
                                      backend_control) {
  if (is_null(bulk_expression) || is_null(bulk_phenotype)) {
    stop("Scissor requires `bulk_expression` and `bulk_phenotype`; a cell metadata phenotype is not a substitute.", call. = FALSE)
  }
  aligned <- .sn_scissor_align_inputs(
    object = object,
    state_by = state_by,
    sample_by = sample_by,
    bulk_expression = bulk_expression,
    bulk_phenotype = bulk_phenotype,
    family = family,
    control = backend_control
  )
  metadata <- object[[]]
  cells <- colnames(object)
  alpha <- backend_control$alpha %||% NULL
  cutoff <- backend_control$cutoff %||% 0.2
  preprocessing_seed <- .sn_scissor_validate_seed(backend_control$preprocessing_seed %||% 123L)
  save_file <- tempfile("shennong-scissor-", fileext = ".RData")
  on.exit(unlink(save_file), add = TRUE)

  backend_warnings <- character()
  fit <- backend_control$result
  if (is_null(fit)) {
    check_installed("Scissor", reason = "to run phenotype-guided Scissor selection.")
    set.seed(preprocessing_seed)
    selected_expression <- aligned$expression$matrix
    placeholder_counts <- Matrix::sparseMatrix(
      i = integer(),
      j = integer(),
      x = numeric(),
      dims = dim(selected_expression),
      dimnames = dimnames(selected_expression)
    )
    working <- SeuratObject::CreateSeuratObject(
      counts = placeholder_counts,
      meta.data = metadata
    )
    suppressWarnings(
      working[["RNA"]] <- SeuratObject::CreateAssayObject(counts = placeholder_counts)
    )
    SeuratObject::LayerData(working[["RNA"]], layer = "data") <- selected_expression

    expression_mean <- Matrix::rowMeans(selected_expression)
    expression_variance <- Matrix::rowMeans(selected_expression ^ 2) - expression_mean ^ 2
    usable_features <- which(is.finite(expression_variance) &
      expression_variance > sqrt(.Machine$double.eps))
    if (length(usable_features) < 2L) {
      stop(
        "The selected Scissor expression layer requires at least two variable features.",
        call. = FALSE
      )
    }
    feature_order <- usable_features[
      order(expression_variance[usable_features], decreasing = TRUE)
    ]
    feature_order <- utils::head(
      feature_order,
      min(backend_control$nfeatures %||% 2000L, length(feature_order))
    )
    variable_features <- rownames(selected_expression)[feature_order]
    SeuratObject::VariableFeatures(working) <- variable_features
    working <- Seurat::ScaleData(
      working,
      features = variable_features,
      verbose = FALSE
    )
    npcs <- min(
      backend_control$npcs %||% 20L,
      length(variable_features),
      ncol(working) - 1L
    )
    if (npcs < 2L) stop("Scissor requires at least two usable PCA dimensions.", call. = FALSE)
    working <- suppressWarnings(Seurat::RunPCA(
      working,
      features = variable_features,
      npcs = npcs,
      verbose = FALSE
    ))
    working <- Seurat::FindNeighbors(working, dims = seq_len(min(10L, npcs)), verbose = FALSE)
    runner <- backend_control$runner %||% Scissor::Scissor
    fit <- withCallingHandlers(
      runner(
        bulk_dataset = aligned$bulk_expression,
        sc_dataset = working,
        phenotype = aligned$phenotype,
        tag = aligned$tag,
        alpha = alpha,
        cutoff = cutoff,
        family = family,
        Save_file = save_file
      ),
      warning = function(warning) {
        backend_warnings <<- c(backend_warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
    )
  }
  if (!is.list(fit) || is_null(fit$Coefs)) {
    stop("The Scissor backend must return a list containing `Coefs`.", call. = FALSE)
  }

  coefficient_names <- names(fit$Coefs)
  coefficients <- as.numeric(fit$Coefs)
  if (!is_null(coefficient_names) && any(nzchar(coefficient_names))) {
    if (any(!nzchar(coefficient_names)) || anyDuplicated(coefficient_names) ||
        !setequal(coefficient_names, cells)) {
      stop("Named Scissor coefficients must match the single-cell barcodes exactly.", call. = FALSE)
    }
    coefficients <- coefficients[match(cells, coefficient_names)]
  }
  if (length(coefficients) != length(cells) || any(!is.finite(coefficients))) {
    stop("Scissor returned finite coefficients that do not align with the single-cell object.", call. = FALSE)
  }
  current <- tibble::tibble(
    cell = cells,
    coefficient = coefficients,
    selection = ifelse(coefficients > 0, "Scissor+", ifelse(coefficients < 0, "Scissor-", "Unselected")),
    state = if (is_null(state_by)) "all_cells" else as.character(metadata[cells, state_by])
  )
  if (!is_null(sample_by)) current$sample <- as.character(metadata[cells, sample_by])

  saved <- new.env(parent = emptyenv())
  if (file.exists(save_file)) load(save_file, envir = saved)
  correlation_matrix <- backend_control$correlation_matrix %||%
    if (exists("X", envir = saved, inherits = FALSE)) get("X", envir = saved) else NULL
  correlation_tables <- .sn_scissor_correlation_tables(
    value = correlation_matrix,
    bulk_samples = colnames(aligned$bulk_expression),
    cell_table = current,
    retain_matrix = backend_control$retain_correlation_matrix %||% FALSE,
    max_rows = as.integer(backend_control$max_correlation_rows %||% 1000000L)
  )
  if (nrow(correlation_tables$summary) > 0L) {
    correlation_columns <- setdiff(
      colnames(correlation_tables$summary),
      c("cell", "coefficient", "selection", "state", "sample")
    )
    current <- dplyr::left_join(
      current,
      correlation_tables$summary[, c("cell", correlation_columns), drop = FALSE],
      by = "cell"
    )
    context_columns <- c("cell", "coefficient", "selection", "state", intersect("sample", colnames(current)))
    correlation_tables$summary <- dplyr::left_join(
      correlation_tables$summary[, c("cell", correlation_columns), drop = FALSE],
      current[, context_columns, drop = FALSE],
      by = "cell"
    )
  }

  ranking <- dplyr::bind_rows(lapply(split(current, current$state), function(group) {
    mean_correlation <- if ("mean_correlation" %in% colnames(group)) {
      mean(group$mean_correlation, na.rm = TRUE)
    } else {
      NA_real_
    }
    if (!is.finite(mean_correlation)) mean_correlation <- NA_real_
    tibble::tibble(
      state = group$state[[1L]],
      priority_score = mean(abs(group$coefficient)),
      phenotype_association = mean(group$coefficient),
      selected_fraction = mean(group$selection != "Unselected"),
      positive_fraction = mean(group$selection == "Scissor+"),
      negative_fraction = mean(group$selection == "Scissor-"),
      mean_correlation = mean_correlation,
      uncertainty = stats::sd(group$coefficient) / sqrt(nrow(group)),
      n_cells = nrow(group)
    )
  }))
  ranking <- ranking[order(ranking$priority_score, decreasing = TRUE), , drop = FALSE]

  samples <- tibble::tibble()
  if (!is_null(sample_by)) {
    samples <- dplyr::bind_rows(lapply(
      split(current, interaction(current$state, current$sample, drop = TRUE)),
      function(group) {
        tibble::tibble(
          sample = group$sample[[1L]],
          state = group$state[[1L]],
          contribution = mean(group$coefficient),
          selected_fraction = mean(group$selection != "Unselected"),
          positive_fraction = mean(group$selection == "Scissor+"),
          negative_fraction = mean(group$selection == "Scissor-"),
          n_cells = nrow(group)
        )
      }
    ))
  }

  reliability <- .sn_scissor_run_reliability(
    fit = fit,
    save_file = save_file,
    correlation_matrix = correlation_matrix,
    aligned = aligned,
    object = object,
    cell_table = current,
    control = backend_control
  )
  fit_parameters <- fit$para %||% list()
  alpha_search_requested <- is_null(alpha)
  selected_alpha <- as.numeric((fit_parameters$alpha %||% alpha %||% NA_real_)[[1L]])
  selected_fraction <- mean(current$selection != "Unselected")
  selection_cutoff_satisfied <- selected_fraction < cutoff
  cutoff_warning <- if (!selection_cutoff_satisfied) {
    sprintf(
      paste0(
        "Scissor selected %.1f%% of cells, which did not fall below the configured ",
        "cutoff of %.1f%%; review the selected alpha and treat the cell set as exploratory."
      ),
      100 * selected_fraction,
      100 * cutoff
    )
  } else {
    character()
  }
  model_table <- tibble::tibble(
    family = family,
    alpha = selected_alpha,
    lambda = as.numeric((fit_parameters$lambda %||% NA_real_)[[1L]]),
    cutoff = cutoff,
    alpha_search_requested = alpha_search_requested,
    selection_cutoff_satisfied = selection_cutoff_satisfied,
    tag = paste(aligned$tag %||% character(), collapse = " | "),
    phenotype_levels = paste(aligned$phenotype_levels, collapse = " | "),
    phenotype_alignment = aligned$phenotype_alignment,
    bulk_samples = ncol(aligned$bulk_expression),
    cells = length(cells),
    shared_features = length(aligned$shared_features),
    preprocessing_seed = preprocessing_seed,
    backend_seed = 123L,
    selected_positive = sum(current$selection == "Scissor+"),
    selected_negative = sum(current$selection == "Scissor-"),
    selected_fraction = selected_fraction
  )

  list(
    ranking = ranking,
    cells = current,
    samples = samples,
    correlations = correlation_tables$summary,
    correlation_matrix = correlation_tables$matrix,
    model_table = model_table,
    reliability = reliability$table,
    null = tibble::tibble(),
    contrast = NULL,
    input = list(
      assay = aligned$assay,
      layer = aligned$layer,
      bulk_samples = colnames(aligned$bulk_expression),
      shared_features = aligned$shared_features,
      phenotype_alignment = aligned$phenotype_alignment,
      phenotype_levels = aligned$phenotype_levels
    ),
    parameters = list(
      family = family,
      alpha = if (alpha_search_requested) NA_real_ else alpha,
      alpha_search_requested = alpha_search_requested,
      cutoff = cutoff,
      tag = aligned$tag,
      reliability_requested = isTRUE(backend_control$reliability),
      retain_correlation_matrix = isTRUE(backend_control$retain_correlation_matrix),
      preprocessing_seed = preprocessing_seed,
      backend_seed = 123L
    ),
    diagnostics = list(
      shared_feature_count = length(aligned$shared_features),
      selected_cell_count = sum(current$selection != "Unselected"),
      positive_cell_count = sum(current$selection == "Scissor+"),
      negative_cell_count = sum(current$selection == "Scissor-"),
      selected_fraction = selected_fraction,
      selection_cutoff_satisfied = selection_cutoff_satisfied,
      selected_alpha = selected_alpha,
      alpha_search_requested = alpha_search_requested,
      inferential_unit = "bulk_sample",
      single_cell_sample_summary = !is_null(sample_by),
      backend_random_seed = 123L,
      correlation_summary_available = nrow(correlation_tables$summary) > 0L,
      reliability_rows = nrow(reliability$table)
    ),
    warnings = unique(c(
      backend_warnings,
      correlation_tables$warnings,
      reliability$warnings,
      cutoff_warning,
      "Scissor does not return a phenotype-permutation null by default; opt into its bootstrap reliability workflow for confirmatory cell-level checks."
    )),
    model = fit_parameters
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
  scissor_control <- if (identical(method, "scissor")) {
    control <- .sn_scissor_backend_control(backend_control, assay = assay, layer = layer)
    control$preprocessing_seed <- .sn_scissor_validate_seed(seed)
    control
  } else {
    NULL
  }
  backend <- switch(
    method,
    augur = .sn_priority_augur(
      object, phenotype, state_by, sample_by, contrast, assay, layer, features,
      max_features, max_cells_per_state, folds, as.integer(permutations), seed
    ),
    rareq = .sn_priority_rareq(object, phenotype, sample_by, contrast, reduction, dims, backend_control$rareq %||% list()),
    scissor = .sn_priority_scissor(
      object,
      state_by,
      sample_by,
      bulk_expression,
      bulk_phenotype,
      family,
      scissor_control
    )
  )
  tables <- list(
    primary = backend$ranking,
    cells = backend$cells,
    sample_contributions = backend$samples,
    permutation_null = backend$null
  )
  if (method == "scissor") {
    tables$states <- backend$ranking
    tables$samples <- backend$samples
    tables$correlations <- backend$correlations %||% tibble::tibble()
    tables$correlation_matrix <- backend$correlation_matrix %||% tibble::tibble()
    tables$model <- backend$model_table %||% tibble::tibble()
    tables$reliability <- backend$reliability %||% tibble::tibble()
  }
  result <- .sn_new_analysis_result(
    analysis_type = "state_priority",
    name = store_name,
    method = method,
    backend = switch(method, augur = "Shennong sample-aware Augur", rareq = "RareQ", scissor = "Scissor"),
    input = utils::modifyList(
      list(cells = ncol(object), phenotype = phenotype, sample_by = sample_by, state_by = state_by),
      backend$input %||% list()
    ),
    parameters = utils::modifyList(
      list(
        contrast = backend$contrast,
        permutations = as.integer(permutations),
        family = if (method == "scissor") family else NULL
      ),
      backend$parameters %||% list()
    ),
    tables = tables,
    embeddings = list(),
    graphs = list(),
    models = list(summary = backend$model %||% list()),
    diagnostics = utils::modifyList(
      list(
        states = nrow(backend$ranking),
        sample_aware = !identical(method, "scissor") && !is_null(sample_by)
      ),
      backend$diagnostics %||% list()
    ),
    warnings = as.character(backend$warnings),
    random_seed = if (identical(method, "scissor")) {
      list(preprocessing = as.integer(seed), scissor_backend = 123L)
    } else {
      seed
    }
  )
  object <- sn_store_result(object, "state_priority", store_name, result)
  object <- .sn_log_seurat_command(object = object, assay = assay, name = "sn_prioritize_states")
  if (isTRUE(return_object)) object else sn_get_result(object, "state_priority", store_name)
}

#' Run phenotype-guided Scissor cell selection
#'
#' `sn_run_scissor()` is the direct Scissor entry point. It returns a cell-first
#' unified Scissor result while sharing the same backend implementation used by
#' [sn_prioritize_states()]. The result also carries state and sample summaries,
#' correlation summaries, model metadata, and optional reliability output.
#'
#' @param object A Seurat object.
#' @param bulk_expression Numeric gene-by-bulk-sample expression matrix with
#'   unique gene row names and sample column names.
#' @param bulk_phenotype A binary or numeric phenotype vector, or a two-column
#'   time/event matrix for `family = "cox"`. Named inputs are aligned strictly
#'   to `bulk_expression`; unnamed inputs retain positional compatibility.
#' @param state_by Optional cell-state metadata column. When omitted, the
#'   direct result remains cell-first and uses one descriptive `all_cells`
#'   state summary.
#' @param sample_by Optional biological-sample metadata column.
#' @param family Scissor response family.
#' @param phenotype Descriptive label stored with the result.
#' @param assay,layer Single-cell expression source.
#' @param store_name Stored result name under `scissor`.
#' @param seed Random seed recorded in provenance.
#' @param backend_control Direct Scissor controls, or a list containing a
#'   `scissor` sub-list. Set `reliability = TRUE` to run bootstrap reliability;
#'   set `retain_correlation_matrix = TRUE` to retain the long sample-cell
#'   correlation table subject to `max_correlation_rows`. Scissor's `cutoff` is
#'   an alpha-search stopping criterion rather than a guaranteed upper bound;
#'   Shennong warns and marks the result exploratory when no candidate alpha
#'   achieves it.
#' @param return_object Return the updated object instead of the result.
#' @return A Seurat object or unified Scissor state-priority result.
#' @examples
#' \dontrun{
#' result <- sn_run_scissor(
#'   object,
#'   bulk_expression = bulk_matrix,
#'   bulk_phenotype = bulk_group,
#'   state_by = "cell_type",
#'   sample_by = "sample",
#'   return_object = FALSE
#' )
#' result$tables$cells
#' result$tables$states
#' }
#' @export
sn_run_scissor <- function(object,
                           bulk_expression,
                           bulk_phenotype,
                           state_by = NULL,
                           sample_by = NULL,
                           family = c("binomial", "gaussian", "cox"),
                           phenotype = "bulk_phenotype",
                           assay = NULL,
                           layer = "data",
                           store_name = "scissor",
                           seed = 717L,
                           backend_control = list(),
                           return_object = TRUE) {
  .sn_validate_result_object(object)
  family <- match.arg(family)
  if (!is_null(state_by) && !state_by %in% colnames(object[[]])) {
    stop("`state_by` must name object metadata for Scissor analysis.", call. = FALSE)
  }
  control <- .sn_scissor_backend_control(backend_control, assay = assay, layer = layer)
  control$preprocessing_seed <- .sn_scissor_validate_seed(seed)
  backend <- .sn_priority_scissor(
    object = object,
    state_by = state_by,
    sample_by = sample_by,
    bulk_expression = bulk_expression,
    bulk_phenotype = bulk_phenotype,
    family = family,
    backend_control = control
  )
  result <- .sn_new_analysis_result(
    analysis_type = "scissor",
    name = store_name,
    method = "scissor",
    backend = "Scissor",
    input = utils::modifyList(
      list(cells = ncol(object), phenotype = phenotype, sample_by = sample_by, state_by = state_by),
      backend$input %||% list()
    ),
    parameters = utils::modifyList(list(family = family), backend$parameters %||% list()),
    tables = list(
      primary = backend$cells,
      cells = backend$cells,
      states = backend$ranking,
      samples = backend$samples,
      sample_contributions = backend$samples,
      correlations = backend$correlations %||% tibble::tibble(),
      correlation_matrix = backend$correlation_matrix %||% tibble::tibble(),
      model = backend$model_table %||% tibble::tibble(),
      reliability = backend$reliability %||% tibble::tibble(),
      permutation_null = backend$null
    ),
    embeddings = list(),
    graphs = list(),
    models = list(summary = backend$model %||% list()),
    diagnostics = utils::modifyList(
      list(states = nrow(backend$ranking), sample_aware = FALSE),
      backend$diagnostics %||% list()
    ),
    warnings = as.character(backend$warnings),
    random_seed = list(
      preprocessing = as.integer(seed),
      scissor_backend = 123L
    )
  )
  object <- sn_store_result(object, "scissor", store_name, result)
  object <- .sn_log_seurat_command(
    object = object,
    assay = backend$input$assay %||% assay,
    name = "sn_run_scissor"
  )
  if (isTRUE(return_object)) object else sn_get_result(object, "scissor", store_name)
}
