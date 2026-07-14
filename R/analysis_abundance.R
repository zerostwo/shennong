.sn_abundance_inputs <- function(object, sample_by, condition_by, cell_type_by, contrast = NULL, extra_columns = NULL) {
  .sn_validate_result_object(object)
  columns <- unique(c(sample_by, condition_by, cell_type_by, extra_columns))
  missing <- setdiff(columns, colnames(object[[]]))
  if (length(missing) > 0L) stop("Missing abundance metadata column(s): ", paste(missing, collapse = ", "), ".", call. = FALSE)
  metadata <- object[[]][, columns, drop = FALSE]
  metadata$cell <- rownames(metadata)
  metadata <- metadata[stats::complete.cases(metadata[, columns, drop = FALSE]), , drop = FALSE]
  if (nrow(metadata) == 0L) stop("No complete cells remain for abundance testing.", call. = FALSE)
  .sn_validate_constant_within_sample(metadata, sample_col = sample_by, group_col = condition_by)

  sample_info <- unique(metadata[, c(sample_by, condition_by), drop = FALSE])
  levels <- unique(as.character(sample_info[[condition_by]]))
  if (is_null(contrast)) {
    if (length(levels) != 2L) stop("`contrast` is required unless `condition_by` has exactly two levels.", call. = FALSE)
    contrast <- levels
  }
  if (!is.character(contrast) || length(contrast) != 2L || anyDuplicated(contrast)) {
    stop("`contrast` must be `c(case, control)` with two distinct condition labels.", call. = FALSE)
  }
  metadata <- metadata[metadata[[condition_by]] %in% contrast, , drop = FALSE]
  sample_info <- unique(metadata[, c(sample_by, condition_by), drop = FALSE])
  group_counts <- table(factor(sample_info[[condition_by]], levels = contrast))
  if (any(group_counts < 2L)) {
    stop("Abundance testing requires at least two biological samples in each contrast group.", call. = FALSE)
  }

  counts <- metadata |>
    dplyr::count(.data[[sample_by]], .data[[condition_by]], .data[[cell_type_by]], name = "count")
  cell_types <- unique(as.character(metadata[[cell_type_by]]))
  grid <- merge(
    sample_info,
    stats::setNames(data.frame(cell_types, stringsAsFactors = FALSE), cell_type_by),
    by = NULL
  )
  sample_data <- merge(grid, counts, by = c(sample_by, condition_by, cell_type_by), all.x = TRUE, sort = FALSE)
  sample_data$count[is.na(sample_data$count)] <- 0L
  totals <- metadata |>
    dplyr::count(.data[[sample_by]], name = "sample_total")
  sample_data <- dplyr::left_join(sample_data, totals, by = sample_by)
  sample_data$proportion <- sample_data$count / sample_data$sample_total
  list(
    metadata = metadata,
    sample_info = sample_info,
    sample_data = tibble::as_tibble(sample_data),
    contrast = contrast,
    group_counts = group_counts
  )
}

.sn_abundance_summary <- function(inputs, sample_by, condition_by, cell_type_by, pseudocount = 1e-06) {
  by_type <- split(inputs$sample_data, as.character(inputs$sample_data[[cell_type_by]]))
  table <- dplyr::bind_rows(lapply(names(by_type), function(feature) {
    current <- by_type[[feature]]
    case <- current$proportion[current[[condition_by]] == inputs$contrast[[1]]]
    control <- current$proportion[current[[condition_by]] == inputs$contrast[[2]]]
    tibble::tibble(
      feature = feature,
      mean_case = mean(case),
      mean_control = mean(control),
      estimate = mean(case) - mean(control),
      log2_ratio = log2((mean(case) + pseudocount) / (mean(control) + pseudocount)),
      n_case = length(case),
      n_control = length(control)
    )
  }))
  contributions <- inputs$sample_data
  n_case <- unname(inputs$group_counts[[1]])
  n_control <- unname(inputs$group_counts[[2]])
  contributions$contribution <- ifelse(
    contributions[[condition_by]] == inputs$contrast[[1]],
    contributions$proportion / n_case,
    -contributions$proportion / n_control
  )
  list(summary = table, contributions = contributions)
}

.sn_propeller_design <- function(inputs, sample_by, condition_by, design, contrast, sample_order) {
  design_data <- inputs$sample_info
  rownames(design_data) <- as.character(design_data[[sample_by]])
  design_data <- design_data[sample_order, , drop = FALSE]
  design_data[[condition_by]] <- factor(as.character(design_data[[condition_by]]), levels = c(contrast[[2]], contrast[[1]]))
  design_data$.sn_condition <- design_data[[condition_by]]
  if (is_null(design)) {
    formula <- stats::as.formula("~ 0 + .sn_condition")
    condition_prefix <- ".sn_condition"
  } else if (inherits(design, "formula")) {
    formula <- design
    condition_prefix <- condition_by
  } else if (is.character(design)) {
    missing <- setdiff(design, colnames(inputs$metadata))
    if (length(missing) > 0L) stop("Design covariate(s) were not found: ", paste(missing, collapse = ", "), ".", call. = FALSE)
    for (column in design) {
      .sn_validate_constant_within_sample(inputs$metadata, sample_col = sample_by, group_col = column)
      values <- tapply(inputs$metadata[[column]], inputs$metadata[[sample_by]], function(value) value[[1]])
      design_data[[column]] <- unname(values[rownames(design_data)])
    }
    formula <- stats::as.formula(paste("~ 0 + .sn_condition", paste(design, collapse = " + "), sep = if (length(design)) " + " else ""))
    condition_prefix <- ".sn_condition"
  } else {
    stop("`design` must be NULL, a formula, or a character vector of sample-level covariates.", call. = FALSE)
  }
  matrix <- stats::model.matrix(formula, data = design_data)
  locate <- function(level) {
    candidates <- c(
      paste0(condition_prefix, level),
      make.names(paste0(condition_prefix, level))
    )
    hit <- which(colnames(matrix) %in% candidates)
    if (length(hit) == 0L) integer() else hit[[1]]
  }
  case_column <- locate(contrast[[1]])
  control_column <- locate(contrast[[2]])
  if (length(case_column) == 0L || length(control_column) == 0L) {
    stop("The abundance design must contain separate no-intercept coefficients for both contrast groups.", call. = FALSE)
  }
  contrast_vector <- rep(0, ncol(matrix))
  contrast_vector[case_column] <- 1
  contrast_vector[control_column] <- -1
  names(contrast_vector) <- colnames(matrix)
  list(matrix = matrix, contrast = contrast_vector, formula = formula, sample_data = design_data)
}

.sn_test_abundance_propeller <- function(inputs,
                                          sample_by,
                                          condition_by,
                                          cell_type_by,
                                          design,
                                          transform,
                                          robust,
                                          trend) {
  check_installed("speckle", reason = "to run Propeller differential abundance.")
  props <- speckle::getTransformedProps(
    clusters = inputs$metadata[[cell_type_by]],
    sample = inputs$metadata[[sample_by]],
    transform = transform
  )
  design_info <- .sn_propeller_design(
    inputs = inputs,
    sample_by = sample_by,
    condition_by = condition_by,
    design = design,
    contrast = inputs$contrast,
    sample_order = colnames(props$Proportions)
  )
  raw <- speckle::propeller.ttest(
    prop.list = props,
    design = design_info$matrix,
    contrasts = design_info$contrast,
    robust = robust,
    trend = trend,
    sort = FALSE
  )
  raw$feature <- rownames(raw)
  rownames(raw) <- NULL
  list(
    raw = tibble::as_tibble(raw),
    statistics = tibble::tibble(
      feature = raw$feature,
      statistic = raw$Tstatistic,
      p_value = raw$P.Value,
      adjusted_p_value = raw$FDR
    ),
    design = design_info
  )
}

.sn_test_abundance_permutation <- function(inputs,
                                            condition_by,
                                            cell_type_by,
                                            permutations,
                                            seed) {
  permutations <- as.integer(permutations)
  if (length(permutations) != 1L || is.na(permutations) || permutations < 99L) {
    stop("`permutations` must be a single integer of at least 99.", call. = FALSE)
  }
  set.seed(seed)
  wide <- stats::reshape(
    as.data.frame(inputs$sample_data[, c(names(inputs$sample_info)[[1]], condition_by, cell_type_by, "proportion")]),
    idvar = c(names(inputs$sample_info)[[1]], condition_by),
    timevar = cell_type_by,
    direction = "wide"
  )
  feature_columns <- grep("^proportion\\.", names(wide), value = TRUE)
  features <- sub("^proportion\\.", "", feature_columns)
  values <- as.matrix(wide[, feature_columns, drop = FALSE])
  labels <- as.character(wide[[condition_by]])
  observed <- colMeans(values[labels == inputs$contrast[[1]], , drop = FALSE]) -
    colMeans(values[labels == inputs$contrast[[2]], , drop = FALSE])
  null <- matrix(NA_real_, nrow = permutations, ncol = length(features), dimnames = list(NULL, features))
  for (index in seq_len(permutations)) {
    shuffled <- sample(labels, replace = FALSE)
    null[index, ] <- colMeans(values[shuffled == inputs$contrast[[1]], , drop = FALSE]) -
      colMeans(values[shuffled == inputs$contrast[[2]], , drop = FALSE])
  }
  p <- (1 + colSums(abs(null) >= rep(abs(observed), each = permutations))) / (permutations + 1)
  null_table <- tibble::tibble(
    permutation = rep(seq_len(permutations), times = length(features)),
    feature = rep(features, each = permutations),
    null_estimate = as.numeric(null)
  )
  list(
    statistics = tibble::tibble(
      feature = features,
      statistic = observed / apply(null, 2, stats::sd),
      p_value = as.numeric(p),
      adjusted_p_value = stats::p.adjust(as.numeric(p), method = "BH"),
      null_sd = apply(null, 2, stats::sd)
    ),
    null = null_table
  )
}

.sn_standardize_milo_abundance <- function(table, cell_type_by) {
  table <- as.data.frame(table)
  table$feature <- if (!is.null(rownames(table)) && all(nzchar(rownames(table)))) rownames(table) else paste0("Nhood", seq_len(nrow(table)))
  rownames(table) <- NULL
  table$estimate <- table$logFC
  table$p_value <- table$PValue
  table$adjusted_p_value <- table$SpatialFDR %||% table$FDR
  if (!is_null(cell_type_by) && cell_type_by %in% names(table)) table$state <- table[[cell_type_by]]
  tibble::as_tibble(table)
}

.sn_standardize_sccoda_abundance <- function(output) {
  table <- output$table %||% output$effects %||% output
  if (!is.data.frame(table)) stop("scCODA adapter output must contain a result table.", call. = FALSE)
  table <- as.data.frame(table, check.names = FALSE)
  feature_col <- intersect(c("feature", "cell_type", "Cell Type", "covariate"), names(table))
  estimate_col <- intersect(c("estimate", "effect", "Final Parameter", "log2_fold_change", "log2_fc"), names(table))
  inclusion_col <- intersect(c("inclusion_probability", "Inclusion probability", "inclusion_prob", "probability"), names(table))
  log2_col <- intersect(c("log2_fold_change", "log2-fold change", "log2_fc"), names(table))
  credible_col <- intersect(c("credible", "is_credible"), names(table))
  p_value_col <- intersect(c("p_value", "p.value", "pvalue"), names(table))
  adjusted_col <- intersect(c("adjusted_p_value", "padj", "FDR", "q_value"), names(table))
  if (length(feature_col) == 0L || length(estimate_col) == 0L) {
    stop("scCODA result table requires feature/cell_type and effect/estimate columns.", call. = FALSE)
  }
  inclusion <- if (length(inclusion_col) > 0L) as.numeric(table[[inclusion_col[[1]]]]) else rep(NA_real_, nrow(table))
  estimate_name <- estimate_col[[1]]
  estimate <- as.numeric(table[[estimate_name]])
  log2_fc <- if (length(log2_col) > 0L) {
    as.numeric(table[[log2_col[[1]]]])
  } else if (estimate_name %in% c("effect", "Final Parameter")) {
    estimate / log(2)
  } else {
    rep(NA_real_, nrow(table))
  }
  supplied_credible <- output$credible_effects %||% NULL
  credible <- if (length(credible_col) > 0L) {
    as.logical(table[[credible_col[[1]]]])
  } else if (is.logical(supplied_credible) && length(supplied_credible) == nrow(table)) {
    as.logical(supplied_credible)
  } else if (identical(estimate_name, "Final Parameter")) {
    estimate != 0
  } else {
    rep(NA, nrow(table))
  }
  p_value <- if (length(p_value_col) > 0L) as.numeric(table[[p_value_col[[1]]]]) else rep(NA_real_, nrow(table))
  adjusted_p_value <- if (length(adjusted_col) > 0L) as.numeric(table[[adjusted_col[[1]]]]) else rep(NA_real_, nrow(table))
  standardized <- tibble::tibble(
    feature = as.character(table[[feature_col[[1]]]]), estimate = estimate,
    log2_fc = log2_fc, inclusion_probability = inclusion,
    credible = as.logical(credible),
    p_value = p_value,
    adjusted_p_value = adjusted_p_value
  )
  list(table = standardized, raw = tibble::as_tibble(table))
}

#' Test differential abundance across biological samples
#'
#' Provides one stable entry point for sample-level Propeller or permutation
#' tests and neighborhood-level Milo testing. Cells are never treated as
#' independent replicates.
#'
#' @param object A Seurat object.
#' @param method Abundance backend.
#' @param sample_by Biological sample column.
#' @param condition_by Sample-level condition column.
#' @param cell_type_by Cell-type/state column.
#' @param design Optional no-intercept formula or sample-level covariate names
#'   for Propeller; covariate names are forwarded to Milo.
#' @param contrast Two condition labels ordered as `c(case, control)`.
#' @param store_name Stored differential-abundance result name.
#' @param transform Propeller proportion transformation.
#' @param permutations Number of sample-label permutations.
#' @param seed Random seed for permutation testing.
#' @param backend_control Named backend argument lists.
#' @param return_object Return the updated object instead of the result.
#'
#' @return A Seurat object or unified differential-abundance result.
#'
#' @examples
#' \dontrun{
#' object <- sn_test_abundance(
#'   object, sample_by = "sample", condition_by = "condition",
#'   cell_type_by = "cell_type", contrast = c("treated", "control")
#' )
#' sn_get_result(object, "differential_abundance", "abundance")
#' }
#' @export
sn_test_abundance <- function(object,
                              method = c("propeller", "milo", "sccoda", "permutation"),
                              sample_by,
                              condition_by,
                              cell_type_by,
                              design = NULL,
                              contrast = NULL,
                              store_name = "abundance",
                              transform = c("logit", "asin"),
                              permutations = 1000L,
                              seed = 717L,
                              backend_control = list(),
                              return_object = TRUE) {
  method <- match.arg(method)
  transform <- match.arg(transform)
  design_columns <- if (inherits(design, "formula")) all.vars(design) else if (is.character(design)) design else NULL
  inputs <- .sn_abundance_inputs(object, sample_by, condition_by, cell_type_by, contrast, extra_columns = design_columns)
  summaries <- .sn_abundance_summary(inputs, sample_by, condition_by, cell_type_by)
  null <- tibble::tibble()
  design_table <- inputs$sample_info

  if (identical(method, "propeller")) {
    backend <- .sn_test_abundance_propeller(
      inputs, sample_by, condition_by, cell_type_by, design, transform,
      robust = backend_control$propeller$robust %||% TRUE,
      trend = backend_control$propeller$trend %||% FALSE
    )
    primary <- dplyr::left_join(summaries$summary, backend$statistics, by = "feature")
    raw <- backend$raw
    design_table <- tibble::as_tibble(backend$design$sample_data, rownames = sample_by)
  } else if (identical(method, "permutation")) {
    backend <- .sn_test_abundance_permutation(inputs, condition_by, cell_type_by, permutations, seed)
    primary <- dplyr::left_join(summaries$summary, backend$statistics, by = "feature")
    raw <- tibble::tibble()
    null <- backend$null
  } else if (identical(method, "sccoda")) {
    output <- if (is.function(backend_control$runner)) {
      backend_control$runner(
        sample_counts = inputs$sample_data, sample_metadata = inputs$sample_info,
        sample_by = sample_by, condition_by = condition_by, cell_type_by = cell_type_by,
        contrast = inputs$contrast, design = design, backend_control = backend_control
      )
    } else if (!is_null(backend_control$result)) {
      backend_control$result
    } else {
      stop("scCODA requires `backend_control$runner` or `backend_control$result` from a scCODA/pertpy run.", call. = FALSE)
    }
    standardized <- .sn_standardize_sccoda_abundance(output)
    observed <- dplyr::rename(summaries$summary, observed_estimate = "estimate")
    primary <- dplyr::left_join(observed, standardized$table, by = "feature")
    raw <- standardized$raw
  } else {
    milo_args <- utils::modifyList(list(
      x = object,
      sample_by = sample_by,
      group_by = condition_by,
      contrast = inputs$contrast,
      covariates = if (is.character(design)) design else NULL,
      annotation_by = cell_type_by,
      store_name = NULL,
      return_object = FALSE,
      return_intermediate = FALSE
    ), backend_control$milo %||% list(), keep.null = TRUE)
    primary <- .sn_standardize_milo_abundance(do.call(sn_run_milo, milo_args), cell_type_by)
    raw <- primary
  }

  result <- list(
    schema_version = "1.0",
    analysis_type = "differential_abundance",
    name = store_name,
    method = method,
    backend = if (identical(method, "milo")) "miloR" else if (identical(method, "propeller")) "speckle" else if (identical(method, "sccoda")) "scCODA" else "Shennong",
    input = list(
      cells = nrow(inputs$metadata),
      samples = nrow(inputs$sample_info),
      sample_by = sample_by,
      condition_by = condition_by,
      cell_type_by = cell_type_by
    ),
    parameters = list(
      contrast = inputs$contrast,
      design = design,
      transform = if (identical(method, "propeller")) transform else NULL,
      permutations = if (identical(method, "permutation")) as.integer(permutations) else NULL
    ),
    tables = list(
      primary = primary,
      sample_proportions = inputs$sample_data,
      sample_contributions = summaries$contributions,
      permutation_null = null,
      backend_raw = raw,
      design = design_table
    ),
    embeddings = list(),
    graphs = list(),
    models = list(),
    diagnostics = list(
      inferential_unit = "sample",
      samples_per_group = as.list(inputs$group_counts),
      zero_count_rows = sum(inputs$sample_data$count == 0L)
    ),
    warnings = character(),
    provenance = .sn_analysis_provenance(random_seed = if (identical(method, "permutation")) seed else NA_integer_)
  )
  sn_validate_result(result)
  object <- sn_store_result(object, "differential_abundance", store_name, result)
  object <- .sn_log_seurat_command(object = object, name = "sn_test_abundance")
  if (isTRUE(return_object)) object else sn_get_result(object, "differential_abundance", store_name)
}
