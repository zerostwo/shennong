# Cell-composition analysis.
#
# Extracted from analysis_metrics.R: composition proportions,
# observed-over-expected enrichment, sample-level composition comparison, and
# their private helpers.

#' Calculate composition proportions
#'
#' Calculates the proportion of different categories within groups from metadata.
#'
#' This function takes either a Seurat object or a data frame (like cell
#' metadata) and computes the percentage composition of a given \code{variable}
#' (for example, cell type) for each category specified by \code{group_by}
#' (for example, sample or condition). When \code{group_by} contains multiple
#' columns, proportions are calculated within each unique combination of those
#' grouping columns.
#'
#' @param x A Seurat object or a data frame.
#' @param group_by Column name or character vector of column names used to
#'   define groups.
#' @param variable Column name whose proportions should be calculated.
#' @param min_cells Minimum number of cells required for a returned
#'   \code{group_by + variable} category to be retained. Defaults to
#'   \code{20}.
#' @param measure One of \code{"proportion"}, \code{"count"}, or
#'   \code{"both"}. Defaults to \code{"proportion"}.
#' @param additional_cols Optional character vector of additional columns to
#'   carry into the output.
#' @param sort_by Optional summary column used to order the primary
#'   \code{group_by} column in the returned table. One of \code{NULL},
#'   \code{"proportion"}, \code{"count"}, or \code{"group_total"}.
#'   Sorting is only applied when \code{group_by} has length 1.
#' @param sort_value Optional level of \code{variable} used for sorting when
#'   \code{sort_by} is supplied. Defaults to the first observed level.
#' @param sort_desc Logical; if \code{TRUE}, sort in descending order.
#' @param object Alias for \code{x}; supply only one of \code{x} and \code{object}.
#' @return A data frame with per-group proportions.
#'
#' @importFrom dplyr count filter pull mutate rename left_join distinct across all_of group_by summarise first select n
#'
#' @examples
#' \dontrun{
#' composition_df <- sn_calculate_composition(
#'   x = seu,
#'   group_by = c("sample_id", "cell_type"),
#'   variable = "Phase",
#'   min_cells = 10,
#'   measure = "both"
#' )
#' composition_df <- sn_calculate_composition(
#'   x = seu,
#'   group_by = "sample_id",
#'   variable = "cell_type",
#'   min_cells = 10
#' )
#' print(composition_df)
#' }
#'
#' @export
sn_calculate_composition <- function(x,
                                     group_by,
                                     variable,
                                     min_cells = 20,
                                     measure = c("proportion", "count", "both"),
                                     additional_cols = NULL,
                                     sort_by = NULL,
                                     sort_value = NULL,
                                     sort_desc = FALSE,
                                     object = NULL) {
  x <- .sn_resolve_object_alias(x, object, missing(x))
  stopifnot(is.character(group_by), length(group_by) >= 1)
  stopifnot(is.character(variable), length(variable) == 1)
  stopifnot(is.numeric(min_cells), length(min_cells) == 1, min_cells >= 0)
  stopifnot(is.null(sort_by) || (is.character(sort_by) && length(sort_by) == 1))
  stopifnot(is.null(sort_value) || length(sort_value) == 1)
  stopifnot(is.logical(sort_desc), length(sort_desc) == 1)
  if (!is.null(additional_cols)) {
    stopifnot(is.character(additional_cols))
  }
  measure <- match.arg(measure)
  if (!is.null(sort_by)) {
    sort_by <- match.arg(sort_by, choices = c("proportion", "count", "group_total"))
    if (length(group_by) != 1L) {
      stop("`sort_by` is only supported when `group_by` contains a single column.", call. = FALSE)
    }
  }

  metadata <- .sn_extract_metric_metadata(x)

  required_cols <- c(group_by, variable)
  missing_cols <- setdiff(required_cols, colnames(metadata))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  if (!is.null(additional_cols)) {
    missing_additional <- setdiff(additional_cols, colnames(metadata))
    if (length(missing_additional) > 0) {
      stop("Missing additional columns: ", paste(missing_additional, collapse = ", "))
    }
  }

  group_levels <- if (is.factor(metadata[[group_by[[1]]]])) levels(metadata[[group_by[[1]]]]) else NULL

  keep <- stats::complete.cases(metadata[, c(group_by, variable), drop = FALSE])
  retained_metadata <- metadata[keep, , drop = FALSE]
  composition_full <- .sn_base_group_count(
    retained_metadata,
    c(group_by, variable),
    name = "count"
  )
  group_totals <- .sn_base_group_count(
    retained_metadata,
    group_by,
    name = "group_total"
  )
  composition_full <- composition_full[
    composition_full$count >= min_cells,
    ,
    drop = FALSE
  ]
  if (nrow(composition_full) > 0L) {
    group_id <- .sn_base_group_id(composition_full, group_by)
    total_group_id <- .sn_base_group_id(group_totals, group_by)
    composition_full$group_total <- group_totals$group_total[
      match(group_id, total_group_id)
    ]
    composition_full$proportion <-
      composition_full$count / composition_full$group_total * 100
  } else {
    composition_full$group_total <- numeric()
    composition_full$proportion <- numeric()
  }
  composition_full <- tibble::as_tibble(composition_full)

  if (nrow(composition_full) == 0) {
    stop("No groups remaining after filtering by `min_cells`. Consider lowering the threshold.")
  }

  if (!is.null(sort_by)) {
    composition_full <- .sn_sort_discrete_levels(
      data = composition_full,
      level_col = group_by[[1]],
      metric_col = sort_by,
      within_col = variable,
      within_value = sort_value,
      decreasing = sort_desc,
      fallback_levels = group_levels
    )
  }

  composition <- switch(
    measure,
    proportion = composition_full |>
      dplyr::select(dplyr::all_of(group_by), dplyr::all_of(variable), "proportion"),
    count = composition_full |>
      dplyr::select(dplyr::all_of(group_by), dplyr::all_of(variable), "count"),
    both = composition_full |>
      dplyr::select(dplyr::all_of(group_by), dplyr::all_of(variable), "count", "group_total", "proportion")
  )

  if (!is.null(additional_cols)) {
    metadata_group <- .sn_base_group_id(
      metadata,
      group_by,
      include_na = TRUE
    )
    metadata_indices <- split(seq_len(nrow(metadata)), metadata_group)
    inconsistent_cols <- additional_cols[vapply(additional_cols, function(col) {
      any(vapply(metadata_indices, function(index) {
        length(unique(metadata[[col]][index])) > 1L
      }, logical(1)))
    }, logical(1))]

    if (length(inconsistent_cols) > 0) {
      warning(
        "Additional columns are not constant within `group_by` groups; ",
        "using the first value for: ",
        paste(inconsistent_cols, collapse = ", "),
        call. = FALSE
      )
    }

    first_indices <- vapply(metadata_indices, `[[`, integer(1), 1L)
    additional_data <- metadata[first_indices, c(group_by, additional_cols), drop = FALSE]
    rownames(additional_data) <- NULL
    composition$.sn_order <- seq_len(nrow(composition))
    composition <- merge(
      as.data.frame(composition),
      additional_data,
      by = group_by,
      all.x = TRUE,
      sort = FALSE
    )
    composition <- composition[order(composition$.sn_order), , drop = FALSE]
    composition$.sn_order <- NULL
    composition <- tibble::as_tibble(composition)
  }

  composition
}

#' Calculate observed-over-expected enrichment
#'
#' Calculates observed-over-expected (RO/E) enrichment for a categorical
#' \code{variable} across metadata groups. RO/E is useful for asking whether a
#' cell type, cluster, state, or annotation is over-represented in a condition
#' relative to the marginal distribution of the full table.
#'
#' @param x A Seurat object or a data frame.
#' @param group_by Column name or character vector of column names used to
#'   define the rows/groups of the contingency table.
#' @param variable Column name defining the categories to test for enrichment.
#' @param pseudocount Numeric scalar added to both observed and expected counts
#'   before computing the ratio. Defaults to \code{0}.
#' @param return_matrix Logical; if \code{TRUE}, return a numeric matrix with
#'   \code{group_by} levels as rows and \code{variable} levels as columns.
#' @param matrix_value Value to place in the matrix when
#'   \code{return_matrix = TRUE}. One of \code{"roe"}, \code{"log2_roe"},
#'   \code{"observed"}, or \code{"expected"}.
#' @param object Alias for \code{x}; supply only one of \code{x} and \code{object}.
#' @return A data frame with \code{observed}, \code{expected}, totals,
#'   \code{roe}, and \code{log2_roe} columns. When \code{return_matrix = TRUE},
#'   a numeric matrix is returned.
#'
#' @examples
#' \dontrun{
#' roe_tbl <- sn_calculate_roe(
#'   seu,
#'   group_by = "condition",
#'   variable = "cell_type"
#' )
#' roe_mat <- sn_calculate_roe(
#'   seu,
#'   group_by = "condition",
#'   variable = "cell_type",
#'   return_matrix = TRUE
#' )
#' }
#'
#' @export
sn_calculate_roe <- function(x,
                             group_by,
                             variable,
                             pseudocount = 0,
                             return_matrix = FALSE,
                             matrix_value = c("roe", "log2_roe", "observed", "expected"),
                             object = NULL) {
  x <- .sn_resolve_object_alias(x, object, missing(x))
  stopifnot(is.character(group_by), length(group_by) >= 1L)
  stopifnot(is.character(variable), length(variable) == 1L)
  stopifnot(is.numeric(pseudocount), length(pseudocount) == 1L, pseudocount >= 0)
  stopifnot(is.logical(return_matrix), length(return_matrix) == 1L)
  matrix_value <- match.arg(matrix_value)

  metadata <- .sn_extract_metric_metadata(x)
  required_cols <- c(group_by, variable)
  missing_cols <- setdiff(required_cols, colnames(metadata))
  if (length(missing_cols) > 0L) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  metadata <- metadata |>
    dplyr::filter(
      dplyr::if_all(dplyr::all_of(group_by), ~ !is.na(.x)),
      !is.na(.data[[variable]])
    )

  if (nrow(metadata) == 0L) {
    stop("No cells remaining after removing missing `group_by` or `variable` values.", call. = FALSE)
  }
  if (!is.factor(metadata[[variable]])) {
    metadata[[variable]] <- as.character(metadata[[variable]])
  }

  group_data <- unique(metadata[, group_by, drop = FALSE])
  rownames(group_data) <- NULL

  variable_levels <- .sn_observed_discrete_levels(metadata[[variable]])
  full_grid <- merge(
    as.data.frame(group_data),
    data.frame(.sn_variable = variable_levels, stringsAsFactors = FALSE),
    by = NULL,
    sort = FALSE
  )
  colnames(full_grid)[colnames(full_grid) == ".sn_variable"] <- variable
  if (is.factor(metadata[[variable]])) {
    full_grid[[variable]] <- factor(full_grid[[variable]], levels = levels(metadata[[variable]]))
  }

  observed <- .sn_base_group_count(
    metadata,
    c(group_by, variable),
    name = "observed"
  )

  full_grid$.sn_order <- seq_len(nrow(full_grid))
  roe_tbl <- merge(
    full_grid,
    observed,
    by = c(group_by, variable),
    all.x = TRUE,
    sort = FALSE
  )
  roe_tbl <- roe_tbl[order(roe_tbl$.sn_order), , drop = FALSE]
  roe_tbl$.sn_order <- NULL
  roe_tbl$observed[is.na(roe_tbl$observed)] <- 0
  grand_total <- sum(roe_tbl$observed)
  row_id <- .sn_base_group_id(roe_tbl, group_by)
  column_id <- factor(roe_tbl[[variable]])
  roe_tbl$row_total <- stats::ave(roe_tbl$observed, row_id, FUN = sum)
  roe_tbl$col_total <- stats::ave(roe_tbl$observed, column_id, FUN = sum)
  roe_tbl$grand_total <- grand_total
  roe_tbl$expected <- roe_tbl$row_total * roe_tbl$col_total / grand_total
  denominator <- roe_tbl$expected + pseudocount
  roe_tbl$roe <- ifelse(
    denominator > 0,
    (roe_tbl$observed + pseudocount) / denominator,
    NA_real_
  )
  roe_tbl$log2_roe <- log2(roe_tbl$roe)
  roe_tbl <- tibble::as_tibble(roe_tbl)

  if (return_matrix) {
    return(.sn_roe_matrix(
      roe_tbl = roe_tbl,
      group_by = group_by,
      variable = variable,
      value = matrix_value
    ))
  }

  roe_tbl
}

.sn_observed_discrete_levels <- function(x) {
  observed <- x[!is.na(x)]
  if (is.factor(x)) {
    observed_chr <- unique(as.character(observed))
    return(levels(x)[levels(x) %in% observed_chr])
  }

  unique(as.character(observed))
}

.sn_roe_matrix <- function(roe_tbl, group_by, variable, value) {
  group_labels <- if (length(group_by) == 1L) {
    as.character(roe_tbl[[group_by]])
  } else {
    do.call(paste, c(as.data.frame(roe_tbl[group_by]), sep = "|"))
  }
  variable_labels <- as.character(roe_tbl[[variable]])

  group_levels <- unique(group_labels)
  variable_levels <- unique(variable_labels)
  mat <- matrix(
    NA_real_,
    nrow = length(group_levels),
    ncol = length(variable_levels),
    dimnames = list(group_levels, variable_levels)
  )
  mat[cbind(match(group_labels, group_levels), match(variable_labels, variable_levels))] <- roe_tbl[[value]]
  mat
}

.sn_complete_sample_composition <- function(composition, sample_info, sample_col, variable) {
  variable_levels <- if (is.factor(composition[[variable]])) {
    levels(composition[[variable]])
  } else {
    unique(as.character(composition[[variable]]))
  }

  full_grid <- merge(
    sample_info,
    data.frame(.sn_variable = variable_levels, stringsAsFactors = FALSE),
    by = NULL
  )
  colnames(full_grid)[colnames(full_grid) == ".sn_variable"] <- variable

  merged <- merge(
    full_grid,
    as.data.frame(composition),
    by = c(sample_col, variable, setdiff(colnames(sample_info), sample_col)),
    all.x = TRUE,
    sort = FALSE
  )

  merged$count[is.na(merged$count)] <- 0
  merged$proportion[is.na(merged$proportion)] <- 0
  sample_totals <- tapply(
    composition$group_total,
    as.character(composition[[sample_col]]),
    function(value) max(value, na.rm = TRUE)
  )
  missing_totals <- !is.finite(merged$group_total)
  merged$group_total[missing_totals] <- unname(sample_totals[
    as.character(merged[[sample_col]][missing_totals])
  ])
  merged
}

.sn_run_composition_test <- function(values_case, values_control, test = c("wilcox", "none")) {
  test <- match.arg(test)
  if (identical(test, "none")) {
    return(NA_real_)
  }

  values_case <- values_case[!is.na(values_case)]
  values_control <- values_control[!is.na(values_control)]

  if (length(values_case) < 2L || length(values_control) < 2L) {
    return(NA_real_)
  }

  if (all(values_case == values_case[[1]]) && all(values_control == values_control[[1]]) &&
      identical(unname(values_case[[1]]), unname(values_control[[1]]))) {
    return(1)
  }

  stats::wilcox.test(values_case, values_control, exact = FALSE)$p.value
}

#' Compare sample-level composition between groups
#'
#' Computes sample-level composition for a categorical variable and then compares
#' the resulting per-sample proportions between two groups. This avoids treating
#' individual cells as independent replicates and is therefore the recommended
#' way to estimate composition fold changes and significance when replicate
#' samples are available.
#'
#' @param x A Seurat object or a data frame containing cell-level metadata.
#' @param sample_by Column defining biological samples.
#' @param group_by Column defining the group or condition to compare between
#'   samples.
#' @param variable Column whose composition should be compared, for example cell
#'   type or cell-cycle phase.
#' @param contrast Character vector of length 2 giving the group levels as
#'   \code{c(case, control)}.
#' @param min_cells Minimum number of cells required per sample before that
#'   sample is retained for comparison. Defaults to \code{20}.
#' @param pseudocount Cell-count correction used only when either comparison
#'   group has zero mean abundance. The correction is converted to a
#'   sample-specific percentage using that sample's total cell count. Defaults
#'   to the Haldane-Anscombe value of \code{0.5} cells.
#' @param test Statistical test to apply to sample-level proportions. One of
#'   \code{"wilcox"} or \code{"none"}. Defaults to \code{"wilcox"}.
#' @param adjust_method Multiple-testing correction method passed to
#'   \code{stats::p.adjust()}. Defaults to \code{"BH"}.
#' @param additional_cols Optional sample-level columns to carry into the
#'   intermediate composition table before comparison.
#' @param return_sample_data Logical; if \code{TRUE}, return both the summary
#'   table and the completed sample-level composition table.
#'
#' @return A data frame with one row per \code{variable} level. When
#'   \code{return_sample_data = TRUE}, a list with \code{summary} and
#'   \code{sample_data} is returned.
#'
#' @examples
#' \dontrun{
#' comparison <- sn_compare_composition(
#'   seu,
#'   sample_by = "sample",
#'   group_by = "condition",
#'   variable = "cell_type",
#'   contrast = c("treated", "control")
#' )
#' }
#'
#' @export
sn_compare_composition <- function(x,
                                   sample_by = NULL,
                                   group_by = NULL,
                                   variable,
                                   contrast,
                                   min_cells = 20,
                                   pseudocount = 0.5,
                                   test = c("wilcox", "none"),
                                   adjust_method = "BH",
                                   additional_cols = NULL,
                                   return_sample_data = FALSE) {
  stopifnot(is.character(sample_by), length(sample_by) == 1L)
  stopifnot(is.character(group_by), length(group_by) == 1L)
  stopifnot(is.character(variable), length(variable) == 1L)
  if (!is.character(contrast) || length(contrast) != 2L || anyNA(contrast) ||
      any(!nzchar(contrast)) || anyDuplicated(contrast)) {
    stop("`contrast` must contain two distinct, non-missing group labels as `c(case, control)`.", call. = FALSE)
  }
  stopifnot(is.numeric(min_cells), length(min_cells) == 1L, min_cells >= 0)
  stopifnot(is.numeric(pseudocount), length(pseudocount) == 1L, pseudocount >= 0)
  stopifnot(is.character(adjust_method), length(adjust_method) == 1L)
  stopifnot(is.logical(return_sample_data), length(return_sample_data) == 1L)
  if (!is.null(additional_cols)) {
    stopifnot(is.character(additional_cols))
  }
  test <- match.arg(test)

  metadata <- .sn_extract_metric_metadata(x)
  required_cols <- c(sample_by, group_by, variable)
  missing_cols <- setdiff(required_cols, colnames(metadata))
  if (length(missing_cols) > 0L) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  .sn_validate_constant_within_sample(metadata, sample_col = sample_by, group_col = group_by)

  sample_info <- metadata |>
    dplyr::filter(!is.na(.data[[sample_by]]), !is.na(.data[[group_by]])) |>
    dplyr::group_by(.data[[sample_by]]) |>
    dplyr::summarise(
      .sn_sample_cells = dplyr::n(),
      dplyr::across(dplyr::all_of(unique(c(group_by, additional_cols))), dplyr::first),
      .groups = "drop"
    ) |>
    dplyr::filter(.data$.sn_sample_cells >= min_cells) |>
    dplyr::filter(.data[[group_by]] %in% contrast) |>
    dplyr::select(-".sn_sample_cells")

  if (nrow(sample_info) == 0L) {
    stop("No samples remaining for the requested contrast.", call. = FALSE)
  }
  observed_contrast <- unique(as.character(sample_info[[group_by]]))
  missing_contrast <- setdiff(contrast, observed_contrast)
  if (length(missing_contrast) > 0L) {
    stop(
      "The requested contrast is missing retained sample(s) for: ",
      paste(missing_contrast, collapse = ", "), ".",
      call. = FALSE
    )
  }

  composition <- sn_calculate_composition(
    x = metadata[metadata[[sample_by]] %in% sample_info[[sample_by]], , drop = FALSE],
    group_by = sample_by,
    variable = variable,
    min_cells = 0,
    measure = "both",
    additional_cols = unique(c(group_by, additional_cols))
  )

  sample_info <- sample_info[sample_info[[sample_by]] %in% unique(composition[[sample_by]]), , drop = FALSE]
  if (nrow(sample_info) == 0L) {
    stop("No samples remaining after filtering by `min_cells`.", call. = FALSE)
  }

  composition <- composition[composition[[sample_by]] %in% sample_info[[sample_by]], , drop = FALSE]
  composition_complete <- .sn_complete_sample_composition(
    composition = composition,
    sample_info = sample_info,
    sample_col = sample_by,
    variable = variable
  )

  comparison_levels <- if (is.factor(metadata[[variable]])) levels(metadata[[variable]]) else unique(as.character(composition_complete[[variable]]))
  comparison_levels <- comparison_levels[comparison_levels %in% unique(as.character(composition_complete[[variable]]))]
  n_comparison_levels <- length(comparison_levels)

  summary_tbl <- lapply(comparison_levels, function(current_level) {
    current_data <- composition_complete[as.character(composition_complete[[variable]]) %in% current_level, , drop = FALSE]
    case_data <- current_data[current_data[[group_by]] %in% contrast[[1]], , drop = FALSE]
    control_data <- current_data[current_data[[group_by]] %in% contrast[[2]], , drop = FALSE]

    mean_case <- mean(case_data$proportion, na.rm = TRUE)
    mean_control <- mean(control_data$proportion, na.rm = TRUE)
    median_case <- stats::median(case_data$proportion, na.rm = TRUE)
    median_control <- stats::median(control_data$proportion, na.rm = TRUE)
    corrected_case <- (case_data$count + pseudocount) /
      (case_data$group_total + pseudocount * n_comparison_levels) * 100
    corrected_control <- (control_data$count + pseudocount) /
      (control_data$group_total + pseudocount * n_comparison_levels) * 100
    mean_case_for_fc <- if (mean_case > 0 && mean_control > 0) {
      mean_case
    } else {
      mean(corrected_case, na.rm = TRUE)
    }
    mean_control_for_fc <- if (mean_case > 0 && mean_control > 0) {
      mean_control
    } else {
      mean(corrected_control, na.rm = TRUE)
    }
    p_value <- .sn_run_composition_test(
      values_case = case_data$proportion,
      values_control = control_data$proportion,
      test = test
    )

    data.frame(
      feature = current_level,
      mean_case = mean_case,
      mean_control = mean_control,
      median_case = median_case,
      median_control = median_control,
      difference = mean_case - mean_control,
      log2_fc = log2(mean_case_for_fc / mean_control_for_fc),
      n_case = nrow(case_data),
      n_control = nrow(control_data),
      p_value = p_value,
      stringsAsFactors = FALSE
    )
  })

  summary_tbl <- .sn_bind_rows(summary_tbl)
  colnames(summary_tbl)[colnames(summary_tbl) == "feature"] <- variable
  summary_tbl$contrast_case <- contrast[[1]]
  summary_tbl$contrast_control <- contrast[[2]]
  summary_tbl$change <- factor(
    ifelse(is.na(summary_tbl$log2_fc), NA_character_, ifelse(summary_tbl$log2_fc >= 0, "Increase", "Decrease")),
    levels = c("Increase", "Decrease")
  )
  summary_tbl$p_adj <- stats::p.adjust(summary_tbl$p_value, method = adjust_method)

  if (return_sample_data) {
    return(list(
      summary = summary_tbl,
      sample_data = composition_complete
    ))
  }

  summary_tbl
}

.sn_base_group_id <- function(data, columns, include_na = FALSE) {
  values <- lapply(data[, columns, drop = FALSE], function(value) {
    if (isTRUE(include_na)) {
      if (is.factor(value)) {
        factor(
          value,
          levels = c(levels(value), NA),
          exclude = NULL,
          ordered = is.ordered(value)
        )
      } else {
        factor(value, exclude = NULL)
      }
    } else if (is.factor(value)) {
      droplevels(value)
    } else {
      factor(value)
    }
  })
  do.call(
    interaction,
    c(values, list(drop = TRUE, lex.order = TRUE))
  )
}

.sn_base_group_count <- function(data, columns, name = "count") {
  if (nrow(data) == 0L) {
    result <- data[FALSE, columns, drop = FALSE]
    result[[name]] <- integer()
    return(result)
  }
  group_id <- .sn_base_group_id(data, columns)
  codes <- as.integer(group_id)
  first <- match(seq_len(nlevels(group_id)), codes)
  result <- data[first, columns, drop = FALSE]
  result[[name]] <- tabulate(codes, nbins = nlevels(group_id))
  rownames(result) <- NULL
  result
}
