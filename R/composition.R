#' Calculate Composition Proportions
#'
#' Calculates the proportion of different categories within groups from metadata.
#'
#' This function takes either a Seurat object or a data frame (like cell metadata)
#' and computes the percentage composition of a given `variable` (e.g., cell type)
#' for each category specified by `group_by` (e.g., sample, condition).
#' It allows filtering out groups with fewer than a minimum number of cells
#' and optionally adding other metadata columns to the results.
#'
#' @param x A Seurat object or a data frame (e.g., metadata).
#' @param group_by A character string specifying the column name in the metadata
#'   to group by (e.g., "sample_id", "treatment").
#' @param variable A character string specifying the column name in the metadata
#'   for which to calculate proportions (e.g., "cell_type", "cluster_id").
#' @param min_cells An integer specifying the minimum number of cells required
#'   for a group (defined by `group_by`) to be included in the analysis.
#'   Defaults to 20. Set to 0 to include all groups.
#' @param additional_cols A character vector specifying the names of additional
#'   metadata columns to include in the output table. These columns should ideally
#'   have values that are constant within each `group_by` category. If values
#'   are not constant, the function will use the first value encountered for each group
#'   and issue a warning. Defaults to NULL.
#'
#' @return A data frame with the following columns:
#'   \item{group_by}{The grouping variable categories.}
#'   \item{variable}{The variable categories whose proportions are calculated.}
#'   \item{proportion}{The calculated proportion (percentage) for each variable category within each group.}
#'   \item{...}{Additional columns specified by `additional_cols`.}
#'
#' @importFrom dplyr count filter pull mutate rename left_join distinct across all_of group_by summarise first select n
#'
#' @export
#' @examples
#' \dontrun{
#' # Example 1: Using a Seurat object
#' # Assuming 'seu' is a Seurat object with metadata columns
#' # 'sample_id', 'cell_type', and 'treatment'
#' composition_df <- sn_calculate_composition(
#'   x = seu,
#'   group_by = "sample_id",
#'   variable = "cell_type",
#'   min_cells = 10
#' )
#' print(composition_df)
#'
#' # Example 2: Adding additional metadata columns
#' composition_with_treatment <- sn_calculate_composition(
#'   x = seu,
#'   group_by = "sample_id",
#'   variable = "cell_type",
#'   min_cells = 10,
#'   additional_cols = "treatment"
#' )
#' print(composition_with_treatment)
#'
#' # Example 3: Using a plain data frame
#' meta_df <- data.frame(
#'   sample = rep(c("A", "B", "C"), each = 50),
#'   ctype = sample(c("Tcell", "Bcell", "Myeloid"), 150, replace = TRUE),
#'   condition = rep(c("Control", "Treated", "Control"), each = 50),
#'   stringsAsFactors = FALSE
#' )
#' composition_from_df <- sn_calculate_composition(
#'   x = meta_df,
#'   group_by = "sample",
#'   variable = "ctype",
#'   min_cells = 5,
#'   additional_cols = "condition"
#' )
#' print(composition_from_df)
#'
#' # Example 4: Handling groups with fewer cells than `min_cells`
#' small_meta_df <- data.frame(
#'   sample = c("A", "A", "B"),
#'   ctype = c("Tcell", "Bcell", "Tcell"),
#'   stringsAsFactors = FALSE
#' )
#' composition_small <- sn_calculate_composition(
#'   x = small_meta_df,
#'   group_by = "sample",
#'   variable = "ctype",
#'   min_cells = 2
#' )
#' print(composition_small)
#' }
sn_calculate_composition <- function(x, group_by, variable, min_cells = 20, additional_cols = NULL) {
  # --- Input Validation ---
  stopifnot(is.character(group_by), length(group_by) == 1)
  stopifnot(is.character(variable), length(variable) == 1)
  stopifnot(is.numeric(min_cells), length(min_cells) == 1, min_cells >= 0)
  if (!is.null(additional_cols)) stopifnot(is.character(additional_cols))

  # --- Extract Metadata ---
  metadata <- if (inherits(x, "Seurat")) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
      stop("Package 'SeuratObject' is required for Seurat objects. Please install it.")
    }
    x@meta.data
  } else if (is.data.frame(x)) {
    x
  } else {
    stop("Input `x` must be a Seurat object or a data frame.")
  }

  # --- Check Required Columns ---
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

  # --- Filter by Minimum Cells ---
  metadata <- metadata %>%
    mutate(across(all_of(c(group_by, variable)), as.character)) %>%
    group_by(across(all_of(group_by))) %>%
    mutate(n_cells_group = n()) %>%
    ungroup() %>%
    filter(n_cells_group >= min_cells)

  if (nrow(metadata) == 0) {
    stop("No groups remaining after filtering by `min_cells`. Consider lowering the threshold.")
  }

  # --- Calculate Proportions ---
  proportions <- metadata %>%
    filter(!is.na(.data[[group_by]]), !is.na(.data[[variable]])) %>%
    count(across(all_of(c(group_by, variable))), name = "count") %>%
    group_by(across(all_of(group_by))) %>%
    mutate(proportion = count / sum(count) * 100) %>%
    ungroup() %>%
    select(all_of(group_by), all_of(variable), proportion)

  # --- Add Additional Columns ---
  if (!is.null(additional_cols)) {
    additional_data <- metadata %>%
      distinct(across(all_of(c(group_by, additional_cols)))) %>%
      group_by(across(all_of(group_by))) %>%
      summarise(across(all_of(additional_cols), first), .groups = "drop")

    proportions <- proportions %>%
      left_join(additional_data, by = group_by)
  }

  return(proportions)
}
