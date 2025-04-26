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
#' # Assuming 'seu' is a Seurat object with metadata columns
#' # 'sample_id', 'cell_type', and 'treatment'
#'
#' # Calculate cell type composition per sample
#' composition_df <- sn_calculate_composition(
#'   x = seu,
#'   group_by = "sample_id",
#'   variable = "cell_type",
#'   min_cells = 10
#' )
#'
#' # Calculate cell type composition per sample, adding the 'treatment' column
#' composition_with_treatment <- sn_calculate_composition(
#'   x = seu,
#'   group_by = "sample_id",
#'   variable = "cell_type",
#'   min_cells = 10,
#'   additional_cols = "treatment"
#' )
#'
#' # Using a plain data frame
#' meta_df <- data.frame(
#'  sample = rep(c("A", "B", "C"), each = 50),
#'  ctype = sample(c("Tcell", "Bcell", "Myeloid"), 150, replace = TRUE),
#'  condition = rep(c("Control", "Treated", "Control"), each = 50),
#'  stringsAsFactors = FALSE
#' )
#'
#' composition_from_df <- sn_calculate_composition(
#'   x = meta_df,
#'   group_by = "sample",
#'   variable = "ctype",
#'   min_cells = 5,
#'   additional_cols = "condition"
#' )
#' print(composition_from_df)
#' }
sn_calculate_composition <- function(x, group_by, variable, min_cells = 20, additional_cols = NULL) {
  # --- Input Validation ---
  if (!is.character(group_by) || length(group_by) != 1) {
    stop("`group_by` must be a single character string specifying a column name.")
  }
  if (!is.character(variable) || length(variable) != 1) {
    stop("`variable` must be a single character string specifying a column name.")
  }
  if (!is.numeric(min_cells) || length(min_cells) != 1 || min_cells < 0) {
    stop("`min_cells` must be a single non-negative number.")
  }
  if (!is.null(additional_cols) && !is.character(additional_cols)) {
     stop("`additional_cols` must be NULL or a character vector of column names.")
  }

  # --- Extract Metadata ---
  if (inherits(x, what = "Seurat")) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
        stop("Package 'SeuratObject' needed for this function to work with Seurat objects. Please install it.", call. = FALSE)
    }
    metadata <- x@meta.data
  } else if (is.data.frame(x)) {
    metadata <- as.data.frame(x) # Ensure it's a basic data frame
  } else {
     stop("Input `x` must be a Seurat object or a data frame.")
  }

  # --- Check required columns exist ---
  required_cols <- c(group_by, variable)
  missing_cols <- setdiff(required_cols, colnames(metadata))
  if (length(missing_cols) > 0) {
    stop("The following required columns are missing in the metadata: ",
         paste(missing_cols, collapse=", "))
  }
  if (!is.null(additional_cols)) {
      missing_additional <- setdiff(additional_cols, colnames(metadata))
       if (length(missing_additional) > 0) {
        stop("The following additional columns are missing in the metadata: ",
             paste(missing_additional, collapse=", "))
      }
  }

  # Ensure factor levels don't cause issues if columns are factors
  metadata[[group_by]] <- as.character(metadata[[group_by]])
  metadata[[variable]] <- as.character(metadata[[variable]])

  # --- Filter by Minimum Cells ---
  if (min_cells > 0) {
    cell_counts <- metadata |>
      dplyr::count(.data[[group_by]], name = "n_cells_group") # Use a unique name

    keep_group_bys <- cell_counts |>
      dplyr::filter(.data$n_cells_group >= min_cells) |>
      dplyr::pull(.data[[group_by]])

    original_groups <- unique(metadata[[group_by]])
    removed_groups <- setdiff(original_groups, keep_group_bys)

    if(length(removed_groups) > 0) {
        message("Groups removed due to having fewer than ", min_cells, " cells: ",
                paste(removed_groups, collapse=", "))
    }
     if(length(keep_group_bys) == 0) {
        stop("No groups remaining after filtering by `min_cells`. Consider lowering the threshold.")
    }

    metadata <- metadata |>
      dplyr::filter(.data[[group_by]] %in% keep_group_bys)
  }

  # --- Calculate Proportions ---
  # Use factors to ensure all levels of `variable` are present for each `group_by`
  # even if count is 0, before calculating proportions. Handle potential NAs.
  # NAs in group_by or variable columns would disrupt table creation, filter them.
  metadata_filtered <- metadata |>
      dplyr::filter(!is.na(.data[[group_by]]) & !is.na(.data[[variable]]))

  if (nrow(metadata_filtered) == 0) {
      warning("No valid rows remaining after filtering NAs in group_by or variable columns.")
      # Return an empty data frame with correct columns
       empty_df <- data.frame(matrix(ncol = 3 + length(additional_cols), nrow = 0))
       colnames(empty_df) <- c(group_by, variable, "proportion", additional_cols)
       return(empty_df)
  }


  # Create the contingency table and calculate proportions
  # Using factors ensures that all levels of `variable` are considered for each `group_by`,
  # potentially resulting in 0 proportions which is desired.
  prop_table <- table(
      factor(metadata_filtered[[group_by]]),
      factor(metadata_filtered[[variable]])
    ) |>
    prop.table(margin = 1) # Calculate row-wise proportions

  adata <- prop_table |>
    as.data.frame() |> # Converts table to Var1, Var2, Freq format
    dplyr::mutate(Freq = .data$Freq * 100) |>
    # stats::na.omit() # prop.table might produce NaN if a group has 0 cells *total* (filtered out by min_cells)
                        # Let's keep rows with 0 proportion, NAs shouldn't occur post min_cell filter unless issues upstream.
    dplyr::rename(
        "{group_by}" := "Var1",
        "{variable}" := "Var2",
        proportion = "Freq"
    ) |>
    # Convert factors back to characters for easier joining later
    dplyr::mutate(dplyr::across(c(.data[[group_by]], .data[[variable]]), as.character))


  # --- Add additional columns to the output if specified ---
  if (!is.null(additional_cols)) {
      # Select unique combinations of group_by and additional_cols
      # Warn if additional_cols are not unique per group_by level
      meta_summary <- metadata |>
          dplyr::select(dplyr::all_of(c(group_by, additional_cols))) |>
          dplyr::distinct()
      print(meta_summary)
      # # Check for non-uniqueness
      # non_unique_check <- meta_summary |>
      #     dplyr::count(.data[[group_by]]) |>
      #     dplyr::filter(dplyr::n() > 1)

      # if(nrow(non_unique_check) > 0) {
      #     warning("Values in `additional_cols` are not unique for the following groups: ",
      #             paste(non_unique_check[[group_by]], collapse=", "),
      #             ". Using the first encountered value for each group.")
      #     # Resolve non-uniqueness explicitly by taking the first value
      #     meta_summary <- metadata |>
      #         dplyr::select(dplyr::all_of(c(group_by, additional_cols))) |>
      #         dplyr::group_by(.data[[group_by]]) |>
      #         # Summarise using 'first' for each additional column
      #         dplyr::summarise(dplyr::across(dplyr::all_of(additional_cols), dplyr::first), .groups = "drop")
      # }

      # Join the additional columns to the proportion table
      adata <- dplyr::left_join(adata, meta_summary, by = group_by)
  }

  return(adata)
}
