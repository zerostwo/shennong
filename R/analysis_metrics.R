#' Calculate LISI score
#'
#' This function calculates the Local Intrinsic Dimensionality-based Outlier score (Lisi score) for a Seurat object.
#'
#' @param x A Seurat object.
#' @param reduction The dimensionality reduction method used to generate the embeddings (default is "pca").
#' @param label The column name in object@meta.data that specifies the sample labels (default is "sample").
#'
#' @return A data frame with the Lisi score for each cell, along with the cell ID.
#'
#' @importFrom tibble rownames_to_column
#'
#' @examples
#' # Calculate Lisi score for a Seurat object
#'
#' @export
sn_calculate_lisi <-
  function(x,
           reduction = "pca",
           label = "sample") {
    check_installed_github(pkg = "lisi", repo = "immunogenomics/lisi")
    check_installed("SeuratObject")
    lisi_score <- lisi::compute_lisi(
      X = SeuratObject::Embeddings(object = x, reduction = reduction),
      meta_data = x@meta.data,
      label_colnames = label
    ) |> rownames_to_column("cell_id")
    return(lisi_score)
  }

#' Calculate ROGUE score for Seurat Object
#'
#' This function calculates ROGUE score based on Seurat object.
#'
#' @param x A Seurat object.
#' @param cluster Column name in metadata specifying cluster labels.
#' @param sample Column name in metadata specifying sample labels.
#' @param span The span parameter for rogue estimation.
#' @param assay Assay used for ROGUE calculation. Defaults to \code{"RNA"}.
#' @param layer Layer used as the input count matrix. Defaults to \code{"counts"}.
#'
#' @return A data.frame of ROGUE score per cluster/sample or per cluster.
#'
#' @examples
#' \dontrun{
#' pbmc <- sn_load_data("pbmc1k")
#' pbmc <- sn_run_cluster(pbmc, normalization_method = "seurat", verbose = FALSE)
#' rogue_tbl <- sn_calculate_rogue(pbmc, cluster = "seurat_clusters")
#' head(rogue_tbl)
#' }
#' @export
sn_calculate_rogue <- function(
  x,
  cluster = NULL,
  sample = NULL,
  span = 0.9,
  assay = "RNA",
  layer = "counts"
) {
  check_installed_github(pkg = "ROGUE", repo = "PaulingLiu/ROGUE")
  if (!inherits(x, "Seurat")) {
    stop("Input x must be a Seurat object.")
  }

  counts <- .sn_get_seurat_layer_data(object = x, assay = assay, layer = layer)
  counts <- Matrix::as.matrix(counts)
  metadata <- x@meta.data

  counts <- ROGUE::matr.filter(counts, min.cells = 10, min.genes = 10)

  message("Calculating entropy...")
  entropy <- ROGUE::SE_fun(counts)

  message("Calculating ROGUE score...")
  rogue_result <- ROGUE::CalculateRogue(entropy, platform = "UMI")

  if (!is_null(cluster) && !is_null(sample)) {
    if (!cluster %in% colnames(metadata)) {
      stop("Specified cluster column not found in metadata.")
    }
    if (!sample %in% colnames(metadata)) {
      stop("Specified sample column not found in metadata.")
    }

    rogue_result <- ROGUE::rogue(
      counts,
      labels = as.character(metadata[[cluster]]),
      samples = as.character(metadata[[sample]]),
      platform = "UMI",
      span = span
    )

    rogue_result <- as.data.frame(rogue_result)
    rogue_result$cluster <- rownames(rogue_result)
    rownames(rogue_result) <- NULL
  }

  return(rogue_result)
}

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
#' composition_df <- sn_calculate_composition(
#'   x = seu,
#'   group_by = "sample_id",
#'   variable = "cell_type",
#'   min_cells = 10
#' )
#' print(composition_df)
#' }
sn_calculate_composition <- function(x, group_by, variable, min_cells = 20, additional_cols = NULL) {
  stopifnot(is.character(group_by), length(group_by) == 1)
  stopifnot(is.character(variable), length(variable) == 1)
  stopifnot(is.numeric(min_cells), length(min_cells) == 1, min_cells >= 0)
  if (!is.null(additional_cols)) stopifnot(is.character(additional_cols))

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

  metadata <- metadata |>
    dplyr::mutate(dplyr::across(dplyr::all_of(c(group_by, variable)), as.character)) |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_by))) |>
    dplyr::mutate(n_cells_group = dplyr::n()) |>
    dplyr::ungroup() |>
    dplyr::filter(.data$n_cells_group >= min_cells)

  if (nrow(metadata) == 0) {
    stop("No groups remaining after filtering by `min_cells`. Consider lowering the threshold.")
  }

  proportions <- metadata |>
    dplyr::filter(!is.na(.data[[group_by]]), !is.na(.data[[variable]])) |>
    dplyr::count(dplyr::across(dplyr::all_of(c(group_by, variable))), name = "count") |>
    dplyr::group_by(dplyr::across(dplyr::all_of(group_by))) |>
    dplyr::mutate(proportion = .data$count / sum(.data$count) * 100) |>
    dplyr::ungroup() |>
    dplyr::select(dplyr::all_of(group_by), dplyr::all_of(variable), dplyr::all_of("proportion"))

  if (!is.null(additional_cols)) {
    inconsistent_cols <- additional_cols[
      vapply(additional_cols, function(col) {
        any(
          metadata |>
            dplyr::group_by(dplyr::across(dplyr::all_of(group_by))) |>
            dplyr::summarise(
              is_inconsistent = length(unique(.data[[col]])) > 1,
              .groups = "drop"
            ) |>
            dplyr::pull(.data$is_inconsistent)
        )
      }, logical(1))
    ]

    if (length(inconsistent_cols) > 0) {
      warning(
        "Additional columns are not constant within `group_by` groups; ",
        "using the first value for: ",
        paste(inconsistent_cols, collapse = ", "),
        call. = FALSE
      )
    }

    additional_data <- metadata |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_by))) |>
      dplyr::summarise(dplyr::across(dplyr::all_of(additional_cols), dplyr::first), .groups = "drop")

    proportions <- proportions |>
      dplyr::left_join(additional_data, by = group_by)
  }

  proportions
}
