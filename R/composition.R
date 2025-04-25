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
#' @importFrom rlang .data
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


#' calculate Startrac.dist (tissue distribution preference)
#' @import data.table
#' @importFrom plyr aaply
#' @importFrom stats chisq.test
#' @param dat.tb data.frame. Each line for a cell, and these columns as required: `majorCluster`, `loc`
#' @param byPatient logical. whether calculate the index for each patient. (default: FALSE)
#' @param colname.cluster character. which column specify the cluster (default: "majorCluster")
#' @param colname.patient character. which column specify the patient  (default: "patient")
#' @param colname.tissue character. which column specify the tissue  (default: "loc")
#' @param method character. method to use, one of "chisq", "fisher", and "freq"  (default: "chisq")
#' @param min.rowSum integer. rows with rowSum <= this value will be filtered out (default: 0)
#' @details calculate Startrac.dist (tissue distribution preference).
#' @return an array full of R_{o/e} (method="chisq") or list with components of OR, p.value etc. from fisher.test (method="fisher")
#' @export
sn_calculate_odds_ratio <- function(dat.tb,
                                    byPatient = F,
                                    colname.cluster = "majorCluster",
                                    colname.patient = "patient",
                                    colname.tissue = "loc",
                                    method = "chisq",
                                    min.rowSum = 0) {
  if (method == "freq") {
    ncount.sampleID <- dat.tb[, .(N.tot = .N), by = c(colname.patient, colname.tissue)]
    ncount.sampleID_mcls <- dat.tb[, .(N.mcls = .N), by = c(colname.patient, colname.tissue, colname.cluster)]
    ncount.sampleID_mcls <- reshape2::melt(
        reshape2::dcast(
            ncount.sampleID_mcls,
            as.formula(paste(colname.patient, "+", colname.tissue, "~", colname.cluster)),
            value.var = "N.mcls",
            fill = 0
        ),
        id.vars = c(colname.patient, colname.tissue),
        variable.name = colname.cluster,
        value.name = "N.mcls"
    )
    ncount.sampleID_mcls <- merge(ncount.sampleID_mcls, ncount.sampleID, by = c(colname.patient, colname.tissue))
    ncount.sampleID_mcls <- data.table::as.data.table(ncount.sampleID_mcls)
    ncount.sampleID_mcls[, freq.mcls := N.mcls / N.tot]
    res <- ncount.sampleID_mcls[,
      {
        loc.vec <- unique(.SD$loc)
        o.tb <- data.table::as.data.table(plyr::ldply(loc.vec, function(xx) {
          freq.x <- .SD[loc == xx, ][["freq.mcls"]]
          freq.y <- .SD[loc != xx, ][["freq.mcls"]]
          if (length(freq.x) >= 3) {
            oo.t <- t.test(freq.x, freq.y)
            oo.w <- wilcox.test(freq.x, freq.y)
            data.table(
              loc = xx,
              logFC = log2(mean(freq.x) / mean(freq.y)),
              p.value.t = oo.t$p.value,
              p.value.w = oo.w$p.value
            )
          } else {
            NULL
          }
        }))
      },
      by = c(colname.cluster)
    ][order(colname.cluster, colname.tissue), ]
    res[, FDR.t := p.adjust(p.value.t, "BH")]
    res[, FDR.w := p.adjust(p.value.w, "BH")]

    res[, char.sig := ""]
    res[FDR.t < 0.01 & p.value.t < 0.05, char.sig := "\U2020"]
    res[FDR.t < 0.05, char.sig := "\U2731"]
    res[FDR.t < 0.01, char.sig := "\U2731\U2731"]
    # res[FDR.t < 0.05,char.sig:="*"]
    # res[FDR.t < 0.01,char.sig:="**"]
    dist.charSig.tb <- dcast(res, majorCluster ~ loc, value.var = "char.sig")
    dist.logFC.tb <- dcast(res, majorCluster ~ loc, value.var = "logFC")
    dist.FDR.tb <- dcast(res, majorCluster ~ loc, value.var = "FDR.t")
    ret <- list(
      "res" = res,
      "dist.logFC.tb" = dist.logFC.tb,
      "dist.FDR.tb" = dist.FDR.tb,
      "dist.charSig.tb" = dist.charSig.tb
    )
  } else {
    if (byPatient == F) {
      N.o <- table(dat.tb[[colname.cluster]], dat.tb[[colname.tissue]])
      if (method == "chisq") {
        res.chisq <- chisq.test(N.o)
        R.oe <- (res.chisq$observed) / (res.chisq$expected)
        ret <- R.oe
      } else if (method == "fisher") {
        count.dist <- N.o[rowSums(N.o) > min.rowSum, , drop = F]
        count.dist.melt.ext.tb <- .table.fisher(count.dist)
        dist.p.tb <- dcast(count.dist.melt.ext.tb, rid ~ cid, value.var = "p.value")
        dist.padj.tb <- dcast(count.dist.melt.ext.tb, rid ~ cid, value.var = "p.adj")
        dist.OR.tb <- dcast(count.dist.melt.ext.tb, rid ~ cid, value.var = "OR")
        ret <- list(
          "count.dist" = count.dist.melt.ext.tb,
          "p.tb" = dist.p.tb,
          "padj.tb" = dist.padj.tb,
          "OR.tb" = dist.OR.tb
        )
      }
    } else {
      N.o.byPatient <- table(
        dat.tb[[colname.patient]],
        dat.tb[[colname.cluster]], dat.tb[[colname.tissue]]
      )
      ret <- aaply(N.o.byPatient, 1, function(x) {
        if (method == "chisq") {
          res.chisq <- chisq.test(x)
          return((res.chisq$observed) / (res.chisq$expected))
        } else {
          res.fisher <- .table.fisher(x)
          return(dcast(res.fisher, rid ~ cid, value.var = "OR"))
        }
      })
    }
  }
  return(ret)
}

.table.fisher <- function(count.dist) {
  sum.col <- colSums(count.dist)
  sum.row <- rowSums(count.dist)
  count.dist.tb <- as.data.frame(unclass(count.dist))
  data.table::setDT(count.dist.tb, keep.rownames = T)
  count.dist.melt.tb <- reshape2::melt(count.dist.tb, id.vars = "rn")
  colnames(count.dist.melt.tb) <- c("rid", "cid", "count")
  count.dist.melt.ext.tb <- data.table::as.data.table(plyr::ldply(seq_len(nrow(count.dist.melt.tb)), function(i) {
    this.row <- count.dist.melt.tb$rid[i]
    this.col <- count.dist.melt.tb$cid[i]
    this.c <- count.dist.melt.tb$count[i]
    other.col.c <- sum.col[this.col] - this.c
    this.m <- matrix(
      c(
        this.c,
        sum.row[this.row] - this.c,
        other.col.c,
        sum(sum.col) - sum.row[this.row] - other.col.c
      ),
      ncol = 2
    )
    res.test <- fisher.test(this.m)
    data.frame(
      rid = this.row,
      cid = this.col,
      p.value = res.test$p.value,
      OR = res.test$estimate
    )
  }))
  count.dist.melt.ext.tb <- merge(count.dist.melt.tb, count.dist.melt.ext.tb,
    by = c("rid", "cid")
  )
  count.dist.melt.ext.tb[, p.adj := p.adjust(p.value, "BH")]
  return(count.dist.melt.ext.tb)
}
