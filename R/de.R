.sn_store_de_result <- function(object, store_name, result) {
  misc_data <- methods::slot(object, "misc")
  misc_data$de_results <- misc_data$de_results %||% list()
  misc_data$de_results[[store_name]] <- result
  methods::slot(object, "misc") <- misc_data
  object
}

.sn_run_seurat_de <- function(object,
                              analysis = c("markers", "contrast"),
                              ident_1 = NULL,
                              ident_2 = NULL,
                              group_by = NULL,
                              features = NULL,
                              assay = "RNA",
                              slot = "data",
                              test_use = "wilcox",
                              only_pos = NULL,
                              logfc_threshold = 0.1,
                              min_pct = 0.25,
                              verbose = TRUE,
                              ...) {
  analysis <- match.arg(analysis)

  if (analysis == "markers") {
    only_pos <- only_pos %||% TRUE
    result <- Seurat::FindAllMarkers(
      object = object,
      assay = assay,
      slot = slot,
      features = features,
      group.by = group_by,
      test.use = test_use,
      only.pos = only_pos,
      min.pct = min_pct,
      logfc.threshold = logfc_threshold,
      verbose = verbose,
      ...
    )
    return(result)
  }

  only_pos <- only_pos %||% FALSE
  Seurat::FindMarkers(
    object = object,
    ident.1 = ident_1,
    ident.2 = ident_2,
    group.by = group_by,
    assay = assay,
    slot = slot,
    features = features,
    test.use = test_use,
    only.pos = only_pos,
    min.pct = min_pct,
    logfc.threshold = logfc_threshold,
    verbose = verbose,
    ...
  )
}

.sn_run_pseudobulk_de <- function(object,
                                  ident_1,
                                  ident_2,
                                  group_by,
                                  sample_col,
                                  subset_by = NULL,
                                  subset_levels = NULL,
                                  assay = "RNA",
                                  layer = "counts",
                                  features = NULL,
                                  pseudobulk_method = c("DESeq2", "edgeR"),
                                  min_cells_per_sample = 10,
                                  verbose = TRUE) {
  pseudobulk_method <- match.arg(pseudobulk_method)
  metadata <- object[[]]

  if (!group_by %in% colnames(metadata)) {
    stop(glue("Column '{group_by}' was not found in metadata."))
  }
  if (!sample_col %in% colnames(metadata)) {
    stop(glue("Column '{sample_col}' was not found in metadata."))
  }
  if (!is_null(subset_by) && !subset_by %in% colnames(metadata)) {
    stop(glue("Column '{subset_by}' was not found in metadata."))
  }

  counts <- .sn_get_seurat_layer_data(object = object, assay = assay, layer = layer)
  if (!is_null(features)) {
    features <- intersect(features, rownames(counts))
    counts <- counts[features, , drop = FALSE]
  }

  subset_values <- subset_levels %||% if (is_null(subset_by)) "all" else unique(as.character(metadata[[subset_by]]))
  results <- vector("list", length(subset_values))
  names(results) <- subset_values

  for (current_subset in subset_values) {
    keep_cells <- if (is_null(subset_by)) {
      rep(TRUE, nrow(metadata))
    } else {
      as.character(metadata[[subset_by]]) == current_subset
    }

    meta_subset <- metadata[keep_cells, , drop = FALSE]
    counts_subset <- counts[, rownames(meta_subset), drop = FALSE]
    meta_subset$pb_group <- as.character(meta_subset[[group_by]])
    meta_subset$pb_sample <- as.character(meta_subset[[sample_col]])
    meta_subset$pb_key <- paste(meta_subset$pb_sample, meta_subset$pb_group, sep = "___")

    cell_totals <- table(meta_subset$pb_key)
    keep_keys <- names(cell_totals)[cell_totals >= min_cells_per_sample]
    if (length(keep_keys) == 0) {
      next
    }

    keep_key_cells <- meta_subset$pb_key %in% keep_keys
    meta_subset <- meta_subset[keep_key_cells, , drop = FALSE]
    counts_subset <- counts_subset[, rownames(meta_subset), drop = FALSE]

    if (!all(c(ident_1, ident_2) %in% unique(meta_subset$pb_group))) {
      next
    }

    aggregated <- rowsum(
      x = t(as.matrix(counts_subset)),
      group = meta_subset$pb_key,
      reorder = FALSE
    )
    aggregated <- t(aggregated)

    pb_meta <- unique(meta_subset[, c("pb_key", "pb_sample", "pb_group"), drop = FALSE])
    pb_meta <- pb_meta[match(colnames(aggregated), pb_meta$pb_key), , drop = FALSE]

    replicate_counts <- table(pb_meta$pb_group)
    if (any(replicate_counts[c(ident_1, ident_2)] < 2)) {
      if (verbose) {
        log_warn(glue(
          "Skipping pseudobulk DE for subset '{current_subset}' because fewer than 2 samples were available in at least one comparison group."
        ))
      }
      next
    }

    if (pseudobulk_method == "DESeq2") {
      check_installed("DESeq2")
      col_data <- data.frame(
        row.names = pb_meta$pb_key,
        sample = pb_meta$pb_sample,
        group = factor(pb_meta$pb_group, levels = c(ident_2, ident_1))
      )
      dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = round(aggregated),
        colData = col_data,
        design = ~ group
      )
      dds <- DESeq2::DESeq(dds, quiet = !verbose)
      result <- as.data.frame(
        DESeq2::results(dds, contrast = c("group", ident_1, ident_2))
      )
    } else {
      check_installed("edgeR")
      group_factor <- factor(pb_meta$pb_group, levels = c(ident_2, ident_1))
      dge <- edgeR::DGEList(counts = round(aggregated), group = group_factor)
      dge <- edgeR::calcNormFactors(dge)
      design <- stats::model.matrix(~ group_factor)
      dge <- edgeR::estimateDisp(dge, design = design)
      fit <- edgeR::glmQLFit(dge, design = design)
      test <- edgeR::glmQLFTest(fit, coef = 2)
      result <- edgeR::topTags(test, n = Inf, sort.by = "none")$table
      result$baseMean <- rowMeans(aggregated)
      result$log2FoldChange <- result$logFC
      result$pvalue <- result$PValue
      result$padj <- result$FDR
    }

    result$gene <- rownames(result)
    result$comparison <- paste(ident_1, "vs", ident_2)
    if (!is_null(subset_by)) {
      result[[subset_by]] <- current_subset
    }
    results[[current_subset]] <- tibble::as_tibble(result)
  }

  dplyr::bind_rows(results)
}

#' Run differential expression analysis on a Seurat object
#'
#' This function provides a single entry point for:
#' \itemize{
#'   \item marker discovery across all groups (`Seurat::FindAllMarkers()`),
#'   \item direct contrasts between two groups (`Seurat::FindMarkers()`), and
#'   \item pseudobulk contrasts aggregated by sample with `DESeq2` or `edgeR`.
#' }
#'
#' Results can be stored back into `object@misc$de_results[[store_name]]`, so
#' downstream helpers such as `sn_plot_dot()` can reuse the same marker table.
#'
#' @param object A \code{Seurat} object.
#' @param analysis One of \code{"markers"}, \code{"contrast"}, or
#'   \code{"pseudobulk"}. If \code{NULL}, the function infers the analysis from
#'   the other arguments.
#' @param ident_1,ident_2 Group labels for direct contrasts. These are required
#'   for \code{"contrast"} and \code{"pseudobulk"} analyses.
#' @param group_by Optional metadata column that defines the groups to compare.
#'   When omitted, Seurat identities are used.
#' @param subset_by Optional metadata column used to repeat a contrast within
#'   each subset, for example per cell type.
#' @param subset_levels Optional character vector of subset values to analyze.
#'   Defaults to all observed values in \code{subset_by}.
#' @param sample_col Metadata column containing sample IDs. Required for
#'   \code{"pseudobulk"} analyses.
#' @param assay Assay used for DE analysis. Defaults to \code{"RNA"}.
#' @param layer Count layer used for pseudobulk aggregation. Defaults to
#'   \code{"counts"}.
#' @param slot Seurat slot used for marker/contrast analyses. Defaults to
#'   \code{"data"}.
#' @param features Optional feature subset to test.
#' @param test_use Statistical test used by Seurat marker/contrast analyses.
#' @param pseudobulk_method One of \code{"DESeq2"} or \code{"edgeR"}.
#' @param only_pos Whether to return only positive markers. Defaults to
#'   \code{TRUE} for \code{"markers"} and \code{FALSE} otherwise.
#' @param logfc_threshold,min_pct Standard Seurat marker filtering arguments.
#' @param p_val_cutoff Adjusted p-value threshold used when storing result
#'   metadata.
#' @param de_logfc Absolute log fold-change threshold used when storing result
#'   metadata.
#' @param min_cells_per_sample Minimum cells required for a sample/group
#'   pseudobulk profile to be retained.
#' @param store_name Key used under \code{object@misc$de_results}.
#' @param return_object If \code{TRUE}, return the updated Seurat object with
#'   stored DE results. Otherwise return the result table.
#' @param verbose Whether to emit progress information.
#' @param ... Additional arguments passed through to Seurat DE functions.
#'
#' @return Either a DE result table or an updated \code{Seurat} object.
#' @export
sn_find_de <- function(
  object,
  analysis = NULL,
  ident_1 = NULL,
  ident_2 = NULL,
  group_by = NULL,
  subset_by = NULL,
  subset_levels = NULL,
  sample_col = NULL,
  assay = "RNA",
  layer = "counts",
  slot = "data",
  features = NULL,
  test_use = "wilcox",
  pseudobulk_method = c("DESeq2", "edgeR"),
  only_pos = NULL,
  logfc_threshold = 0.1,
  min_pct = 0.25,
  p_val_cutoff = 0.05,
  de_logfc = 0.25,
  min_cells_per_sample = 10,
  store_name = "default",
  return_object = TRUE,
  verbose = TRUE,
  ...
) {
  if (!inherits(object, "Seurat")) {
    stop("Input must be a Seurat object.")
  }

  if (is_null(analysis)) {
    analysis <- if (!is_null(sample_col)) {
      "pseudobulk"
    } else if (is_null(ident_1) && is_null(ident_2)) {
      "markers"
    } else {
      "contrast"
    }
  }
  analysis <- match.arg(analysis, c("markers", "contrast", "pseudobulk"))
  pseudobulk_method <- match.arg(pseudobulk_method)

  if (analysis %in% c("contrast", "pseudobulk") && (is_null(ident_1) || is_null(ident_2))) {
    stop("`ident_1` and `ident_2` are required for contrast and pseudobulk analyses.")
  }

  if (!is_null(subset_by) && !subset_by %in% colnames(object[[]])) {
    stop(glue("Column '{subset_by}' was not found in metadata."))
  }

  result <- NULL
  prepared <- NULL
  restore_context <- NULL

  if (analysis == "pseudobulk") {
    if (is_null(sample_col)) {
      stop("`sample_col` is required for pseudobulk analyses.")
    }
    pseudobulk_group_by <- group_by
    if (is_null(pseudobulk_group_by)) {
      pseudobulk_group_by <- ".shennong_ident"
      object[[pseudobulk_group_by]] <- as.character(Seurat::Idents(object))
    }
    result <- .sn_run_pseudobulk_de(
      object = object,
      ident_1 = ident_1,
      ident_2 = ident_2,
      group_by = pseudobulk_group_by,
      sample_col = sample_col,
      subset_by = subset_by,
      subset_levels = subset_levels,
      assay = assay,
      layer = layer,
      features = features,
      pseudobulk_method = pseudobulk_method,
      min_cells_per_sample = min_cells_per_sample,
      verbose = verbose
    )
  } else {
    analysis_object <- object
    if (identical(slot, "counts")) {
      prepared <- .sn_prepare_seurat_analysis_input(
        object = object,
        assay = assay,
        layer = layer
      )
      analysis_object <- prepared$object
      restore_context <- prepared$context
    } else {
      SeuratObject::DefaultAssay(analysis_object) <- assay
    }

    if (is_null(subset_by)) {
      result <- .sn_run_seurat_de(
        object = analysis_object,
        analysis = analysis,
        ident_1 = ident_1,
        ident_2 = ident_2,
        group_by = group_by,
        features = features,
        assay = assay,
        slot = slot,
        test_use = test_use,
        only_pos = only_pos,
        logfc_threshold = logfc_threshold,
        min_pct = min_pct,
        verbose = verbose,
        ...
      )
    } else {
      subset_values <- subset_levels %||% unique(as.character(analysis_object[[subset_by]][, 1]))
      subset_results <- vector("list", length(subset_values))
      names(subset_results) <- subset_values

      for (current_subset in subset_values) {
        subset_cells <- rownames(analysis_object[[]])[as.character(analysis_object[[subset_by]][, 1]) == current_subset]
        subset_object <- analysis_object[, subset_cells]
        current_result <- .sn_run_seurat_de(
          object = subset_object,
          analysis = analysis,
          ident_1 = ident_1,
          ident_2 = ident_2,
          group_by = group_by,
          features = features,
          assay = assay,
          slot = slot,
          test_use = test_use,
          only_pos = only_pos,
          logfc_threshold = logfc_threshold,
          min_pct = min_pct,
          verbose = verbose,
          ...
        )

        current_result <- tibble::rownames_to_column(
          as.data.frame(current_result),
          var = "gene"
        )
        current_result[[subset_by]] <- current_subset
        current_result$comparison <- if (analysis == "contrast") {
          paste(ident_1, "vs", ident_2)
        } else {
          NA_character_
        }
        subset_results[[current_subset]] <- tibble::as_tibble(current_result)
      }

      result <- dplyr::bind_rows(subset_results)
    }

    if (!is_null(restore_context)) {
      object <- .sn_restore_seurat_analysis_input(object = analysis_object, context = restore_context)
    } else {
      SeuratObject::DefaultAssay(object) <- SeuratObject::DefaultAssay(analysis_object)
    }
  }

  if (!is.data.frame(result)) {
    result <- tibble::rownames_to_column(as.data.frame(result), var = "gene")
  } else if (!"gene" %in% colnames(result)) {
    result <- tibble::rownames_to_column(as.data.frame(result), var = "gene")
  }

  group_col <- NULL
  if ("cluster" %in% colnames(result)) {
    group_col <- "cluster"
  } else if (!is_null(subset_by) && subset_by %in% colnames(result)) {
    group_col <- subset_by
  } else if (!is_null(group_by) && group_by %in% colnames(result)) {
    group_col <- group_by
  }

  rank_col <- c("avg_log2FC", "avg_logFC", "log2FoldChange", "logFC", "stat")
  rank_col <- rank_col[rank_col %in% colnames(result)][1] %||% NA_character_
  p_col <- c("p_val_adj", "padj", "FDR")
  p_col <- p_col[p_col %in% colnames(result)][1] %||% NA_character_

  stored_result <- list(
    schema_version = "1.0.0",
    package_version = as.character(utils::packageVersion("Shennong")),
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    table = tibble::as_tibble(result),
    analysis = analysis,
    group_by = group_by,
    group_col = group_col,
    ident_1 = ident_1,
    ident_2 = ident_2,
    subset_by = subset_by,
    sample_col = sample_col,
    assay = assay,
    layer = layer,
    slot = slot,
    rank_col = rank_col,
    p_col = p_col,
    p_val_cutoff = p_val_cutoff,
    de_logfc = de_logfc,
    min_pct = min_pct,
    logfc_threshold = logfc_threshold,
    n_genes = nrow(result),
    test_use = if (analysis == "pseudobulk") pseudobulk_method else test_use
  )

  object <- .sn_store_de_result(
    object = object,
    store_name = store_name,
    result = stored_result
  )

  if (return_object) {
    return(.sn_log_seurat_command(object = object, assay = assay, name = "sn_find_de"))
  }

  stored_result$table
}
