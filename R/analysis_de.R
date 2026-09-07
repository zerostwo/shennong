.sn_normalize_de_method <- function(analysis, method = NULL) {
  if (is_null(method)) {
    if (analysis == "pseudobulk") {
      return("DESeq2")
    }
    return("wilcox")
  }

  method_key <- tolower(method)

  if (analysis %in% c("markers", "contrast")) {
    if (method_key %in% c("cosgr", "cosg")) {
      return("COSGR")
    }

    return(method)
  }

  pseudobulk_methods <- c(
    deseq2 = "DESeq2",
    edger = "edgeR",
    limma = "limma"
  )

  if (!method_key %in% names(pseudobulk_methods)) {
    stop("For pseudobulk analyses, `method` must be one of 'DESeq2', 'edgeR', or 'limma'.")
  }

  pseudobulk_methods[[method_key]]
}

.sn_resolve_de_feature_universe <- function(object,
                                            assay,
                                            layer,
                                            features = NULL,
                                            method = NULL,
                                            verbose = TRUE) {
  .sn_validate_seurat_assay_layer(object = object, assay = assay, layer = layer)
  matched_layers <- .sn_match_seurat_layers(
    object = object,
    assay = assay,
    layer = layer
  )
  available <- unique(unlist(lapply(matched_layers, function(current_layer) {
    rownames(SeuratObject::LayerData(
      object = object,
      assay = assay,
      layer = current_layer
    ))
  }), use.names = FALSE))
  if (length(available) == 0L || anyNA(available) || any(!nzchar(available))) {
    stop(
      "The selected assay layer does not expose a valid feature universe for differential expression.",
      call. = FALSE
    )
  }

  requested <- NULL
  if (!is_null(features)) {
    if (!is.character(features)) {
      stop("`features` must be NULL or a character vector of feature names.", call. = FALSE)
    }
    supplied <- unique(features[!is.na(features) & nzchar(features)])
    requested <- intersect(supplied, available)
    missing_features <- setdiff(supplied, requested)
    if (length(missing_features) > 0L && isTRUE(verbose)) {
      .sn_log_warn(
        "`features` contains {length(missing_features)} feature(s) not present in assay '{assay}' layer '{layer}'; ",
        "ignoring examples: {paste(utils::head(missing_features, 5), collapse = ', ')}."
      )
    }
    if (length(requested) == 0L) {
      stop(
        "`features` did not contain any features present in assay '", assay,
        "' layer '", layer, "'.",
        call. = FALSE
      )
    }
  }

  cosgr_full_layer_test <- identical(method, "COSGR")
  tested <- if (is_null(requested) || cosgr_full_layer_test) available else requested
  list(
    requested_features = requested,
    tested_features = tested,
    tested_features_source = if (cosgr_full_layer_test && !is_null(requested)) {
      "assay_layer_full_test_then_result_filter"
    } else if (is_null(requested)) {
      "assay_layer"
    } else {
      "requested_features"
    }
  )
}

.sn_run_seurat_de <- function(object,
                              analysis = c("markers", "contrast"),
                              ident_1 = NULL,
                              ident_2 = NULL,
                              group_by = NULL,
                              features = NULL,
                              assay = "RNA",
                              layer = "data",
                              method = "wilcox",
                              only_pos = NULL,
                              logfc_threshold = 0.1,
                              min_pct = 0.25,
                              verbose = TRUE,
                              ...) {
  analysis <- match.arg(analysis)
  method <- .sn_normalize_de_method(analysis = analysis, method = method)
  target_layer <- .sn_guess_seurat_target_layer(layer)

  if (identical(method, "COSGR")) {
    if (!identical(analysis, "markers")) {
      stop("`method = 'COSGR'` is only supported when `analysis = 'markers'`.")
    }

    check_installed_github(pkg = "COSG", repo = "genecell/COSGR")
    prepared <- .sn_prepare_seurat_layer_alias(
      object = object,
      assay = assay,
      source_layer = layer,
      target_layer = target_layer
    )
    analysis_object <- prepared$object
    original_idents <- Seurat::Idents(analysis_object)

    on.exit({
      Seurat::Idents(analysis_object) <- original_idents
    }, add = TRUE)

    if (!is_null(group_by)) {
      Seurat::Idents(analysis_object) <- analysis_object[[group_by, drop = TRUE]]
    }

    cosg_result <- COSG::cosg(
      object = analysis_object,
      groups = "all",
      assay = assay,
      slot = target_layer,
      ...
    )

    result <- dplyr::bind_rows(lapply(colnames(cosg_result$names), function(current_group) {
      tibble::tibble(
        cluster = current_group,
        gene = as.character(cosg_result$names[[current_group]]),
        cosg_score = as.numeric(cosg_result$scores[[current_group]]),
        rank = seq_len(nrow(cosg_result$names))
      )
    }))

    if (!is_null(features)) {
      result <- dplyr::filter(result, .data$gene %in% features)
    }

    return(result)
  }

  if (analysis == "markers") {
    only_pos <- only_pos %||% TRUE
    result <- .sn_with_default_seurat_acceleration(
      Seurat::FindAllMarkers(
        object = object,
        assay = assay,
        slot = target_layer,
        features = features,
        group.by = group_by,
        test.use = method,
        only.pos = only_pos,
        min.pct = min_pct,
        logfc.threshold = logfc_threshold,
        verbose = verbose,
        ...
      ),
      object = object,
      assay = assay
    )
    return(result)
  }

  only_pos <- only_pos %||% FALSE
  .sn_with_default_seurat_acceleration(
    Seurat::FindMarkers(
      object = object,
      ident.1 = ident_1,
      ident.2 = ident_2,
      group.by = group_by,
      assay = assay,
      slot = target_layer,
      features = features,
      test.use = method,
      only.pos = only_pos,
      min.pct = min_pct,
      logfc.threshold = logfc_threshold,
      verbose = verbose,
      ...
    ),
    object = object,
    assay = assay
  )
}

.sn_validate_pseudobulk_count_layer <- function(counts,
                                                layer,
                                                require_integer = TRUE) {
  if (!.sn_name_declares_count_scale(layer)) {
    stop(
      "Pseudobulk differential expression requires a raw or corrected count layer; normalized expression layers are not supported.",
      call. = FALSE
    )
  }
  counts <- .sn_as_sparse_matrix(counts)
  stored_values <- if (inherits(counts, "sparseMatrix")) counts@x else as.numeric(counts)
  if (anyNA(stored_values) || any(!is.finite(stored_values))) {
    stop("The pseudobulk count layer contains missing or non-finite values.", call. = FALSE)
  }
  if (any(stored_values < 0)) {
    stop("The pseudobulk count layer contains negative values.", call. = FALSE)
  }
  if (isTRUE(require_integer) &&
      any(abs(stored_values - round(stored_values)) > 1e-8)) {
    stop(
      "DESeq2 pseudobulk analysis requires integer-valued raw or corrected counts; values will not be rounded silently.",
      call. = FALSE
    )
  }
  counts
}

.sn_validate_de_subset_levels <- function(metadata, subset_by, subset_levels) {
  if (is_null(subset_by)) {
    if (!is_null(subset_levels)) {
      stop("`subset_levels` can only be supplied together with `subset_by`.", call. = FALSE)
    }
    return(invisible(NULL))
  }
  if (!is.character(subset_by) || length(subset_by) != 1L ||
      is.na(subset_by) || !nzchar(subset_by)) {
    stop("`subset_by` must be one non-empty metadata column name.", call. = FALSE)
  }
  if (!subset_by %in% colnames(metadata)) {
    stop(glue("Column '{subset_by}' was not found in metadata."), call. = FALSE)
  }
  if (is_null(subset_levels)) {
    return(invisible(NULL))
  }
  if (!is.character(subset_levels) || length(subset_levels) == 0L ||
      anyNA(subset_levels) || any(!nzchar(subset_levels)) ||
      anyDuplicated(subset_levels)) {
    stop("`subset_levels` must be a non-empty character vector of distinct, non-missing subset values.", call. = FALSE)
  }
  observed <- unique(as.character(metadata[[subset_by]]))
  missing_levels <- setdiff(subset_levels, observed)
  if (length(missing_levels) > 0L) {
    stop(
      "`subset_levels` value(s) were not observed in `", subset_by, "`: ",
      paste(missing_levels, collapse = ", "), ".",
      call. = FALSE
    )
  }
  invisible(NULL)
}

.sn_pseudobulk_profile_keys <- function(sample, group) {
  profiles <- data.frame(
    sample = as.character(sample),
    group = as.character(group),
    stringsAsFactors = FALSE
  )
  unique_profiles <- profiles[!duplicated(profiles), , drop = FALSE]
  profile_index <- integer(nrow(profiles))
  for (index in seq_len(nrow(unique_profiles))) {
    profile_index[
      profiles$sample == unique_profiles$sample[[index]] &
        profiles$group == unique_profiles$group[[index]]
    ] <- index
  }
  if (any(profile_index == 0L)) {
    stop("Internal error while constructing pseudobulk profile identifiers.", call. = FALSE)
  }
  paste0("profile_", profile_index)
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
                                  method = c("DESeq2", "edgeR", "limma"),
                                  min_cells_per_sample = 10,
                                  design = NULL,
                                  contrast = NULL,
                                  verbose = TRUE) {
  method <- match.arg(method)
  metadata <- object[[]]
  if (!is.character(ident_1) || length(ident_1) != 1L ||
      !is.character(ident_2) || length(ident_2) != 1L ||
      is.na(ident_1) || is.na(ident_2) || identical(ident_1, ident_2)) {
    stop("`ident_1` and `ident_2` must be distinct, non-missing group labels.", call. = FALSE)
  }
  if (!is.numeric(min_cells_per_sample) || length(min_cells_per_sample) != 1L ||
      is.na(min_cells_per_sample) || min_cells_per_sample < 1L ||
      min_cells_per_sample != as.integer(min_cells_per_sample)) {
    stop("`min_cells_per_sample` must be a positive integer.", call. = FALSE)
  }
  min_cells_per_sample <- as.integer(min_cells_per_sample)
  if (!is_null(design) && !inherits(design, "formula")) {
    stop("For pseudobulk analyses, `design` must be a formula when supplied.", call. = FALSE)
  }
  if (!is_null(design) && .sn_bulk_has_random_effect(design)) {
    stop("Pseudobulk DESeq2, edgeR, and limma designs cannot contain random effects.", call. = FALSE)
  }

  design_variables <- if (inherits(design, "formula")) all.vars(design) else character()
  required_columns <- unique(c(group_by, sample_col, subset_by, design_variables))
  missing_columns <- setdiff(required_columns, colnames(metadata))
  if (length(missing_columns) > 0L) {
    stop("Pseudobulk metadata column(s) missing: ", paste(missing_columns, collapse = ", "), ".", call. = FALSE)
  }

  counts <- .sn_get_seurat_layer_data(object = object, assay = assay, layer = layer)
  counts <- .sn_validate_pseudobulk_count_layer(
    counts,
    layer = layer,
    require_integer = identical(method, "DESeq2")
  )
  metadata <- metadata[colnames(counts), , drop = FALSE]
  incomplete <- !stats::complete.cases(metadata[, required_columns, drop = FALSE])
  empty_identifiers <- !nzchar(as.character(metadata[[group_by]])) |
    !nzchar(as.character(metadata[[sample_col]]))
  if (any(incomplete | empty_identifiers)) {
    stop(
      "Pseudobulk sample, group, subset, and design metadata must be complete and non-empty.",
      call. = FALSE
    )
  }

  if (!is_null(features)) {
    features <- intersect(features, rownames(counts))
    if (length(features) == 0L) {
      stop("None of the requested pseudobulk features were found in the selected assay layer.", call. = FALSE)
    }
    counts <- counts[features, , drop = FALSE]
  }

  contrast <- contrast %||% c(group_by, ident_1, ident_2)
  if (!is.character(contrast) || length(contrast) != 3L ||
      !identical(contrast[[1]], group_by) ||
      !identical(contrast[[2]], ident_1) ||
      !identical(contrast[[3]], ident_2)) {
    stop(
      "For pseudobulk analyses, `contrast` must be c(group_by, ident_1, ident_2).",
      call. = FALSE
    )
  }
  if (!is_null(design) && !group_by %in% design_variables) {
    stop("The pseudobulk `design` formula must include the `group_by` variable.", call. = FALSE)
  }

  subset_values <- subset_levels %||% if (is_null(subset_by)) {
    "all"
  } else {
    observed_subsets <- unique(as.character(metadata[[subset_by]]))
    observed_subsets[!is.na(observed_subsets) & nzchar(observed_subsets)]
  }
  results <- vector("list", length(subset_values))
  names(results) <- subset_values
  contrast_estimands <- vector("list", length(subset_values))
  names(contrast_estimands) <- subset_values

  for (current_subset in subset_values) {
    keep_cells <- if (is_null(subset_by)) {
      rep(TRUE, nrow(metadata))
    } else {
      as.character(metadata[[subset_by]]) == current_subset
    }
    keep_cells <- keep_cells & as.character(metadata[[group_by]]) %in% c(ident_1, ident_2)

    meta_subset <- metadata[keep_cells, , drop = FALSE]
    counts_subset <- counts[, rownames(meta_subset), drop = FALSE]
    meta_subset$.pb_group <- as.character(meta_subset[[group_by]])
    meta_subset$.pb_sample <- as.character(meta_subset[[sample_col]])
    meta_subset$.pb_key <- .sn_pseudobulk_profile_keys(
      meta_subset$.pb_sample,
      meta_subset$.pb_group
    )

    cell_totals <- table(meta_subset$.pb_key)
    keep_keys <- names(cell_totals)[cell_totals >= min_cells_per_sample]
    keep_key_cells <- meta_subset$.pb_key %in% keep_keys
    meta_subset <- meta_subset[keep_key_cells, , drop = FALSE]
    counts_subset <- counts_subset[, rownames(meta_subset), drop = FALSE]
    if (nrow(meta_subset) == 0L ||
        !all(c(ident_1, ident_2) %in% unique(meta_subset$.pb_group))) {
      next
    }

    profiles <- unique(meta_subset[, c(".pb_sample", ".pb_group"), drop = FALSE])
    groups_by_sample <- split(profiles$.pb_group, profiles$.pb_sample)
    groups_by_sample <- lapply(groups_by_sample, unique)
    repeated_samples <- vapply(groups_by_sample, length, integer(1)) > 1L
    auto_design <- is_null(design)
    if (auto_design) {
      independent <- all(!repeated_samples)
      paired <- all(vapply(
        groups_by_sample,
        function(groups) setequal(groups, c(ident_1, ident_2)),
        logical(1)
      ))
      if (!independent && !paired) {
        stop(
          "Pseudobulk samples have a mixture of paired and unpaired group profiles. Supply an explicit full-rank design or use unique biological sample IDs.",
          call. = FALSE
        )
      }
      if (paired && length(groups_by_sample) < 2L) {
        next
      }
      design_formula <- if (paired) ~.pb_sample + .pb_group else ~.pb_group
      current_contrast <- c(".pb_group", ident_1, ident_2)
    } else {
      if (any(repeated_samples) && !sample_col %in% design_variables) {
        stop(
          "The pseudobulk `design` must include `sample_by` when the same biological unit contributes both contrast groups.",
          call. = FALSE
        )
      }
      design_formula <- design
      current_contrast <- contrast
    }

    profile_columns <- unique(c(sample_col, group_by, design_variables))
    for (column in profile_columns) {
      values_per_profile <- split(meta_subset[[column]], meta_subset$.pb_key)
      nonconstant <- vapply(values_per_profile, function(values) {
        length(unique(values)) != 1L
      }, logical(1))
      if (any(nonconstant)) {
        stop(
          "Pseudobulk design variable '", column,
          "' is not constant within aggregated sample/group profiles.",
          call. = FALSE
        )
      }
    }

    aggregated <- .sn_aggregate_columns_by_group(
      x = counts_subset,
      groups = meta_subset$.pb_key
    )
    pb_meta <- meta_subset[!duplicated(meta_subset$.pb_key), , drop = FALSE]
    pb_meta <- pb_meta[match(colnames(aggregated), pb_meta$.pb_key), , drop = FALSE]
    rownames(pb_meta) <- pb_meta$.pb_key

    replicate_counts <- table(factor(pb_meta$.pb_group, levels = c(ident_1, ident_2)))
    if (any(replicate_counts < 2L)) {
      if (verbose) {
        .sn_log_warn(
          "Skipping pseudobulk DE for subset '{current_subset}' because fewer than 2 biological samples were available in at least one comparison group."
        )
      }
      next
    }

    formula_variables <- all.vars(design_formula)
    for (column in formula_variables) {
      if (is.character(pb_meta[[column]]) || is.logical(pb_meta[[column]])) {
        pb_meta[[column]] <- factor(pb_meta[[column]])
      }
    }
    contrast_variable <- current_contrast[[1]]
    pb_meta[[contrast_variable]] <- factor(
      as.character(pb_meta[[contrast_variable]]),
      levels = c(current_contrast[[3]], current_contrast[[2]])
    )
    .sn_bulk_validate_contrast(pb_meta, current_contrast, design = design_formula)
    design_info <- .sn_bulk_design(pb_meta, design_formula)
    if (.sn_bulk_has_random_effect(design_formula)) {
      stop(
        "Pseudobulk DE does not support random-effect designs with DESeq2, edgeR, or limma; ",
        "encode a supported fixed effect or use `sn_find_bulk_de(method = \"dream\")` ",
        "on an explicit bulk matrix.",
        call. = FALSE
      )
    }
    contrast_estimand <- .sn_bulk_contrast_estimand(
      pb_meta,
      design_formula,
      current_contrast
    )
    contrast_vector <- contrast_estimand$vector
    contrast_estimands[[as.character(current_subset)]] <- contrast_estimand

    if (identical(method, "DESeq2")) {
      check_installed("DESeq2")
      dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = aggregated,
        colData = pb_meta,
        design = design_formula
      )
      dds <- DESeq2::DESeq(dds, quiet = !verbose)
      coefficient_names <- DESeq2::resultsNames(dds)
      if (length(contrast_vector) != length(coefficient_names)) {
        stop(
          "DESeq2 coefficient count does not match the validated pseudobulk contrast vector.",
          call. = FALSE
        )
      }
      result <- as.data.frame(DESeq2::results(
        dds,
        contrast = unname(contrast_vector)
      ))
    } else if (identical(method, "edgeR")) {
      check_installed("edgeR")
      dge <- edgeR::DGEList(counts = aggregated)
      dge <- edgeR::calcNormFactors(dge)
      dge <- edgeR::estimateDisp(dge, design = design_info$matrix)
      fit <- edgeR::glmQLFit(dge, design = design_info$matrix)
      test <- edgeR::glmQLFTest(fit, contrast = contrast_vector)
      result <- edgeR::topTags(test, n = Inf, sort.by = "none")$table
      result$baseMean <- Matrix::rowMeans(aggregated)
      result$log2FoldChange <- result$logFC
      result$pvalue <- result$PValue
      result$padj <- result$FDR
    } else {
      check_installed(c("limma", "edgeR"))
      dge <- edgeR::DGEList(counts = aggregated)
      dge <- edgeR::calcNormFactors(dge)
      transformed <- limma::voom(dge, design = design_info$matrix, plot = FALSE)
      fit <- limma::lmFit(transformed, design = design_info$matrix)
      fit <- limma::contrasts.fit(
        fit,
        contrasts = matrix(contrast_vector, ncol = 1L)
      )
      fit <- limma::eBayes(fit)
      result <- limma::topTable(fit, coef = 1L, number = Inf, sort.by = "none")
      result$baseMean <- Matrix::rowMeans(aggregated)
      result$log2FoldChange <- result$logFC
      result$pvalue <- result$P.Value
      result$padj <- result$adj.P.Val
    }

    result$gene <- rownames(result)
    result$comparison <- paste(ident_1, "vs", ident_2)
    if (!is_null(subset_by)) {
      result[[subset_by]] <- current_subset
    }
    results[[current_subset]] <- tibble::as_tibble(result)
  }

  output <- dplyr::bind_rows(results)
  if (nrow(output) == 0L) {
    stop(
      "No pseudobulk comparison retained at least two biological replicates per contrast group.",
      call. = FALSE
    )
  }
  attr(output, "shennong_contrast_estimands") <- Filter(Negate(is.null), contrast_estimands)
  output
}

#' Run differential expression analysis on a Seurat object
#'
#' This function provides a single entry point for:
#' \itemize{
#'   \item marker discovery across all groups,
#'   \item direct contrasts between two groups, and
#'   \item pseudobulk contrasts aggregated by sample.
#' }
#'
#' Results are stored in the canonical Shennong result registry, so
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
#' @param sample_by Metadata column containing sample IDs. Required for
#'   \code{"pseudobulk"} analyses.
#' @param assay Assay used for DE analysis. Defaults to \code{"RNA"}.
#' @param layer Assay layer used for DE analysis. Defaults to \code{"data"} for
#'   marker and contrast analyses and to \code{"counts"} for pseudobulk
#'   analyses.
#' @param features Optional feature subset to test. The resolved feature set is
#'   retained as \code{input$tested_features} in stored single-cell DE results,
#'   so downstream stored-DE ORA can reuse the actual hypothesis universe.
#' @param modality Input modality. \code{"auto"} selects single-cell analysis
#'   for Seurat objects and bulk analysis for matrices, lists, and
#'   \code{SummarizedExperiment} objects.
#' @param metadata Optional sample metadata for bulk analysis.
#' @param design A fixed- or mixed-effects formula for bulk analysis. For
#'   pseudobulk analysis, an optional fixed-effects formula over cell metadata;
#'   include \code{sample_by} when the same biological unit contributes both
#'   contrast groups. When omitted, independent and completely paired designs
#'   are detected automatically.
#' @param contrast Character triple giving the contrast as variable, numerator,
#'   and denominator. For pseudobulk analysis this must agree with
#'   \code{c(group_by, ident_1, ident_2)}.
#' @param backend_control Bulk backend controls or a custom bulk
#'   \code{runner}/precomputed \code{result}.
#' @param method Statistical method. For \code{"markers"} and
#'   \code{"contrast"}, this can be any Seurat \code{test.use} value or
#'   \code{"COSGR"} for marker discovery. For \code{"pseudobulk"}, choose one
#'   of \code{"DESeq2"}, \code{"edgeR"}, or \code{"limma"}. For bulk input,
#'   choose \code{"auto"}, \code{"edger"}, \code{"deseq2"}, \code{"limma"},
#'   or \code{"dream"}.
#' @param only_pos Whether to return only positive markers. Defaults to
#'   \code{TRUE} for \code{"markers"} and \code{FALSE} otherwise.
#' @param logfc_threshold,min_pct Standard Seurat marker filtering arguments.
#' @param p_val_cutoff Adjusted p-value threshold used when storing result
#'   metadata.
#' @param de_logfc Absolute log fold-change threshold used when storing result
#'   metadata.
#' @param min_cells_per_sample Minimum cells required for a sample/group
#'   pseudobulk profile to be retained.
#' @param result_id Stable identifier used as the stored-result key and canonical
#'   result name.
#' @param return_object If \code{TRUE}, return the updated Seurat object with
#'   stored DE results. Otherwise return the result table.
#' @param verbose Whether to emit progress information.
#' @param ... Additional arguments passed through to the selected DE method.
#'
#' @return For single-cell input, either a DE result table or an updated
#'   \code{Seurat} object. For bulk input, a validated Shennong bulk-DE result.
#'
#' @examples
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'   counts <- matrix(rpois(20 * 24, lambda = 1), nrow = 20, ncol = 24)
#'   rownames(counts) <- c(
#'     paste0("GENE", 1:14),
#'     "CD3D", "CD3E", "TRAC", "MS4A1", "CD79A", "HLA-DRA"
#'   )
#'   colnames(counts) <- paste0("cell", 1:24)
#'   counts[c("CD3D", "CD3E", "TRAC"), 1:12] <-
#'     counts[c("CD3D", "CD3E", "TRAC"), 1:12] + 20
#'   counts[c("MS4A1", "CD79A", "HLA-DRA"), 13:24] <-
#'     counts[c("MS4A1", "CD79A", "HLA-DRA"), 13:24] + 20
#'
#'   obj <- sn_initialize_seurat_object(counts, species = "human")
#'   obj$cell_type <- rep(c("Tcell", "Bcell"), each = 12)
#'   Seurat::Idents(obj) <- obj$cell_type
#'   obj <- Seurat::NormalizeData(obj, verbose = FALSE)
#'
#'   marker_tbl <- sn_find_de(
#'     obj,
#'     analysis = "markers",
#'     group_by = "cell_type",
#'     layer = "data",
#'     min_pct = 0,
#'     logfc_threshold = 0,
#'     return_object = FALSE,
#'     verbose = FALSE
#'   )
#'   head(marker_tbl)
#'
#'   obj <- sn_find_de(
#'     obj,
#'     analysis = "markers",
#'     group_by = "cell_type",
#'     layer = "data",
#'     min_pct = 0,
#'     logfc_threshold = 0,
#'     result_id = "celltype_markers",
#'     return_object = TRUE,
#'     verbose = FALSE
#'   )
#'   sn_list_results(obj, type = "de")
#'   sn_get_result(obj, type = "de", result_id = "celltype_markers")
#' }
#'
#' bulk_counts <- matrix(
#'   stats::rpois(40 * 6, 20), nrow = 40,
#'   dimnames = list(paste0("gene_", 1:40), paste0("sample_", 1:6))
#' )
#' bulk_metadata <- data.frame(
#'   condition = factor(rep(c("control", "treated"), each = 3)),
#'   row.names = colnames(bulk_counts)
#' )
#' bulk_de <- sn_find_de(
#'   bulk_counts,
#'   metadata = bulk_metadata,
#'   design = ~condition,
#'   contrast = c("condition", "treated", "control"),
#'   backend_control = list(result = data.frame(
#'     gene = rownames(bulk_counts), logFC = 0, PValue = 1, FDR = 1
#'   ))
#' )
#' @export
sn_find_de <- function(
  object,
  analysis = NULL,
  ident_1 = NULL,
  ident_2 = NULL,
  group_by = NULL,
  subset_by = NULL,
  subset_levels = NULL,
  sample_by = NULL,
  assay = "RNA",
  layer = NULL,
  features = NULL,
  method = NULL,
  only_pos = NULL,
  logfc_threshold = 0.1,
  min_pct = 0.25,
  p_val_cutoff = 0.05,
  de_logfc = 0.25,
  min_cells_per_sample = 10,
  result_id = "default",
  return_object = TRUE,
  verbose = TRUE,
  modality = c("auto", "single_cell", "bulk"),
  metadata = NULL,
  design = ~condition,
  contrast = NULL,
  backend_control = list(),
  ...
) {
  result_id <- .sn_validate_result_id(result_id)
  assay_missing <- missing(assay)
  design_missing <- missing(design)
  modality <- match.arg(modality)
  if (identical(modality, "auto")) {
    modality <- if (inherits(object, "Seurat")) "single_cell" else "bulk"
  }

  if (identical(modality, "bulk")) {
    if (is_null(contrast)) {
      stop("`contrast` is required for bulk differential expression.", call. = FALSE)
    }
    bulk_method <- tolower(method %||% "auto")
    return(.sn_find_bulk_de(
      object = object,
      metadata = metadata,
      design = design,
      contrast = contrast,
      method = bulk_method,
      assay = if (assay_missing) NULL else assay,
      result_id = if (identical(result_id, "default")) "bulk_de" else result_id,
      backend_control = backend_control
    ))
  }

  if (!inherits(object, "Seurat")) {
    stop("`modality = 'single_cell'` requires a Seurat object.", call. = FALSE)
  }

  if (is_null(analysis)) {
    analysis <- if (!is_null(sample_by)) {
      "pseudobulk"
    } else if (is_null(ident_1) && is_null(ident_2)) {
      "markers"
    } else {
      "contrast"
    }
  }
  analysis <- match.arg(analysis, c("markers", "contrast", "pseudobulk"))
  layer <- layer %||% if (analysis == "pseudobulk") "counts" else "data"
  method <- .sn_normalize_de_method(analysis = analysis, method = method)
  de_feature_universe <- .sn_resolve_de_feature_universe(
    object = object,
    assay = assay,
    layer = layer,
    features = features,
    method = method,
    verbose = verbose
  )
  features <- de_feature_universe$requested_features

  if (analysis %in% c("contrast", "pseudobulk") && (is_null(ident_1) || is_null(ident_2))) {
    stop("`ident_1` and `ident_2` are required for contrast and pseudobulk analyses.")
  }

  .sn_validate_de_subset_levels(
    metadata = object[[]],
    subset_by = subset_by,
    subset_levels = subset_levels
  )

  de_acceleration_patches <- if (
    analysis %in% c("markers", "contrast") && !identical(method, "COSGR")
  ) {
    "seurat"
  } else {
    character(0)
  }
  .sn_with_acceleration_provenance_context({
  result <- NULL
  pseudobulk_estimands <- NULL

  if (analysis == "pseudobulk") {
    if (is_null(sample_by)) {
      stop("`sample_by` is required for pseudobulk analyses.")
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
      sample_col = sample_by,
      subset_by = subset_by,
      subset_levels = subset_levels,
      assay = assay,
      layer = layer,
      features = features,
      method = method,
      min_cells_per_sample = min_cells_per_sample,
      design = if (design_missing) NULL else design,
      contrast = contrast,
      verbose = verbose
    )
    pseudobulk_estimands <- attr(result, "shennong_contrast_estimands", exact = TRUE)
  } else {
    prepared <- .sn_prepare_seurat_layer_alias(
      object = object,
      assay = assay,
      source_layer = layer,
      target_layer = .sn_guess_seurat_target_layer(layer)
    )
    analysis_object <- prepared$object
    restore_context <- prepared$context

    if (is_null(subset_by)) {
      result <- .sn_run_seurat_de(
        object = analysis_object,
        analysis = analysis,
        ident_1 = ident_1,
        ident_2 = ident_2,
        group_by = group_by,
        features = features,
        assay = assay,
        layer = layer,
        method = method,
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
          layer = layer,
          method = method,
          only_pos = only_pos,
          logfc_threshold = logfc_threshold,
          min_pct = min_pct,
          verbose = verbose,
          ...
        )

        current_result <- tibble::rownames_to_column(as.data.frame(current_result), var = "gene")
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

    object <- .sn_restore_seurat_layer_alias(object = analysis_object, context = restore_context)
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

  rank_col <- c("cosg_score", "avg_log2FC", "avg_logFC", "log2FoldChange", "logFC", "stat", "rank")
  rank_col <- rank_col[rank_col %in% colnames(result)][1] %||% NA_character_
  p_col <- c("p_val_adj", "padj", "FDR")
  p_col <- p_col[p_col %in% colnames(result)][1] %||% NA_character_

  stored_result <- list(
    schema_version = .sn_analysis_result_schema_version(),
    package_version = as.character(utils::packageVersion("Shennong")),
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    table = tibble::as_tibble(result),
    analysis = analysis,
    method = method,
    group_by = group_by,
    group_col = group_col,
    ident_1 = ident_1,
    ident_2 = ident_2,
    subset_by = subset_by,
    sample_col = sample_by,
    design = if (analysis == "pseudobulk") {
      if (design_missing) "automatic independent-or-paired design" else paste(deparse(design), collapse = "")
    } else {
      NULL
    },
    contrast = if (analysis == "pseudobulk") contrast %||% c(group_by, ident_1, ident_2) else NULL,
    parameters = if (analysis == "pseudobulk") {
      list(contrast_estimands = pseudobulk_estimands %||% list())
    } else {
      list()
    },
    input = list(
      assay = assay,
      layer = layer,
      requested_features = de_feature_universe$requested_features,
      tested_features = de_feature_universe$tested_features,
      tested_features_count = length(de_feature_universe$tested_features),
      tested_features_source = de_feature_universe$tested_features_source
    ),
    assay = assay,
    layer = layer,
    rank_col = rank_col,
    p_col = p_col,
    p_val_cutoff = p_val_cutoff,
    de_logfc = de_logfc,
    min_pct = min_pct,
    logfc_threshold = logfc_threshold,
    n_genes = nrow(result),
    provenance = .sn_contextual_analysis_provenance()
  )

  object <- sn_store_result(
    object = object,
    type = "de",
    result_id = result_id,
    result = stored_result
  )

  if (return_object) {
    .sn_log_seurat_command(object = object, assay = assay, name = "sn_find_de")
  } else {
    tibble::as_tibble(result)
  }
  }, patches = de_acceleration_patches)
}
