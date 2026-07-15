.sn_bulk_input <- function(object, metadata = NULL, assay = NULL) {
  source <- "matrix"
  if (inherits(object, "SummarizedExperiment")) {
    check_installed("SummarizedExperiment", reason = "to read a SummarizedExperiment bulk input.")
    assay_names <- SummarizedExperiment::assayNames(object)
    assay <- assay %||% assay_names[[1]]
    if (!assay %in% assay_names) stop("Bulk assay '", assay, "' was not found.", call. = FALSE)
    expression <- SummarizedExperiment::assay(object, assay)
    metadata <- metadata %||% as.data.frame(SummarizedExperiment::colData(object))
    source <- "SummarizedExperiment"
  } else if (is.list(object) && !is.data.frame(object)) {
    expression <- object$counts %||% object$expression %||% object$matrix
    metadata <- metadata %||% object$metadata %||% object$col_data %||% object$samples
    source <- "list"
  } else {
    expression <- object
  }
  if (inherits(expression, "Matrix")) expression <- as.matrix(expression)
  if (is.data.frame(expression)) expression <- as.matrix(expression)
  if (!is.matrix(expression) || !is.numeric(expression)) {
    stop("Bulk input must resolve to a numeric feature-by-sample matrix.", call. = FALSE)
  }
  if (is.null(rownames(expression)) || any(!nzchar(rownames(expression))) || anyDuplicated(rownames(expression))) {
    stop("Bulk expression requires unique, non-empty feature row names.", call. = FALSE)
  }
  if (is.null(colnames(expression)) || any(!nzchar(colnames(expression))) || anyDuplicated(colnames(expression))) {
    stop("Bulk expression requires unique, non-empty sample column names.", call. = FALSE)
  }
  if (any(!is.finite(expression), na.rm = TRUE)) {
    stop("Bulk expression contains non-finite values.", call. = FALSE)
  }
  storage.mode(expression) <- "double"
  if (is_null(metadata)) {
    metadata <- data.frame(sample = colnames(expression), row.names = colnames(expression), check.names = FALSE)
  } else {
    metadata <- as.data.frame(metadata, check.names = FALSE)
    if (is.null(rownames(metadata)) || any(!colnames(expression) %in% rownames(metadata))) {
      sample_col <- intersect(c("sample", "sample_id", "Sample", "SampleID"), colnames(metadata))
      if (length(sample_col) == 0L) {
        stop("Sample metadata row names (or a sample ID column) must match expression columns.", call. = FALSE)
      }
      rownames(metadata) <- as.character(metadata[[sample_col[[1]]]])
    }
    missing_samples <- setdiff(colnames(expression), rownames(metadata))
    if (length(missing_samples) > 0L) {
      stop("Sample metadata is missing: ", paste(missing_samples, collapse = ", "), ".", call. = FALSE)
    }
    metadata <- metadata[colnames(expression), , drop = FALSE]
  }
  metadata$.sn_sample <- rownames(metadata)
  is_counts <- all(expression >= 0) && all(abs(expression - round(expression)) < 1e-8)
  list(matrix = expression, metadata = metadata, source = source, assay = assay, is_counts = is_counts)
}

.sn_bulk_log_expression <- function(input, prior_count = 0.5) {
  expression <- input$matrix
  if (!isTRUE(input$is_counts)) return(expression)
  library_size <- colSums(expression)
  scale <- pmax(library_size, 1) / 1e6
  log2(sweep(expression, 2L, scale, "/") + prior_count)
}

.sn_bulk_result <- function(type, name, method, input, parameters, tables,
                            embeddings = list(), models = list(), diagnostics = list(),
                            warnings = character(), backend = method, seed = NA_integer_) {
  result <- list(
    schema_version = "1.0", analysis_type = type, name = name,
    method = method, backend = backend,
    input = input, parameters = parameters, tables = tables,
    embeddings = embeddings, graphs = list(), models = models,
    diagnostics = diagnostics, warnings = as.character(warnings),
    provenance = .sn_analysis_provenance(random_seed = seed)
  )
  sn_validate_result(result)
  result
}

#' Assess bulk transcriptomics sample quality
#'
#' Computes library size, detected features, expression distributions, sample
#' PCA, sample correlations, and robust multivariate outlier flags.
#'
#' @param object A feature-by-sample matrix, a `SummarizedExperiment`, or a list
#'   containing `counts`/`expression` and optional `metadata`.
#' @param metadata Optional sample metadata with rows matching sample names.
#' @param assay Assay name for `SummarizedExperiment` input.
#' @param top_variable Number of variable features used for PCA/correlation.
#' @param outlier_z Robust z-score threshold used for sample flags.
#' @param store_name Result name.
#' @return A validated Shennong bulk-QC result.
#' @export
sn_assess_bulk_qc <- function(object, metadata = NULL, assay = NULL,
                              top_variable = 2000L, outlier_z = 3.5,
                              store_name = "bulk_qc") {
  input <- .sn_bulk_input(object, metadata, assay)
  expression <- .sn_bulk_log_expression(input)
  variance <- apply(expression, 1L, stats::var)
  keep <- names(sort(variance, decreasing = TRUE))[seq_len(min(length(variance), as.integer(top_variable)))]
  pca <- stats::prcomp(t(expression[keep, , drop = FALSE]), center = TRUE, scale. = FALSE)
  pca_table <- tibble::as_tibble(pca$x[, seq_len(min(10L, ncol(pca$x))), drop = FALSE], rownames = "sample")
  correlation <- stats::cor(expression[keep, , drop = FALSE], use = "pairwise.complete.obs")
  sample_table <- tibble::tibble(
    sample = colnames(input$matrix),
    library_size = colSums(input$matrix),
    detected_features = colSums(input$matrix > 0),
    median_expression = apply(expression, 2L, stats::median),
    mean_correlation = (colSums(correlation) - 1) / pmax(1, ncol(correlation) - 1)
  )
  robust_z <- function(x) {
    spread <- stats::mad(x, constant = 1, na.rm = TRUE)
    if (!is.finite(spread) || spread == 0) return(rep(0, length(x)))
    (x - stats::median(x, na.rm = TRUE)) / spread
  }
  sample_table$library_z <- robust_z(log1p(sample_table$library_size))
  sample_table$detected_z <- robust_z(sample_table$detected_features)
  sample_table$correlation_z <- robust_z(sample_table$mean_correlation)
  sample_table$outlier <- abs(sample_table$library_z) > outlier_z |
    abs(sample_table$detected_z) > outlier_z | sample_table$correlation_z < -outlier_z
  distribution <- dplyr::bind_rows(lapply(seq_len(ncol(expression)), function(index) {
    values <- stats::quantile(expression[, index], probs = c(0, .25, .5, .75, 1), na.rm = TRUE)
    tibble::tibble(sample = colnames(expression)[index], quantile = names(values), expression = as.numeric(values))
  }))
  variance_explained <- tibble::tibble(
    component = paste0("PC", seq_along(pca$sdev)),
    variance_explained = pca$sdev^2 / sum(pca$sdev^2)
  )
  .sn_bulk_result(
    "bulk_qc", store_name, "robust_qc",
    list(source = input$source, assay = input$assay, samples = ncol(input$matrix), features = nrow(input$matrix), counts = input$is_counts),
    list(top_variable = length(keep), outlier_z = outlier_z),
    list(primary = sample_table, samples = sample_table, distribution = distribution,
         variance_explained = variance_explained, correlation = as.data.frame(correlation)),
    embeddings = list(pca = pca$x), models = list(pca = pca),
    diagnostics = list(outliers = sum(sample_table$outlier), variable_features = length(keep))
  )
}

.sn_bulk_has_random_effect <- function(design) grepl("\\|", paste(deparse(design), collapse = ""), fixed = FALSE)

.sn_bulk_design <- function(metadata, design) {
  if (!inherits(design, "formula")) stop("`design` must be a formula.", call. = FALSE)
  variables <- all.vars(design)
  missing <- setdiff(variables, colnames(metadata))
  if (length(missing) > 0L) stop("Design variable(s) missing from metadata: ", paste(missing, collapse = ", "), ".", call. = FALSE)
  if (.sn_bulk_has_random_effect(design)) return(list(formula = design, matrix = NULL, rank = NA_integer_))
  design_matrix <- stats::model.matrix(design, data = metadata)
  rank <- qr(design_matrix)$rank
  if (rank < ncol(design_matrix)) stop("The bulk design matrix is not full rank.", call. = FALSE)
  list(formula = design, matrix = design_matrix, rank = rank)
}

.sn_bulk_validate_contrast <- function(metadata, contrast) {
  if (!is.character(contrast) || length(contrast) != 3L) {
    stop("`contrast` must be c(variable, numerator, denominator).", call. = FALSE)
  }
  variable <- contrast[[1]]
  if (!variable %in% colnames(metadata)) stop("Contrast variable '", variable, "' was not found.", call. = FALSE)
  values <- as.character(metadata[[variable]])
  missing <- setdiff(contrast[2:3], unique(values))
  if (length(missing) > 0L) stop("Contrast level(s) missing: ", paste(missing, collapse = ", "), ".", call. = FALSE)
  contrast
}

.sn_bulk_contrast_vector <- function(metadata, design, contrast) {
  metadata <- as.data.frame(metadata)
  for (column in colnames(metadata)) if (is.character(metadata[[column]])) metadata[[column]] <- factor(metadata[[column]])
  reference <- metadata[1, , drop = FALSE]
  for (column in colnames(metadata)) {
    if (is.numeric(metadata[[column]])) reference[[column]] <- stats::median(metadata[[column]], na.rm = TRUE)
    if (is.factor(metadata[[column]])) reference[[column]] <- factor(levels(metadata[[column]])[[1]], levels = levels(metadata[[column]]))
  }
  numerator <- denominator <- reference
  variable <- contrast[[1]]
  if (is.factor(metadata[[variable]])) {
    numerator[[variable]] <- factor(contrast[[2]], levels = levels(metadata[[variable]]))
    denominator[[variable]] <- factor(contrast[[3]], levels = levels(metadata[[variable]]))
  } else {
    numerator[[variable]] <- contrast[[2]]
    denominator[[variable]] <- contrast[[3]]
  }
  drop(stats::model.matrix(design, numerator) - stats::model.matrix(design, denominator))
}

.sn_bulk_standardize_de <- function(table, method) {
  table <- as.data.frame(table, check.names = FALSE)
  if (!"gene" %in% colnames(table)) table <- tibble::rownames_to_column(table, "gene")
  choose <- function(candidates, default = NA_real_) {
    column <- intersect(candidates, colnames(table))
    if (length(column) == 0L) rep(default, nrow(table)) else as.numeric(table[[column[[1]]]])
  }
  tibble::tibble(
    gene = as.character(table$gene),
    log2_fold_change = choose(c("log2FoldChange", "logFC", "coef", "estimate")),
    statistic = choose(c("stat", "F", "t", "LR", "z")),
    p_value = choose(c("pvalue", "PValue", "P.Value", "p_value")),
    adjusted_p_value = choose(c("padj", "FDR", "adj.P.Val", "adjusted_p_value")),
    base_mean = choose(c("baseMean", "logCPM", "AveExpr", "base_mean")),
    method = method
  ) |>
    dplyr::arrange(.data$adjusted_p_value, .data$p_value)
}

.sn_edger_norm_lib_sizes <- function(object) {
  if (exists("normLibSizes", envir = asNamespace("edgeR"), inherits = FALSE)) {
    return(edgeR::normLibSizes(object))
  }
  edgeR::calcNormFactors(object)
}

.sn_bulk_de_edger <- function(input, design_info, contrast, control) {
  check_installed("edgeR", reason = "to run edgeR bulk differential expression.")
  y <- edgeR::DGEList(counts = round(input$matrix))
  keep <- edgeR::filterByExpr(y, design = design_info$matrix)
  y <- .sn_edger_norm_lib_sizes(y[keep, , keep.lib.sizes = FALSE])
  y <- edgeR::estimateDisp(y, design_info$matrix, robust = control$robust %||% TRUE)
  fit <- edgeR::glmQLFit(y, design_info$matrix, robust = control$robust %||% TRUE)
  vector <- .sn_bulk_contrast_vector(input$metadata, design_info$formula, contrast)
  test <- edgeR::glmQLFTest(fit, contrast = vector)
  list(table = edgeR::topTags(test, n = Inf, sort.by = "none")$table, model = fit,
       diagnostics = list(retained_features = sum(keep), filtered_features = sum(!keep), contrast_vector = vector))
}

.sn_bulk_de_deseq2 <- function(input, design_info, contrast, control) {
  check_installed("DESeq2", reason = "to run DESeq2 bulk differential expression.")
  dataset <- DESeq2::DESeqDataSetFromMatrix(round(input$matrix), input$metadata, design_info$formula)
  keep <- rowSums(DESeq2::counts(dataset) >= (control$min_count %||% 10)) >= (control$min_samples %||% 2L)
  dataset <- dataset[keep, ]
  dataset <- DESeq2::DESeq(dataset, quiet = TRUE)
  result <- DESeq2::results(dataset, contrast = contrast, independentFiltering = control$independent_filtering %||% TRUE)
  shrink_applied <- FALSE
  if (isTRUE(control$shrink %||% TRUE)) {
    coefficient <- grep(paste0("^", contrast[[1]], "_", contrast[[2]], "_vs_", contrast[[3]], "$"), DESeq2::resultsNames(dataset), value = TRUE)
    shrunken <- tryCatch({
      if (length(coefficient) == 1L && requireNamespace("apeglm", quietly = TRUE)) {
        DESeq2::lfcShrink(dataset, coef = coefficient, res = result, type = "apeglm")
      } else {
        DESeq2::lfcShrink(dataset, contrast = contrast, res = result, type = "normal")
      }
    }, error = function(error) NULL)
    if (!is_null(shrunken)) {
      result <- shrunken
      shrink_applied <- TRUE
    }
  }
  list(table = as.data.frame(result), model = dataset,
       diagnostics = list(retained_features = sum(keep), filtered_features = sum(!keep),
                          size_factors = DESeq2::sizeFactors(dataset), shrink_applied = shrink_applied))
}

.sn_bulk_de_limma <- function(input, design_info, contrast, control) {
  check_installed("limma", reason = "to run limma bulk differential expression.")
  if (isTRUE(input$is_counts)) {
    check_installed("edgeR", reason = "to run limma-voom on bulk counts.")
    y <- .sn_edger_norm_lib_sizes(edgeR::DGEList(counts = round(input$matrix)))
    transformed <- limma::voom(y, design_info$matrix, plot = FALSE)
    fit <- limma::lmFit(transformed, design_info$matrix)
  } else {
    transformed <- input$matrix
    fit <- limma::lmFit(transformed, design_info$matrix)
  }
  vector <- .sn_bulk_contrast_vector(input$metadata, design_info$formula, contrast)
  fit <- limma::contrasts.fit(fit, contrasts = matrix(vector, ncol = 1L))
  fit <- limma::eBayes(fit, robust = control$robust %||% FALSE)
  list(table = limma::topTable(fit, number = Inf, sort.by = "none"), model = fit,
       diagnostics = list(retained_features = nrow(input$matrix), filtered_features = 0L, contrast_vector = vector), transformed = transformed)
}

.sn_bulk_de_dream <- function(input, design_info, contrast, control) {
  check_installed("variancePartition", reason = "to run dream repeated-measures differential expression.")
  check_installed("edgeR", reason = "to create dream precision weights.")
  expression <- if (isTRUE(input$is_counts)) {
    variancePartition::voomWithDreamWeights(edgeR::DGEList(counts = round(input$matrix)), design_info$formula, input$metadata, plot = FALSE)
  } else input$matrix
  fit <- variancePartition::dream(expression, design_info$formula, input$metadata)
  fit <- limma::eBayes(fit)
  coefficient <- grep(paste0("^", contrast[[1]], contrast[[2]], "$"), colnames(fit$coefficients), value = TRUE)
  if (length(coefficient) != 1L) {
    stop("dream could not identify a single coefficient for the requested contrast; relevel the denominator explicitly.", call. = FALSE)
  }
  list(table = limma::topTable(fit, coef = coefficient, number = Inf, sort.by = "none"), model = fit,
       diagnostics = list(retained_features = nrow(input$matrix), coefficient = coefficient), transformed = expression)
}

#' Find differential expression in bulk transcriptomics data
#'
#' @param object Bulk input accepted by `sn_assess_bulk_qc()`.
#' @param metadata Optional sample metadata.
#' @param design A fixed- or mixed-effects model formula.
#' @param contrast Character triple: variable, numerator, denominator.
#' @param method One of `auto`, `edger`, `deseq2`, `limma`, or `dream`.
#' @param assay Assay name for `SummarizedExperiment` input.
#' @param store_name Result name.
#' @param backend_control Backend options or a custom `runner`/precomputed `result`.
#' @return A validated Shennong bulk-DE result.
#' @export
sn_find_bulk_de <- function(object, metadata = NULL, design = ~condition,
                            contrast, method = c("auto", "edger", "deseq2", "limma", "dream"),
                            assay = NULL, store_name = "bulk_de", backend_control = list()) {
  .sn_find_bulk_de(
    object = object,
    metadata = metadata,
    design = design,
    contrast = contrast,
    method = method,
    assay = assay,
    store_name = store_name,
    backend_control = backend_control
  )
}

.sn_find_bulk_de <- function(object, metadata = NULL, design = ~condition,
                             contrast, method = c("auto", "edger", "deseq2", "limma", "dream"),
                             assay = NULL, store_name = "bulk_de", backend_control = list()) {
  method <- match.arg(method)
  input <- .sn_bulk_input(object, metadata, assay)
  contrast <- .sn_bulk_validate_contrast(input$metadata, contrast)
  design_info <- .sn_bulk_design(input$metadata, design)
  selected <- method
  if (identical(method, "auto")) {
    selected <- if (.sn_bulk_has_random_effect(design)) "dream" else if (input$is_counts) "edger" else "limma"
  }
  if (!input$is_counts && selected %in% c("edger", "deseq2")) stop(selected, " requires non-negative integer counts.", call. = FALSE)
  output <- if (is.function(backend_control$runner)) {
    backend_control$runner(input = input, design = design, contrast = contrast, method = selected, backend_control = backend_control)
  } else if (!is_null(backend_control$result)) {
    backend_control$result
  } else {
    switch(selected,
      edger = .sn_bulk_de_edger(input, design_info, contrast, backend_control),
      deseq2 = .sn_bulk_de_deseq2(input, design_info, contrast, backend_control),
      limma = .sn_bulk_de_limma(input, design_info, contrast, backend_control),
      dream = .sn_bulk_de_dream(input, design_info, contrast, backend_control)
    )
  }
  table <- .sn_bulk_standardize_de(output$table %||% output, selected)
  .sn_bulk_result(
    "bulk_de", store_name, selected,
    list(source = input$source, assay = input$assay, samples = ncol(input$matrix), features = nrow(input$matrix), counts = input$is_counts),
    list(design = paste(deparse(design), collapse = ""), contrast = contrast, requested_method = method,
         independent_filtering = backend_control$independent_filtering %||% TRUE, shrink = backend_control$shrink %||% TRUE),
    list(primary = table, differential_expression = table),
    models = list(fit = output$model %||% NULL),
    diagnostics = c(output$diagnostics %||% list(), list(design_rank = design_info$rank, tested_features = nrow(table)))
  )
}

#' Score pathways in bulk expression samples
#'
#' @param object Bulk input accepted by `sn_assess_bulk_qc()`.
#' @param signatures Named list or two-column data frame of gene sets.
#' @param method Scoring method: mean, GSVA, or ssGSEA.
#' @param metadata Optional sample metadata.
#' @param assay Assay name for `SummarizedExperiment` input.
#' @param min_genes Minimum matched genes per pathway.
#' @param store_name Result name.
#' @param backend_control Backend-specific controls.
#' @return A validated bulk pathway-score result.
#' @export
sn_score_bulk_pathways <- function(object, signatures, method = c("mean", "gsva", "ssgsea"),
                                   metadata = NULL, assay = NULL, min_genes = 2L,
                                   store_name = "bulk_pathways", backend_control = list()) {
  method <- match.arg(method)
  input <- .sn_bulk_input(object, metadata, assay)
  expression <- .sn_bulk_log_expression(input)
  signatures <- .sn_normalize_program_signatures(signatures)
  matched <- .sn_match_program_features(signatures, rownames(expression), min_genes)
  scores <- if (identical(method, "mean")) {
    .sn_score_programs_mean(expression, matched$signatures)
  } else {
    .sn_score_programs_gsva(expression, matched$signatures, method, backend_control)
  }
  score_table <- .sn_program_score_table(scores, level = "sample")
  names(score_table)[match(c("entity", "program"), names(score_table))] <- c("sample", "pathway")
  .sn_bulk_result(
    "bulk_pathway", store_name, method,
    list(source = input$source, samples = ncol(input$matrix), features = nrow(input$matrix)),
    list(min_genes = min_genes),
    list(primary = score_table, scores = score_table, coverage = matched$coverage),
    embeddings = list(scores = t(scores)), diagnostics = list(pathways = nrow(scores), samples = ncol(scores))
  )
}

#' Run weighted gene co-expression network analysis
#'
#' @param object Bulk input accepted by `sn_assess_bulk_qc()`.
#' @param metadata Optional sample metadata.
#' @param traits Optional metadata columns tested against module eigengenes.
#' @param power Soft-thresholding power. If `NULL`, it is selected from `powers`.
#' @param powers Candidate powers for automatic selection.
#' @param min_module_size Minimum module size.
#' @param merge_cut_height Module merge threshold.
#' @param assay Assay name for `SummarizedExperiment` input.
#' @param store_name Result name.
#' @param backend_control Additional `blockwiseModules` arguments or custom output.
#' @return A validated WGCNA result.
#' @export
sn_run_wgcna <- function(object, metadata = NULL, traits = NULL, power = NULL,
                         powers = c(1:10, 12, 14, 16, 18, 20), min_module_size = 30L,
                         merge_cut_height = 0.25, assay = NULL,
                         store_name = "wgcna", backend_control = list()) {
  check_installed("WGCNA", reason = "to run weighted co-expression network analysis.")
  wgcna_attached <- "package:WGCNA" %in% search()
  if (!wgcna_attached) {
    suppressPackageStartupMessages(base::attachNamespace("WGCNA"))
    on.exit(detach("package:WGCNA"), add = TRUE)
  }
  input <- .sn_bulk_input(object, metadata, assay)
  expression <- .sn_bulk_log_expression(input)
  variable <- apply(expression, 1L, stats::var)
  expression <- expression[is.finite(variable) & variable > 0, , drop = FALSE]
  dat_expr <- t(expression)
  quality <- WGCNA::goodSamplesGenes(dat_expr, verbose = 0)
  dat_expr <- dat_expr[quality$goodSamples, quality$goodGenes, drop = FALSE]
  selected_power <- power
  power_table <- tibble::tibble()
  if (is_null(selected_power)) {
    candidate <- WGCNA::pickSoftThreshold(dat_expr, powerVector = powers, verbose = 0)
    power_table <- tibble::as_tibble(candidate$fitIndices)
    acceptable <- which(candidate$fitIndices[, "SFT.R.sq"] >= (backend_control$scale_free_r2 %||% 0.8))
    selected_power <- if (length(acceptable) > 0L) candidate$fitIndices[acceptable[[1]], "Power"] else powers[[which.max(candidate$fitIndices[, "SFT.R.sq"])]]
  }
  output <- if (is.function(backend_control$runner)) {
    backend_control$runner(dat_expr = dat_expr, power = selected_power, metadata = input$metadata, backend_control = backend_control)
  } else if (!is_null(backend_control$result)) {
    backend_control$result
  } else {
    defaults <- list(datExpr = dat_expr, power = selected_power, TOMType = "unsigned",
                     minModuleSize = as.integer(min_module_size), mergeCutHeight = merge_cut_height,
                     numericLabels = FALSE, pamRespectsDendro = FALSE, verbose = 0)
    fit <- do.call(WGCNA::blockwiseModules, utils::modifyList(defaults, backend_control$blockwise %||% list()))
    list(model = fit, colors = fit$colors, eigengenes = WGCNA::orderMEs(fit$MEs))
  }
  colors <- as.character(output$colors %||% output$model$colors)
  genes <- colnames(dat_expr)
  if (length(colors) != length(genes)) stop("WGCNA output colors must match retained genes.", call. = FALSE)
  eigengenes <- as.data.frame(output$eigengenes %||% output$model$MEs, check.names = FALSE)
  eigengenes$sample <- rownames(eigengenes)
  modules <- tibble::tibble(gene = genes, module = colors)
  traits <- traits %||% setdiff(colnames(input$metadata), ".sn_sample")
  associations <- tibble::tibble()
  if (length(traits) > 0L && ncol(eigengenes) > 1L) {
    trait_data <- input$metadata[rownames(eigengenes), traits, drop = FALSE]
    encoded <- stats::model.matrix(~ . - 1, data = trait_data)
    module_matrix <- as.matrix(eigengenes[, setdiff(colnames(eigengenes), "sample"), drop = FALSE])
    correlations <- stats::cor(module_matrix, encoded, use = "pairwise.complete.obs")
    p_values <- WGCNA::corPvalueStudent(correlations, nrow(module_matrix))
    associations <- dplyr::bind_rows(lapply(seq_len(nrow(correlations)), function(i) {
      tibble::tibble(module = rownames(correlations)[i], trait = colnames(correlations),
                     correlation = as.numeric(correlations[i, ]), p_value = as.numeric(p_values[i, ]))
    }))
    associations$adjusted_p_value <- stats::p.adjust(associations$p_value, method = "BH")
  }
  .sn_bulk_result(
    "bulk_network", store_name, "wgcna",
    list(source = input$source, samples = nrow(dat_expr), features = ncol(dat_expr)),
    list(power = selected_power, min_module_size = min_module_size, merge_cut_height = merge_cut_height, traits = traits),
    list(primary = modules, modules = modules, eigengenes = tibble::as_tibble(eigengenes),
         trait_associations = associations, power_selection = power_table),
    models = list(fit = output$model %||% NULL),
    diagnostics = list(modules = length(setdiff(unique(colors), "grey")), excluded_samples = sum(!quality$goodSamples), excluded_genes = sum(!quality$goodGenes))
  )
}

.sn_bulk_feature_data <- function(input, features) {
  metadata <- input$metadata
  expression <- .sn_bulk_log_expression(input)
  missing <- setdiff(features, c(colnames(metadata), rownames(expression)))
  if (length(missing) > 0L) stop("Feature(s) not found in metadata or expression: ", paste(missing, collapse = ", "), ".", call. = FALSE)
  values <- lapply(features, function(feature) {
    if (feature %in% colnames(metadata)) metadata[[feature]] else as.numeric(expression[feature, rownames(metadata)])
  })
  names(values) <- features
  as.data.frame(values, check.names = FALSE, row.names = rownames(metadata))
}

#' Run Cox proportional-hazards models for bulk features
#'
#' @param object Bulk input accepted by `sn_assess_bulk_qc()`.
#' @param time,event Metadata columns containing follow-up time and event status.
#' @param features Expression features or numeric metadata columns.
#' @param covariates Optional adjustment variables.
#' @param metadata Optional sample metadata.
#' @param assay Assay name for `SummarizedExperiment` input.
#' @param store_name Result name.
#' @return A validated survival result with one adjusted Cox model per feature.
#' @export
sn_run_survival <- function(object, time, event, features, covariates = NULL,
                            metadata = NULL, assay = NULL, store_name = "bulk_survival") {
  check_installed("survival", reason = "to fit Cox proportional-hazards models.")
  input <- .sn_bulk_input(object, metadata, assay)
  required <- c(time, event, covariates)
  missing <- setdiff(required, colnames(input$metadata))
  if (length(missing) > 0L) stop("Survival metadata column(s) missing: ", paste(missing, collapse = ", "), ".", call. = FALSE)
  feature_data <- .sn_bulk_feature_data(input, features)
  base <- cbind(input$metadata[, required, drop = FALSE], feature_data)
  rows <- vector("list", length(features))
  models <- stats::setNames(vector("list", length(features)), features)
  for (i in seq_along(features)) {
    feature <- features[[i]]
    model_data <- base[, c(time, event, feature, covariates), drop = FALSE]
    colnames(model_data)[1:3] <- c(".time", ".event", ".feature")
    formula <- stats::as.formula(paste0("survival::Surv(.time, .event) ~ .feature", if (length(covariates)) paste0(" + ", paste(covariates, collapse = " + ")) else ""))
    fit <- survival::coxph(formula, data = model_data, x = TRUE, y = TRUE)
    estimate <- summary(fit)$coefficients[".feature", ]
    interval <- summary(fit)$conf.int[".feature", ]
    ph <- tryCatch(survival::cox.zph(fit)$table[".feature", "p"], error = function(e) NA_real_)
    rows[[i]] <- tibble::tibble(feature = feature, hazard_ratio = unname(interval[["exp(coef)"]]),
      conf_low = unname(interval[["lower .95"]]), conf_high = unname(interval[["upper .95"]]),
      coefficient = unname(estimate[["coef"]]), statistic = unname(estimate[["z"]]),
      p_value = unname(estimate[["Pr(>|z|)"]]), ph_p_value = ph, n = fit$n, events = fit$nevent)
    models[[feature]] <- fit
  }
  table <- dplyr::bind_rows(rows)
  table$adjusted_p_value <- stats::p.adjust(table$p_value, method = "BH")
  .sn_bulk_result(
    "bulk_survival", store_name, "cox",
    list(source = input$source, samples = ncol(input$matrix), features = features),
    list(time = time, event = event, covariates = covariates),
    list(primary = table, survival = table), models = models,
    diagnostics = list(models = length(models), proportional_hazards_warnings = sum(table$ph_p_value < 0.05, na.rm = TRUE))
  )
}

#' Associate bulk features with clinical variables
#'
#' Numeric clinical variables use linear models; two-level categorical variables
#' use Welch tests and multi-level variables use ANOVA.
#'
#' @param object Bulk input accepted by `sn_assess_bulk_qc()`.
#' @param features Expression features or numeric metadata columns.
#' @param clinical_vars Clinical metadata variables to test.
#' @param covariates Optional adjustment variables for linear models.
#' @param metadata Optional sample metadata.
#' @param assay Assay name for `SummarizedExperiment` input.
#' @param store_name Result name.
#' @return A validated clinical-association result.
#' @export
sn_run_clinical_association <- function(object, features, clinical_vars, covariates = NULL,
                                        metadata = NULL, assay = NULL,
                                        store_name = "bulk_clinical") {
  input <- .sn_bulk_input(object, metadata, assay)
  missing <- setdiff(c(clinical_vars, covariates), colnames(input$metadata))
  if (length(missing) > 0L) stop("Clinical metadata column(s) missing: ", paste(missing, collapse = ", "), ".", call. = FALSE)
  feature_data <- .sn_bulk_feature_data(input, features)
  data <- cbind(input$metadata, feature_data)
  rows <- list()
  index <- 0L
  for (feature in features) for (clinical in clinical_vars) {
    index <- index + 1L
    subset <- data[, c(feature, clinical, covariates), drop = FALSE]
    colnames(subset)[1:2] <- c(".feature", ".clinical")
    if (is.numeric(subset$.clinical)) {
      formula <- stats::as.formula(paste0(".feature ~ .clinical", if (length(covariates)) paste0(" + ", paste(covariates, collapse = " + ")) else ""))
      fit <- stats::lm(formula, data = subset)
      coefficient <- summary(fit)$coefficients[".clinical", ]
      rows[[index]] <- tibble::tibble(feature = feature, clinical_variable = clinical, test = "linear_model",
        estimate = unname(coefficient[["Estimate"]]), statistic = unname(coefficient[["t value"]]),
        p_value = unname(coefficient[["Pr(>|t|)"]]), n = stats::nobs(fit))
    } else if (nlevels(factor(subset$.clinical)) == 2L && length(covariates) == 0L) {
      subset$.clinical <- factor(subset$.clinical)
      fit <- stats::t.test(.feature ~ .clinical, data = subset)
      estimate <- unname(diff(fit$estimate))
      rows[[index]] <- tibble::tibble(feature = feature, clinical_variable = clinical,
        test = "welch_t", estimate = estimate, statistic = unname(fit$statistic),
        p_value = fit$p.value, n = sum(stats::complete.cases(subset)))
    } else {
      subset$.clinical <- factor(subset$.clinical)
      formula <- stats::as.formula(paste0(".feature ~ .clinical", if (length(covariates)) paste0(" + ", paste(covariates, collapse = " + ")) else ""))
      fit <- stats::lm(formula, data = subset)
      anova_table <- stats::anova(fit)[1, , drop = FALSE]
      coefficient <- if (nlevels(subset$.clinical) == 2L) stats::coef(fit)[[2]] else NA_real_
      rows[[index]] <- tibble::tibble(feature = feature, clinical_variable = clinical,
        test = "adjusted_anova", estimate = unname(coefficient), statistic = unname(anova_table[["F value"]][[1]]),
        p_value = unname(anova_table[["Pr(>F)"]][[1]]), n = stats::nobs(fit))
    }
  }
  table <- dplyr::bind_rows(rows)
  table$adjusted_p_value <- stats::p.adjust(table$p_value, method = "BH")
  .sn_bulk_result(
    "bulk_clinical", store_name, "model_association",
    list(source = input$source, samples = ncol(input$matrix), features = features),
    list(clinical_variables = clinical_vars, covariates = covariates),
    list(primary = table, associations = table),
    diagnostics = list(tests = nrow(table), significant = sum(table$adjusted_p_value < 0.05, na.rm = TRUE))
  )
}

#' Run a bulk transcriptomics workflow
#'
#' @param object Bulk input accepted by `sn_assess_bulk_qc()`.
#' @param workflow Workflow to run.
#' @param ... Arguments passed to the selected workflow.
#' @return A validated Shennong analysis result.
#' @export
sn_run_bulk <- function(object, workflow = c("qc", "de", "pathway", "network", "survival"), ...) {
  workflow <- match.arg(workflow)
  switch(workflow,
    qc = sn_assess_bulk_qc(object, ...),
    de = sn_find_bulk_de(object, ...),
    pathway = sn_score_bulk_pathways(object, ...),
    network = sn_run_wgcna(object, ...),
    survival = sn_run_survival(object, ...)
  )
}
