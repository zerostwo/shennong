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
                            warnings = character(), backend = method, seed = NA_integer_,
                            schema_version = "1.0.0") {
  result <- list(
    schema_version = schema_version, analysis_type = type, name = name,
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

.sn_survival_event <- function(event) {
  if (is.logical(event)) {
    event <- as.integer(event)
  } else if (is.factor(event)) {
    event <- as.character(event)
  }
  if (is.character(event)) {
    observed <- unique(event[!is.na(event)])
    if (!all(observed %in% c("0", "1"))) {
      stop("Survival event status must contain only 0 (censored) and 1 (event).", call. = FALSE)
    }
    event <- as.numeric(event)
  }
  if (!is.numeric(event)) {
    stop("Survival event status must be logical, numeric 0/1, or character/factor '0'/'1'.", call. = FALSE)
  }
  invalid <- !is.na(event) & (!is.finite(event) | !event %in% c(0, 1))
  if (any(invalid)) {
    stop("Survival event status must contain only 0 (censored) and 1 (event).", call. = FALSE)
  }
  as.integer(event)
}

.sn_survival_capture <- function(expr) {
  warnings <- character()
  error <- NULL
  value <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(warning) {
        warnings <<- c(warnings, conditionMessage(warning))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(condition) {
      error <<- conditionMessage(condition)
      NULL
    }
  )
  list(value = value, warnings = unique(warnings), error = error)
}

.sn_survival_association_row <- function(feature, n, events, status = "error", error = NA_character_) {
  tibble::tibble(
    feature = as.character(feature),
    hazard_ratio = NA_real_,
    conf_low = NA_real_,
    conf_high = NA_real_,
    coefficient = NA_real_,
    standard_error = NA_real_,
    statistic = NA_real_,
    p_value = NA_real_,
    ph_p_value = NA_real_,
    ph_global_p_value = NA_real_,
    concordance = NA_real_,
    concordance_se = NA_real_,
    aic = NA_real_,
    n = as.integer(n),
    events = as.integer(events),
    n_dropped = NA_integer_,
    status = as.character(status),
    error = as.character(error)
  )
}

.sn_survival_performance_row <- function(feature, n, events, status = "error", error = NA_character_) {
  tibble::tibble(
    feature = as.character(feature),
    n = as.integer(n),
    events = as.integer(events),
    concordance = NA_real_,
    concordance_se = NA_real_,
    aic = NA_real_,
    log_likelihood = NA_real_,
    likelihood_ratio_statistic = NA_real_,
    likelihood_ratio_df = NA_real_,
    likelihood_ratio_p_value = NA_real_,
    wald_statistic = NA_real_,
    wald_df = NA_real_,
    wald_p_value = NA_real_,
    score_statistic = NA_real_,
    score_df = NA_real_,
    score_p_value = NA_real_,
    status = as.character(status),
    error = as.character(error)
  )
}

.sn_survival_ph_error <- function(feature, error) {
  tibble::tibble(
    feature = as.character(feature),
    term = c(as.character(feature), "GLOBAL"),
    scope = c("feature", "global"),
    rho = NA_real_,
    chi_square = NA_real_,
    p_value = NA_real_,
    status = "error",
    error = as.character(error)
  )
}

.sn_survival_empty_ph_residuals <- function() {
  tibble::tibble(
    feature = character(),
    term = character(),
    scope = character(),
    time = numeric(),
    transformed_time = numeric(),
    scaled_schoenfeld = numeric(),
    transform = character()
  )
}

.sn_survival_ph_residuals <- function(ph_fit, feature, covariates, covariate_aliases) {
  if (is_null(ph_fit) || is_null(ph_fit$y) || length(ph_fit$y) == 0L) {
    return(.sn_survival_empty_ph_residuals())
  }
  residuals <- as.matrix(ph_fit$y)
  raw_terms <- colnames(residuals)
  if (is_null(raw_terms)) raw_terms <- paste0("term_", seq_len(ncol(residuals)))
  display_terms <- raw_terms
  display_terms[raw_terms == ".feature"] <- feature
  for (i in seq_along(covariate_aliases)) {
    display_terms[raw_terms == covariate_aliases[[i]]] <- covariates[[i]]
  }
  event_time <- as.numeric(ph_fit$time %||% rep(NA_real_, nrow(residuals)))
  transformed_time <- as.numeric(ph_fit$x %||% rep(NA_real_, nrow(residuals)))
  dplyr::bind_rows(lapply(seq_len(ncol(residuals)), function(index) {
    tibble::tibble(
      feature = feature,
      term = display_terms[[index]],
      scope = if (identical(raw_terms[[index]], ".feature")) "feature" else "covariate",
      time = event_time,
      transformed_time = transformed_time,
      scaled_schoenfeld = as.numeric(residuals[, index]),
      transform = as.character(ph_fit$transform %||% NA_character_)
    )
  }))
}

.sn_survival_test_value <- function(test, name, position) {
  if (is_null(test) || length(test) < position) return(NA_real_)
  if (!is_null(names(test)) && name %in% names(test)) return(as.numeric(test[[name]]))
  as.numeric(test[[position]])
}

.sn_survival_fit_cox <- function(model_data, feature, covariates, ties) {
  complete <- stats::complete.cases(model_data)
  if (is.numeric(model_data$.feature)) complete <- complete & is.finite(model_data$.feature)
  fit_data <- model_data[complete, , drop = FALSE]
  n <- nrow(fit_data)
  events <- if (n > 0L) sum(fit_data$.event == 1L) else 0L
  association <- .sn_survival_association_row(feature, n, events)
  association$n_dropped <- nrow(model_data) - n
  performance <- .sn_survival_performance_row(feature, n, events)

  fail <- function(message, warnings = character()) {
    association$error <- message
    performance$error <- message
    list(
      association = association,
      performance = performance,
      proportional_hazards = .sn_survival_ph_error(feature, message),
      proportional_hazards_residuals = .sn_survival_empty_ph_residuals(),
      fit = NULL,
      ph_fit = NULL,
      warnings = warnings
    )
  }
  if (!is.numeric(fit_data$.feature)) return(fail("Survival features must be numeric."))
  if (n < 3L) return(fail("At least three complete samples are required for a Cox model."))
  if (events < 1L) return(fail("At least one event is required for a Cox model."))
  if (length(unique(fit_data$.feature)) < 2L || !is.finite(stats::sd(fit_data$.feature)) || stats::sd(fit_data$.feature) == 0) {
    return(fail("The feature is constant among complete survival samples."))
  }

  covariate_aliases <- if (length(covariates)) paste0(".covariate_", seq_along(covariates)) else character()
  formula <- stats::reformulate(
    c(".feature", covariate_aliases),
    response = "survival::Surv(.time, .event)"
  )
  captured <- .sn_survival_capture(
    survival::coxph(
      formula,
      data = fit_data,
      ties = ties,
      x = TRUE,
      y = TRUE,
      singular.ok = FALSE
    )
  )
  if (is_null(captured$value)) return(fail(captured$error, captured$warnings))
  fit <- captured$value
  summary_fit <- summary(fit)
  coefficient_row <- match(".feature", rownames(summary_fit$coefficients))
  interval_row <- match(".feature", rownames(summary_fit$conf.int))
  if (is.na(coefficient_row) || is.na(interval_row)) {
    return(fail("The fitted Cox model did not return a feature coefficient.", captured$warnings))
  }
  estimate <- summary_fit$coefficients[coefficient_row, , drop = FALSE]
  interval <- summary_fit$conf.int[interval_row, , drop = FALSE]
  finite_model <- all(is.finite(c(
    estimate[1, "coef"], estimate[1, "se(coef)"],
    interval[1, "exp(coef)"], interval[1, "lower .95"], interval[1, "upper .95"]
  )))
  if (!finite_model) return(fail("The Cox model returned a non-finite feature estimate.", captured$warnings))

  concordance <- summary_fit$concordance %||% numeric()
  concordance_value <- .sn_survival_test_value(concordance, "Concordance", 1L)
  concordance_se <- .sn_survival_test_value(concordance, "se(Concordance)", 2L)
  association$hazard_ratio <- unname(interval[1, "exp(coef)"])
  association$conf_low <- unname(interval[1, "lower .95"])
  association$conf_high <- unname(interval[1, "upper .95"])
  association$coefficient <- unname(estimate[1, "coef"])
  association$standard_error <- unname(estimate[1, "se(coef)"])
  association$statistic <- unname(estimate[1, "z"])
  association$p_value <- unname(estimate[1, "Pr(>|z|)"])
  association$concordance <- concordance_value
  association$concordance_se <- concordance_se
  association$aic <- tryCatch(stats::AIC(fit), error = function(e) NA_real_)
  association$status <- "ok"
  association$error <- NA_character_

  performance$concordance <- concordance_value
  performance$concordance_se <- concordance_se
  performance$aic <- association$aic
  performance$log_likelihood <- unname(utils::tail(fit$loglik, 1L))
  performance$likelihood_ratio_statistic <- .sn_survival_test_value(summary_fit$logtest, "test", 1L)
  performance$likelihood_ratio_df <- .sn_survival_test_value(summary_fit$logtest, "df", 2L)
  performance$likelihood_ratio_p_value <- .sn_survival_test_value(summary_fit$logtest, "pvalue", 3L)
  performance$wald_statistic <- .sn_survival_test_value(summary_fit$waldtest, "test", 1L)
  performance$wald_df <- .sn_survival_test_value(summary_fit$waldtest, "df", 2L)
  performance$wald_p_value <- .sn_survival_test_value(summary_fit$waldtest, "pvalue", 3L)
  performance$score_statistic <- .sn_survival_test_value(summary_fit$sctest, "test", 1L)
  performance$score_df <- .sn_survival_test_value(summary_fit$sctest, "df", 2L)
  performance$score_p_value <- .sn_survival_test_value(summary_fit$sctest, "pvalue", 3L)
  performance$status <- "ok"
  performance$error <- NA_character_

  captured_ph <- .sn_survival_capture(survival::cox.zph(fit, terms = TRUE))
  if (is_null(captured_ph$value)) {
    ph_table <- .sn_survival_ph_error(feature, captured_ph$error)
    ph_residuals <- .sn_survival_empty_ph_residuals()
  } else {
    raw_ph <- as.data.frame(captured_ph$value$table, check.names = FALSE)
    raw_terms <- rownames(raw_ph)
    display_terms <- raw_terms
    display_terms[raw_terms == ".feature"] <- feature
    for (i in seq_along(covariate_aliases)) {
      display_terms[raw_terms == covariate_aliases[[i]]] <- covariates[[i]]
    }
    ph_table <- tibble::tibble(
      feature = feature,
      term = display_terms,
      scope = ifelse(raw_terms == "GLOBAL", "global", ifelse(raw_terms == ".feature", "feature", "covariate")),
      rho = if ("rho" %in% colnames(raw_ph)) as.numeric(raw_ph[["rho"]]) else rep(NA_real_, nrow(raw_ph)),
      chi_square = as.numeric(raw_ph[["chisq"]]),
      p_value = as.numeric(raw_ph[["p"]]),
      status = "ok",
      error = NA_character_
    )
    feature_ph <- ph_table$p_value[ph_table$scope == "feature"]
    global_ph <- ph_table$p_value[ph_table$scope == "global"]
    association$ph_p_value <- if (length(feature_ph)) feature_ph[[1]] else NA_real_
    association$ph_global_p_value <- if (length(global_ph)) global_ph[[1]] else NA_real_
    ph_residuals <- .sn_survival_ph_residuals(
      captured_ph$value, feature, covariates, covariate_aliases
    )
  }

  list(
    association = association,
    performance = performance,
    proportional_hazards = ph_table,
    proportional_hazards_residuals = ph_residuals,
    fit = fit,
    ph_fit = captured_ph$value,
    warnings = unique(c(captured$warnings, captured_ph$warnings))
  )
}

.sn_survival_fixed_cutpoint <- function(group_cutpoint, feature) {
  if (length(group_cutpoint) == 1L) return(as.numeric(group_cutpoint))
  if (is_null(names(group_cutpoint)) || !feature %in% names(group_cutpoint)) {
    stop("For multiple fixed cutpoints, `group_cutpoint` must be named by feature.", call. = FALSE)
  }
  as.numeric(group_cutpoint[[feature]])
}

.sn_survival_group_feature <- function(values, feature, method, probability, fixed, labels,
                                       eligible = rep(TRUE, length(values))) {
  if (!is.numeric(values)) stop("Survival features must be numeric before grouping.", call. = FALSE)
  observed <- values[eligible & is.finite(values)]
  if (length(observed) < 2L || length(unique(observed)) < 2L) {
    stop("The feature must have at least two finite values for survival grouping.", call. = FALSE)
  }
  cutpoint <- switch(
    method,
    median = stats::median(observed),
    quantile = as.numeric(stats::quantile(observed, probs = probability, names = FALSE, type = 7)),
    fixed = .sn_survival_fixed_cutpoint(fixed, feature)
  )
  if (length(cutpoint) != 1L || !is.finite(cutpoint)) stop("The survival grouping cutpoint must be finite.", call. = FALSE)
  group <- ifelse(!is.finite(values), NA_character_, ifelse(values <= cutpoint, labels[[1]], labels[[2]]))
  list(group = factor(group, levels = labels), cutpoint = cutpoint)
}

.sn_survival_empty_curve <- function() {
  tibble::tibble(
    feature = character(), group = character(), time = numeric(), n_total = integer(),
    n_risk = integer(), n_event = integer(), n_censor = integer(), survival = numeric(),
    standard_error = numeric(), conf_low = numeric(), conf_high = numeric(),
    cutpoint = numeric(), group_method = character()
  )
}

.sn_survival_empty_groups <- function() {
  tibble::tibble(
    sample = character(), feature = character(), value = numeric(), group = character(),
    cutpoint = numeric(), group_method = character(), included = logical()
  )
}

.sn_survival_empty_risk <- function() {
  tibble::tibble(
    feature = character(), group = character(), time = numeric(), n_risk = integer(),
    n_event = integer(), n_censor = integer(), cutpoint = numeric(),
    group_method = character()
  )
}

.sn_survival_risk_grid <- function(time, risk_times = NULL) {
  if (!is_null(risk_times)) return(sort(unique(as.numeric(risk_times))))
  maximum <- max(time)
  grid <- pretty(c(0, maximum), n = 5L)
  grid <- grid[grid >= 0 & grid <= maximum]
  grid <- sort(unique(c(0, grid, maximum)))
  if (length(grid) > 8L) {
    indices <- unique(round(seq(1, length(grid), length.out = 8L)))
    grid <- grid[indices]
  }
  grid
}

.sn_survival_risk_table <- function(fit_data, feature, cutpoint, method, labels,
                                    risk_times = NULL) {
  grid <- .sn_survival_risk_grid(fit_data$.time, risk_times)
  dplyr::bind_rows(lapply(labels, function(label) {
    group_data <- fit_data[as.character(fit_data$.group) == label, , drop = FALSE]
    tibble::tibble(
      feature = feature,
      group = label,
      time = grid,
      n_risk = as.integer(vapply(grid, function(value) sum(group_data$.time >= value), integer(1))),
      n_event = as.integer(vapply(grid, function(value) {
        sum(group_data$.time == value & group_data$.event == 1L)
      }, integer(1))),
      n_censor = as.integer(vapply(grid, function(value) {
        sum(group_data$.time == value & group_data$.event == 0L)
      }, integer(1))),
      cutpoint = cutpoint,
      group_method = method
    )
  }))
}

.sn_survival_km_error <- function(feature, method, error) {
  list(
    curves = .sn_survival_empty_curve(),
    groups = .sn_survival_empty_groups(),
    risk = .sn_survival_empty_risk(),
    log_rank = tibble::tibble(
      feature = feature, cutpoint = NA_real_, group_method = method,
      chi_square = NA_real_, df = NA_real_, p_value = NA_real_, n = 0L,
      events = 0L, status = "error", error = as.character(error)
    ),
    fit = NULL,
    warnings = character(),
    status = "error",
    error = as.character(error)
  )
}

.sn_survival_fit_km <- function(values, time, event, samples, feature, method,
                                probability, fixed, labels, risk_times = NULL) {
  endpoint_complete <- !is.na(time) & !is.na(event)
  feature_finite <- if (is.numeric(values)) is.finite(values) else rep(FALSE, length(values))
  eligible <- endpoint_complete & feature_finite
  grouped <- tryCatch(
    .sn_survival_group_feature(
      values, feature, method, probability, fixed, labels, eligible = eligible
    ),
    error = function(condition) condition
  )
  if (inherits(grouped, "condition")) return(.sn_survival_km_error(feature, method, conditionMessage(grouped)))
  data <- data.frame(
    .time = time,
    .event = event,
    .feature = values,
    .group = grouped$group,
    row.names = samples,
    check.names = FALSE
  )
  complete <- eligible & !is.na(data$.group)
  assignments <- tibble::tibble(
    sample = samples,
    feature = feature,
    value = as.numeric(values),
    group = as.character(grouped$group),
    cutpoint = grouped$cutpoint,
    group_method = method,
    included = complete
  )
  fit_data <- data[complete, , drop = FALSE]
  if (nrow(fit_data) < 3L) {
    out <- .sn_survival_km_error(feature, method, "At least three complete samples are required for Kaplan-Meier analysis.")
    out$groups <- assignments
    return(out)
  }
  observed_groups <- unique(as.character(fit_data$.group))
  if (length(observed_groups) != 2L) {
    out <- .sn_survival_km_error(feature, method, "The grouping cutpoint did not produce two non-empty groups.")
    out$groups <- assignments
    return(out)
  }
  if (sum(fit_data$.event == 1L) < 1L) {
    out <- .sn_survival_km_error(feature, method, "At least one event is required for Kaplan-Meier analysis.")
    out$groups <- assignments
    return(out)
  }

  formula <- survival::Surv(.time, .event) ~ .group
  captured <- .sn_survival_capture(survival::survfit(formula, data = fit_data, conf.type = "log"))
  if (is_null(captured$value)) {
    out <- .sn_survival_km_error(feature, method, captured$error)
    out$groups <- assignments
    out$warnings <- captured$warnings
    return(out)
  }
  fit <- captured$value
  summary_fit <- summary(fit, censored = TRUE)
  strata <- sub("^\\.group=", "", as.character(summary_fit$strata))
  curves <- tibble::tibble(
    feature = feature,
    group = strata,
    time = as.numeric(summary_fit$time),
    n_total = as.integer(unname(table(fit_data$.group)[strata])),
    n_risk = as.integer(summary_fit$n.risk),
    n_event = as.integer(summary_fit$n.event),
    n_censor = as.integer(summary_fit$n.censor),
    survival = as.numeric(summary_fit$surv),
    standard_error = as.numeric(summary_fit$std.err),
    conf_low = as.numeric(summary_fit$lower),
    conf_high = as.numeric(summary_fit$upper),
    cutpoint = grouped$cutpoint,
    group_method = method
  )
  group_sizes <- table(fit_data$.group)
  baseline <- tibble::tibble(
    feature = feature,
    group = labels,
    time = 0,
    n_total = as.integer(group_sizes[labels]),
    n_risk = as.integer(group_sizes[labels]),
    n_event = 0L,
    n_censor = 0L,
    survival = 1,
    standard_error = 0,
    conf_low = 1,
    conf_high = 1,
    cutpoint = grouped$cutpoint,
    group_method = method
  )
  curves <- dplyr::bind_rows(baseline, curves) |>
    dplyr::arrange(.data$group, .data$time)
  risk <- .sn_survival_risk_table(
    fit_data, feature, grouped$cutpoint, method, labels, risk_times
  )

  captured_log_rank <- .sn_survival_capture(survival::survdiff(formula, data = fit_data))
  if (is_null(captured_log_rank$value)) {
    log_rank <- tibble::tibble(
      feature = feature, cutpoint = grouped$cutpoint, group_method = method,
      chi_square = NA_real_, df = NA_real_, p_value = NA_real_, n = nrow(fit_data),
      events = sum(fit_data$.event == 1L), status = "error", error = captured_log_rank$error
    )
  } else {
    comparison <- captured_log_rank$value
    df <- max(1L, length(comparison$n) - 1L)
    log_rank <- tibble::tibble(
      feature = feature, cutpoint = grouped$cutpoint, group_method = method,
      chi_square = unname(comparison$chisq), df = as.numeric(df),
      p_value = stats::pchisq(comparison$chisq, df = df, lower.tail = FALSE),
      n = nrow(fit_data), events = sum(fit_data$.event == 1L),
      status = "ok", error = NA_character_
    )
  }
  list(
    curves = curves,
    groups = assignments,
    risk = risk,
    log_rank = log_rank,
    fit = fit,
    warnings = unique(c(captured$warnings, captured_log_rank$warnings)),
    status = if (identical(log_rank$status[[1]], "ok")) "ok" else "partial",
    error = log_rank$error[[1]]
  )
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
#' @param group_method Feature grouping used for Kaplan-Meier analysis.
#' @param group_quantile Quantile used when `group_method = "quantile"`.
#' @param group_cutpoint Fixed scalar or feature-named cutpoints.
#' @param group_labels Labels ordered as lower/equal and higher than the cutpoint.
#' @param ties Cox partial-likelihood tie method.
#' @param risk_times Optional non-negative times shown in the risk table. The
#'   default uses at most eight deterministic pretty breaks per feature.
#' @return A validated survival result with one adjusted Cox model per feature.
#' @export
sn_run_survival <- function(object, time, event, features, covariates = NULL,
                            metadata = NULL, assay = NULL, store_name = "bulk_survival",
                            group_method = c("median", "quantile", "fixed", "none"),
                            group_quantile = 0.5, group_cutpoint = NULL,
                            group_labels = c("Low", "High"),
                            ties = c("efron", "breslow", "exact"),
                            risk_times = NULL) {
  check_installed("survival", reason = "to fit Cox proportional-hazards models.")
  input <- .sn_bulk_input(object, metadata, assay)
  if (!is.character(time) || length(time) != 1L || is.na(time) || !nzchar(time)) {
    stop("`time` must name one survival-time metadata column.", call. = FALSE)
  }
  if (!is.character(event) || length(event) != 1L || is.na(event) || !nzchar(event)) {
    stop("`event` must name one event-status metadata column.", call. = FALSE)
  }
  if (identical(time, event)) stop("`time` and `event` must name different columns.", call. = FALSE)
  if (!is.character(features) || length(features) == 0L || anyNA(features) || any(!nzchar(features)) || anyDuplicated(features)) {
    stop("`features` must be a unique, non-empty character vector.", call. = FALSE)
  }
  if (!is_null(covariates) && (!is.character(covariates) || anyNA(covariates) || any(!nzchar(covariates)) || anyDuplicated(covariates))) {
    stop("`covariates` must be NULL or a unique character vector.", call. = FALSE)
  }
  endpoint_overlap <- intersect(c(features, covariates), c(time, event))
  if (length(endpoint_overlap) > 0L) {
    stop(
      "Survival endpoints cannot also be used as features or covariates: ",
      paste(endpoint_overlap, collapse = ", "), ".",
      call. = FALSE
    )
  }
  feature_covariate_overlap <- intersect(features, covariates)
  if (length(feature_covariate_overlap) > 0L) {
    stop(
      "Survival features and covariates must be distinct: ",
      paste(feature_covariate_overlap, collapse = ", "), ".",
      call. = FALSE
    )
  }
  group_method <- match.arg(group_method)
  ties <- match.arg(ties)
  if (!is.character(group_labels) || length(group_labels) != 2L || anyNA(group_labels) || any(!nzchar(group_labels)) || anyDuplicated(group_labels)) {
    stop("`group_labels` must contain two distinct, non-empty labels.", call. = FALSE)
  }
  if (identical(group_method, "quantile") &&
      (!is.numeric(group_quantile) || length(group_quantile) != 1L || !is.finite(group_quantile) || group_quantile <= 0 || group_quantile >= 1)) {
    stop("`group_quantile` must be one finite number strictly between 0 and 1.", call. = FALSE)
  }
  if (identical(group_method, "fixed")) {
    if (!is.numeric(group_cutpoint) || length(group_cutpoint) == 0L || any(!is.finite(group_cutpoint))) {
      stop("`group_cutpoint` must supply finite numeric cutpoints for fixed grouping.", call. = FALSE)
    }
    if (length(group_cutpoint) > 1L && (is_null(names(group_cutpoint)) || any(!features %in% names(group_cutpoint)))) {
      stop("Multiple fixed cutpoints must be named for every requested feature.", call. = FALSE)
    }
  }
  if (!is_null(risk_times)) {
    if (!is.numeric(risk_times) || length(risk_times) == 0L || anyNA(risk_times) ||
        any(!is.finite(risk_times)) || any(risk_times < 0)) {
      stop("`risk_times` must be NULL or finite, non-negative numeric values.", call. = FALSE)
    }
    risk_times <- sort(unique(as.numeric(risk_times)))
  }
  required <- c(time, event, covariates)
  missing <- setdiff(required, colnames(input$metadata))
  if (length(missing) > 0L) stop("Survival metadata column(s) missing: ", paste(missing, collapse = ", "), ".", call. = FALSE)
  time_values <- input$metadata[[time]]
  if (!is.numeric(time_values)) stop("Survival time must be numeric.", call. = FALSE)
  if (any(is.nan(time_values) | is.infinite(time_values))) stop("Survival time must contain only finite values or NA.", call. = FALSE)
  if (any(time_values < 0, na.rm = TRUE)) stop("Survival time cannot be negative.", call. = FALSE)
  event_values <- .sn_survival_event(input$metadata[[event]])
  endpoint_complete <- !is.na(time_values) & !is.na(event_values)
  if (sum(endpoint_complete) < 3L) stop("At least three samples require complete survival time and event status.", call. = FALSE)
  if (sum(event_values[endpoint_complete] == 1L) < 1L) stop("At least one observed event is required.", call. = FALSE)
  feature_data <- .sn_bulk_feature_data(input, features)
  rows <- vector("list", length(features))
  performance_rows <- vector("list", length(features))
  ph_rows <- vector("list", length(features))
  ph_residual_rows <- vector("list", length(features))
  curve_rows <- vector("list", length(features))
  risk_rows <- vector("list", length(features))
  group_rows <- vector("list", length(features))
  log_rank_rows <- vector("list", length(features))
  models <- stats::setNames(vector("list", length(features)), features)
  km_models <- stats::setNames(vector("list", length(features)), features)
  ph_models <- stats::setNames(vector("list", length(features)), features)
  backend_warnings <- character()
  covariate_aliases <- if (length(covariates)) paste0(".covariate_", seq_along(covariates)) else character()
  for (i in seq_along(features)) {
    feature <- features[[i]]
    model_data <- data.frame(
      .time = time_values,
      .event = event_values,
      .feature = feature_data[[feature]],
      check.names = FALSE
    )
    for (j in seq_along(covariates)) model_data[[covariate_aliases[[j]]]] <- input$metadata[[covariates[[j]]]]
    cox <- .sn_survival_fit_cox(model_data, feature, covariates, ties)
    km <- if (identical(group_method, "none")) {
      list(
        curves = .sn_survival_empty_curve(), groups = .sn_survival_empty_groups(),
        risk = .sn_survival_empty_risk(),
        log_rank = tibble::tibble(
          feature = feature, cutpoint = NA_real_, group_method = "none",
          chi_square = NA_real_, df = NA_real_, p_value = NA_real_, n = 0L,
          events = 0L, status = "not_requested", error = NA_character_
        ),
        fit = NULL, warnings = character(), status = "not_requested", error = NA_character_
      )
    } else {
      .sn_survival_fit_km(
        values = feature_data[[feature]], time = time_values, event = event_values,
        samples = rownames(input$metadata), feature = feature, method = group_method,
        probability = group_quantile, fixed = group_cutpoint, labels = group_labels,
        risk_times = risk_times
      )
    }
    cox$association$km_status <- km$status
    cox$association$km_error <- km$error
    rows[[i]] <- cox$association
    performance_rows[[i]] <- cox$performance
    ph_rows[[i]] <- cox$proportional_hazards
    ph_residual_rows[[i]] <- cox$proportional_hazards_residuals
    curve_rows[[i]] <- km$curves
    risk_rows[[i]] <- km$risk
    group_rows[[i]] <- km$groups
    log_rank_rows[[i]] <- km$log_rank
    models[i] <- list(cox$fit)
    km_models[i] <- list(km$fit)
    ph_models[i] <- list(cox$ph_fit)
    if (length(cox$warnings)) {
      backend_warnings <- c(backend_warnings, paste0("Cox feature '", feature, "': ", cox$warnings))
    }
    if (length(km$warnings)) {
      backend_warnings <- c(backend_warnings, paste0("Kaplan-Meier feature '", feature, "': ", km$warnings))
    }
  }
  table <- dplyr::bind_rows(rows)
  table$adjusted_p_value <- stats::p.adjust(table$p_value, method = "BH")
  log_rank <- dplyr::bind_rows(log_rank_rows)
  log_rank$adjusted_p_value <- stats::p.adjust(log_rank$p_value, method = "BH")
  curves <- dplyr::bind_rows(curve_rows)
  risk <- dplyr::bind_rows(risk_rows)
  if (nrow(risk) == 0L) risk <- .sn_survival_empty_risk()
  cumulative_hazard <- if (nrow(curves) == 0L) {
    tibble::tibble(
      feature = character(), group = character(), time = numeric(), cumulative_hazard = numeric(),
      conf_low = numeric(), conf_high = numeric(), cutpoint = numeric(), group_method = character()
    )
  } else {
    dplyr::transmute(
      curves,
      feature = .data$feature,
      group = .data$group,
      time = .data$time,
      cumulative_hazard = -log(pmax(.data$survival, .Machine$double.xmin)),
      conf_low = -log(pmax(.data$conf_high, .Machine$double.xmin)),
      conf_high = -log(pmax(.data$conf_low, .Machine$double.xmin)),
      cutpoint = .data$cutpoint,
      group_method = .data$group_method
    )
  }
  models$.kaplan_meier <- km_models
  models$.proportional_hazards <- ph_models
  result <- .sn_bulk_result(
    "bulk_survival", store_name, "cox",
    list(source = input$source, samples = ncol(input$matrix), features = features),
    list(
      time = time, event = event, covariates = covariates, ties = ties,
      group_method = group_method, group_quantile = group_quantile,
      group_cutpoint = group_cutpoint, group_labels = group_labels,
      risk_times = risk_times,
      p_adjust_method = "BH"
    ),
    list(
      primary = table,
      survival = table,
      curves = curves,
      risk = risk,
      cumulative_hazard = cumulative_hazard,
      groups = dplyr::bind_rows(group_rows),
      log_rank = log_rank,
      proportional_hazards = dplyr::bind_rows(ph_rows),
      proportional_hazards_residuals = dplyr::bind_rows(ph_residual_rows),
      performance = dplyr::bind_rows(performance_rows)
    ),
    models = models,
    diagnostics = list(
      requested_models = length(features),
      models = sum(table$status == "ok"),
      failed_models = sum(table$status == "error"),
      kaplan_meier_models = sum(table$km_status == "ok"),
      proportional_hazards_warnings = sum(table$ph_p_value < 0.05, na.rm = TRUE),
      global_proportional_hazards_warnings = sum(table$ph_global_p_value < 0.05, na.rm = TRUE),
      endpoint_complete_samples = sum(endpoint_complete)
    ),
    warnings = unique(backend_warnings[nzchar(backend_warnings)]),
    schema_version = "1.0.0"
  )
  sn_validate_result(result)
  result
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
