.sn_core_metabolic_signatures <- function() {
  list(
    glycolysis = c("HK1", "HK2", "GPI", "PFKM", "PFKP", "ALDOA", "GAPDH", "PGK1", "PGAM1", "ENO1", "PKM", "LDHA"),
    tca_cycle = c("CS", "ACO2", "IDH2", "OGDH", "DLST", "SUCLG1", "SDHA", "SDHB", "FH", "MDH2"),
    oxidative_phosphorylation = c("NDUFS1", "NDUFS2", "SDHA", "UQCRC1", "UQCRC2", "COX4I1", "COX5A", "ATP5F1A", "ATP5F1B"),
    fatty_acid_oxidation = c("CPT1A", "CPT2", "ACADM", "ACADL", "ACADVL", "HADHA", "HADHB", "ETFA", "ETFB"),
    fatty_acid_synthesis = c("ACLY", "ACACA", "FASN", "SCD", "ELOVL6", "HACD3", "HSD17B12"),
    glutaminolysis = c("SLC1A5", "SLC38A1", "GLS", "GLUD1", "GOT1", "GOT2", "GPT2", "SLC25A11"),
    pentose_phosphate = c("G6PD", "PGLS", "PGD", "RPIA", "RPE", "TKT", "TALDO1"),
    one_carbon = c("SHMT1", "SHMT2", "MTHFD1", "MTHFD2", "MTHFD1L", "MTHFR", "TYMS", "DHFR"),
    nucleotide_synthesis = c("PPAT", "GART", "PFAS", "PAICS", "ADSL", "ATIC", "CAD", "DHODH", "UMPS"),
    cholesterol_biosynthesis = c("HMGCS1", "HMGCR", "MVK", "MVD", "FDPS", "FDFT1", "SQLE", "LSS", "DHCR7")
  )
}

#' Retrieve curated core metabolic signatures
#'
#' @param species Species label. Matching is case-insensitive, so the same
#'   symbols support common human and mouse feature conventions.
#' @return A named list of pathway gene vectors.
#' @export
sn_metabolic_signatures <- function(species = c("human", "mouse")) {
  species <- match.arg(species)
  .sn_core_metabolic_signatures()
}

.sn_metabolism_score_matrix <- function(object,
                                        signatures,
                                        scoring_method,
                                        assay,
                                        layer,
                                        min_genes,
                                        backend_control) {
  expression <- .sn_annotation_expression(object, assay = assay, layer = layer)
  matched <- .sn_match_program_features(signatures, rownames(expression$matrix), min_genes = min_genes)
  matrix <- expression$matrix
  if (scoring_method %in% c("gsva", "ssgsea") && inherits(matrix, "sparseMatrix") && ncol(matrix) > 5000L) {
    stop("Per-cell GSVA/ssGSEA above 5,000 cells is disabled; use UCell or aggregate the object first.", call. = FALSE)
  }
  scores <- switch(
    scoring_method,
    mean = .sn_score_programs_mean(matrix, matched$signatures),
    ucell = .sn_score_programs_ucell(matrix, matched$signatures, backend_control$ucell %||% list()),
    gsva = .sn_score_programs_gsva(matrix, matched$signatures, "gsva", backend_control$gsva %||% list()),
    ssgsea = .sn_score_programs_gsva(matrix, matched$signatures, "ssgsea", backend_control$ssgsea %||% list())
  )
  scores <- as.matrix(scores)
  scores <- scores[names(matched$signatures), colnames(object), drop = FALSE]
  list(scores = scores, coverage = matched$coverage, expression = expression)
}

.sn_extract_scmetabolism_scores <- function(object) {
  assay <- object[["METABOLISM"]]
  candidates <- list(
    tryCatch(methods::slot(assay, "score"), error = function(e) NULL),
    tryCatch(assay@misc$score, error = function(e) NULL),
    tryCatch(SeuratObject::LayerData(assay, layer = "score"), error = function(e) NULL),
    tryCatch(SeuratObject::LayerData(assay, layer = "data"), error = function(e) NULL)
  )
  hit <- candidates[!vapply(candidates, is_null, logical(1))]
  if (length(hit) == 0L) stop("Could not extract a metabolism score matrix from the scMetabolism result.", call. = FALSE)
  as.matrix(hit[[1]])
}

.sn_standardize_metabolism_output <- function(output, cells, method) {
  if (is.character(output) && length(output) == 1L && file.exists(output)) {
    output <- utils::read.delim(output, row.names = 1, check.names = FALSE)
  }
  if (is.data.frame(output) && all(c("cell", "pathway", "score") %in% names(output))) {
    table <- tibble::as_tibble(output)
    table$method <- table$method %||% method
    return(list(table = table, matrix = NULL))
  }
  output <- as.matrix(output)
  if (all(cells %in% rownames(output)) && !all(cells %in% colnames(output))) output <- t(output)
  if (!all(cells %in% colnames(output))) {
    stop("Metabolism backend output must have cells in columns or provide cell/pathway/score columns.", call. = FALSE)
  }
  output <- output[, cells, drop = FALSE]
  table <- dplyr::bind_rows(lapply(seq_len(nrow(output)), function(index) {
    tibble::tibble(
      cell = cells, pathway = rownames(output)[[index]],
      score = as.numeric(output[index, ]), method = method
    )
  }))
  list(table = table, matrix = output)
}

.sn_run_metabolism_backend <- function(object,
                                       method,
                                       assay,
                                       layer,
                                       backend_control) {
  if (is.function(backend_control$runner)) {
    output <- backend_control$runner(object = object, assay = assay, layer = layer, method = method)
    return(.sn_standardize_metabolism_output(output, colnames(object), method))
  }
  if (!is_null(backend_control$result)) {
    return(.sn_standardize_metabolism_output(backend_control$result, colnames(object), method))
  }
  if (identical(method, "scmetabolism")) {
    check_installed("scMetabolism", reason = "to run the scMetabolism backend.")
    controls <- backend_control
    controls$runner <- NULL
    controls$result <- NULL
    defaults <- list(
      obj = object, method = controls$scoring_method %||% "AUCell",
      imputation = FALSE, ncores = 1, metabolism.type = controls$collection %||% "KEGG"
    )
    controls$scoring_method <- NULL
    controls$collection <- NULL
    scored <- do.call(scMetabolism::sc.metabolism.Seurat, utils::modifyList(defaults, controls, keep.null = TRUE))
    return(.sn_standardize_metabolism_output(.sn_extract_scmetabolism_scores(scored), colnames(object), method))
  }
  stop(
    "The ", method, " adapter requires `backend_control$runner` or `backend_control$result` from the optional external backend.",
    call. = FALSE
  )
}

.sn_metabolism_sample_scores <- function(scores, metadata, sample_by, condition_by, group_by) {
  scores$sample <- if (is_null(sample_by)) scores$cell else as.character(metadata[scores$cell, sample_by])
  scores$condition <- if (is_null(condition_by)) "all" else as.character(metadata[scores$cell, condition_by])
  scores$group <- if (is_null(group_by)) "all" else as.character(metadata[scores$cell, group_by])
  key <- interaction(scores$sample, scores$condition, scores$group, scores$pathway, drop = TRUE, lex.order = TRUE)
  dplyr::bind_rows(lapply(split(scores, key), function(table) {
    tibble::tibble(
      sample = table$sample[[1]], condition = table$condition[[1]],
      group = table$group[[1]], pathway = table$pathway[[1]],
      score = mean(table$score, na.rm = TRUE), n_cells = nrow(table)
    )
  }))
}

.sn_metabolism_differential <- function(sample_scores, contrast = NULL) {
  conditions <- unique(sample_scores$condition[!is.na(sample_scores$condition)])
  if (length(conditions) < 2L) return(tibble::tibble())
  contrast <- contrast %||% conditions[seq_len(min(2L, length(conditions)))]
  if (length(contrast) != 2L || any(!contrast %in% conditions)) {
    stop("`contrast` must contain two observed condition levels.", call. = FALSE)
  }
  groups <- split(sample_scores, interaction(sample_scores$group, sample_scores$pathway, drop = TRUE, lex.order = TRUE))
  result <- dplyr::bind_rows(lapply(groups, function(table) {
    x <- table$score[table$condition == contrast[[1]]]
    y <- table$score[table$condition == contrast[[2]]]
    p_value <- if (length(x) > 0L && length(y) > 0L) {
      suppressWarnings(stats::wilcox.test(x, y, exact = FALSE)$p.value)
    } else {
      NA_real_
    }
    tibble::tibble(
      group = table$group[[1]], pathway = table$pathway[[1]],
      comparison = paste(contrast[[1]], "vs", contrast[[2]]),
      estimate = mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE),
      mean_1 = mean(x, na.rm = TRUE), mean_2 = mean(y, na.rm = TRUE),
      n_1 = length(x), n_2 = length(y), p_value = p_value
    )
  }))
  result$adjusted_p_value <- stats::p.adjust(result$p_value, method = "BH")
  result
}

.sn_add_metabolism_metadata <- function(object, scores, store_name) {
  pathways <- unique(scores$pathway)
  metadata <- data.frame(row.names = colnames(object))
  prefix <- gsub("[^[:alnum:]_]+", "_", store_name)
  for (pathway in pathways) {
    current <- scores[scores$pathway == pathway, , drop = FALSE]
    values <- stats::setNames(current$score, current$cell)
    suffix <- gsub("[^[:alnum:]_]+", "_", pathway)
    metadata[[paste(prefix, suffix, sep = "_")]] <- as.numeric(values[colnames(object)])
  }
  SeuratObject::AddMetaData(object, metadata = metadata)
}

#' Run unified single-cell metabolic activity analysis
#'
#' The default gene-set workflow uses curated signatures with UCell, GSVA,
#' ssGSEA, or a sparse-aware mean score. Optional scMetabolism, scFEA, and
#' Compass outputs are standardized through the same result contract.
#'
#' @param object A Seurat object.
#' @param method Metabolism backend.
#' @param signatures Optional named pathway gene sets. Defaults to curated core pathways.
#' @param scoring_method Gene-set scoring method for `method = "geneset"`.
#' @param assay,layer Expression assay and layer.
#' @param sample_by Sample/patient metadata column used as the inferential unit.
#' @param condition_by Optional condition metadata column.
#' @param group_by Optional cell-type/state column for stratified comparisons.
#' @param contrast Optional two condition levels.
#' @param species Species used for the curated signatures.
#' @param min_genes Minimum matched genes per pathway.
#' @param store_name Stored result name.
#' @param backend_control Backend options. For scFEA/Compass, provide a
#'   `runner` function or parsed `result`; this keeps heavyweight
#'   runtimes outside the default R dependency set.
#' @param return_object Return the modified object or unified result.
#' @return A Seurat object or unified metabolism result.
#' @examples
#' \dontrun{
#' object <- sn_run_metabolism(
#'   object, scoring_method = "ucell", sample_by = "patient",
#'   condition_by = "condition", group_by = "cell_type"
#' )
#' metabolism <- sn_get_result(object, "metabolism", "metabolism")
#' metabolism$tables$differential
#' }
#' @export
sn_run_metabolism <- function(object,
                              method = c("geneset", "scmetabolism", "scfea", "compass"),
                              signatures = NULL,
                              scoring_method = c("ucell", "gsva", "ssgsea", "mean"),
                              assay = NULL,
                              layer = "data",
                              sample_by = NULL,
                              condition_by = NULL,
                              group_by = NULL,
                              contrast = NULL,
                              species = NULL,
                              min_genes = 3L,
                              store_name = "metabolism",
                              backend_control = list(),
                              return_object = TRUE) {
  .sn_validate_result_object(object)
  method <- match.arg(method)
  scoring_method <- match.arg(scoring_method)
  metadata <- object[[]]
  required <- c(sample_by, condition_by, group_by)
  required <- required[!is.na(required) & nzchar(required)]
  missing <- setdiff(required, colnames(metadata))
  if (length(missing) > 0L) stop("Metadata column(s) not found: ", paste(missing, collapse = ", "), call. = FALSE)
  if (!is_null(sample_by) && !is_null(condition_by)) {
    .sn_validate_constant_within_sample(metadata, sample_col = sample_by, group_col = condition_by)
  }
  species <- species %||% tryCatch(sn_get_species(object), error = function(e) "human")
  coverage <- tibble::tibble()
  expression_info <- list(assay = assay %||% SeuratObject::DefaultAssay(object), layer = layer)
  if (identical(method, "geneset")) {
    signatures <- signatures %||% sn_metabolic_signatures(if (tolower(species) %in% c("mouse", "mm10", "mus_musculus")) "mouse" else "human")
    scored <- .sn_metabolism_score_matrix(
      object, signatures, scoring_method, assay, layer, min_genes, backend_control
    )
    standardized <- .sn_standardize_metabolism_output(scored$scores, colnames(object), method)
    scores <- standardized$table
    coverage <- scored$coverage
    expression_info <- scored$expression
    backend <- scoring_method
  } else {
    standardized <- .sn_run_metabolism_backend(object, method, assay, layer, backend_control)
    scores <- standardized$table
    backend <- method
  }
  scores <- scores[scores$cell %in% colnames(object), , drop = FALSE]
  if (nrow(scores) == 0L) stop("Metabolism backend returned no scores for object cells.", call. = FALSE)
  object <- .sn_add_metabolism_metadata(object, scores, store_name)
  sample_scores <- .sn_metabolism_sample_scores(scores, metadata, sample_by, condition_by, group_by)
  differential <- .sn_metabolism_differential(sample_scores, contrast)
  warnings <- if (is_null(sample_by) && !is_null(condition_by)) {
    "Cells were used as exploratory units because `sample_by` was not supplied."
  } else {
    character()
  }
  result <- list(
    schema_version = "1.0", analysis_type = "metabolism", name = store_name,
    method = method, backend = backend,
    input = list(
      assay = expression_info$assay, layer = expression_info$layer,
      sample_by = sample_by, condition_by = condition_by, group_by = group_by,
      species = species
    ),
    parameters = list(scoring_method = scoring_method, min_genes = min_genes, contrast = contrast),
    tables = list(
      primary = scores, coverage = coverage,
      sample_scores = sample_scores, differential = differential
    ),
    embeddings = list(), graphs = list(), models = list(),
    diagnostics = list(
      pathways = length(unique(scores$pathway)), cells = length(unique(scores$cell)),
      samples = length(unique(sample_scores$sample)),
      inferential_unit = if (is_null(sample_by)) "cell" else "sample"
    ),
    warnings = warnings, provenance = .sn_analysis_provenance()
  )
  sn_validate_result(result)
  object <- sn_store_result(object, "metabolism", store_name, result)
  object <- .sn_log_seurat_command(object, assay = expression_info$assay, name = "sn_run_metabolism")
  if (isTRUE(return_object)) object else sn_get_result(object, "metabolism", store_name)
}
