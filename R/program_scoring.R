.sn_normalize_program_signatures <- function(signatures, species = NULL) {
  if (is.data.frame(signatures)) {
    program_col <- intersect(c("program", "signature", "set", "name"), colnames(signatures))[[1]] %||% NULL
    gene_col <- intersect(c("gene", "feature", "symbol"), colnames(signatures))[[1]] %||% NULL
    if (is_null(program_col) || is_null(gene_col)) {
      stop("Signature data frames require program/signature and gene/feature columns.", call. = FALSE)
    }
    signatures <- split(as.character(signatures[[gene_col]]), as.character(signatures[[program_col]]))
  } else if (is.character(signatures) && !is.null(names(signatures)) && any(nzchar(names(signatures)))) {
    signatures <- split(as.character(signatures), names(signatures))
  } else if (is.character(signatures)) {
    if (is_null(species)) {
      stop("`species` is required when character `signatures` are bundled signature queries.", call. = FALSE)
    }
    queries <- unique(as.character(signatures))
    signatures <- stats::setNames(lapply(queries, function(query) {
      sn_get_signatures(species = species, category = query)
    }), queries)
  }
  if (!is.list(signatures) || length(signatures) == 0L || is.null(names(signatures)) || any(!nzchar(names(signatures)))) {
    stop("`signatures` must resolve to a non-empty named list of gene vectors.", call. = FALSE)
  }
  signatures <- lapply(signatures, function(genes) {
    genes <- unique(as.character(unlist(genes, use.names = FALSE)))
    genes[!is.na(genes) & nzchar(genes)]
  })
  signatures[lengths(signatures) > 0L]
}

.sn_match_program_features <- function(signatures, features, min_genes = 1L) {
  feature_keys <- toupper(sub("\\.[0-9]+$", "", features))
  feature_lookup <- stats::setNames(features, feature_keys)
  matched <- lapply(signatures, function(genes) {
    unique(unname(feature_lookup[toupper(sub("\\.[0-9]+$", "", genes))])) |>
      stats::na.omit() |>
      as.character()
  })
  coverage <- tibble::tibble(
    program = names(signatures),
    n_genes = unname(lengths(signatures)),
    n_matched = unname(lengths(matched)),
    coverage = unname(lengths(matched) / pmax(1L, lengths(signatures))),
    matched_genes = vapply(matched, paste, collapse = ";", character(1)),
    missing_genes = vapply(names(signatures), function(program) {
      keys <- toupper(sub("\\.[0-9]+$", "", signatures[[program]]))
      paste(signatures[[program]][!keys %in% feature_keys], collapse = ";")
    }, character(1))
  )
  keep <- lengths(matched) >= as.integer(min_genes)
  if (!any(keep)) {
    stop("No signature retained at least `min_genes` matched features.", call. = FALSE)
  }
  if (any(!keep)) {
    .sn_log_warn(
      "Dropping {sum(!keep)} signature(s) with fewer than {min_genes} matched features: ",
      "{paste(names(matched)[!keep], collapse = ', ')}."
    )
  }
  list(signatures = matched[keep], coverage = coverage, retained = names(matched)[keep])
}

.sn_score_programs_mean <- function(matrix, signatures) {
  scores <- vapply(signatures, function(features) {
    Matrix::colMeans(matrix[features, , drop = FALSE])
  }, numeric(ncol(matrix)))
  scores <- t(scores)
  rownames(scores) <- names(signatures)
  colnames(scores) <- colnames(matrix)
  scores
}

.sn_score_programs_ucell <- function(matrix, signatures, control) {
  check_installed("UCell", reason = "to run UCell program scoring.")
  defaults <- list(matrix = matrix, features = signatures, name = "", ncores = 1)
  scores <- do.call(UCell::ScoreSignatures_UCell, utils::modifyList(defaults, control, keep.null = TRUE))
  scores <- as.matrix(scores)
  scores <- t(scores)
  rownames(scores) <- names(signatures)
  scores
}

.sn_score_programs_aucell <- function(matrix, signatures, control) {
  check_installed("AUCell", reason = "to run AUCell program scoring.")
  ranking_defaults <- list(exprMat = matrix, plotStats = FALSE, verbose = FALSE)
  rankings <- do.call(
    AUCell::AUCell_buildRankings,
    utils::modifyList(ranking_defaults, control$rankings %||% list(), keep.null = TRUE)
  )
  auc_defaults <- list(
    geneSets = signatures,
    rankings = rankings,
    aucMaxRank = max(1L, min(nrow(matrix), max(ceiling(0.05 * nrow(matrix)), max(lengths(signatures))))),
    verbose = FALSE
  )
  auc <- do.call(
    AUCell::AUCell_calcAUC,
    utils::modifyList(auc_defaults, control$auc %||% list(), keep.null = TRUE)
  )
  AUCell::getAUC(auc)
}

.sn_score_programs_gsva <- function(matrix, signatures, method, control) {
  check_installed("GSVA", reason = paste0("to run ", method, " program scoring."))
  namespace <- asNamespace("GSVA")
  if (exists(if (method == "gsva") "gsvaParam" else "ssgseaParam", envir = namespace, inherits = FALSE)) {
    constructor <- get(if (method == "gsva") "gsvaParam" else "ssgseaParam", envir = namespace)
    defaults <- list(exprData = matrix, geneSets = signatures, minSize = 1L, verbose = FALSE)
    parameter_args <- utils::modifyList(defaults, control$parameters %||% list(), keep.null = TRUE)
    parameter_args <- parameter_args[names(parameter_args) %in% names(formals(constructor))]
    parameter <- do.call(constructor, parameter_args)
    return(do.call(GSVA::gsva, c(list(param = parameter), control$run %||% list())))
  }
  defaults <- list(expr = matrix, gset.idx.list = signatures, method = method, min.sz = 1L, verbose = FALSE)
  do.call(GSVA::gsva, utils::modifyList(defaults, control$legacy %||% list(), keep.null = TRUE))
}

.sn_program_score_table <- function(scores, level, group_by = NULL) {
  rows <- lapply(seq_len(nrow(scores)), function(index) {
    tibble::tibble(
      entity = colnames(scores),
      program = rownames(scores)[[index]],
      score = as.numeric(scores[index, ]),
      level = level,
      group_by = group_by %||% NA_character_
    )
  })
  dplyr::bind_rows(rows)
}

.sn_add_program_metadata <- function(object, scores, name) {
  metadata <- data.frame(row.names = colnames(object))
  prefix <- gsub("[^[:alnum:]_]+", "_", name)
  for (program in rownames(scores)) {
    suffix <- gsub("[^[:alnum:]_]+", "_", program)
    metadata[[paste(prefix, suffix, sep = "_")]] <- as.numeric(scores[program, colnames(object)])
  }
  SeuratObject::AddMetaData(object, metadata = metadata)
}

#' Score gene programs in cells or aggregated samples
#'
#' Provides one stable interface for sparse-aware UCell, AUCell, GSVA, ssGSEA,
#' and dependency-free mean-expression scoring. Per-cell scores are added to
#' Seurat metadata; aggregated scores remain in the stored result.
#'
#' @param object A \code{Seurat} object.
#' @param signatures Named list, program/gene data frame, named gene vector, or
#'   bundled signature query vector.
#' @param method Scoring backend. UCell is the default for per-cell data.
#' @param assay,layer Expression source.
#' @param name Stored-result and metadata prefix. A stable method-derived name
#'   is used when omitted.
#' @param group_by Optional metadata column used to average expression before
#'   GSVA/ssGSEA or other group-level scoring.
#' @param species Species required for bundled signature queries.
#' @param min_genes Minimum matched features required per signature.
#' @param backend_control Named backend-specific control list.
#' @param return_object Return the updated object or the unified result.
#'
#' @return A Seurat object or unified program-scoring result.
#'
#' @examples
#' \dontrun{
#' object <- sn_score_programs(
#'   object,
#'   signatures = list(T_cell = c("CD3D", "CD3E")),
#'   method = "ucell",
#'   name = "immune_programs"
#' )
#' sn_get_result(object, "program_scoring", "immune_programs")
#' }
#'
#' @export
sn_score_programs <- function(object,
                              signatures,
                              method = c("ucell", "aucell", "gsva", "ssgsea", "mean"),
                              assay = NULL,
                              layer = "data",
                              name = NULL,
                              group_by = NULL,
                              species = NULL,
                              min_genes = 1L,
                              backend_control = list(),
                              return_object = TRUE) {
  .sn_validate_result_object(object)
  method <- match.arg(method)
  name <- name %||% paste0("programs_", method)
  species <- species %||% tryCatch(sn_get_species(object), error = function(e) NULL)
  signatures <- .sn_normalize_program_signatures(signatures, species = species)
  expression <- .sn_annotation_expression(object, assay = assay, layer = layer)
  matched <- .sn_match_program_features(signatures, rownames(expression$matrix), min_genes = min_genes)
  matrix <- expression$matrix
  level <- "cell"
  if (!is_null(group_by)) {
    if (!group_by %in% colnames(object[[]])) {
      stop("`group_by` column '", group_by, "' was not found in object metadata.", call. = FALSE)
    }
    groups <- as.character(object[[group_by, drop = TRUE]])
    names(groups) <- colnames(object)
    matrix <- .sn_group_average_matrix(matrix, groups[colnames(matrix)])
    level <- "group"
  }
  if (method %in% c("gsva", "ssgsea") && is_null(group_by) && inherits(matrix, "sparseMatrix") && ncol(matrix) > 5000L) {
    stop(
      "Per-cell ", method, " scoring of more than 5,000 sparse columns is disabled to avoid accidental densification; ",
      "supply `group_by` or use UCell.",
      call. = FALSE
    )
  }

  scores <- switch(
    method,
    mean = .sn_score_programs_mean(matrix, matched$signatures),
    ucell = .sn_score_programs_ucell(matrix, matched$signatures, backend_control$ucell %||% list()),
    aucell = .sn_score_programs_aucell(matrix, matched$signatures, backend_control$aucell %||% list()),
    gsva = .sn_score_programs_gsva(matrix, matched$signatures, "gsva", backend_control$gsva %||% list()),
    ssgsea = .sn_score_programs_gsva(matrix, matched$signatures, "ssgsea", backend_control$ssgsea %||% list())
  )
  scores <- as.matrix(scores)
  scores <- scores[names(matched$signatures), , drop = FALSE]
  if (is.null(colnames(scores))) colnames(scores) <- colnames(matrix)
  score_table <- .sn_program_score_table(scores, level = level, group_by = group_by)
  if (identical(level, "cell")) {
    object <- .sn_add_program_metadata(object, scores, name = name)
  }

  result <- list(
    schema_version = "1.0",
    analysis_type = "program_scoring",
    name = name,
    method = method,
    backend = method,
    input = list(
      assay = expression$assay,
      layer = expression$layer,
      cells = ncol(object),
      features = nrow(expression$matrix),
      group_by = group_by,
      species = species
    ),
    parameters = list(min_genes = min_genes, level = level),
    tables = list(scores = score_table, coverage = matched$coverage),
    embeddings = list(),
    graphs = list(),
    models = list(),
    diagnostics = list(
      retained_programs = names(matched$signatures),
      dropped_programs = setdiff(names(signatures), names(matched$signatures)),
      minimum_coverage = min(matched$coverage$coverage[matched$coverage$program %in% names(matched$signatures)])
    ),
    warnings = character(),
    provenance = .sn_analysis_provenance(random_seed = backend_control$seed %||% NA_integer_)
  )
  sn_validate_result(result)
  object <- sn_store_result(object, "program_scoring", name, result)
  object <- .sn_log_seurat_command(object = object, assay = expression$assay, name = "sn_score_programs")
  if (isTRUE(return_object)) object else sn_get_result(object, "program_scoring", name)
}

#' Test program activity between conditions
#'
#' Aggregates cell-level scores to the sample level before inference when
#' \code{sample_by} is supplied. This prevents cells from being treated as
#' independent biological replicates.
#'
#' @param object A Seurat object containing a stored program-scoring result.
#' @param score_name Stored scoring result name.
#' @param condition_by Condition metadata column.
#' @param sample_by Sample/patient metadata column. Strongly recommended for
#'   inference.
#' @param group_by Optional cell-type or state column used for stratified tests.
#' @param contrast Optional two condition levels; defaults to the first two.
#' @param method \code{"wilcox"} or \code{"limma"}.
#' @param store_name Result name.
#' @param return_object Return the object or unified result.
#'
#' @return A Seurat object or unified program-comparison result.
#'
#' @examples
#' \dontrun{
#' object <- sn_test_programs(
#'   object, "immune_programs", condition_by = "condition",
#'   sample_by = "patient", group_by = "cell_type"
#' )
#' }
#'
#' @export
sn_test_programs <- function(object,
                             score_name,
                             condition_by,
                             sample_by = NULL,
                             group_by = NULL,
                             contrast = NULL,
                             method = c("wilcox", "limma"),
                             store_name = NULL,
                             return_object = TRUE) {
  .sn_validate_result_object(object)
  method <- match.arg(method)
  metadata_columns <- c(condition_by, sample_by, group_by)
  missing <- setdiff(metadata_columns[!is.na(metadata_columns) & nzchar(metadata_columns)], colnames(object[[]]))
  if (length(missing) > 0L) {
    stop("Metadata column(s) not found: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  scored <- sn_get_result(object, "program_scoring", score_name)
  scores <- scored$tables$scores
  if (!identical(unique(scores$level), "cell")) {
    stop("`sn_test_programs()` currently requires cell-level stored scores.", call. = FALSE)
  }
  metadata <- object[[]]
  indices <- match(scores$entity, rownames(metadata))
  scores$condition <- as.character(metadata[[condition_by]][indices])
  scores$sample <- if (is_null(sample_by)) scores$entity else as.character(metadata[[sample_by]][indices])
  scores$group <- if (is_null(group_by)) "all" else as.character(metadata[[group_by]][indices])
  if (is_null(sample_by)) {
    .sn_log_warn("`sample_by` was not supplied; cells are exploratory units, not biological replicates.")
  }
  aggregate_key <- interaction(scores$sample, scores$condition, scores$group, scores$program, drop = TRUE, lex.order = TRUE)
  aggregated <- dplyr::bind_rows(lapply(split(scores, aggregate_key), function(table) {
    tibble::tibble(
      sample = table$sample[[1]],
      condition = table$condition[[1]],
      group = table$group[[1]],
      program = table$program[[1]],
      score = mean(table$score, na.rm = TRUE),
      n_cells = nrow(table)
    )
  }))
  conditions <- unique(aggregated$condition)
  contrast <- contrast %||% conditions[seq_len(min(2L, length(conditions)))]
  if (length(contrast) != 2L || any(!contrast %in% conditions)) {
    stop("`contrast` must contain two observed condition levels.", call. = FALSE)
  }
  comparison <- paste(contrast[[1]], "vs", contrast[[2]])
  test_groups <- split(aggregated, interaction(aggregated$group, aggregated$program, drop = TRUE, lex.order = TRUE))
  tests <- dplyr::bind_rows(lapply(test_groups, function(table) {
    table <- table[table$condition %in% contrast, , drop = FALSE]
    x <- table$score[table$condition == contrast[[1]]]
    y <- table$score[table$condition == contrast[[2]]]
    estimate <- mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)
    if (length(x) == 0L || length(y) == 0L) {
      p_value <- NA_real_
    } else if (method == "wilcox") {
      p_value <- suppressWarnings(stats::wilcox.test(x, y, exact = FALSE)$p.value)
    } else {
      check_installed("limma", reason = "to test program scores with limma.")
      response <- c(x, y)
      design <- stats::model.matrix(~0 + factor(c(rep(contrast[[1]], length(x)), rep(contrast[[2]], length(y))), levels = contrast))
      fit <- limma::lmFit(matrix(response, nrow = 1), design)
      contrast_matrix <- matrix(c(1, -1), ncol = 1)
      fit <- limma::eBayes(limma::contrasts.fit(fit, contrast_matrix))
      p_value <- fit$p.value[[1]]
    }
    tibble::tibble(
      group = table$group[[1]],
      program = table$program[[1]],
      comparison = comparison,
      estimate = estimate,
      mean_1 = mean(x, na.rm = TRUE),
      mean_2 = mean(y, na.rm = TRUE),
      n_1 = length(x),
      n_2 = length(y),
      p_value = p_value
    )
  }))
  tests$adjusted_p_value <- stats::p.adjust(tests$p_value, method = "BH")
  store_name <- store_name %||% paste0(score_name, "_", condition_by)
  result <- list(
    schema_version = "1.0", analysis_type = "program_comparison", name = store_name,
    method = method, backend = method,
    input = list(score_name = score_name, condition_by = condition_by, sample_by = sample_by, group_by = group_by),
    parameters = list(contrast = contrast),
    tables = list(primary = tests, sample_scores = aggregated),
    embeddings = list(), graphs = list(), models = list(),
    diagnostics = list(inferential_unit = if (is_null(sample_by)) "cell" else "sample"),
    warnings = if (is_null(sample_by)) "Cells were used as exploratory units because sample_by was not supplied." else character(),
    provenance = .sn_analysis_provenance()
  )
  object <- sn_store_result(object, "program_comparison", store_name, result)
  if (isTRUE(return_object)) object else sn_get_result(object, "program_comparison", store_name)
}
