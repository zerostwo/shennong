.sn_standardize_grn_edges <- function(x, method) {
  if (!is.data.frame(x)) {
    stop("GRN edges must be provided as a data frame.", call. = FALSE)
  }
  x <- tibble::as_tibble(x)
  find_column <- function(candidates) {
    hit <- intersect(candidates, names(x))
    if (length(hit) == 0L) NULL else hit[[1]]
  }
  source_column <- find_column(c("source", "regulator", "tf", "regulatoryGene", "regulatory_gene"))
  target_column <- find_column(c("target", "target_gene", "gene"))
  weight_column <- find_column(c("weight", "importance", "score"))
  if (is_null(source_column) || is_null(target_column) || is_null(weight_column)) {
    stop(
      "GRN edges require source/regulator, target, and weight/importance columns.",
      call. = FALSE
    )
  }
  out <- tibble::tibble(
    source = as.character(x[[source_column]]),
    target = as.character(x[[target_column]]),
    weight = suppressWarnings(as.numeric(x[[weight_column]])),
    method = as.character(if ("method" %in% names(x)) x$method else method)
  )
  out <- out[
    !is.na(out$source) & nzchar(out$source) &
      !is.na(out$target) & nzchar(out$target) & is.finite(out$weight),
    , drop = FALSE
  ]
  if (nrow(out) == 0L) stop("GRN edges contain no complete finite interactions.", call. = FALSE)
  out <- out[order(abs(out$weight), decreasing = TRUE), , drop = FALSE]
  out[!duplicated(paste(out$source, out$target, sep = "\r")), , drop = FALSE]
}

.sn_grn_regulons <- function(edges, top_targets = 50L, min_weight = 0) {
  edges <- edges[abs(edges$weight) >= min_weight, , drop = FALSE]
  groups <- split(edges, edges$source)
  regulons <- dplyr::bind_rows(lapply(names(groups), function(regulator) {
    current <- groups[[regulator]]
    current <- utils::head(current[order(abs(current$weight), decreasing = TRUE), , drop = FALSE], as.integer(top_targets))
    tibble::tibble(
      regulon = regulator,
      regulator = regulator,
      target = current$target,
      weight = current$weight,
      rank = seq_len(nrow(current))
    )
  }))
  if (nrow(regulons) == 0L) stop("No regulon retained the requested target filters.", call. = FALSE)
  regulons
}

.sn_grn_activity_table <- function(activity, object, regulons, expression, method, control) {
  if (!is_null(activity)) {
    activity <- tibble::as_tibble(activity)
    cell_column <- intersect(c("cell", "entity", "cell_id"), names(activity))[[1]] %||% NULL
    regulon_column <- intersect(c("regulon", "source", "tf", "program"), names(activity))[[1]] %||% NULL
    score_column <- intersect(c("score", "activity", "auc", "value"), names(activity))[[1]] %||% NULL
    if (is_null(cell_column) || is_null(regulon_column) || is_null(score_column)) {
      stop("GRN activity requires cell, regulon, and score/activity columns.", call. = FALSE)
    }
    out <- tibble::tibble(
      cell = as.character(activity[[cell_column]]),
      regulon = as.character(activity[[regulon_column]]),
      score = suppressWarnings(as.numeric(activity[[score_column]])),
      method = as.character(if ("method" %in% names(activity)) activity$method else method)
    )
    out <- out[out$cell %in% colnames(object) & is.finite(out$score), , drop = FALSE]
    if (nrow(out) == 0L) stop("GRN activity contains no scores for cells in `object`.", call. = FALSE)
    return(out)
  }

  signatures <- split(regulons$target, regulons$regulon)
  matched <- .sn_match_program_features(signatures, rownames(expression$matrix), min_genes = control$min_targets %||% 1L)
  score_method <- match.arg(control$activity_method %||% "mean", c("mean", "ucell"))
  scores <- if (identical(score_method, "ucell")) {
    .sn_score_programs_ucell(expression$matrix, matched$signatures, control$activity_control %||% list())
  } else {
    .sn_score_programs_mean(expression$matrix, matched$signatures)
  }
  dplyr::bind_rows(lapply(rownames(scores), function(regulon) {
    cells <- colnames(scores)
    values <- as.numeric(scores[regulon, ])
    tibble::tibble(
      cell = cells, regulon = regulon,
      score = values, method = paste(method, score_method, sep = ":")
    )
  }))
}

.sn_grn_specificity <- function(activity, object, group_by) {
  if (is_null(group_by)) {
    return(tibble::tibble(
      regulon = unique(activity$regulon), group = "all",
      mean_activity = vapply(unique(activity$regulon), function(regulon) {
        mean(activity$score[activity$regulon == regulon], na.rm = TRUE)
      }, numeric(1)),
      specificity_score = NA_real_, metric = "not_computed"
    ))
  }
  if (!group_by %in% colnames(object[[]])) stop("`group_by` was not found in object metadata.", call. = FALSE)
  groups <- stats::setNames(as.character(object[[group_by, drop = TRUE]]), colnames(object))
  activity$group <- unname(groups[activity$cell])
  activity <- activity[!is.na(activity$group) & nzchar(activity$group), , drop = FALSE]
  dplyr::bind_rows(lapply(split(activity, activity$regulon), function(current) {
    overall <- mean(current$score, na.rm = TRUE)
    total_ss <- sum((current$score - overall) ^ 2, na.rm = TRUE)
    group_means <- tapply(current$score, current$group, mean, na.rm = TRUE)
    group_sizes <- table(current$group)
    between_ss <- sum(as.numeric(group_sizes[names(group_means)]) * (group_means - overall) ^ 2)
    specificity <- if (is.finite(total_ss) && total_ss > 0) between_ss / total_ss else 0
    tibble::tibble(
      regulon = current$regulon[[1]], group = names(group_means),
      mean_activity = as.numeric(group_means), specificity_score = specificity,
      metric = "eta_squared"
    )
  }))
}

.sn_run_grn_genie3 <- function(object, regulators, assay, layer, control) {
  check_installed("GENIE3", reason = "to infer a GRN with the GENIE3 backend.")
  expression <- .sn_annotation_expression(object, assay = assay, layer = layer)
  matrix <- expression$matrix
  means <- Matrix::rowMeans(matrix)
  variances <- Matrix::rowMeans(matrix ^ 2) - means ^ 2
  features <- names(utils::head(sort(variances, decreasing = TRUE), control$n_features %||% min(2000L, nrow(matrix))))
  regulators <- intersect(regulators %||% features, features)
  if (length(regulators) == 0L) stop("No requested `regulators` were found among retained expression features.", call. = FALSE)
  current <- matrix[features, , drop = FALSE]
  if (length(current) > (control$max_dense_values %||% 5e6)) {
    stop("GENIE3 input exceeds `backend_control$max_dense_values`; reduce features/cells.", call. = FALSE)
  }
  set.seed(as.integer(control$seed %||% 717L))
  weight_matrix <- GENIE3::GENIE3(
    exprMatrix = as.matrix(current), regulators = regulators,
    nTrees = as.integer(control$n_trees %||% 1000L),
    nCores = as.integer(control$n_cores %||% 1L), verbose = FALSE
  )
  edges <- GENIE3::getLinkList(weight_matrix, reportMax = as.integer(control$max_edges %||% 10000L))
  list(edges = edges, expression = expression, model = weight_matrix)
}

#' Infer and summarize a gene regulatory network
#'
#' @param object A Seurat object.
#' @param method GRN backend. GENIE3 runs in R. SCENIC, pySCENIC, and GRNBoost2
#'   accept explicit runner/result adapters so their external motif databases and
#'   runtimes remain visible.
#' @param name Stored result name.
#' @param assay,layer Expression assay and layer.
#' @param regulators Optional regulator genes. Strongly recommended for GENIE3.
#' @param group_by Optional metadata column used to quantify regulon specificity.
#' @param backend_control Backend controls or an external `runner`/`result`.
#' @param return_object Return the modified object or unified result.
#' @return A Seurat object or unified GRN result with edge, regulon, activity,
#'   and specificity tables.
#' @examples
#' \dontrun{
#' object <- sn_run_grn(object, method = "genie3", regulators = c("STAT1", "IRF1"),
#'                      group_by = "cell_type")
#' grn <- sn_get_result(object, "grn", "grn_genie3")
#' grn$tables$regulons
#' }
#' @export
sn_run_grn <- function(object,
                       method = c("genie3", "pyscenic", "scenic", "grnboost2"),
                       name = NULL,
                       assay = NULL,
                       layer = "data",
                       regulators = NULL,
                       group_by = NULL,
                       backend_control = list(),
                       return_object = TRUE) {
  .sn_validate_result_object(object)
  method <- match.arg(method)
  name <- name %||% paste0("grn_", method)
  if (is.function(backend_control$runner)) {
    output <- backend_control$runner(
      object = object, method = method, assay = assay, layer = layer,
      regulators = regulators, group_by = group_by, backend_control = backend_control
    )
  } else if (!is_null(backend_control$result)) {
    output <- backend_control$result
  } else if (identical(method, "genie3")) {
    output <- .sn_run_grn_genie3(object, regulators, assay, layer, backend_control)
  } else {
    stop("The ", method, " adapter requires `backend_control$runner` or `backend_control$result`.", call. = FALSE)
  }
  if (!is.list(output)) stop("GRN backend output must be a list.", call. = FALSE)
  edges <- .sn_standardize_grn_edges(output$edges %||% output$network %||% output$links, method)
  regulons <- .sn_grn_regulons(
    edges, top_targets = backend_control$top_targets %||% 50L,
    min_weight = backend_control$min_weight %||% 0
  )
  expression <- output$expression %||% .sn_annotation_expression(object, assay = assay, layer = layer)
  activity <- .sn_grn_activity_table(output$activity %||% output$auc, object, regulons, expression, method, backend_control)
  specificity <- if (!is_null(output$specificity)) {
    tibble::as_tibble(output$specificity)
  } else {
    .sn_grn_specificity(activity, object, group_by)
  }

  prefix <- gsub("[^[:alnum:]_]+", "_", name)
  metadata <- data.frame(row.names = colnames(object))
  for (regulon in unique(activity$regulon)) {
    values <- stats::setNames(activity$score[activity$regulon == regulon], activity$cell[activity$regulon == regulon])
    metadata[[paste(prefix, gsub("[^[:alnum:]_]+", "_", regulon), sep = "_")]] <- as.numeric(values[colnames(object)])
  }
  object <- SeuratObject::AddMetaData(object, metadata = metadata)
  result <- list(
    schema_version = "1.0", analysis_type = "grn", name = name,
    method = method, backend = method,
    input = list(
      assay = expression$assay %||% assay %||% SeuratObject::DefaultAssay(object),
      layer = expression$layer %||% layer, regulators = regulators,
      group_by = group_by, cells = ncol(object)
    ),
    parameters = list(
      top_targets = backend_control$top_targets %||% 50L,
      activity_method = backend_control$activity_method %||% "mean"
    ),
    tables = list(
      primary = edges, edges = edges, regulons = regulons,
      activity = activity, specificity = specificity
    ),
    embeddings = list(), graphs = list(network = edges),
    models = list(backend = output$model %||% NULL),
    diagnostics = list(
      edges = nrow(edges), regulons = length(unique(regulons$regulon)),
      activity_cells = length(unique(activity$cell)),
      specificity_metric = unique(as.character(specificity$metric %||% NA_character_))
    ),
    warnings = as.character(output$warnings %||% character()),
    provenance = .sn_analysis_provenance(random_seed = backend_control$seed %||% NA_integer_)
  )
  sn_validate_result(result)
  object <- sn_store_result(object, "grn", name, result)
  object <- .sn_log_seurat_command(object, assay = result$input$assay, name = "sn_run_grn")
  if (isTRUE(return_object)) object else sn_get_result(object, "grn", name)
}
