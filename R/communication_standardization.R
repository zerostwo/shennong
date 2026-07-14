.sn_communication_column <- function(table, candidates) {
  hits <- intersect(candidates, names(table))
  if (length(hits) == 0L) NULL else hits[[1]]
}

.sn_empty_communication_table <- function() {
  tibble::tibble(
    source = character(),
    target = character(),
    ligand = character(),
    receptor = character(),
    score = double(),
    p_value = double(),
    q_value = double(),
    rank = double(),
    method = character(),
    condition = character(),
    sample = character(),
    pathway = character(),
    target_genes = character(),
    evidence_source = character(),
    spatial_distance = double()
  )
}

.sn_communication_values <- function(table, candidates, default = NA) {
  column <- .sn_communication_column(table, candidates)
  if (is_null(column)) rep(default, nrow(table)) else table[[column]]
}

.sn_standardize_communication <- function(table,
                                          method,
                                          sender = NULL,
                                          receiver = NULL,
                                          condition = NULL,
                                          sample = NULL,
                                          evidence_source = method) {
  table <- tibble::as_tibble(table)
  if (nrow(table) == 0L) return(.sn_empty_communication_table())
  source <- .sn_communication_values(table, c("source", "sender", "source_cell", "cell_type1"), default = sender %||% NA_character_)
  target <- .sn_communication_values(table, c("target", "receiver", "target_cell", "cell_type2"), default = receiver %||% NA_character_)
  ligand <- .sn_communication_values(table, c("ligand", "ligand_complex", "test_ligand", "gene_a", "partner_a"), default = NA_character_)
  receptor <- .sn_communication_values(table, c("receptor", "receptor_complex", "gene_b", "partner_b"), default = NA_character_)
  pathway <- .sn_communication_values(table, c("pathway", "pathway_name", "pathway_name_1"), default = NA_character_)
  target_genes <- .sn_communication_values(table, c("target_genes", "targets", "geneset"), default = NA_character_)
  p_value <- suppressWarnings(as.numeric(.sn_communication_values(table, c("p_value", "pval", "p.value", "PValue", "pvalue"))))
  q_value <- suppressWarnings(as.numeric(.sn_communication_values(table, c("q_value", "p_val_adj", "FDR", "adjusted_p_value"))))
  score_column <- .sn_communication_column(table, c(
    "score", "prob", "magnitude", "lr_means", "mean", "pearson", "aupr_corrected",
    "scaled_weight", "weight", "interaction_score", "prioritization_score"
  ))
  rank_column <- .sn_communication_column(table, c("aggregate_rank", "magnitude_rank", "rank"))
  score <- if (!is_null(score_column)) {
    suppressWarnings(as.numeric(table[[score_column]]))
  } else if (!is_null(rank_column)) {
    1 - suppressWarnings(as.numeric(table[[rank_column]]))
  } else if (any(is.finite(p_value))) {
    -log10(pmax(p_value, .Machine$double.xmin))
  } else {
    rep(NA_real_, nrow(table))
  }
  if (all(!is.finite(q_value)) && any(is.finite(p_value))) q_value <- stats::p.adjust(p_value, method = "BH")
  interaction_rank <- rank(-score, ties.method = "average", na.last = "keep")

  standardized <- tibble::tibble(
    source = as.character(source),
    target = as.character(target),
    ligand = as.character(ligand),
    receptor = as.character(receptor),
    score = score,
    p_value = p_value,
    q_value = q_value,
    rank = as.numeric(interaction_rank),
    method = method,
    condition = as.character(.sn_communication_values(table, c("condition", "group", "contrast"), default = condition %||% NA_character_)),
    sample = as.character(.sn_communication_values(table, c("sample"), default = sample %||% NA_character_)),
    pathway = as.character(pathway),
    target_genes = as.character(target_genes),
    evidence_source = evidence_source,
    spatial_distance = suppressWarnings(as.numeric(.sn_communication_values(table, c("spatial_distance"))))
  )
  raw_extra <- table[, setdiff(names(table), names(standardized)), drop = FALSE]
  dplyr::bind_cols(standardized, raw_extra)
}

.sn_communication_edge_key <- function(table) {
  keys <- c("source", "target", "ligand", "receptor")
  values <- lapply(table[keys], function(value) {
    value <- as.character(value)
    value[is.na(value) | !nzchar(value)] <- "<missing>"
    value
  })
  do.call(paste, c(values, sep = "\r"))
}

.sn_finite_min <- function(value) {
  value <- value[is.finite(value)]
  if (length(value) == 0L) NA_real_ else min(value)
}

.sn_finite_max <- function(value) {
  value <- value[is.finite(value)]
  if (length(value) == 0L) NA_real_ else max(value)
}

.sn_communication_consensus <- function(table) {
  if (nrow(table) == 0L) return(.sn_empty_communication_table())
  method_sizes <- table(table$method)
  table$.sn_rank_score <- vapply(seq_len(nrow(table)), function(index) {
    current_rank <- table$rank[[index]]
    current_size <- unname(method_sizes[[table$method[[index]]]])
    if (!is.finite(current_rank)) return(NA_real_)
    1 - (current_rank - 1) / max(1, current_size - 1)
  }, numeric(1))
  table$.sn_edge <- .sn_communication_edge_key(table)
  groups <- split(table, table$.sn_edge)
  consensus <- dplyr::bind_rows(lapply(groups, function(current) {
    methods <- unique(current$method)
    per_method <- tibble::tibble(
      method = methods,
      x = vapply(methods, function(method) {
        .sn_finite_max(current$.sn_rank_score[current$method == method])
      }, numeric(1))
    )
    normalized <- per_method$x
    normalized[!is.finite(normalized)] <- NA_real_
    tibble::tibble(
      source = current$source[[1]],
      target = current$target[[1]],
      ligand = current$ligand[[1]],
      receptor = current$receptor[[1]],
      score = mean(normalized, na.rm = TRUE),
      p_value = .sn_finite_min(current$p_value),
      q_value = .sn_finite_min(current$q_value),
      rank = NA_real_,
      method = "consensus",
      condition = paste(unique(stats::na.omit(current$condition)), collapse = ";"),
      sample = paste(unique(stats::na.omit(current$sample)), collapse = ";"),
      pathway = paste(unique(stats::na.omit(current$pathway)), collapse = ";"),
      target_genes = paste(unique(stats::na.omit(current$target_genes)), collapse = ";"),
      evidence_source = paste(sort(unique(current$method)), collapse = ";"),
      spatial_distance = if (all(!is.finite(current$spatial_distance))) NA_real_ else mean(current$spatial_distance, na.rm = TRUE),
      n_methods = nrow(per_method),
      method_concordance = if (sum(is.finite(normalized)) < 2L) NA_real_ else 1 / (1 + stats::sd(normalized, na.rm = TRUE))
    )
  }))
  consensus$p_value[!is.finite(consensus$p_value)] <- NA_real_
  consensus$q_value[!is.finite(consensus$q_value)] <- NA_real_
  consensus$score[is.nan(consensus$score)] <- NA_real_
  consensus$rank <- rank(-consensus$score, ties.method = "average", na.last = "keep")
  consensus[order(consensus$rank), , drop = FALSE]
}

.sn_communication_concordance <- function(table) {
  methods <- unique(table$method)
  if (length(methods) < 2L) return(tibble::tibble())
  table$.sn_edge <- .sn_communication_edge_key(table)
  ranks <- lapply(split(table, table$method), function(current) {
    edges <- unique(current$.sn_edge)
    aggregated <- tibble::tibble(
      edge = edges,
      x = vapply(edges, function(edge) {
        .sn_finite_min(current$rank[current$.sn_edge == edge])
      }, numeric(1))
    )
    stats::setNames(
      aggregated$x,
      aggregated$edge
    )
  })
  pairs <- utils::combn(methods, 2, simplify = FALSE)
  dplyr::bind_rows(lapply(pairs, function(pair) {
    shared <- intersect(names(ranks[[pair[[1]]]]), names(ranks[[pair[[2]]]]))
    tibble::tibble(
      method_1 = pair[[1]],
      method_2 = pair[[2]],
      shared_edges = length(shared),
      rank_correlation = if (length(shared) < 3L) NA_real_ else stats::cor(
        ranks[[pair[[1]]]][shared], ranks[[pair[[2]]]][shared], method = "spearman", use = "complete.obs"
      )
    )
  }))
}

.sn_lr_components <- function(value) {
  value <- as.character(value)
  value <- sub("^(simple|complex):", "", value)
  unique(unlist(strsplit(value, "[+_&]", perl = TRUE), use.names = FALSE))
}

.sn_communication_sample_evidence <- function(object,
                                               interactions,
                                               group_by,
                                               sample_by,
                                               condition_by,
                                               assay,
                                               layer) {
  if (is_null(sample_by)) return(tibble::tibble())
  metadata <- object[[]]
  required <- c(group_by, sample_by, condition_by)
  missing <- setdiff(required, colnames(metadata))
  if (length(missing) > 0L) stop("Missing communication sample metadata: ", paste(missing, collapse = ", "), ".", call. = FALSE)
  if (!is_null(condition_by)) .sn_validate_constant_within_sample(metadata, sample_col = sample_by, group_col = condition_by)
  expression <- .sn_annotation_expression(object, assay = assay, layer = layer)
  samples <- unique(as.character(metadata[[sample_by]]))
  interactions <- unique(interactions[, c("source", "target", "ligand", "receptor"), drop = FALSE])
  rows <- list()
  counter <- 0L
  for (sample in samples) {
    sample_cells <- rownames(metadata)[as.character(metadata[[sample_by]]) == sample]
    condition <- if (is_null(condition_by)) NA_character_ else unique(as.character(metadata[[condition_by]][metadata[[sample_by]] == sample]))[[1]]
    groups <- split(sample_cells, as.character(metadata[sample_cells, group_by]))
    for (index in seq_len(nrow(interactions))) {
      current <- interactions[index, , drop = FALSE]
      source <- current$source[[1]]
      target <- current$target[[1]]
      ligand <- current$ligand[[1]]
      receptor <- current$receptor[[1]]
      source_cells <- groups[[source]] %||% character()
      target_cells <- groups[[target]] %||% character()
      ligand_features <- intersect(.sn_lr_components(ligand), rownames(expression$matrix))
      receptor_features <- intersect(.sn_lr_components(receptor), rownames(expression$matrix))
      ligand_expression <- if (length(source_cells) == 0L || length(ligand_features) == 0L) NA_real_ else mean(Matrix::rowMeans(expression$matrix[ligand_features, source_cells, drop = FALSE]))
      receptor_expression <- if (length(target_cells) == 0L || length(receptor_features) == 0L) NA_real_ else mean(Matrix::rowMeans(expression$matrix[receptor_features, target_cells, drop = FALSE]))
      counter <- counter + 1L
      rows[[counter]] <- tibble::tibble(
        source = source,
        target = target,
        ligand = ligand,
        receptor = receptor,
        sample = sample,
        condition = condition,
        ligand_expression = ligand_expression,
        receptor_expression = receptor_expression,
        score = sqrt(pmax(ligand_expression, 0) * pmax(receptor_expression, 0))
      )
    }
  }
  dplyr::bind_rows(rows)
}

.sn_compare_communication_samples <- function(sample_evidence, contrast = NULL) {
  if (nrow(sample_evidence) == 0L || all(is.na(sample_evidence$condition))) return(tibble::tibble())
  levels <- unique(sample_evidence$condition[!is.na(sample_evidence$condition)])
  if (is_null(contrast)) {
    if (length(levels) != 2L) return(tibble::tibble())
    contrast <- levels
  }
  if (length(contrast) != 2L || any(!contrast %in% levels)) {
    stop("`contrast` must contain two conditions present in the sample-level evidence.", call. = FALSE)
  }
  keys <- .sn_communication_edge_key(sample_evidence)
  groups <- split(sample_evidence, keys)
  table <- dplyr::bind_rows(lapply(groups, function(current) {
    case <- current$score[current$condition == contrast[[1]]]
    control <- current$score[current$condition == contrast[[2]]]
    p <- if (sum(is.finite(case)) >= 2L && sum(is.finite(control)) >= 2L) stats::wilcox.test(case, control, exact = FALSE)$p.value else NA_real_
    mean_case <- if (any(is.finite(case))) mean(case, na.rm = TRUE) else NA_real_
    mean_control <- if (any(is.finite(control))) mean(control, na.rm = TRUE) else NA_real_
    tibble::tibble(
      source = current$source[[1]], target = current$target[[1]],
      ligand = current$ligand[[1]], receptor = current$receptor[[1]],
      mean_case = mean_case, mean_control = mean_control,
      estimate = mean_case - mean_control,
      p_value = p, n_case = sum(is.finite(case)), n_control = sum(is.finite(control)),
      contrast_case = contrast[[1]], contrast_control = contrast[[2]]
    )
  }))
  table$adjusted_p_value <- stats::p.adjust(table$p_value, method = "BH")
  table
}

.sn_nichenet_target_links <- function(ligand_target_matrix, ligands, geneset, top_n = 100L) {
  if (is_null(ligand_target_matrix)) return(tibble::tibble())
  ligands <- intersect(ligands, colnames(ligand_target_matrix))
  targets <- intersect(geneset, rownames(ligand_target_matrix))
  links <- dplyr::bind_rows(lapply(ligands, function(ligand) {
    weights <- ligand_target_matrix[targets, ligand]
    keep <- order(weights, decreasing = TRUE)[seq_len(min(length(weights), top_n))]
    tibble::tibble(ligand = ligand, target_gene = targets[keep], weight = as.numeric(weights[keep]))
  }))
  links[is.finite(links$weight), , drop = FALSE]
}

.sn_read_cellphonedb_output <- function(output_dir) {
  p_path <- file.path(output_dir, "pvalues.txt")
  m_path <- file.path(output_dir, "means.txt")
  if (!file.exists(p_path) || !file.exists(m_path)) stop("CellPhoneDB output is missing `pvalues.txt` or `means.txt`.", call. = FALSE)
  p_values <- utils::read.delim(p_path, check.names = FALSE)
  means <- utils::read.delim(m_path, check.names = FALSE)
  id_candidates <- intersect(c("id_cp_interaction", "interacting_pair"), names(p_values))
  if (length(id_candidates) == 0L) stop("CellPhoneDB output lacks an interaction identifier.", call. = FALSE)
  id <- id_candidates[[1]]
  pair_columns <- intersect(grep("\\|", names(p_values), value = TRUE), names(means))
  if (length(pair_columns) == 0L) stop("CellPhoneDB output contains no shared sender-receiver columns.", call. = FALSE)
  metadata_columns <- setdiff(names(p_values), pair_columns)
  ligand_column <- .sn_communication_column(p_values, c("gene_a", "partner_a"))
  receptor_column <- .sn_communication_column(p_values, c("gene_b", "partner_b"))
  if (is_null(ligand_column) || is_null(receptor_column)) stop("CellPhoneDB output lacks ligand or receptor columns.", call. = FALSE)
  rows <- lapply(pair_columns, function(pair) {
    pair_parts <- strsplit(pair, "\\|", fixed = FALSE)[[1]]
    if (length(pair_parts) != 2L) stop("CellPhoneDB sender-receiver column is malformed: ", pair, call. = FALSE)
    tibble::tibble(
      interaction_id = p_values[[id]],
      source = pair_parts[[1]], target = pair_parts[[2]],
      ligand = sub("^(simple|complex):", "", as.character(p_values[[ligand_column]])),
      receptor = sub("^(simple|complex):", "", as.character(p_values[[receptor_column]])),
      score = as.numeric(means[[pair]]), p_value = as.numeric(p_values[[pair]])
    )
  })
  list(table = dplyr::bind_rows(rows), raw = list(pvalues = p_values[, metadata_columns, drop = FALSE], means = means))
}
