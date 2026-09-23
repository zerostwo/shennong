.sn_enrichment_cache_env <- local({
  env <- new.env(parent = emptyenv())
  env$msigdb_terms <- new.env(parent = emptyenv())
  env$symbol_to_entrez <- new.env(parent = emptyenv())
  env
})

.sn_enrichment_with_rng_preserved <- function(expr) {
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  force(expr)
}

.sn_enrich_cache_key <- function(...) {
  paste(vapply(list(...), as.character, character(1)), collapse = "::")
}

.sn_enrich_parse_msigdb_database <- function(database,
                                             collection = NULL,
                                             subcollection = NULL) {
  database <- toupper(database)
  collection <- if (!is_null(collection)) toupper(collection) else NULL
  subcollection <- if (!is_null(subcollection)) toupper(subcollection) else NULL

  if (identical(database, "MSIGDB")) {
    if (is_null(collection) || !nzchar(collection)) {
      stop("When `database = \"MSIGDB\"`, `collection` must be supplied.", call. = FALSE)
    }
    return(list(collection = collection, subcollection = subcollection))
  }

  parts <- strsplit(database, ":", fixed = TRUE)[[1]]
  parsed_collection <- parts[[1]]
  known_collections <- c("H", paste0("C", 1:8))
  if (!parsed_collection %in% known_collections) {
    return(NULL)
  }

  parsed_subcollection <- if (length(parts) > 1) {
    paste(parts[-1], collapse = ":")
  } else {
    subcollection
  }

  list(
    collection = parsed_collection,
    subcollection = parsed_subcollection
  )
}

.sn_enrich_get_msigdb_terms <- function(species, collection, subcollection = NULL) {
  cache_key <- .sn_enrich_cache_key(
    species,
    as.character(utils::packageVersion("msigdbr")),
    toupper(collection),
    toupper(subcollection %||% "")
  )
  if (exists(cache_key, envir = .sn_enrichment_cache_env$msigdb_terms, inherits = FALSE)) {
    return(get(cache_key, envir = .sn_enrichment_cache_env$msigdb_terms, inherits = FALSE))
  }

  msig_args <- list(
    species = species,
    collection = collection
  )
  if (!is_null(subcollection) && nzchar(subcollection)) {
    msig_args$subcollection <- subcollection
  }

  msigdbr_tbl <- do.call(msigdbr::msigdbr, msig_args)
  if (nrow(msigdbr_tbl) == 0) {
    sub_msg <- if (!is_null(subcollection) && nzchar(subcollection)) paste0(":", subcollection) else ""
    stop(
      glue("No MSigDB terms were returned for collection '{collection}{sub_msg}' and species '{species}'."),
      call. = FALSE
    )
  }

  out <- msigdbr_tbl |>
    dplyr::transmute(
      term = .data$gs_name,
      description = .data$gs_description,
      gene = .data$gene_symbol
    ) |>
    dplyr::distinct()
  assign(cache_key, out, envir = .sn_enrichment_cache_env$msigdb_terms)
  out
}

.sn_enrich_normalize_database_labels <- function(database) {
  unique(toupper(as.character(database)))
}

.sn_enrich_output_label <- function(database) {
  gsub("[^[:alnum:]_.-]+", "_", toupper(database))
}

.sn_enrich_id_label <- function(value) {
  gsub("[^[:alnum:]_.-]+", "_", as.character(value))
}

.sn_enrich_parse_mapping <- function(mapping) {
  if (is_null(mapping)) {
    return(NULL)
  }
  if (!inherits(mapping, "formula") || length(mapping) != 3L) {
    stop(
      "`mapping` must be a two-sided formula such as `gene ~ cluster`, ",
      "`gene ~ score`, or `gene ~ score | group`.",
      call. = FALSE
    )
  }

  gene_vars <- all.vars(mapping[[2L]])
  if (length(gene_vars) != 1L) {
    stop("The left-hand side of `mapping` must contain exactly one gene column.", call. = FALSE)
  }

  rhs <- mapping[[3L]]
  grouped_gsea <- is.call(rhs) && identical(rhs[[1L]], as.name("|"))
  if (grouped_gsea) {
    score_vars <- all.vars(rhs[[2L]])
    group_vars <- all.vars(rhs[[3L]])
    if (length(score_vars) != 1L || length(group_vars) == 0L) {
      stop(
        "Grouped GSEA mapping must use `gene ~ score | group` with one score ",
        "column and one or more grouping columns.",
        call. = FALSE
      )
    }
    return(list(
      gene_col = gene_vars[[1L]],
      value_col = score_vars[[1L]],
      group_cols = group_vars,
      grouped_gsea = TRUE,
      formula = mapping
    ))
  }

  rhs_vars <- all.vars(rhs)
  if (length(rhs_vars) == 0L) {
    stop("The right-hand side of `mapping` must contain at least one column.", call. = FALSE)
  }
  list(
    gene_col = gene_vars[[1L]],
    value_col = if (length(rhs_vars) == 1L) rhs_vars[[1L]] else NULL,
    group_cols = rhs_vars,
    grouped_gsea = FALSE,
    formula = mapping
  )
}

.sn_enrich_resolve_mapping <- function(mapping = NULL,
                                       gene_clusters = NULL) {
  if (!is_null(mapping) && !is_null(gene_clusters)) {
    stop("Supply only one of `mapping` and compatibility alias `gene_clusters`.", call. = FALSE)
  }
  .sn_enrich_parse_mapping(mapping %||% gene_clusters)
}

.sn_enrich_resolve_input <- function(x,
                                     source_de_result_id = NULL) {
  if (inherits(x, "Seurat")) {
    source_de_result_id <- .sn_resolve_stored_result_id(
      object = x,
      type = "de",
      result_id = source_de_result_id,
      preferred_analysis = "markers"
    )
    de_result <- sn_get_result(x, type = "de", result_id = source_de_result_id)
    return(list(
      input = de_result$tables$primary,
      object = x,
      source_de_result_id = source_de_result_id,
      de_result = de_result
    ))
  }

  list(
    input = x,
    object = NULL,
    source_de_result_id = source_de_result_id,
    de_result = NULL
  )
}

.sn_enrich_stored_de_universe <- function(de_result, object) {
  if (!is.list(de_result)) {
    stop("The stored DE result must be a list.", call. = FALSE)
  }
  de_input <- de_result$input
  if (isTRUE(de_input$tested_features_source %in% c("assay_layer", "requested_features", "assay_layer_full_test_then_result_filter"))) {
    stop("This stored DE result records candidate features rather than verified backend tests; rerun DE or supply `universe` explicitly.", call. = FALSE)
  }
  has_input_universe <- is.list(de_input) &&
    "tested_features" %in% names(de_input)
  has_top_level_universe <- "tested_features" %in% names(de_result)

  if (has_input_universe || has_top_level_universe) {
    tested_features <- if (has_input_universe) {
      de_input$tested_features
    } else {
      de_result$tested_features
    }
    if (!is.character(tested_features) || length(tested_features) == 0L ||
        anyNA(tested_features) || any(!nzchar(trimws(tested_features))) ||
        anyDuplicated(tested_features)) {
      stop(
        "The stored DE `tested_features` universe must be a non-empty character vector of distinct, non-missing feature names.",
        call. = FALSE
      )
    }
    return(list(
      universe = tested_features,
      source = "stored_de_tested_features"
    ))
  }

  universe_assay <- de_result$assay %||%
    de_result$input$assay %||%
    SeuratObject::DefaultAssay(object)
  if (!universe_assay %in% names(object@assays)) {
    stop(
      "The legacy stored DE assay '", universe_assay,
      "' is unavailable, so an ORA background cannot be reconstructed; supply `universe` explicitly.",
      call. = FALSE
    )
  }
  list(
    universe = rownames(object[[universe_assay]]),
    source = paste0("stored_de_assay_fallback:", universe_assay)
  )
}

.sn_enrich_comparison_backgrounds <- function(input, de_result, mapping = NULL) {
  backgrounds <- de_result$input$tested_features_by_comparison
  columns <- de_result$input$tested_feature_groups
  if (is.null(backgrounds) || length(columns) == 0L) return(NULL)
  if (!all(c("gene", columns) %in% names(backgrounds)) || !all(columns %in% names(input))) {
    stop("Stored DE comparison backgrounds are incomplete; rerun DE or supply `universe` explicitly.", call. = FALSE)
  }
  comparisons <- unique(as.data.frame(input[, columns, drop = FALSE]))
  if (!is.null(mapping)) {
    # A shared background cannot silently pool distinct hypothesis families.
    varying <- columns[vapply(comparisons, function(x) length(unique(x)) > 1L, logical(1))]
    mapped <- if (!is.null(mapping$formula)) all.vars(mapping$formula)[-1L] else mapping$value_col
    if (!all(varying %in% mapped)) {
      stop("The ORA mapping must retain all DE comparison columns: ",
           paste(varying, collapse = ", "), "; or supply `universe` explicitly.", call. = FALSE)
    }
  }
  key_columns <- unique(c(columns, if (!is.null(mapping)) mapped else character()))
  if (!all(key_columns %in% names(input))) stop("ORA grouping columns are absent from the DE table.", call. = FALSE)
  comparisons <- unique(as.data.frame(input[, key_columns, drop = FALSE]))
  lapply(seq_len(nrow(comparisons)), function(i) {
    key <- comparisons[i, , drop = FALSE]
    matches <- function(table, keys = columns) {
      Reduce(`&`, lapply(keys, function(column) as.character(table[[column]]) == as.character(key[[column]])))
    }
    universe <- unique(as.character(backgrounds$gene[matches(backgrounds)]))
    if (!length(universe)) {
      stop("A stored DE comparison has no verified background; rerun DE or supply `universe` explicitly.", call. = FALSE)
    }
    list(key = key, input = input[matches(input, key_columns), , drop = FALSE], universe = universe)
  })
}

.sn_enrich_run_comparisons <- function(run_one, database, comparisons, gene_col, grouped = FALSE) {
  results <- lapply(comparisons, function(comparison) {
    # Reuse the single-comparison backend with its exact hypothesis universe.
    run <- run_one
    environment(run) <- list2env(list(input = comparison$input,
      universe = comparison$universe, mapping = NULL), parent = environment(run_one))
    run(database)
  })
  if (length(results) == 1L && !grouped) {
    return(results[[1L]] %||% data.frame(ID = character(), Description = character()))
  }
  labels <- make.unique(vapply(comparisons, function(x) paste(x$key[1, ], collapse = "."), character(1)))
  tables <- lapply(seq_along(results), function(i) {
    table <- as.data.frame(.sn_as_enrichment_table(results[[i]]))
    if (ncol(table) == 0L) table <- data.frame(ID = character(), Description = character())
    table$Cluster <- rep(labels[[i]], nrow(table))
    for (column in names(comparisons[[i]]$key)) {
      table[[column]] <- rep(comparisons[[i]]$key[[column]], nrow(table))
    }
    table
  })
  table <- as.data.frame(dplyr::bind_rows(tables))
  is_go <- database %in% c("GO", "GOBP", "GOMF", "GOCC")
  methods::new("compareClusterResult", compareClusterResult = table,
    geneClusters = stats::setNames(lapply(seq_along(comparisons), function(i) {
      if (identical(database, "KEGG") && methods::is(results[[i]], "enrichResult")) return(results[[i]]@gene)
      unique(comparisons[[i]]$input[[gene_col]])
    }), labels),
    fun = if (is_go) "enrichGO" else if (identical(database, "KEGG")) "enrichKEGG" else "enricher",
    keytype = if (is_go) "SYMBOL" else "UNKNOWN", readable = FALSE, .call = match.call())
}

.sn_enrich_filter_stored_de_ora <- function(input,
                                             de_result,
                                             p_adjusted_cutoff = 0.05,
                                             direction = c("up", "down", "both"),
                                             logfc_threshold = 0) {
  direction <- match.arg(direction)
  input <- as.data.frame(input, check.names = FALSE)
  adjusted_candidates <- unique(c(
    de_result$p_col %||% character(),
    "adjusted_p_value", "p_val_adj", "padj", "FDR", "q_value"
  ))
  adjusted_candidates <- adjusted_candidates[
    !is.na(adjusted_candidates) & adjusted_candidates %in% colnames(input)
  ]
  # Ranking scores such as COSG specificity are not signed fold changes and
  # therefore cannot determine up/down direction for ORA.
  signed_effect_columns <- c("avg_log2FC", "avg_logFC", "log2FoldChange", "logFC", "effect")
  preferred_effect <- de_result$rank_col %||% character()
  preferred_effect <- preferred_effect[!is.na(preferred_effect) & preferred_effect %in% signed_effect_columns]
  effect_candidates <- unique(c(preferred_effect, signed_effect_columns))
  effect_candidates <- effect_candidates[
    !is.na(effect_candidates) & effect_candidates %in% colnames(input)
  ]
  if (length(adjusted_candidates) == 0L) {
    stop(
      "ORA from a stored DE result requires an adjusted-p-value column ",
      "(`p_val_adj`, `padj`, `FDR`, or `adjusted_p_value`).",
      call. = FALSE
    )
  }
  if (length(effect_candidates) == 0L) {
    stop(
      "ORA from a stored DE result requires a signed effect column so gene ",
      "direction can be selected explicitly.",
      call. = FALSE
    )
  }
  adjusted <- suppressWarnings(as.numeric(input[[adjusted_candidates[[1]]]]))
  effect <- suppressWarnings(as.numeric(input[[effect_candidates[[1]]]]))
  keep_direction <- switch(
    direction,
    up = effect > logfc_threshold,
    down = effect < -logfc_threshold,
    both = abs(effect) > logfc_threshold
  )
  keep <- is.finite(adjusted) & adjusted <= p_adjusted_cutoff &
    is.finite(effect) & keep_direction
  filtered <- input[keep, , drop = FALSE]
  if (nrow(filtered) == 0L) {
    stop(
      "No stored DE genes passed `de_p_adjusted_cutoff`, ",
      "`de_logfc_threshold`, and `de_direction`.",
      call. = FALSE
    )
  }
  attr(filtered, "sn_de_ora_selection") <- list(
    adjusted_p_column = adjusted_candidates[[1]],
    effect_column = effect_candidates[[1]],
    p_adjusted_cutoff = p_adjusted_cutoff,
    logfc_threshold = logfc_threshold,
    direction = direction,
    input_rows = nrow(input),
    retained_rows = nrow(filtered)
  )
  filtered
}

.sn_enrich_resolve_analysis <- function(input,
                                        mapping = NULL,
                                        analysis = NULL) {
  if (!is_null(analysis)) {
    analysis <- match.arg(analysis, c("ora", "gsea"))
    if (isTRUE(mapping$grouped_gsea) && identical(analysis, "ora")) {
      stop("A `gene ~ score | group` mapping requires `analysis = \"gsea\"`.", call. = FALSE)
    }
    return(analysis)
  }

  if (is.numeric(input) && !is.null(names(input))) {
    return("gsea")
  }

  if (is.character(input)) {
    return("ora")
  }

  if (is.data.frame(input)) {
    if (!is_null(mapping)) {
      if (isTRUE(mapping$grouped_gsea)) {
        return("gsea")
      }
      if (length(mapping$group_cols) > 1L) {
        return("ora")
      }
      value <- input[[mapping$value_col]]
      if (is.numeric(value)) {
        return("gsea")
      }
      return("ora")
    }
  }

  "ora"
}

.sn_enrich_resolve_gene_vector <- function(input, gene_col = "gene") {
  if (is.character(input)) {
    genes <- input
  } else if (is.data.frame(input)) {
    if (!gene_col %in% colnames(input)) {
      stop(glue("Column '{gene_col}' was not found in `x`."), call. = FALSE)
    }
    genes <- input[[gene_col]]
  } else {
    stop("ORA input must be a character vector or a data frame with a gene column.", call. = FALSE)
  }

  genes <- trimws(as.character(genes))
  invalid <- is.na(genes) | !nzchar(genes)
  if (any(invalid)) {
    stop("ORA gene identifiers cannot be missing or empty.", call. = FALSE)
  }
  unique(genes)
}

.sn_enrich_collapse_gene_list <- function(
    gene_list,
    duplicate_gene_method = c("error", "max_abs", "max", "mean")) {
  duplicate_gene_method <- match.arg(duplicate_gene_method)
  if (!is.numeric(gene_list) || is.null(names(gene_list))) {
    stop("GSEA input must be a named numeric vector.", call. = FALSE)
  }

  gene_ids <- trimws(as.character(names(gene_list)))
  if (anyNA(gene_ids) || any(!nzchar(gene_ids))) {
    stop("GSEA gene identifiers cannot be missing or empty.", call. = FALSE)
  }
  if (anyNA(gene_list) || any(!is.finite(gene_list))) {
    stop("GSEA ranking statistics must all be finite and non-missing.", call. = FALSE)
  }
  names(gene_list) <- gene_ids

  duplicated_ids <- unique(gene_ids[duplicated(gene_ids) | duplicated(gene_ids, fromLast = TRUE)])
  if (length(duplicated_ids) > 0L) {
    if (identical(duplicate_gene_method, "error")) {
      stop(
        "GSEA gene identifiers must be unique. Duplicated identifier(s): ",
        paste(utils::head(duplicated_ids, 5L), collapse = ", "),
        if (length(duplicated_ids) > 5L) "..." else "",
        ". Set `duplicate_gene_method` explicitly to collapse duplicates.",
        call. = FALSE
      )
    }
    grouped <- split(seq_along(gene_list), gene_ids)
    collapsed <- vapply(grouped, function(index) {
      values <- gene_list[index]
      switch(
        duplicate_gene_method,
        max_abs = values[[which.max(abs(values))]],
        max = max(values),
        mean = mean(values)
      )
    }, numeric(1))
    gene_list <- stats::setNames(as.numeric(collapsed), names(collapsed))
  }

  sort(gene_list, decreasing = TRUE)
}

.sn_enrich_resolve_gene_list <- function(input,
                                         mapping = NULL,
                                         duplicate_gene_method = c("error", "max_abs", "max", "mean")) {
  duplicate_gene_method <- match.arg(duplicate_gene_method)
  if (is.numeric(input) && !is.null(names(input))) {
    gene_list <- stats::setNames(as.numeric(input), names(input))
    return(.sn_enrich_collapse_gene_list(gene_list, duplicate_gene_method))
  }

  if (!is.data.frame(input)) {
    stop("For GSEA, `x` must be a named numeric vector or a data frame containing gene and ranking columns.", call. = FALSE)
  }

  if (is_null(mapping)) {
    stop(
      "For GSEA with data-frame input, `mapping` must be supplied as ",
      "`gene ~ ranking_column` or `gene ~ ranking_column | group`.",
      call. = FALSE
    )
  }

  gene_col <- mapping$gene_col
  score_col <- mapping$value_col
  if (!gene_col %in% colnames(input)) {
    stop(glue("Column '{gene_col}' was not found in `x`."), call. = FALSE)
  }
  if (!score_col %in% colnames(input)) {
    stop(glue("Column '{score_col}' was not found in `x`."), call. = FALSE)
  }
  if (!is.numeric(input[[score_col]])) {
    stop(glue("Column '{score_col}' must be numeric for GSEA."), call. = FALSE)
  }

  gene_list <- input[[score_col]]
  names(gene_list) <- as.character(input[[gene_col]])
  .sn_enrich_collapse_gene_list(gene_list, duplicate_gene_method)
}

.sn_enrich_resolve_grouped_gene_lists <- function(
    input,
    mapping,
    duplicate_gene_method = c("error", "max_abs", "max", "mean")) {
  duplicate_gene_method <- match.arg(duplicate_gene_method)
  if (!is.data.frame(input) || !isTRUE(mapping$grouped_gsea)) {
    stop("Grouped GSEA requires data-frame input and `gene ~ score | group` mapping.", call. = FALSE)
  }

  required <- c(mapping$gene_col, mapping$value_col, mapping$group_cols)
  missing <- setdiff(required, colnames(input))
  if (length(missing) > 0L) {
    stop("Column(s) not found in `x`: ", paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  invalid_group <- vapply(mapping$group_cols, function(column) {
    values <- input[[column]]
    anyNA(values) || any(!nzchar(trimws(as.character(values))))
  }, logical(1))
  if (any(invalid_group)) {
    stop(
      "Grouped GSEA grouping columns cannot contain missing or empty values: ",
      paste(mapping$group_cols[invalid_group], collapse = ", "), ".",
      call. = FALSE
    )
  }

  group_values <- lapply(
    mapping$group_cols,
    function(column) as.character(input[[column]])
  )
  groups <- do.call(
    interaction,
    c(group_values, list(drop = TRUE, lex.order = TRUE, sep = "."))
  )
  indices <- split(seq_len(nrow(input)), groups, drop = TRUE)
  gene_lists <- lapply(indices, function(index) {
    .sn_enrich_resolve_gene_list(
      input = input[index, , drop = FALSE],
      mapping = mapping,
      duplicate_gene_method = duplicate_gene_method
    )
  })
  if (length(gene_lists) == 0L) {
    stop("Grouped GSEA mapping produced no non-empty groups.", call. = FALSE)
  }
  gene_lists
}

.sn_enrich_muffle_empty_warning <- function(expr) {
  withCallingHandlers(
    expr,
    warning = function(w) {
      message <- conditionMessage(w)
      benign_patterns <- c(
        "No enrichment found"
      )
      if (any(vapply(benign_patterns, grepl, logical(1), x = message, fixed = TRUE))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

.sn_enrich_result_ids <- function(result_id,
                                  databases,
                                  source_de_result_id = NULL,
                                  analysis = NULL) {
  if (is_null(result_id)) {
    source_label <- source_de_result_id %||% "enrichment"
    analysis_label <- analysis %||% "analysis"
    generated <- paste(
      .sn_enrich_id_label(source_label),
      tolower(.sn_enrich_id_label(analysis_label)),
      .sn_enrich_output_label(databases),
      sep = "."
    )
    return(stats::setNames(generated, databases))
  }

  if (!is.character(result_id) || length(result_id) == 0L) {
    stop("`result_id` must be NULL or contain one or more result names.", call. = FALSE)
  }
  result_id <- vapply(result_id, .sn_validate_result_id, character(1))
  if (length(databases) == 1) {
    return(result_id[[1]])
  }

  if (length(result_id) == 1) {
    return(stats::setNames(
      paste(result_id[[1]], .sn_enrich_output_label(databases), sep = "."),
      databases
    ))
  }

  if (length(result_id) != length(databases)) {
    stop("`result_id` must have length 1 or match the length of `database`.", call. = FALSE)
  }

  stats::setNames(as.character(result_id), databases)
}

.sn_enrich_symbol_to_entrez <- function(genes, org_db) {
  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (length(genes) == 0) {
    return(data.frame(SYMBOL = character(0), ENTREZID = character(0), stringsAsFactors = FALSE))
  }
  genes_sorted <- sort(genes)

  cache_key <- .sn_enrich_cache_key(
    org_db,
    as.character(utils::packageVersion(org_db)),
    paste(genes_sorted, collapse = "|")
  )
  if (exists(cache_key, envir = .sn_enrichment_cache_env$symbol_to_entrez, inherits = FALSE)) {
    return(get(cache_key, envir = .sn_enrichment_cache_env$symbol_to_entrez, inherits = FALSE))
  }

  out <- clusterProfiler::bitr(
    geneID = genes_sorted,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org_db
  )
  if (nrow(out) > 0L && "SYMBOL" %in% colnames(out)) {
    out <- out[order(match(out$SYMBOL, genes_sorted)), , drop = FALSE]
  }
  assign(cache_key, out, envir = .sn_enrichment_cache_env$symbol_to_entrez)
  out
}

#' Run gene set enrichment analysis
#'
#' Runs GO, KEGG, or MSigDB enrichment using \pkg{clusterProfiler}. It supports
#' both over-representation analysis (ORA) and ranked-list GSEA. The enrichment
#' input can be a gene vector, a ranked numeric vector, a data frame, or a
#' Seurat object paired with \code{source_de_result_id} to reuse stored DE results.
#'
#' @param x A character vector of genes, a named numeric vector for GSEA, a
#'   data frame, or a \code{Seurat} object when enriching a stored DE result.
#' @param mapping Optional formula mapping input columns to enrichment roles.
#'   Use \code{gene ~ group} for grouped ORA, \code{gene ~ score} for one
#'   global GSEA ranking, and \code{gene ~ score | group} for grouped GSEA.
#'   Multiple ORA or GSEA grouping columns can be joined with "+".
#' @param analysis Optional explicit analysis mode. Named numeric vectors and
#'   numeric formula RHS values infer GSEA, while character/categorical inputs
#'   infer ORA. Supply \code{analysis = "ora"} explicitly when numeric formula
#'   RHS values are intended as group codes rather than ranking statistics.
#' @param species One of \code{"human"} or \code{"mouse"}.
#' @param database One or more databases. Supported values include GO/KEGG
#'   databases such as \code{"GOBP"} and MSigDB collections such as
#'   \code{"H"}, \code{"C2"}, or \code{"C2:CP:REACTOME"}.
#' @param collection Optional MSigDB collection used when
#'   \code{database = "MSIGDB"}.
#' @param subcollection Optional MSigDB subcollection used when
#'   \code{database = "MSIGDB"} or when you want to override the parsed
#'   subcollection for a collection-level request such as \code{"C2"}.
#' @param pvalue_cutoff Cutoff passed unchanged to clusterProfiler. ORA applies
#'   the upstream raw-p, adjusted-p, and q-value reporting rules. GSEA cutoff
#'   behavior is defined by the validated clusterProfiler/enrichit version and
#'   is therefore recorded in the conformance contract.
#' @param p_adjust_method Multiple-testing adjustment method passed as
#'   `pAdjustMethod`.
#' @param qvalue_cutoff ORA q-value cutoff passed as `qvalueCutoff`. It is not
#'   used by GSEA.
#' @param universe Optional ORA background gene universe in the same symbol
#'   namespace as `x`. Stored-DE ORA defaults to the source result's recorded
#'   \code{input$tested_features}; only legacy results without that field fall
#'   back to the source assay feature space. For KEGG the universe is converted
#'   to ENTREZID together with the query genes. Supplying a universe for GSEA is
#'   an error.
#' @param min_gs_size,max_gs_size Minimum and maximum tested gene-set sizes.
#' @param gsea_exponent GSEA running-score exponent. For reproducible stochastic
#'   GSEA results with clusterProfiler 4.20, call `set.seed()` immediately before
#'   `sn_run_enrichment()`, as for the direct upstream call.
#' @param duplicate_gene_method Policy for duplicate identifiers in a GSEA
#'   ranked list. The default, `"error"`, avoids silent changes. Explicit
#'   alternatives are `"max_abs"`, `"max"`, and `"mean"`.
#' @param result_id Optional name used when storing the enrichment result on a
#'   Seurat object. When omitted, Shennong combines the source DE result,
#'   analysis mode, and database, for example \code{bulk.gsea.H}. When multiple
#'   databases are requested, one name is generated per database. Automatic
#'   IDs receive a numeric suffix when a previous result already uses the name.
#'   An explicit scalar name receives database suffixes for a multi-database
#'   request; a vector can instead name every result directly.
#' @param source_de_result_id Optional stored DE-result name associated with the
#'   enrichment input.
#' @param de_p_adjusted_cutoff,de_logfc_threshold Significance and absolute
#'   effect thresholds applied when ORA consumes a stored DE result.
#' @param de_direction Direction retained from a stored DE result for ORA:
#'   `"up"`, `"down"`, or an explicit `"both"`.
#' @param return_object Return the updated Seurat object when \code{TRUE}; this
#'   requires Seurat input. Otherwise return one unified analysis result.
#' @param prefix Optional filename prefix when writing results.
#' @param outdir Optional output directory. If supplied, each enrichment result
#'   is saved as an `.rds` file.
#' @param object Alias for \code{x}; supply only one of \code{x} and \code{object}.
#' @param gene_clusters Compatibility alias for \code{mapping}. New code should
#'   use \code{mapping}; supply only one of the two arguments.
#'
#' @return A Seurat object or one unified enrichment result. For multiple
#'   databases, \code{tables$primary} includes a \code{database} column and
#'   \code{models$database_results} contains per-database envelopes. Native
#'   clusterProfiler objects are available under \code{models$backend_results}.
#'   Use \code{sn_get_enrichment_result()} to select a table.
#'
#' @examples
#' \dontrun{
#' sn_run_enrichment(
#'   x = c("CD3D", "IL7R", "LTB"),
#'   species = "human",
#'   database = c("GOBP", "H")
#' )
#' sn_run_enrichment(
#'   x = marker_table,
#'   mapping = gene ~ avg_log2FC | cell_type,
#'   species = "human",
#'   database = "H"
#' )
#' }
#'
#' @export
sn_run_enrichment <- function(
  x,
  mapping = NULL,
  analysis = NULL,
  species = NULL,
  database = "GOBP",
  collection = NULL,
  subcollection = NULL,
  pvalue_cutoff = 0.05,
  p_adjust_method = "BH",
  qvalue_cutoff = 0.2,
  universe = NULL,
  min_gs_size = 10,
  max_gs_size = 500,
  gsea_exponent = 1,
  duplicate_gene_method = c("error", "max_abs", "max", "mean"),
  result_id = NULL,
  source_de_result_id = NULL,
  return_object = inherits(x, "Seurat"),
  prefix = NULL,
  outdir = NULL,
  object = NULL,
  de_p_adjusted_cutoff = 0.05,
  de_logfc_threshold = 0,
  de_direction = c("up", "down", "both"),
  gene_clusters = NULL) {
  x <- .sn_resolve_object_alias(x, object, missing(x))
  resolved <- .sn_enrich_resolve_input(
    x = x,
    source_de_result_id = source_de_result_id
  )
  input <- resolved$input
  object <- resolved$object
  if (isTRUE(return_object) && is.null(object)) {
    stop("`return_object = TRUE` requires a Seurat input.", call. = FALSE)
  }
  source_de_result_id <- resolved$source_de_result_id %||% source_de_result_id

  species <- species %||% if (!is_null(object)) tryCatch(sn_get_species(object), error = function(...) NULL) else NULL
  species <- species %||% "human"
  if (!is.character(species) || length(species) != 1L || is.na(species)) {
    stop("`species` must be one of \"human\" or \"mouse\".", call. = FALSE)
  }
  species <- match.arg(tolower(species), c("human", "mouse"))
  check_installed(pkg = "clusterProfiler")

  databases <- .sn_enrich_normalize_database_labels(database)
  msigdb_cfgs <- stats::setNames(
    lapply(databases, .sn_enrich_parse_msigdb_database, collection = collection, subcollection = subcollection),
    databases
  )
  if (any(vapply(msigdb_cfgs, Negate(is.null), logical(1)))) {
    check_installed(pkg = "msigdbr")
  }
  if (any(databases %in% c("GO", "GOBP", "GOMF", "GOCC", "KEGG"))) {
    check_installed(pkg = if (species == "human") "org.Hs.eg.db" else "org.Mm.eg.db")
  }

  org_db <- switch(EXPR = species,
    "human" = "org.Hs.eg.db",
    "mouse" = "org.Mm.eg.db"
  )
  organism <- switch(EXPR = species,
    "human" = "hsa",
    "mouse" = "mmu"
  )

  mapping <- .sn_enrich_resolve_mapping(
    mapping = mapping,
    gene_clusters = gene_clusters
  )
  gene_col <- mapping$gene_col %||% "gene"
  analysis <- .sn_enrich_resolve_analysis(
    input = input,
    mapping = mapping,
    analysis = analysis
  )

  validate_probability <- function(value, name) {
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value) || value < 0 || value > 1) {
      stop("`", name, "` must be one finite number between 0 and 1.", call. = FALSE)
    }
    as.numeric(value)
  }
  pvalue_cutoff <- validate_probability(pvalue_cutoff, "pvalue_cutoff")
  qvalue_cutoff <- validate_probability(qvalue_cutoff, "qvalue_cutoff")
  p_adjust_method <- match.arg(p_adjust_method, stats::p.adjust.methods)
  duplicate_gene_method <- match.arg(duplicate_gene_method)
  de_direction <- match.arg(de_direction)
  de_p_adjusted_cutoff <- validate_probability(de_p_adjusted_cutoff, "de_p_adjusted_cutoff")
  if (!is.numeric(de_logfc_threshold) || length(de_logfc_threshold) != 1L ||
      is.na(de_logfc_threshold) || !is.finite(de_logfc_threshold) || de_logfc_threshold < 0) {
    stop("`de_logfc_threshold` must be one finite non-negative number.", call. = FALSE)
  }
  if (!is.numeric(min_gs_size) || length(min_gs_size) != 1L ||
      is.na(min_gs_size) || !is.finite(min_gs_size) || min_gs_size < 1L ||
      min_gs_size > .Machine$integer.max || min_gs_size != as.integer(min_gs_size)) {
    stop("`min_gs_size` must be one positive integer.", call. = FALSE)
  }
  if (!is.numeric(max_gs_size) || length(max_gs_size) != 1L ||
      is.na(max_gs_size) || !is.finite(max_gs_size) ||
      max_gs_size > .Machine$integer.max || max_gs_size != as.integer(max_gs_size) ||
      max_gs_size < min_gs_size) {
    stop("`max_gs_size` must be one integer greater than or equal to `min_gs_size`.", call. = FALSE)
  }
  min_gs_size <- as.integer(min_gs_size)
  max_gs_size <- as.integer(max_gs_size)
  if (!is.numeric(gsea_exponent) || length(gsea_exponent) != 1L ||
      is.na(gsea_exponent) || !is.finite(gsea_exponent) || gsea_exponent < 0) {
    stop("`gsea_exponent` must be one finite non-negative number.", call. = FALSE)
  }
  universe_source <- if (is.null(universe)) NULL else "user"
  automatic_background <- identical(analysis, "ora") && is.null(universe)
  comparison_backgrounds <- NULL
  if (identical(analysis, "ora") && is.null(universe) &&
      !is.null(resolved$de_result) && !is.null(object)) {
    stored_universe <- .sn_enrich_stored_de_universe(
      de_result = resolved$de_result,
      object = object
    )
    universe <- stored_universe$universe
    universe_source <- stored_universe$source
  }
  if (!is.null(universe)) {
    if (identical(analysis, "gsea")) {
      stop("`universe` is an ORA parameter and cannot be supplied for GSEA.", call. = FALSE)
    }
    if (!is.character(universe)) {
      stop("`universe` must be NULL or a character vector of gene identifiers.", call. = FALSE)
    }
    universe <- .sn_enrich_resolve_gene_vector(universe)
  }

  de_ora_selection <- NULL
  if (identical(analysis, "ora") && !is_null(resolved$de_result)) {
    input <- .sn_enrich_filter_stored_de_ora(
      input = input,
      de_result = resolved$de_result,
      p_adjusted_cutoff = de_p_adjusted_cutoff,
      direction = de_direction,
      logfc_threshold = de_logfc_threshold
    )
    de_ora_selection <- attr(input, "sn_de_ora_selection")
    if (automatic_background) {
      comparison_backgrounds <- .sn_enrich_comparison_backgrounds(input, resolved$de_result, mapping)
    }
    if (is_null(mapping)) {
      group_candidates <- unique(c(
        resolved$de_result$group_col %||% character(),
        "cluster", "comparison"
      ))
      group_candidates <- group_candidates[
        !is.na(group_candidates) & group_candidates %in% colnames(input)
      ]
      group_hits <- group_candidates[vapply(group_candidates, function(column) {
        values <- as.character(input[[column]])
        length(unique(values[!is.na(values) & nzchar(values)])) > 1L
      }, logical(1))]
      group_col <- if (length(group_hits) > 0L) group_hits[[1]] else NULL
      if (!is.null(group_col)) {
        mapping <- .sn_enrich_parse_mapping(
          stats::reformulate(group_col, response = gene_col)
        )
        de_ora_selection$group_column <- group_col
      }
    }
    de_ora_selection$universe_source <- universe_source
    de_ora_selection$universe_size <- length(universe)
  }

  grouped_gsea <- identical(analysis, "gsea") && isTRUE(mapping$grouped_gsea)
  if (grouped_gsea) {
    gene_lists <- .sn_enrich_resolve_grouped_gene_lists(
      input = input,
      mapping = mapping,
      duplicate_gene_method = duplicate_gene_method
    )
    gene_list <- NULL
  } else if (identical(analysis, "gsea")) {
    gene_list <- .sn_enrich_resolve_gene_list(
      input = input,
      mapping = mapping,
      duplicate_gene_method = duplicate_gene_method
    )
    gene_lists <- NULL
  } else {
    gene_list <- NULL
    gene_lists <- NULL
  }

  run_one <- function(current_database) {
    current_cfg <- msigdb_cfgs[[current_database]]
    with_enrichment_acceleration <- function(expr) {
      if (current_database %in% c("GO", "GOBP", "GOMF", "GOCC")) {
        return(.sn_enrichment_with_rng_preserved(.sn_with_default_acceleration(
          expr,
          patches = "clusterprofiler"
        )))
      }
      .sn_enrichment_with_rng_preserved(.sn_with_acceleration_disabled(expr))
    }
    .sn_log_info("Running {toupper(analysis)} analysis for the {current_database} database.")

    if (current_database %in% c("GO", "GOBP", "GOMF", "GOCC")) {
      ont <- switch(EXPR = current_database,
        "GO" = "ALL",
        "GOBP" = "BP",
        "GOMF" = "MF",
        "GOCC" = "CC"
      )

      if (grouped_gsea) {
        result <- .sn_enrich_muffle_empty_warning(
          with_enrichment_acceleration(
            clusterProfiler::compareCluster(
              geneClusters = gene_lists,
              fun = "gseGO",
              ont = ont,
              OrgDb = org_db,
              keyType = "SYMBOL",
              exponent = gsea_exponent,
              minGSSize = min_gs_size,
              maxGSSize = max_gs_size,
              pvalueCutoff = pvalue_cutoff,
              pAdjustMethod = p_adjust_method
            )
          )
        )
      } else if (identical(analysis, "gsea")) {
        result <- .sn_enrich_muffle_empty_warning(
          with_enrichment_acceleration(
            clusterProfiler::gseGO(
              geneList = gene_list,
              ont = ont,
              OrgDb = org_db,
              keyType = "SYMBOL",
              exponent = gsea_exponent,
              minGSSize = min_gs_size,
              maxGSSize = max_gs_size,
              pvalueCutoff = pvalue_cutoff,
              pAdjustMethod = p_adjust_method
            )
          )
        )
      } else if (is_null(mapping)) {
        result <- .sn_enrich_muffle_empty_warning(
          with_enrichment_acceleration(
            clusterProfiler::enrichGO(
              gene = .sn_enrich_resolve_gene_vector(input, gene_col = gene_col),
              ont = ont,
              OrgDb = org_db,
              keyType = "SYMBOL",
              pvalueCutoff = pvalue_cutoff,
              pAdjustMethod = p_adjust_method,
              universe = universe,
              qvalueCutoff = qvalue_cutoff,
              minGSSize = min_gs_size,
              maxGSSize = max_gs_size
            )
          )
        )
      } else {
        result <- .sn_enrich_muffle_empty_warning(
          with_enrichment_acceleration(
            clusterProfiler::compareCluster(
              geneClusters = mapping$formula,
              fun = "enrichGO",
              ont = ont,
              data = input,
              OrgDb = org_db,
              keyType = "SYMBOL",
              pvalueCutoff = pvalue_cutoff,
              pAdjustMethod = p_adjust_method,
              universe = universe,
              qvalueCutoff = qvalue_cutoff,
              minGSSize = min_gs_size,
              maxGSSize = max_gs_size
            )
          )
        )
      }

      return(result)
    }

    if (identical(current_database, "KEGG")) {
      if (grouped_gsea) {
        kegg_gene_lists <- lapply(gene_lists, function(current_gene_list) {
          gid <- .sn_enrich_symbol_to_entrez(
            genes = names(current_gene_list),
            org_db = org_db
          )
          converted <- current_gene_list[gid$SYMBOL]
          names(converted) <- gid$ENTREZID
          .sn_enrich_collapse_gene_list(
            converted,
            duplicate_gene_method = duplicate_gene_method
          )
        })
        result <- .sn_enrich_muffle_empty_warning(
          with_enrichment_acceleration(
            clusterProfiler::compareCluster(
              geneClusters = kegg_gene_lists,
              fun = "gseKEGG",
              organism = organism,
              exponent = gsea_exponent,
              minGSSize = min_gs_size,
              maxGSSize = max_gs_size,
              pvalueCutoff = pvalue_cutoff,
              pAdjustMethod = p_adjust_method
            )
          )
        )
      } else if (identical(analysis, "gsea")) {
        gid <- .sn_enrich_symbol_to_entrez(
          genes = names(gene_list),
          org_db = org_db
        )
        kegg_gene_list <- gene_list[gid$SYMBOL]
        names(kegg_gene_list) <- gid$ENTREZID
        kegg_gene_list <- .sn_enrich_collapse_gene_list(
          kegg_gene_list,
          duplicate_gene_method = duplicate_gene_method
        )
        result <- .sn_enrich_muffle_empty_warning(
          with_enrichment_acceleration(
            clusterProfiler::gseKEGG(
              geneList = kegg_gene_list,
              organism = organism,
              exponent = gsea_exponent,
              minGSSize = min_gs_size,
              maxGSSize = max_gs_size,
              pvalueCutoff = pvalue_cutoff,
              pAdjustMethod = p_adjust_method
            )
          )
        )
      } else if (is_null(mapping)) {
        gid <- .sn_enrich_symbol_to_entrez(
          genes = .sn_enrich_resolve_gene_vector(input, gene_col = gene_col),
          org_db = org_db
        )
        kegg_universe <- if (is.null(universe)) {
          NULL
        } else {
          .sn_enrich_symbol_to_entrez(universe, org_db = org_db)$ENTREZID
        }
        result <- .sn_enrich_muffle_empty_warning(
          with_enrichment_acceleration(
            clusterProfiler::enrichKEGG(
              gene = gid$ENTREZID,
              organism = organism,
              pvalueCutoff = pvalue_cutoff,
              pAdjustMethod = p_adjust_method,
              universe = unique(kegg_universe),
              minGSSize = min_gs_size,
              maxGSSize = max_gs_size,
              qvalueCutoff = qvalue_cutoff
            )
          )
        )
      } else {
        gene_ids <- unique(as.character(input[[mapping$gene_col]]))
        gid <- .sn_enrich_symbol_to_entrez(
          genes = gene_ids,
          org_db = org_db
        )
        kegg_input <- dplyr::inner_join(
          x = as.data.frame(input),
          y = gid,
          by = stats::setNames("SYMBOL", mapping$gene_col)
        ) |>
          dplyr::select(-dplyr::all_of(mapping$gene_col)) |>
          dplyr::rename(!!mapping$gene_col := dplyr::all_of("ENTREZID"))
        kegg_universe <- if (is.null(universe)) {
          NULL
        } else {
          .sn_enrich_symbol_to_entrez(universe, org_db = org_db)$ENTREZID
        }

        result <- .sn_enrich_muffle_empty_warning(
          with_enrichment_acceleration(
            clusterProfiler::compareCluster(
              geneClusters = mapping$formula,
              data = kegg_input,
              fun = "enrichKEGG",
              pvalueCutoff = pvalue_cutoff,
              pAdjustMethod = p_adjust_method,
              universe = unique(kegg_universe),
              minGSSize = min_gs_size,
              maxGSSize = max_gs_size,
              qvalueCutoff = qvalue_cutoff,
              organism = organism
            )
          )
        )
      }

      if (!is.null(result) && !inherits(result, "gseaResult")) {
        result <- clusterProfiler::setReadable(result,
          OrgDb = org_db,
          keyType = "ENTREZID"
        )
      }
      return(result)
    }

    if (!is_null(current_cfg)) {
      msigdb_tbl <- .sn_enrich_get_msigdb_terms(
        species = species,
        collection = current_cfg$collection,
        subcollection = current_cfg$subcollection
      )
      term2gene <- dplyr::select(msigdb_tbl, "term", "gene")
      term2name <- msigdb_tbl |>
        dplyr::select("term", "description") |>
        dplyr::distinct()

      if (grouped_gsea) {
        result <- .sn_enrich_muffle_empty_warning(
          with_enrichment_acceleration(
            clusterProfiler::compareCluster(
              geneClusters = gene_lists,
              fun = clusterProfiler::GSEA,
              exponent = gsea_exponent,
              minGSSize = min_gs_size,
              maxGSSize = max_gs_size,
              pvalueCutoff = pvalue_cutoff,
              pAdjustMethod = p_adjust_method,
              TERM2GENE = term2gene,
              TERM2NAME = term2name
            )
          )
        )
      } else if (identical(analysis, "gsea")) {
        result <- .sn_enrich_muffle_empty_warning(
          with_enrichment_acceleration(
            clusterProfiler::GSEA(
              geneList = gene_list,
              exponent = gsea_exponent,
              minGSSize = min_gs_size,
              maxGSSize = max_gs_size,
              pvalueCutoff = pvalue_cutoff,
              pAdjustMethod = p_adjust_method,
              TERM2GENE = term2gene,
              TERM2NAME = term2name
            )
          )
        )
      } else if (is_null(mapping)) {
        result <- .sn_enrich_muffle_empty_warning(
          with_enrichment_acceleration(
            clusterProfiler::enricher(
              gene = .sn_enrich_resolve_gene_vector(input, gene_col = gene_col),
              TERM2GENE = term2gene,
              TERM2NAME = term2name,
              pvalueCutoff = pvalue_cutoff,
              pAdjustMethod = p_adjust_method,
              universe = universe,
              minGSSize = min_gs_size,
              maxGSSize = max_gs_size,
              qvalueCutoff = qvalue_cutoff
            )
          )
        )
      } else {
        result <- .sn_enrich_muffle_empty_warning(
          with_enrichment_acceleration(
            clusterProfiler::compareCluster(
              geneClusters = mapping$formula,
              data = input,
              fun = clusterProfiler::enricher,
              TERM2GENE = term2gene,
              TERM2NAME = term2name,
              pvalueCutoff = pvalue_cutoff,
              pAdjustMethod = p_adjust_method,
              universe = universe,
              minGSSize = min_gs_size,
              maxGSSize = max_gs_size,
              qvalueCutoff = qvalue_cutoff
            )
          )
        )
      }

      return(result)
    }

    stop(glue("Unsupported database '{current_database}'."), call. = FALSE)
  }

  .sn_with_acceleration_provenance_context({
    results <- stats::setNames(lapply(databases, function(database) {
      if (is.null(comparison_backgrounds)) return(run_one(database))
      .sn_enrich_run_comparisons(run_one, database, comparison_backgrounds, gene_col, grouped = !is.null(mapping))
    }), databases)

    if (!is_null(outdir)) {
      outdir <- sn_set_path(path = outdir)
      prefix_part <- if (!is_null(prefix) && nzchar(prefix)) paste0(prefix, ".") else ""
      for (current_database in names(results)) {
        filename <- glue("{prefix_part}enrichment.{.sn_enrich_output_label(current_database)}.rds")
        saveRDS(results[[current_database]], file = file.path(outdir, filename))
      }
    }

    {
      result_ids <- .sn_enrich_result_ids(
        result_id = result_id,
        databases = names(results),
        source_de_result_id = source_de_result_id,
        analysis = analysis
      )
      stored_results <- stats::setNames(vector("list", length(results)), names(results))
      for (current_database in names(results)) {
        current_result_id <- if (length(databases) == 1) {
          unname(result_ids[[1L]])
        } else {
          unname(result_ids[[current_database]])
        }
        current_result_id <- .sn_resolve_new_result_id(
          object, "enrichment", if (is.null(result_id)) NULL else current_result_id,
          current_result_id
        )
        stored_results[[current_database]] <- .sn_build_enrichment_result(
          result = results[[current_database]],
          result_id = current_result_id,
          analysis = analysis,
          database = current_database,
          species = species,
          source_de_result_id = source_de_result_id,
          gene_col = gene_col,
          score_col = if (identical(analysis, "gsea") && !is_null(mapping)) mapping$value_col else NULL,
          parameters = list(
            backend_versions = {
              backend_packages <- c(
                "clusterProfiler",
                if (requireNamespace("enrichit", quietly = TRUE)) "enrichit" else character(),
                if (!is.null(msigdb_cfgs[[current_database]])) "msigdbr" else character(),
                if (current_database %in% c("GO", "GOBP", "GOMF", "GOCC", "KEGG")) org_db else character()
              )
              stats::setNames(
                vapply(backend_packages, function(package) {
                  as.character(utils::packageVersion(package))
                }, character(1)),
                backend_packages
              )
            },
            pvalue_cutoff = pvalue_cutoff,
            p_adjust_method = p_adjust_method,
            qvalue_cutoff = if (identical(analysis, "ora")) qvalue_cutoff else NULL,
            universe = if (identical(analysis, "ora")) universe else NULL,
            universe_source = if (identical(analysis, "ora")) universe_source else NULL,
            universe_size = if (identical(analysis, "ora")) length(universe) else NULL,
            comparison_backgrounds = if (is.null(comparison_backgrounds)) NULL else lapply(
              comparison_backgrounds, function(x) x[c("key", "universe")]
            ),
            min_gs_size = min_gs_size,
            max_gs_size = max_gs_size,
            gsea_exponent = if (identical(analysis, "gsea")) gsea_exponent else NULL,
            duplicate_gene_method = if (identical(analysis, "gsea")) duplicate_gene_method else NULL,
            group_columns = if (grouped_gsea) mapping$group_cols else NULL,
            de_ora_selection = if (identical(analysis, "ora")) de_ora_selection else NULL
          )
        )
        if (isTRUE(return_object)) {
          object <- sn_store_result(object, "enrichment", current_result_id, stored_results[[current_database]])
        }
      }

      if (isTRUE(return_object)) {
        return(object)
      }
    }

    if (length(stored_results) == 1L) return(stored_results[[1L]])
    combined <- dplyr::bind_rows(lapply(names(stored_results), function(database) {
      table <- stored_results[[database]]$tables$primary
      table$database <- rep(database, nrow(table))
      table
    }))
    ids <- vapply(stored_results, function(result) result$result_id, character(1))
    .sn_prepare_result(list(
      table = combined, analysis = analysis, method = analysis, backend = "clusterProfiler",
      database = names(results), species = species, source_de_result_id = source_de_result_id,
      parameters = list(by_database = lapply(stored_results, function(result) result$parameters)),
      models = list(backend_results = results, database_results = stored_results),
      diagnostics = list(result_ids = ids), provenance = .sn_contextual_analysis_provenance()
    ), type = "enrichment", result_id = paste(unname(ids), collapse = "+"))
  }, patches = "clusterprofiler")
}
#' Store an enrichment result on a Seurat object
#'
#' This helper stores enrichment output inside
#' the canonical Shennong result registry so interpretation and writing
#' helpers can reuse it later.
#'
#' @param object A \code{Seurat} object.
#' @param result An enrichment result object or data frame coercible with
#'   \code{as.data.frame()}.
#' @param result_id Stable identifier for the stored enrichment result.
#' @param analysis One of \code{"ora"} or \code{"gsea"}.
#' @param database Database used for enrichment, for example \code{"GOBP"}.
#' @param species Species label used in the enrichment run.
#' @param source_de_result_id Optional stored DE result name that produced the input
#'   ranked gene list or gene set.
#' @param gene_col Column containing gene symbols when the enrichment input came
#'   from a data frame.
#' @param score_col Column containing ranking scores for GSEA inputs.
#' @param parameters Named list of effective enrichment parameters retained for
#'   discovery and reproducibility.
#' @param return_object If \code{TRUE}, return the updated Seurat object.
#'
#' @return A \code{Seurat} object or a stored-result list.
#'
#' @examples
#' if (requireNamespace("Seurat", quietly = TRUE)) {
#'   counts <- matrix(rpois(10 * 12, lambda = 1), nrow = 10, ncol = 12)
#'   rownames(counts) <- c(
#'     "CD3D", "CD3E", "TRAC", "LTB", "MS4A1",
#'     "CD79A", "HLA-DRA", "LYZ", "ACTB", "MALAT1"
#'   )
#'   colnames(counts) <- paste0("cell", 1:12)
#'   obj <- sn_initialize_seurat_object(counts, species = "human")
#'   enrich_tbl <- tibble::tibble(
#'     ID = c("GO:0001", "GO:0002"),
#'     Description = c("immune response", "lymphocyte activation"),
#'     NES = c(2.1, 1.7),
#'     p.adjust = c(0.01, 0.03)
#'   )
#'   obj <- sn_store_enrichment(obj, enrich_tbl, result_id = "demo_gsea")
#'   sn_list_results(obj, type = "enrichment")
#'   sn_get_result(obj, type = "enrichment", result_id = "demo_gsea")
#' }
#' @export
sn_store_enrichment <- function(object,
                                result,
                                result_id = "default",
                                analysis = c("ora", "gsea"),
                                database = "GOBP",
                                species = NULL,
                                source_de_result_id = NULL,
                                gene_col = "gene",
                                score_col = NULL,
                                parameters = list(),
                                return_object = TRUE) {
  result_id <- .sn_validate_result_id(result_id)
  .sn_validate_seurat_object(object)

  analysis <- match.arg(analysis)
  if (!is.list(parameters) ||
      (length(parameters) > 0L &&
        (is.null(names(parameters)) || any(!nzchar(names(parameters)))))) {
    stop("`parameters` must be a named list.", call. = FALSE)
  }
  stored_result <- .sn_build_enrichment_result(
    result, result_id, analysis, database, species, source_de_result_id,
    gene_col, score_col, parameters
  )

  object <- sn_store_result(
    object = object,
    type = "enrichment",
    result_id = result_id,
    result = stored_result
  )

  if (return_object) {
    return(.sn_log_seurat_command(object = object, name = "sn_store_enrichment"))
  }

  sn_get_result(
    object = object,
    type = "enrichment",
    result_id = result_id
  )
}

.sn_build_enrichment_result <- function(result, result_id, analysis, database, species,
                                        source_de_result_id, gene_col, score_col, parameters) {
  stored_result <- list(
    schema_version = .sn_analysis_result_schema_version(),
    package_version = as.character(utils::packageVersion("Shennong")),
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    table = .sn_as_enrichment_table(result),
    method = analysis,
    backend = "clusterProfiler",
    models = list(backend_results = stats::setNames(list(result), database)),
    analysis = analysis,
    database = database,
    species = species,
    source_de_result_id = source_de_result_id,
    gene_col = gene_col,
    score_col = score_col,
    parameters = parameters,
    provenance = .sn_contextual_analysis_provenance()
  )
  .sn_prepare_result(stored_result, type = "enrichment", result_id = result_id)
}

#' Deprecated alias of `sn_run_enrichment()`
#'
#' `sn_enrich()` is a deprecated compatibility alias. Use [sn_run_enrichment()] directly;
#' the alias will be removed in a future release.
#'
#' @param ... Named arguments passed on to [sn_run_enrichment()].
#'
#' @return Result of `sn_run_enrichment(...)`.
#'
#' @export
sn_enrich <- function(...) {
  .Deprecated("sn_run_enrichment", package = "Shennong")
  sn_run_enrichment(...)
}
