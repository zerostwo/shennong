.sn_enrichment_cache_env <- local({
  env <- new.env(parent = emptyenv())
  env$msigdb_terms <- new.env(parent = emptyenv())
  env$symbol_to_entrez <- new.env(parent = emptyenv())
  env
})

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

.sn_enrich_parse_formula <- function(gene_clusters) {
  if (is_null(gene_clusters)) {
    return(NULL)
  }
  if (!inherits(gene_clusters, "formula")) {
    stop("`gene_clusters` must be a two-sided formula such as `gene ~ cluster` or `gene ~ log2fc`.", call. = FALSE)
  }

  vars <- all.vars(gene_clusters)
  if (length(vars) != 2) {
    stop("`gene_clusters` must be a two-sided formula with one gene column and one grouping or ranking column.", call. = FALSE)
  }

  list(
    gene_col = vars[[1]],
    value_col = vars[[2]]
  )
}

.sn_enrich_resolve_input <- function(x,
                                     source_de_name = NULL) {
  if (inherits(x, "Seurat")) {
    de_results <- x@misc$de_results %||% list()
    if (length(de_results) == 0L) {
      stop("When `x` is a Seurat object, no stored DE results were found in `x@misc$de_results`.", call. = FALSE)
    }

    if (is_null(source_de_name) || !nzchar(source_de_name)) {
      available_names <- names(de_results)
      marker_names <- available_names[vapply(
        de_results,
        function(entry) identical(entry$analysis %||% NULL, "markers"),
        logical(1)
      )]
      latest_name <- function(candidates) {
        if (length(candidates) == 0L) {
          return(NULL)
        }
        created_at <- vapply(
          candidates,
          function(candidate) de_results[[candidate]]$created_at %||% "",
          character(1)
        )
        candidates[[order(created_at, decreasing = TRUE, na.last = TRUE)[[1]]]]
      }

      source_de_name <- if ("default" %in% available_names) {
        "default"
      } else if (length(available_names) == 1L) {
        available_names[[1]]
      } else {
        latest_name(marker_names) %||% latest_name(available_names)
      }

      .sn_log_info("`source_de_name` was not supplied; using stored DE result '{source_de_name}'.")
    }

    if (!source_de_name %in% names(de_results)) {
      stop(glue("Stored DE result '{source_de_name}' was not found in `x@misc$de_results`."), call. = FALSE)
    }
    return(list(
      input = de_results[[source_de_name]]$table,
      object = x,
      source_de_name = source_de_name
    ))
  }

  list(
    input = x,
    object = NULL,
    source_de_name = source_de_name
  )
}

.sn_enrich_resolve_analysis <- function(input,
                                        mapping = NULL,
                                        analysis = NULL) {
  if (!is_null(analysis)) {
    return(match.arg(analysis, c("ora", "gsea")))
  }

  if (is.numeric(input) && !is.null(names(input))) {
    return("gsea")
  }

  if (is.character(input)) {
    return("ora")
  }

  if (is.data.frame(input)) {
    if (!is_null(mapping)) {
      value <- input[[mapping$value_col]]
      if (is.numeric(value)) {
        stop(
          "`analysis` must be supplied when the formula RHS is numeric; ",
          "numeric group codes and GSEA ranking statistics are ambiguous.",
          call. = FALSE
        )
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
      "For GSEA with data-frame input, `gene_clusters` must be supplied as `gene ~ ranking_column`.",
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

.sn_enrich_store_names <- function(store_name, databases) {
  if (length(databases) == 1) {
    return(store_name[[1]])
  }

  if (length(store_name) == 1) {
    return(stats::setNames(
      paste(store_name[[1]], .sn_enrich_output_label(databases), sep = "."),
      databases
    ))
  }

  if (length(store_name) != length(databases)) {
    stop("`store_name` must have length 1 or match the length of `database`.", call. = FALSE)
  }

  stats::setNames(as.character(store_name), databases)
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
#' Seurat object paired with \code{source_de_name} to reuse stored DE results.
#'
#' @param x A character vector of genes, a named numeric vector for GSEA, a
#'   data frame, or a \code{Seurat} object when enriching a stored DE result.
#' @param gene_clusters Optional two-sided formula describing the gene column
#'   and the grouping/ranking column. Examples include \code{gene ~ cluster} for
#'   grouped ORA and \code{gene ~ log2fc} for one global GSEA ranking. The
#'   clusterProfiler grouped-GSEA formula \code{gene | score ~ group} is not yet
#'   supported and fails explicitly rather than being treated as global GSEA.
#' @param analysis Optional explicit analysis mode. Named numeric vectors infer
#'   GSEA and character/categorical inputs infer ORA. Supply `analysis`
#'   explicitly for a numeric formula RHS because numeric group codes and GSEA
#'   ranking statistics are otherwise ambiguous.
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
#'   namespace as `x`. For KEGG it is converted to ENTREZID together with the
#'   query genes. Supplying a universe for GSEA is an error.
#' @param min_gs_size,max_gs_size Minimum and maximum tested gene-set sizes.
#' @param gsea_exponent GSEA running-score exponent. For reproducible stochastic
#'   GSEA results with clusterProfiler 4.20, call `set.seed()` immediately before
#'   `sn_enrich()`, as for the direct upstream call.
#' @param duplicate_gene_method Policy for duplicate identifiers in a GSEA
#'   ranked list. The default, `"error"`, avoids silent changes. Explicit
#'   alternatives are `"max_abs"`, `"max"`, and `"mean"`.
#' @param store_name Name used when storing the enrichment result on a Seurat
#'   object. When multiple databases are requested, the database label is
#'   appended automatically unless a vector of names is supplied.
#' @param source_de_name Optional stored DE-result name associated with the
#'   enrichment input.
#' @param return_object Logical; when \code{TRUE} and a Seurat object is
#'   available, return the updated Seurat object instead of raw enrichment
#'   results.
#' @param prefix Optional filename prefix when writing results.
#' @param outdir Optional output directory. If supplied, each enrichment result
#'   is saved as an `.rds` file.
#' @param object Alias for \code{x}; supply only one of \code{x} and \code{object}.
#'
#' @return A single `clusterProfiler` result, a named list of results when
#'   multiple databases are requested, or a \code{Seurat} object when
#'   \code{return_object = TRUE}.
#'
#' @examples
#' \dontrun{
#' sn_enrich(
#'   x = c("CD3D", "IL7R", "LTB"),
#'   species = "human",
#'   database = c("GOBP", "H")
#' )
#' }
#'
#' @export
sn_enrich <- function(
  x,
  gene_clusters = NULL,
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
  store_name = "default",
  source_de_name = NULL,
  return_object = inherits(x, "Seurat"),
  prefix = NULL,
  outdir = NULL,
  object = NULL
) {
  x <- .sn_resolve_object_alias(x, object, missing(x))
  resolved <- .sn_enrich_resolve_input(
    x = x,
    source_de_name = source_de_name
  )
  input <- resolved$input
  object <- resolved$object
  source_de_name <- resolved$source_de_name %||% source_de_name

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

  mapping <- .sn_enrich_parse_formula(gene_clusters)
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
  if (!is.null(universe)) {
    if (identical(analysis, "gsea")) {
      stop("`universe` is an ORA parameter and cannot be supplied for GSEA.", call. = FALSE)
    }
    if (!is.character(universe)) {
      stop("`universe` must be NULL or a character vector of gene identifiers.", call. = FALSE)
    }
    universe <- .sn_enrich_resolve_gene_vector(universe)
  }

  if (identical(analysis, "gsea")) {
    gene_list <- .sn_enrich_resolve_gene_list(
      input = input,
      mapping = mapping,
      duplicate_gene_method = duplicate_gene_method
    )
  } else {
    gene_list <- NULL
  }

  run_one <- function(current_database) {
    current_cfg <- msigdb_cfgs[[current_database]]
    with_enrichment_acceleration <- function(expr) {
      if (current_database %in% c("GO", "GOBP", "GOMF", "GOCC")) {
        return(.sn_with_default_acceleration(
          expr,
          patches = "clusterprofiler"
        ))
      }
      .sn_with_acceleration_disabled(expr)
    }
    .sn_log_info("Running {toupper(analysis)} analysis for the {current_database} database.")

    if (current_database %in% c("GO", "GOBP", "GOMF", "GOCC")) {
      ont <- switch(EXPR = current_database,
        "GO" = "ALL",
        "GOBP" = "BP",
        "GOMF" = "MF",
        "GOCC" = "CC"
      )

      if (identical(analysis, "gsea")) {
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
              geneClusters = gene_clusters,
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
      if (identical(analysis, "gsea")) {
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
              geneClusters = stats::as.formula(glue("{mapping$gene_col} ~ {mapping$value_col}")),
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

      if (identical(analysis, "gsea")) {
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
              geneClusters = gene_clusters,
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
    results <- stats::setNames(lapply(databases, run_one), databases)

    if (!is_null(outdir)) {
      outdir <- sn_set_path(path = outdir)
      prefix_part <- if (!is_null(prefix) && nzchar(prefix)) paste0(prefix, ".") else ""
      for (current_database in names(results)) {
        filename <- glue("{prefix_part}enrichment.{.sn_enrich_output_label(current_database)}.rds")
        saveRDS(results[[current_database]], file = file.path(outdir, filename))
      }
    }

    if (!is_null(object)) {
      store_names <- .sn_enrich_store_names(store_name = store_name, databases = names(results))
      for (current_database in names(results)) {
        current_store <- if (length(databases) == 1) {
          store_names
        } else {
          store_names[[current_database]]
        }
        object <- sn_store_enrichment(
          object = object,
          result = results[[current_database]],
          store_name = current_store,
          analysis = analysis,
          database = current_database,
          species = species,
          source_de_name = source_de_name,
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
            min_gs_size = min_gs_size,
            max_gs_size = max_gs_size,
            gsea_exponent = if (identical(analysis, "gsea")) gsea_exponent else NULL,
            duplicate_gene_method = if (identical(analysis, "gsea")) duplicate_gene_method else NULL
          ),
          return_object = TRUE
        )
      }

      if (isTRUE(return_object)) {
        return(object)
      }
    }

    if (length(results) == 1) {
      return(results[[1]])
    }

    results
  }, patches = "clusterprofiler")
}
#' Store an enrichment result on a Seurat object
#'
#' This helper stores enrichment output inside
#' `object@misc$enrichment_results[[store_name]]` so interpretation and writing
#' helpers can reuse it later.
#'
#' @param object A \code{Seurat} object.
#' @param result An enrichment result object or data frame coercible with
#'   \code{as.data.frame()}.
#' @param store_name Name used under \code{object@misc$enrichment_results}.
#' @param analysis One of \code{"ora"} or \code{"gsea"}.
#' @param database Database used for enrichment, for example \code{"GOBP"}.
#' @param species Species label used in the enrichment run.
#' @param source_de_name Optional stored DE result name that produced the input
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
#'   obj <- sn_store_enrichment(obj, enrich_tbl, store_name = "demo_gsea")
#'   names(obj@misc$enrichment_results)
#' }
#' @export
sn_store_enrichment <- function(object,
                                result,
                                store_name = "default",
                                analysis = c("ora", "gsea"),
                                database = "GOBP",
                                species = NULL,
                                source_de_name = NULL,
                                gene_col = "gene",
                                score_col = NULL,
                                parameters = list(),
                                return_object = TRUE) {
  .sn_validate_seurat_object(object)

  analysis <- match.arg(analysis)
  if (!is.list(parameters) ||
      (length(parameters) > 0L &&
        (is.null(names(parameters)) || any(!nzchar(names(parameters)))))) {
    stop("`parameters` must be a named list.", call. = FALSE)
  }
  stored_result <- list(
    schema_version = "1.0.0",
    package_version = as.character(utils::packageVersion("Shennong")),
    created_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
    table = .sn_as_enrichment_table(result),
    analysis = analysis,
    database = database,
    species = species,
    source_de_name = source_de_name,
    gene_col = gene_col,
    score_col = score_col,
    parameters = parameters,
    provenance = .sn_contextual_analysis_provenance()
  )

  object <- .sn_store_misc_result(
    object = object,
    collection = "enrichment_results",
    store_name = store_name,
    result = stored_result
  )

  if (return_object) {
    return(.sn_log_seurat_command(object = object, name = "sn_store_enrichment"))
  }

  .sn_get_misc_result(
    object = object,
    collection = "enrichment_results",
    store_name = store_name
  )
}
