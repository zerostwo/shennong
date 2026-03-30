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
      return(if (is.numeric(value)) "gsea" else "ora")
    }
  }

  "ora"
}

.sn_enrich_resolve_gene_vector <- function(input, gene_col = "gene") {
  if (is.character(input)) {
    return(unique(as.character(input)))
  }

  if (is.data.frame(input)) {
    if (!gene_col %in% colnames(input)) {
      stop(glue("Column '{gene_col}' was not found in `x`."), call. = FALSE)
    }
    return(unique(as.character(input[[gene_col]])))
  }

  stop("ORA input must be a character vector or a data frame with a gene column.", call. = FALSE)
}

.sn_enrich_resolve_gene_list <- function(input,
                                         mapping = NULL) {
  if (is.numeric(input) && !is.null(names(input))) {
    gene_list <- stats::setNames(as.numeric(input), names(input))
    return(sort(gene_list, decreasing = TRUE))
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
  gene_list <- tapply(gene_list, names(gene_list), max)
  gene_list <- stats::setNames(as.numeric(gene_list), names(gene_list))
  sort(gene_list, decreasing = TRUE)
}

.sn_enrich_filter_by_pvalue <- function(result, pvalue_cutoff) {
  slot_name <- if (inherits(result, "compareClusterResult")) {
    "compareClusterResult"
  } else if (inherits(result, c("enrichResult", "gseaResult"))) {
    "result"
  } else {
    return(result)
  }

  result_table <- methods::slot(result, slot_name)
  p_col <- c("pvalue", "p.value")[c("pvalue", "p.value") %in% colnames(result_table)][1] %||% NA_character_
  if (is.na(p_col)) {
    return(result)
  }

  filtered <- result_table[!is.na(result_table[[p_col]]) & result_table[[p_col]] <= pvalue_cutoff, , drop = FALSE]
  methods::slot(result, slot_name) <- filtered
  result
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
#'   grouped ORA and \code{gene ~ log2fc} for GSEA.
#' @param analysis Optional explicit analysis mode. If omitted, Shennong
#'   infers ORA versus GSEA from the input type or the formula RHS column type.
#' @param species One of \code{"human"} or \code{"mouse"}.
#' @param database One or more databases. Supported values include GO/KEGG
#'   databases such as \code{"GOBP"} and MSigDB collections such as
#'   \code{"H"}, \code{"C2"}, or \code{"C2:CP:REACTOME"}.
#' @param collection Optional MSigDB collection used when
#'   \code{database = "MSIGDB"}.
#' @param subcollection Optional MSigDB subcollection used when
#'   \code{database = "MSIGDB"} or when you want to override the parsed
#'   subcollection for a collection-level request such as \code{"C2"}.
#' @param pvalue_cutoff Raw p-value cutoff used to filter returned enrichment
#'   tables after the underlying enrichment call completes.
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
  store_name = "default",
  source_de_name = NULL,
  return_object = inherits(x, "Seurat"),
  prefix = NULL,
  outdir = NULL
) {
  resolved <- .sn_enrich_resolve_input(
    x = x,
    source_de_name = source_de_name
  )
  input <- resolved$input
  object <- resolved$object
  source_de_name <- resolved$source_de_name %||% source_de_name

  species <- species %||% if (!is_null(object)) tryCatch(sn_get_species(object), error = function(...) NULL) else NULL
  species <- species %||% "human"

  if (species == "human") {
    check_installed(pkg = c("clusterProfiler", "org.Hs.eg.db"))
  }
  if (species == "mouse") {
    check_installed(pkg = c("clusterProfiler", "org.Mm.eg.db"))
  }

  databases <- .sn_enrich_normalize_database_labels(database)
  msigdb_cfgs <- stats::setNames(
    lapply(databases, .sn_enrich_parse_msigdb_database, collection = collection, subcollection = subcollection),
    databases
  )
  if (any(vapply(msigdb_cfgs, Negate(is.null), logical(1)))) {
    check_installed(pkg = "msigdbr")
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

  if (identical(analysis, "gsea")) {
    gene_list <- .sn_enrich_resolve_gene_list(
      input = input,
      mapping = mapping
    )
  } else {
    gene_list <- NULL
  }

  run_one <- function(current_database) {
    current_cfg <- msigdb_cfgs[[current_database]]
    .sn_log_info("Running {toupper(analysis)} analysis for the {current_database} database.")

    if (current_database %in% c("GO", "GOBP", "GOMF", "GOCC")) {
      ont <- switch(EXPR = current_database,
        "GO" = "ALL",
        "GOBP" = "BP",
        "GOMF" = "MF",
        "GOCC" = "CC"
      )

      if (identical(analysis, "gsea")) {
        result <- clusterProfiler::gseGO(
          geneList = gene_list,
          ont = ont,
          OrgDb = org_db,
          keyType = "SYMBOL",
          pvalueCutoff = 1
        )
      } else if (is_null(mapping)) {
        result <- clusterProfiler::enrichGO(
          gene = .sn_enrich_resolve_gene_vector(input, gene_col = gene_col),
          ont = ont,
          OrgDb = org_db,
          keyType = "SYMBOL",
          pvalueCutoff = 1,
          qvalueCutoff = 1
        )
      } else {
        result <- clusterProfiler::compareCluster(
          geneClusters = gene_clusters,
          fun = "enrichGO",
          ont = ont,
          data = input,
          OrgDb = org_db,
          keyType = "SYMBOL",
          pvalueCutoff = 1,
          qvalueCutoff = 1
        )
      }

      return(.sn_enrich_filter_by_pvalue(result, pvalue_cutoff))
    }

    if (identical(current_database, "KEGG")) {
      if (identical(analysis, "gsea")) {
        gid <- .sn_enrich_symbol_to_entrez(
          genes = names(gene_list),
          org_db = org_db
        )
        kegg_gene_list <- gene_list[gid$SYMBOL]
        names(kegg_gene_list) <- gid$ENTREZID
        kegg_gene_list <- sort(kegg_gene_list, decreasing = TRUE)
        result <- clusterProfiler::gseKEGG(
          geneList = kegg_gene_list,
          organism = organism,
          pvalueCutoff = 1
        )
      } else if (is_null(mapping)) {
        gid <- .sn_enrich_symbol_to_entrez(
          genes = .sn_enrich_resolve_gene_vector(input, gene_col = gene_col),
          org_db = org_db
        )
        result <- clusterProfiler::enrichKEGG(
          gene = gid$ENTREZID,
          organism = organism,
          pvalueCutoff = 1,
          qvalueCutoff = 1
        )
      } else {
        gene_ids <- unique(as.character(input[[mapping$gene_col]]))
        gid <- .sn_enrich_symbol_to_entrez(
          genes = gene_ids,
          org_db = org_db
        )
        kegg_input <- dplyr::full_join(
          x = as.data.frame(input),
          y = gid,
          by = stats::setNames("SYMBOL", mapping$gene_col)
        ) |>
          dplyr::select(-dplyr::all_of(mapping$gene_col)) |>
          dplyr::rename(!!mapping$gene_col := dplyr::all_of("ENTREZID"))

        result <- clusterProfiler::compareCluster(
          geneClusters = stats::as.formula(glue("{mapping$gene_col} ~ {mapping$value_col}")),
          data = kegg_input,
          fun = "enrichKEGG",
          pvalueCutoff = 1,
          qvalueCutoff = 1,
          organism = organism
        )
      }

      if (!inherits(result, "gseaResult")) {
        result <- clusterProfiler::setReadable(result,
          OrgDb = org_db,
          keyType = "ENTREZID"
        )
      }
      return(.sn_enrich_filter_by_pvalue(result, pvalue_cutoff))
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
        result <- clusterProfiler::GSEA(
          geneList = gene_list,
          TERM2GENE = term2gene,
          TERM2NAME = term2name,
          pvalueCutoff = 1
        )
      } else if (is_null(mapping)) {
        result <- clusterProfiler::enricher(
          gene = .sn_enrich_resolve_gene_vector(input, gene_col = gene_col),
          TERM2GENE = term2gene,
          TERM2NAME = term2name,
          pvalueCutoff = 1,
          qvalueCutoff = 1
        )
      } else {
        result <- clusterProfiler::compareCluster(
          geneClusters = gene_clusters,
          data = input,
          fun = clusterProfiler::enricher,
          TERM2GENE = term2gene,
          TERM2NAME = term2name,
          pvalueCutoff = 1,
          qvalueCutoff = 1
        )
      }

      return(.sn_enrich_filter_by_pvalue(result, pvalue_cutoff))
    }

    stop(glue("Unsupported database '{current_database}'."), call. = FALSE)
  }

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
}
