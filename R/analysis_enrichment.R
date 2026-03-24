#' Run gene set enrichment analysis
#'
#' Runs GO or KEGG enrichment for a gene vector or grouped gene table using
#' `clusterProfiler`. It supports both over-representation analysis (ORA) and
#' ranked-list GSEA.
#'
#' @param x A character vector of gene symbols, or a data frame used together
#'   with `gene_clusters` / ranked-GSEA columns.
#' @param object Optional \code{Seurat} object used to store the enrichment
#'   result in \code{object@misc$enrichment_results[[store_name]]}. When
#'   supplied, \code{return_object} defaults to \code{TRUE}.
#' @param gene_clusters An optional grouping formula passed to
#'   `clusterProfiler::compareCluster()`.
#' @param analysis One of `"ora"` or `"gsea"`.
#' @param species One of `"human"` or `"mouse"`.
#' @param database One of `"GO"`, `"GOBP"`, `"GOMF"`, `"GOCC"`, `"KEGG"`, or
#'   an MSigDB collection handled through \pkg{msigdbr}. Supported MSigDB forms
#'   include `"H"`, `"C1"` through `"C8"`, and collection plus subcollection
#'   strings such as `"C2:CP:REACTOME"` or `"C5:GO:BP"`.
#' @param gene_col Column containing gene symbols when `x` is a data frame.
#' @param score_col Column containing ranking scores for `analysis = "gsea"`.
#' @param msigdb_subcollection Optional MSigDB subcollection used when
#'   \code{database} names an MSigDB collection such as \code{"C2"} or
#'   \code{"C5"}. Ignored for GO and KEGG.
#' @param pvalue_cutoff Numeric p-value cutoff used by the enrichment method.
#' @param store_name Name used when storing the enrichment result on
#'   \code{object}.
#' @param source_de_name Optional stored DE-result name associated with the
#'   enrichment input.
#' @param return_object Logical; when \code{TRUE} and \code{object} is supplied,
#'   return the updated Seurat object instead of the raw enrichment result.
#' @param prefix Optional filename prefix when writing results.
#' @param outdir Optional output directory. If supplied, the enrichment result
#'   is also saved as an `.rds` file.
#'
#' @return A `clusterProfiler` enrichment result object.
#'
#' @examples
#' \dontrun{
#' sn_enrich(
#'   x = c("CD3D", "IL7R", "LTB"),
#'   species = "human",
#'   database = "GOBP"
#' )
#' }
#'
#' @export
sn_enrich <- function(
  x,
  object = NULL,
  gene_clusters = NULL,
  analysis = c("ora", "gsea"),
  species = NULL,
  database = "GOBP",
  gene_col = "gene",
  score_col = NULL,
  msigdb_subcollection = NULL,
  pvalue_cutoff = 0.05,
  store_name = "default",
  source_de_name = NULL,
  return_object = !is.null(object),
  prefix = NULL,
  outdir = NULL
) {
  analysis <- match.arg(analysis)
  database <- toupper(database)
  if (!is_null(object) && !inherits(object, "Seurat")) {
    stop("`object` must be a Seurat object when supplied.")
  }
  species <- species %||% if (!is_null(object)) tryCatch(sn_get_species(object), error = function(...) NULL) else NULL
  species <- species %||% "human"

  # Check if required packages are installed
  if (species == "human") {
    check_installed(pkg = c("clusterProfiler", "org.Hs.eg.db"))
  }
  if (species == "mouse") {
    check_installed(pkg = c("clusterProfiler", "org.Mm.eg.db"))
  }
  parse_msigdb_database <- function(database, msigdb_subcollection = NULL) {
    parts <- strsplit(database, ":", fixed = TRUE)[[1]]
    collection <- parts[[1]]
    known_collections <- c("H", paste0("C", 1:8))
    if (!collection %in% known_collections) {
      return(NULL)
    }

    subcollection <- if (length(parts) > 1) {
      paste(parts[-1], collapse = ":")
    } else {
      msigdb_subcollection
    }

    list(
      collection = collection,
      subcollection = subcollection
    )
  }

  get_msigdb_terms <- function(species, collection, subcollection = NULL) {
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

    msigdbr_tbl |>
      dplyr::transmute(
        term = .data$gs_name,
        description = .data$gs_description,
        gene = .data$gene_symbol
      ) |>
      dplyr::distinct()
  }

  msigdb_cfg <- parse_msigdb_database(database, msigdb_subcollection = msigdb_subcollection)
  if (!is_null(msigdb_cfg)) {
    check_installed(pkg = "msigdbr")
  }

  # Determine the organism database based on the species input
  org_db <- switch(EXPR = species,
    "human" = "org.Hs.eg.db",
    "mouse" = "org.Mm.eg.db"
  )

  # Determine the organism based on the species input
  organism <- switch(EXPR = species,
    "human" = "hsa",
    "mouse" = "mmu"
  )

  # Print a message indicating which database is being used for gene enrichment
  # analysis
  log_info(glue("{database} database {toupper(analysis)} analysis..."))

  if (analysis == "gsea") {
    gene_list <- NULL
    if (is.numeric(x) && !is.null(names(x))) {
      gene_list <- sort(x, decreasing = TRUE)
    } else if (is.data.frame(x)) {
      if (!gene_col %in% colnames(x)) {
        stop(glue("Column '{gene_col}' was not found in `x`."))
      }
      inferred_score_col <- score_col %||% c("avg_log2FC", "avg_logFC", "log2FoldChange", "logFC", "stat")[c("avg_log2FC", "avg_logFC", "log2FoldChange", "logFC", "stat") %in% colnames(x)][1]
      if (is_null(inferred_score_col)) {
        stop("`score_col` is required for GSEA when `x` is a data frame.")
      }
      score_col <- inferred_score_col
      gene_list <- x[[score_col]]
      names(gene_list) <- x[[gene_col]]
      gene_list <- tapply(gene_list, names(gene_list), max)
      gene_list <- stats::setNames(as.numeric(gene_list), names(gene_list))
      gene_list <- sort(gene_list, decreasing = TRUE)
    } else {
      stop("For GSEA, `x` must be a named numeric vector or a data frame containing gene and score columns.")
    }
  }

  # Perform gene enrichment analysis for Gene Ontology (GO) terms
  if (database %in% c("GO", "GOBP", "GOMF", "GOCC")) {
    ont <- switch(EXPR = database,
      "GO" = "ALL",
      "GOBP" = "BP",
      "GOMF" = "MF",
      "GOCC" = "CC"
    )
    if (analysis == "gsea") {
      results <- clusterProfiler::gseGO(
        geneList = gene_list,
        ont = ont,
        OrgDb = org_db,
        keyType = "SYMBOL",
        pvalueCutoff = pvalue_cutoff
      )
    } else if (is_null(x = gene_clusters)) {
      results <- clusterProfiler::enrichGO(
        gene = x,
        ont = ont,
        OrgDb = org_db,
        keyType = "SYMBOL",
        pvalueCutoff = pvalue_cutoff
      )
    } else {
      # If gene_clusters is not NULL, perform enrichment analysis using
      # compareCluster
      results <- clusterProfiler::compareCluster(
        geneClusters = gene_clusters,
        fun = "enrichGO",
        ont = ont,
        data = x,
        OrgDb = org_db,
        keyType = "SYMBOL",
        pvalueCutoff = pvalue_cutoff
      )
    }
  }

  # Perform gene enrichment analysis for KEGG pathways
  if (database == "KEGG") {
    if (analysis == "gsea") {
      gid <- clusterProfiler::bitr(
        geneID = names(gene_list),
        fromType = "SYMBOL",
        toType = "ENTREZID",
        OrgDb = org_db
      )
      gene_list <- gene_list[gid$SYMBOL]
      names(gene_list) <- gid$ENTREZID
      gene_list <- sort(gene_list, decreasing = TRUE)
      results <- clusterProfiler::gseKEGG(
        geneList = gene_list,
        organism = organism,
        pvalueCutoff = pvalue_cutoff
      )
    } else if (is_null(x = gene_clusters)) {
      gid <-
        clusterProfiler::bitr(
          geneID = x,
          "SYMBOL",
          "ENTREZID",
          OrgDb = org_db
        )
      results <- clusterProfiler::enrichKEGG(
        gene = gid$ENTREZID,
        organism = organism,
        pvalueCutoff = pvalue_cutoff
      )
    } else {
      # If gene_clusters is not NULL, perform enrichment analysis using
      # compareCluster
      x <- as.data.frame(x = x)
      gene <- all.vars(gene_clusters)[1]
      gene_id <- unique(x = x[, gene])
      gid <-
        clusterProfiler::bitr(
          geneID = gene_id,
          fromType = "SYMBOL",
          toType = "ENTREZID",
          OrgDb = org_db
        )
      x <-
        dplyr::full_join(
          x = x,
          y = gid,
          by = c(gene = "SYMBOL")
        ) |>
        dplyr::select(-gene) |>
        dplyr::rename(gene = dplyr::all_of("ENTREZID"))
      results <-
        clusterProfiler::compareCluster(
          geneClusters = gene_clusters,
          data = x,
          fun = "enrichKEGG",
          pvalueCutoff = pvalue_cutoff,
          organism = organism
        )
    }
    if (!inherits(results, "gseaResult")) {
      results <- clusterProfiler::setReadable(results,
        OrgDb = org_db,
        keyType = "ENTREZID"
      )
    }
  }

  if (!is_null(msigdb_cfg)) {
    msigdb_tbl <- get_msigdb_terms(
      species = species,
      collection = msigdb_cfg$collection,
      subcollection = msigdb_cfg$subcollection
    )

    term2gene <- dplyr::select(msigdb_tbl, "term", "gene")
    term2name <- msigdb_tbl |>
      dplyr::select("term", "description") |>
      dplyr::distinct()

    if (analysis == "gsea") {
      results <- clusterProfiler::GSEA(
        geneList = gene_list,
        TERM2GENE = term2gene,
        TERM2NAME = term2name,
        pvalueCutoff = pvalue_cutoff
      )
    } else if (is_null(x = gene_clusters)) {
      genes <- if (is.data.frame(x)) {
        if (!gene_col %in% colnames(x)) {
          stop(glue("Column '{gene_col}' was not found in `x`."))
        }
        unique(as.character(x[[gene_col]]))
      } else {
        unique(as.character(x))
      }

      results <- clusterProfiler::enricher(
        gene = genes,
        TERM2GENE = term2gene,
        TERM2NAME = term2name,
        pvalueCutoff = pvalue_cutoff
      )
    } else {
      results <- clusterProfiler::compareCluster(
        geneClusters = gene_clusters,
        data = x,
        fun = clusterProfiler::enricher,
        TERM2GENE = term2gene,
        TERM2NAME = term2name,
        pvalueCutoff = pvalue_cutoff
      )
    }
  }

  # Save the results as an RDS file
  if (!is_null(x = outdir)) {
    outdir <- sn_set_path(path = outdir)
    if (!is_null(x = prefix)) {
      prefix <- paste0(prefix, ".")
    }
    file <- glue("{outdir}/{prefix}enrichment.{database}.rds")
    saveRDS(object = results, file = file)
    # sn_write_xlsx(
    #   x = results@compareClusterResult,
    #   filename = glue("{outdir}/{prefix}enrichment.{database}.xlsx"),
    #   sheet_names = databases
    # )
  }

  if (!is_null(object)) {
    object <- sn_store_enrichment(
      object = object,
      result = results,
      store_name = store_name,
      analysis = analysis,
      database = database,
      species = species,
      source_de_name = source_de_name,
      gene_col = gene_col,
      score_col = score_col,
      return_object = TRUE
    )

    if (isTRUE(return_object)) {
      return(object)
    }
  }

  results
}
