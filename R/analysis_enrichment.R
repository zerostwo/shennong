#' Run gene set enrichment analysis
#'
#' Runs GO or KEGG enrichment for a gene vector or grouped gene table using
#' `clusterProfiler`. It supports both over-representation analysis (ORA) and
#' ranked-list GSEA.
#'
#' @param x A character vector of gene symbols, or a data frame used together
#'   with `gene_clusters` / ranked-GSEA columns.
#' @param gene_clusters An optional grouping formula passed to
#'   `clusterProfiler::compareCluster()`.
#' @param analysis One of `"ora"` or `"gsea"`.
#' @param species One of `"human"` or `"mouse"`.
#' @param database One of `"GO"`, `"GOBP"`, `"GOMF"`, `"GOCC"`, or `"KEGG"`.
#' @param gene_col Column containing gene symbols when `x` is a data frame.
#' @param score_col Column containing ranking scores for `analysis = "gsea"`.
#' @param pvalue_cutoff Numeric p-value cutoff used by the enrichment method.
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
  gene_clusters = NULL,
  analysis = c("ora", "gsea"),
  species = "human",
  database = "GOBP",
  gene_col = "gene",
  score_col = NULL,
  pvalue_cutoff = 0.05,
  prefix = NULL,
  outdir = NULL
) {
  analysis <- match.arg(analysis)
  # Check if required packages are installed
  if (species == "human") {
    check_installed(pkg = c("clusterProfiler", "org.Hs.eg.db"))
  }
  if (species == "mouse") {
    check_installed(pkg = c("clusterProfiler", "org.Mm.eg.db"))
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

  # Return the enriched terms or pathways
  return(results)
}
