#' @importFrom logger log_info
#' @export
sn_enrich <- function(
    x,
    gene_clusters = NULL,
    species = "human",
    database = "GOBP",
    pvalue_cutoff = 0.05,
    prefix = NULL,
    outdir = NULL) {
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
  log_info("{database} database gene enrichment analysis...")

  # Perform gene enrichment analysis for Gene Ontology (GO) terms
  if (database %in% c("GO", "GOBP", "GOMF", "GOCC")) {
    ont <- switch(EXPR = database,
      "GO" = "ALL",
      "GOBP" = "BP",
      "GOMF" = "MF",
      "GOCC" = "CC"
    )
    # If gene_clusters is NULL, perform enrichment analysis using enrichGO
    if (is.null(x = gene_clusters)) {
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
    # If gene_clusters is NULL, perform enrichment analysis using enrichKEGG
    if (is.null(x = gene_clusters)) {
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
    # Convert the gene IDs to readable format
    results <- clusterProfiler::setReadable(results,
      OrgDb = org_db,
      keyType = "ENTREZID"
    )
  }

  # Save the results as an RDS file
  if (!is.null(x = outdir)) {
    outdir <- sn_set_path(path = outdir)
    if (!is.null(x = prefix)) {
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
