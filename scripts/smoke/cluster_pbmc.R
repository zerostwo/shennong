suppressPackageStartupMessages({
  stopifnot(requireNamespace("devtools", quietly = TRUE))
  stopifnot(requireNamespace("Seurat", quietly = TRUE))
  stopifnot(requireNamespace("harmony", quietly = TRUE))
  stopifnot(requireNamespace("clusterProfiler", quietly = TRUE))
  stopifnot(requireNamespace("org.Hs.eg.db", quietly = TRUE))
  stopifnot(requireNamespace("lisi", quietly = TRUE))
})

devtools::load_all(".", quiet = TRUE)

prepare_sample <- function(dataset) {
  object <- sn_load_data(dataset = dataset)
  object$sample <- dataset
  object$study <- "pbmc_demo"
  object <- suppressWarnings(
    sn_filter_cells(
      object,
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
      plot = FALSE
    )
  )
  suppressWarnings(sn_filter_genes(object, min_cells = 3, plot = FALSE))
}

run_single_cluster <- function(object) {
  suppressWarnings(
    sn_run_cluster(
      object = object,
      normalization_method = "seurat",
      species = "human",
      nfeatures = 800,
      block_genes = NULL,
      npcs = 30,
      dims = 1:20,
      resolution = 0.5,
      verbose = FALSE
    )
  )
}

find_enrichment_result <- function(markers) {
  for (current_cluster in unique(as.character(markers$cluster))) {
    marker_subset <- markers |>
      dplyr::filter(.data$cluster == current_cluster) |>
      dplyr::arrange(dplyr::desc(.data$avg_log2FC))

    enrichment <- suppressWarnings(
      sn_enrich(
        marker_subset,
        analysis = "gsea",
        species = "human",
        database = "GOBP",
        gene_col = "gene",
        score_col = "avg_log2FC",
        pvalue_cutoff = 0.1
      )
    )

    enrichment_tbl <- as.data.frame(enrichment)
    if (nrow(enrichment_tbl) > 0) {
      return(enrichment_tbl)
    }
  }

  stop("No non-empty enrichment result was produced for the PBMC marker sets.")
}

pbmc1k <- prepare_sample("pbmc1k")
pbmc3k <- prepare_sample("pbmc3k")
pbmc1k_clustered <- run_single_cluster(pbmc1k)

stopifnot("seurat_clusters" %in% colnames(pbmc1k_clustered[[]]))
stopifnot("umap" %in% names(pbmc1k_clustered@reductions))

merged <- merge(
  x = pbmc1k,
  y = pbmc3k,
  add.cell.ids = c("pbmc1k", "pbmc3k")
)
merged$sample <- ifelse(grepl("^pbmc1k_", colnames(merged)), "pbmc1k", "pbmc3k")
merged$study <- "pbmc_demo"

integrated <- suppressWarnings(
  sn_run_cluster(
    object = merged,
    batch = "sample",
    normalization_method = "seurat",
    hvg_group_by = "sample",
    species = "human",
    nfeatures = 800,
    block_genes = NULL,
    npcs = 30,
    dims = 1:20,
    resolution = 0.5,
    verbose = FALSE
  )
)

stopifnot("harmony" %in% names(integrated@reductions))
stopifnot("seurat_clusters" %in% colnames(integrated[[]]))

composition <- sn_calculate_composition(
  integrated,
  group_by = "sample",
  variable = "seurat_clusters",
  min_cells = 0
)
stopifnot(nrow(composition) > 0)

lisi_tbl <- sn_calculate_lisi(
  integrated,
  reduction = "harmony",
  label = "sample"
)
stopifnot(nrow(lisi_tbl) == ncol(integrated))

integrated <- suppressWarnings(
  sn_find_de(
    integrated,
    analysis = "markers",
    group_by = "seurat_clusters",
    assay = "RNA",
    slot = "data",
    logfc_threshold = 0.25,
    min_pct = 0.2,
    store_name = "cluster_markers",
    return_object = TRUE,
    verbose = FALSE
  )
)
markers <- integrated@misc$de_results$cluster_markers$table
stopifnot(nrow(markers) > 0)

enrichment_tbl <- find_enrichment_result(markers)
stopifnot(nrow(enrichment_tbl) > 0)

cat("PBMC workflow smoke validation passed.\n")
