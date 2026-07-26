suppressPackageStartupMessages({
  stopifnot(requireNamespace("devtools", quietly = TRUE))
  stopifnot(requireNamespace("Seurat", quietly = TRUE))
  stopifnot(requireNamespace("harmony", quietly = TRUE))
  stopifnot(requireNamespace("clusterProfiler", quietly = TRUE))
  stopifnot(requireNamespace("org.Hs.eg.db", quietly = TRUE))
  stopifnot(requireNamespace("lisi", quietly = TRUE))
  stopifnot(requireNamespace("qs2", quietly = TRUE))
})

devtools::load_all(".", quiet = TRUE)

prepare_sample <- function(object, sample_id) {
  object <- subset(object, cells = colnames(object)[as.character(object$real_sample) == sample_id])
  object$sample <- as.character(object$real_sample)
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

fixture_path <- file.path(
  Sys.getenv("SHENNONG_REAL_DATA_DIR"),
  "single-cell",
  "kotliarov_pbmc.qs2"
)
stopifnot(file.exists(fixture_path))
fixture <- qs2::qs_read(fixture_path)
stopifnot(inherits(fixture, "Seurat"), "real_sample" %in% colnames(fixture[[]]))
sample_ids <- unique(as.character(fixture$real_sample))
stopifnot(length(sample_ids) >= 2L)

sample_a <- prepare_sample(fixture, sample_ids[[1]])
sample_b <- prepare_sample(fixture, sample_ids[[2]])
sample_a_clustered <- run_single_cluster(sample_a)

stopifnot("seurat_clusters" %in% colnames(sample_a_clustered[[]]))
stopifnot("umap" %in% names(sample_a_clustered@reductions))

merged <- merge(
  x = sample_a,
  y = sample_b,
  add.cell.ids = c("sample_a", "sample_b")
)
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
