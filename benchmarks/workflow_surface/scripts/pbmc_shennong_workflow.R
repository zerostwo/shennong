library(Shennong)

data("pbmc_small", package = "Shennong")

pbmc <- pbmc_small

pbmc <- sn_filter_cells(
  pbmc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  plot = FALSE
)

pbmc <- sn_filter_genes(
  pbmc,
  min_cells = 3,
  plot = FALSE
)

pbmc <- sn_run_cluster(
  pbmc,
  normalization_method = "seurat",
  resolution = 0.6
)

pbmc <- sn_find_de(
  pbmc,
  analysis = "markers",
  group_by = "seurat_clusters",
  layer = "data",
  store_name = "cluster_markers",
  return_object = TRUE,
  verbose = FALSE
)

result_index <- sn_list_results(pbmc)

saveRDS(
  list(
    object = pbmc,
    result_index = result_index
  ),
  file = "shennong_pbmc_results.rds"
)
