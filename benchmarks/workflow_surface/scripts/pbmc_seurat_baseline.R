library(Seurat)

pbmc <- qs2::qs_read(file.path(
  Sys.getenv("SHENNONG_REAL_DATA_DIR"),
  "single-cell",
  "kotliarov_pbmc.qs2"
))

pbmc <- subset(
  pbmc,
  subset = nFeature_RNA > 200 & nCount_RNA > 500 & percent.mt < 10
)

pbmc <- NormalizeData(
  pbmc,
  verbose = FALSE
)

pbmc <- FindVariableFeatures(
  pbmc,
  verbose = FALSE
)

pbmc <- ScaleData(
  pbmc,
  verbose = FALSE
)

pbmc <- RunPCA(
  pbmc,
  verbose = FALSE
)

pbmc <- FindNeighbors(
  pbmc,
  dims = 1:20,
  verbose = FALSE
)

pbmc <- FindClusters(
  pbmc,
  resolution = 0.6,
  verbose = FALSE
)

pbmc <- RunUMAP(
  pbmc,
  dims = 1:20,
  verbose = FALSE
)

markers <- FindAllMarkers(
  pbmc,
  only.pos = TRUE,
  verbose = FALSE
)

saveRDS(
  list(
    object = pbmc,
    markers = markers
  ),
  file = "baseline_pbmc_results.rds"
)
