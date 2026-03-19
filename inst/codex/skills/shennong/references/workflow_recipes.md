# Workflow Recipes

## 1. Quick Example Workflow

```r
library(Shennong)

object <- sn_load_data("pbmc3k")
object <- sn_filter_cells(
  object,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  plot = FALSE
)
object <- sn_filter_genes(object, plot = FALSE)
object <- sn_run_cluster(object, normalization_method = "seurat")
object <- sn_find_de(
  object,
  analysis = "markers",
  group_by = "seurat_clusters",
  store_name = "cluster_markers",
  return_object = TRUE
)
sn_plot_dim(object, reduction = "umap", group_by = "seurat_clusters")
sn_plot_dot(object, features = "top_markers", de_name = "cluster_markers", n = 3)
```

## 2. Compare Two Conditions Within Each Cell Type

```r
object <- sn_find_de(
  object,
  analysis = "contrast",
  ident_1 = "treated",
  ident_2 = "control",
  group_by = "condition",
  subset_by = "cell_type",
  return_object = FALSE
)
```

Use this when the biological question is:
"within each cell type, which genes differ between two conditions?"

## 3. Run Pseudobulk Differential Expression

```r
pb_result <- sn_find_de(
  object,
  analysis = "pseudobulk",
  ident_1 = "treated",
  ident_2 = "control",
  group_by = "condition",
  subset_by = "cell_type",
  sample_col = "sample",
  pseudobulk_method = "edgeR",
  return_object = FALSE
)
```

Use pseudobulk when sample-level replication matters and cell-level tests are
not appropriate.

## 4. Compare Integration Before and After Harmony

```r
merged <- merge(obj1, y = obj2, add.cell.ids = c("s1", "s2"))
merged$sample <- ifelse(grepl("^s1_", colnames(merged)), "s1", "s2")

before <- sn_run_cluster(merged, normalization_method = "seurat")
after <- sn_run_cluster(
  merged,
  batch = "sample",
  normalization_method = "seurat",
  hvg_group_by = "sample"
)

sn_calculate_lisi(before, reduction = "pca", label = "sample")
sn_calculate_lisi(after, reduction = "harmony", label = "sample")
```

## 5. ORA and GSEA

ORA from a simple gene set:

```r
sn_enrich(
  c("CD3D", "CD3E", "TRAC", "LCK"),
  analysis = "ora",
  species = "human",
  database = "GOBP"
)
```

GSEA from a ranked marker table:

```r
gsea_result <- sn_enrich(
  marker_table,
  analysis = "gsea",
  species = "human",
  database = "GOBP",
  gene_col = "gene",
  score_col = "avg_log2FC"
)
```

## 6. Troubleshooting Patterns

- If a workflow uses a custom counts layer, pass both `assay` and `layer`.
- If Harmony integration fails, confirm that the `batch` column exists in
  `object@meta.data`.
- If `sn_plot_dot(features = "top_markers")` fails, confirm that
  `sn_find_de(..., return_object = TRUE)` was run first.
- If `sn_enrich(analysis = "gsea")` fails, ensure the ranking input is a named
  numeric vector or a data frame with both gene and score columns.
