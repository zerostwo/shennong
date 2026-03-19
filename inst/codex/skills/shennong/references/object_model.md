# Object Model

## Primary Analysis Object

Shennong primarily operates on Seurat objects.

Important places to look:

- `object@meta.data`: cell-level metadata such as `sample`, `study`,
  `seurat_clusters`, QC columns, and annotations
- assays and layers: count-like analyses may use explicit `assay` and `layer`
- `object@reductions`: PCA, UMAP, Harmony, and other embeddings
- `object@graphs`: neighbor graphs used by clustering
- `object@commands`: recorded workflow steps for reproducibility
- `object@misc`: package-specific stored results such as DE outputs

## Stored Differential Expression Results

`sn_find_de(..., return_object = TRUE)` stores results under:

`object@misc$de_results[[store_name]]`

Each stored entry includes:

- `table`: the result table
- `analysis`: `markers`, `contrast`, or `pseudobulk`
- `group_by` and `group_col`: grouping metadata
- `ident_1`, `ident_2`: comparison labels for direct contrasts
- `subset_by`: optional stratifying metadata column
- `rank_col`: ranking column used for top markers or GSEA
- `p_col`: adjusted p-value column when available

This is the structure that `sn_plot_dot(features = "top_markers")` reuses.

## Count-Layer Handling

Several Shennong workflows accept explicit `assay` and `layer` arguments. Use
them when counts live in a non-default layer such as
`decontaminated_counts`.

Examples:

- `sn_run_cluster(assay = "RNA", layer = "decontaminated_counts")`
- `sn_find_doublets(assay = "RNA", layer = "decontaminated_counts")`
- `sn_calculate_rogue(assay = "RNA", layer = "decontaminated_counts")`

## Example Data

`sn_load_data()` returns:

- a Seurat object for filtered example data
- a sparse count matrix for raw example data

The packaged example registry currently focuses on PBMC data sets such as
`pbmc1k` and `pbmc3k`.
