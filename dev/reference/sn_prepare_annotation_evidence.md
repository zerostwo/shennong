# Prepare cluster-annotation evidence from a Seurat object

Prepare cluster-annotation evidence from a Seurat object

## Usage

``` r
sn_prepare_annotation_evidence(
  object,
  de_result_id = NULL,
  cluster_by = NULL,
  n_markers = 10,
  marker_selection = c("specific", "top"),
  enrichment_result_id = NULL,
  n_terms = 5,
  enrichment_selection = c("specific", "top"),
  include_qc = TRUE,
  reduction = "umap",
  n_neighbor_clusters = 3
)
```

## Arguments

- object:

  A `Seurat` object.

- de_result_id:

  Optional stored marker-result name in the canonical Shennong result
  registry. When omitted, Shennong prefers `"default"`, then a single
  available result, and otherwise the most recent marker result.

- cluster_by:

  Metadata column containing cluster labels.

- n_markers:

  Number of top markers to retain per cluster.

- marker_selection:

  How to choose marker genes for annotation evidence: `"specific"`
  prefers genes that are relatively unique to one cluster, while `"top"`
  keeps the raw top-ranked genes.

- enrichment_result_id:

  Optional stored enrichment result used to add cluster-level functional
  evidence to the annotation prompt.

- n_terms:

  Number of enrichment terms to retain per cluster when
  `enrichment_result_id` is supplied.

- enrichment_selection:

  How to choose pathway/function terms for annotation evidence:
  `"specific"` prefers terms concentrated in fewer clusters, while
  `"top"` keeps the raw top-ranked terms.

- include_qc:

  Logical; whether to attach cluster-level QC summaries such as
  mitochondrial burden, failed-QC fractions, and doublet fractions when
  available in metadata.

- reduction:

  Optional dimensional reduction name used to summarize cluster
  neighborhood geometry, for example `"umap"`. Use `NULL` to disable
  geometry evidence.

- n_neighbor_clusters:

  Number of nearest clusters to report from the reduction centroid
  distances. Canonical lineage heuristic hints derived from known marker
  programs are included automatically when the required genes are
  present.

## Value

A structured list ready for prompt construction.

## Examples

``` r
if (requireNamespace("Seurat", quietly = TRUE)) {
  counts <- matrix(rpois(20 * 24, lambda = 1), nrow = 20, ncol = 24)
  rownames(counts) <- c(
    paste0("GENE", 1:14),
    "CD3D", "CD3E", "TRAC", "MS4A1", "CD79A", "HLA-DRA"
  )
  colnames(counts) <- paste0("cell", 1:24)
  counts[c("CD3D", "CD3E", "TRAC"), 1:12] <-
    counts[c("CD3D", "CD3E", "TRAC"), 1:12] + 20
  counts[c("MS4A1", "CD79A", "HLA-DRA"), 13:24] <-
    counts[c("MS4A1", "CD79A", "HLA-DRA"), 13:24] + 20
  obj <- sn_initialize_seurat_object(counts, species = "human")
  obj$cell_type <- rep(c("Tcell", "Bcell"), each = 12)
  Seurat::Idents(obj) <- obj$cell_type
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- sn_find_de(obj, analysis = "markers", group_by = "cell_type",
    layer = "data", min_pct = 0, logfc_threshold = 0,
    result_id = "celltype_markers", return_object = TRUE, verbose = FALSE
  )
  evidence <- sn_prepare_annotation_evidence(
    obj,
    de_result_id = "celltype_markers",
    cluster_by = "cell_type"
  )
  names(evidence)
}
```
