# Input initialization, QC, feature identity and layer storage

Use only the stages needed by the requested analysis. Each section is an API
reference, not a requirement to run all listed methods.

## 1

Start with data access and preprocessing:
   initialize or load the object, infer species if needed, and run QC with
   `sn_filter_cells()` and `sn_filter_genes()`.
   `sn_filter_genes()` can combine `min_cells` with bundled GENCODE-based
   `gene_class` or exact `gene_type` filtering for human and mouse workflows.
   Preserve BPCells `IterableMatrix` counts by passing them directly to
   `sn_initialize_seurat_object()`. When finding doublets on a BPCells-backed
   object, require a donor/capture `group_by` column and default to
   `n_workers = 1`, which materializes one sample-sized sparse matrix at a time.
   For in-memory `dgCMatrix` input, scDblFinder native clustering is the default
   and supported default-call shapes use the ShennongOpt `scdblfinder` fast
   path; grouped and non-default calls remain on the upstream path. Use
   `cluster_backend = "shennong"` only when the caller explicitly wants
   `sn_run_cluster()` assignments. `sn_score_programs()` uses UCell by default
   for per-cell programs, and LISI metrics use the ShennongOpt patch when its
   contract is met.
   Use `sn_set_layer_backend()` when selected Seurat layers must move between
   BPCells and in-memory `dgCMatrix` storage; materialize only the layers needed
   by an in-memory-only operation.
