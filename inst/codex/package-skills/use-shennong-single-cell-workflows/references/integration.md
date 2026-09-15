# Clustering, integration, CITE-seq and integration assessment

Use only the stages needed by the requested analysis. Each section is an API
reference, not a requirement to run all listed methods.

## 2

Run clustering or batch integration with `sn_run_cluster()`. Use
   `hvg_features = c(...)` when the user has known marker genes that should be
   forced into the backend feature set, use
   `rare_feature_method = "gini"` or `"local_markers"` when Shennong should
   automatically add rare-aware genes, use `block_genes = c(...)` to exclude
   bundled signature queries such as `cellCycle.G2M`, `ribo`, or `mito` and/or
   custom gene symbols from internally selected HVGs in log-normalization or
   SCTransform workflows, and use
	   `integration_method = "harmony"`, `"coralysis"`, `"seurat_cca"`,
	   `"seurat_rpca"`, `"scvi"`, `"scanvi"`, `"scpoli"`, or `"bbknn"` when a
	   specific batch-integration backend is requested.
	   Pass a vector such as `c("unintegrated", "harmony", "coralysis")` when
	   methods should share preprocessing and coexist in one object. Key
	   `integration_control` by method, then inspect
	   `object@misc$integration_comparison$grid` and `$results` for the stored
	   preprocessing/embedding IDs, reduction, graph, cluster, UMAP, and optional
	   t-SNE names. Multiple values for scalar controls such as `nfeatures`,
	   `npcs`, `resolution`, clustering algorithm, rare-feature count, or Harmony
	   theta form a conditional Cartesian grid; natural vectors such as `dims`
	   remain intact. UMAP is the default and t-SNE requires
	   `run_tsne = TRUE`.
	   For large grids, set `checkpoint_dir`; matching calls resume completed
	   combinations by default. Read `$performance` for time and memory
	   provenance, and call `sn_get_integration_control_template()` for the complete
	   control template of each backend.
	   scVI/scANVI/scPoli honor
	   the requested `assay`/`layer`; BBKNN derives PCA from the same selected
	   layer and supplies the graph used for clustering and UMAP. Python runners
	   retain complete expression matrices as sparse objects; only bounded
	   neural-network minibatches and low-dimensional outputs may be dense. For scVI/scANVI,
	   Shennong manages a
	   shared pixi scverse project under `~/.shennong/pixi/scvi/`, uses a unique
	   package-owned temporary run directory unless an explicit directory is
	   supplied, and imports the latent reduction back
	   into Seurat; scANVI requires `integration_control = list(label_by = ...)`.
	   For object-level Python wrappers, use the explicit `keep_run_dir` and
	   `max_artifact_import_gb` arguments to control retained diagnostics and the
	   bounded validation/materialization budget. With cleanup requested, a
	   supplied output path is only a parent for a marked child and is preserved.
   scPoli uses the shared `scarches` pixi family and accepts optional
   `integration_control = list(label_by = ...)` for prototype supervision.
   Coralysis stores the trained reference SingleCellExperiment by default for
   label transfer; use `integration_control = list(store_sce = FALSE)` only for
   clustering-only runs.
   Use `sn_get_pixi_paths()` when users ask where Python environments live, use
   `sn_list_pixi_environments()` and `sn_get_pixi_config_path()` to inspect bundled
   configs under `inst/pixi/`, and pass
   `integration_control = list(accelerator = "auto", mirror = "auto")` when
   GPU/CPU selection and China-friendly mirror configuration should be handled
   by Shennong.
   Default to tested exact Pixi `0.69.0`; overrides must name an immutable
   release, mutable `latest` is rejected, and custom download URLs require an
   explicit SHA-256. Temporary successful runs are cleaned and failure metadata
   is sanitized unless retention is explicit.
   Re-running `sn_run_cluster()` on its own output reuses stages only when the
   selected layer content and relevant metadata digests still match; use
   `rerun_from = "integration"` or `reuse = FALSE` when a stage
	   must be forced to recompute. Leiden clustering auto-installs `leidenbase` by
	   default unless `auto_install = FALSE`.
	   Use `umap_control = list(n.neighbors = ..., min.dist = ..., spread = ...)`
	   with `rerun_from = "umap"` when only the two-dimensional embedding needs
	   retuning.
	   Use `normalization_method = "sctransform"` with `batch = ...` only when an
	   SCTransform-normalized Harmony integration is requested.
	   For CITE-seq objects with paired RNA and ADT assays, use
	   `modality = "cite_seq"` plus `multimodal_method = "wnn"` for Seurat
	   weighted nearest-neighbor clustering on `weighted.nn` / `wsnn`,
	   `multimodal_method = "totalvi"` for scvi-tools totalVI RNA+ADT latent
	   integration, or `multimodal_method = "coralysis"` to run native
	   Coralysis on the ADT protein assay. Use `multimodal_method = "mmochi"`
	   when ADT alignment should be driven by MMoCHi landmark registration
	   before protein-only clustering; it can run with `batch = NULL` for a
	   single CITE-seq sample.

## 3

Assess integration quality or cluster_by structure with
   `sn_assess_integration()`, `sn_calculate_lisi()`,
   `sn_compare_integrations()` for a scib-metrics benchmark of the native
   reductions retained in a multi-method/grid object (one score per unique
   embedding, with a matching unintegrated baseline per preprocessing group),
   `sn_calculate_variance_explained()`,
   `sn_calculate_isolated_label_score()`, or
   `sn_identify_challenging_groups()` when sample mixing, rare groups,
   isolated labels, or difficult-to-separate populations matter.
