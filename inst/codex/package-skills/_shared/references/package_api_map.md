# Shennong Package API Map

This reference summarizes the exported Shennong surface so an agent can choose
the right `sn_*` entry point without scanning the whole package.

Core object rule:
- Current workflows are centered on Seurat objects.
- Prefer Shennong APIs over raw Seurat calls when the package already exposes
  the needed behavior.
- Prefer stored-result workflows over ad hoc `object@misc` access.

## Data and Import

- `sn_load_data()`: load packaged or remote example datasets; vectorized datasets such as `c("pbmc1k", "pbmc3k")` return a merged Seurat object for filtered matrices or a named list for raw matrices
- `sn_download_zenodo()`: download reusable files from public Zenodo records without a token, or from restricted/private records when a token is supplied
- `sn_list_10x_paths()`: discover 10x `outs/`, filtered/raw matrices, H5, or metrics paths
- `sn_read()`: import tabular, serialized (`qs`/`qs2`), and bioinformatics formats
- `sn_write()`: export tabular, serialized (`qs`/`qs2`), and supported omics formats
- `sn_upload_zenodo()`: upload reusable data files to Zenodo through `zen4R`, with a Shennong checksum/version manifest for later reuse
- `sn_add_data_from_anndata()`: add exported AnnData metadata and embeddings

Datasets:
- `pbmc_small`
- `pbmc_small_raw`
- `marker_genes`
- `hom_genes`
- `shennong_gene_annotations`
- `shennong_signature_catalog`

## Preprocessing and QC

- `sn_get_species()`: infer or retrieve species
- `sn_initialize_seurat_object()`: initialize a Seurat object, including single-path or multi-path 10x import from `sn_list_10x_paths()`
- `sn_standardize_gene_symbols()`: standardize gene symbols
- `sn_normalize_data()`: normalize with supported workflows
- `sn_score_cell_cycle()`: cell-cycle scoring
- `sn_filter_genes()`: gene filtering
- `sn_filter_cells()`: cell-level QC filtering
- `sn_assess_qc()`: summarize QC outcomes and before/after status
- `sn_find_doublets()`: doublet detection via `scDblFinder`
- `sn_remove_ambient_contamination()`: ambient RNA correction

## Clustering and Integration

- `sn_run_cluster()`: single-dataset clustering or batch integration; supports Seurat log-normalization, SCTransform, Harmony, Coralysis, Seurat CCA/RPCA, and pixi-managed scVI/scANVI integration through `integration_method`. SCTransform integration currently uses Harmony. For `integration_method = "scvi"` or `"scanvi"`, Shennong writes selected counts/metadata under `~/.shennong/runs/`, manages the shared scVI-family pixi project under `~/.shennong/pixi/scvi/`, imports the learned latent reduction, then continues Seurat neighbors/clustering/UMAP; scANVI requires `integration_control = list(label_by = ...)`. Use `integration_control = list(accelerator = "auto", mirror = "auto")` for CUDA/CPU auto-selection and Shennong-level mirror configuration. Rare-aware feature augmentation can combine `gini` and `local_markers`, with advanced thresholds kept in `rare_feature_control`. Use `hvg_features` to merge user-supplied marker genes into the final ScaleData/PCA feature set. Re-running on the returned object reuses matching stages by default; use `rerun_from` or `reuse = FALSE` for forced recompute. Leiden clustering auto-installs `leidenbase` unless `auto_install = FALSE`.
- `sn_run_scvi()` / `sn_run_scanvi()`: explicit wrappers for the corresponding `sn_run_cluster()` integration methods.
- `sn_transfer_labels()`: query-first reference `label_by` transfer wrapper. Defaults to Seurat anchors, can use `method = "coralysis"` for Coralysis `ReferenceMapping()` when the reference stores a trained Coralysis SingleCellExperiment, and supports semi-supervised scVI-family transfer with `method = "scanvi"` or `method = "scarches"`.
- `sn_simulate()`: method-based simulation entry point; currently supports `method = "scdesign3"` for Seurat or SingleCellExperiment inputs and returns Seurat, SingleCellExperiment, sparse counts, or the raw scDesign3 result.
- `sn_plot_heatmap()`: focused heatmap for user-selected genes, with cell-level and group-averaged modes, optional grouping/splitting, and selected-feature scaling.

## Python Runtime Helpers

- `sn_check_pixi()`: check whether a pixi executable is available.
- `sn_install_pixi()` / `sn_ensure_pixi()`: install or ensure the standalone pixi binary when Python backends need it.
- `sn_pixi_paths()`: inspect the `~/.shennong/pixi/` layout for scVI/scANVI and other Python method families.
- `sn_list_pixi_environments()` / `sn_pixi_config_path()`: discover bundled pixi configs under `inst/pixi/`.
- `sn_prepare_pixi_environment()` / `sn_call_pixi_environment()`: materialize a bundled config into `~/.shennong/pixi/<family>/` and run commands inside it.
- `sn_call_scvi()`, `sn_call_scanvi()`, `sn_call_scarches()`, `sn_call_scpoli()`, `sn_call_infercnvpy()`, `sn_call_cellphonedb()`, `sn_call_cell2location()`, `sn_call_tangram()`, `sn_call_squidpy()`, `sn_call_spatialdata()`, `sn_call_stlearn()`: environment-specific command-call helpers. `scanvi` shares the `scvi` environment; `scpoli` shares the `scarches` environment.
- `sn_run_scarches(object = ...)`, `sn_run_scpoli(object = ...)`, `sn_run_infercnvpy(object = ...)`, `sn_run_cellphonedb(object = ...)`, `sn_run_cell2location(object = ...)`, `sn_run_tangram(object = ...)`, `sn_run_squidpy(object = ...)`, `sn_run_spatialdata(object = ...)`, `sn_run_stlearn(object = ...)`: object-level Python wrappers. They export Seurat input under `~/.shennong/runs/`, run family-local scripts from `inst/pixi/<family>/scripts/`, import cell-level metadata/reductions when produced, and record manifests under `object@misc`.
- The same `sn_run_*()` wrappers still support command mode when `object` is omitted, for example `sn_run_tangram(args = c("-c", "import tangram"))`.
- `sn_detect_accelerator()`: detect CUDA-capable NVIDIA GPUs and report CPU fallback status.
- `sn_configure_pixi_mirror()`: write Shennong-level pixi mirror configuration for default, China, TUNA, USTC, or BFSU sources.

## Diagnostics and Benchmarking

- `sn_assess_integration()`: aggregate integration quality summary
- `sn_sweep_cluster_resolution()`: compare candidate Seurat resolutions with cluster-quality diagnostics and a recommended resolution
- `sn_detect_rare_cells()`: rare-cell diagnostics
- `sn_calculate_lisi()`: batch_by mixing and label_by mixing scores
- `sn_calculate_silhouette()`: silhouette widths
- `sn_calculate_graph_connectivity()`: graph connectivity
- `sn_calculate_pcr_batch()`: principal-component batch_by effect score
- `sn_calculate_variance_explained()`: rank metadata variables by embedding variance explained
- `sn_calculate_clustering_agreement()`: clustering versus reference agreement
- `sn_calculate_isolated_label_score()`: isolated-label preservation
- `sn_calculate_cluster_purity()`: cluster purity
- `sn_calculate_cluster_entropy()`: cluster_by entropy
- `sn_identify_challenging_groups()`: rare or difficult groups
- `sn_calculate_rogue()`: cluster purity / heterogeneity score

## Annotation, Markers, and Pathways

- `sn_run_celltypist()`: external CellTypist-based annotation
- `sn_find_de()`: markers, contrasts, and pseudobulk DE
- `sn_enrich()`: ORA or GSEA from vectors, tables, or stored DE
- `sn_list_signatures()`: list bundled signatures
- `sn_get_signatures()`: retrieve signatures by path or category
- `sn_add_signature()`: add a signature to the editable registry
- `sn_update_signature()`: update a signature in the editable registry
- `sn_delete_signature()`: delete a signature from the editable registry

## Composition and Comparative Analysis

- `sn_calculate_composition()`: grouped counts and proportions
- `sn_calculate_roe()`: observed-over-expected enrichment for categorical composition tables
- `sn_compare_composition()`: compare sample-level composition between groups
- `sn_run_milo()`: neighborhood differential abundance with miloR
- `sn_store_milo()`: persist milo results
- `sn_get_milo_result()`: retrieve milo results
- `sn_deconvolve_bulk()`: bulk deconvolution from single-cell reference data
- `sn_set_cibersortx_credentials()`: store CIBERSORTx credentials
- `sn_store_deconvolution()`: persist deconvolution results
- `sn_get_deconvolution_result()`: retrieve deconvolution results
- `sn_run_cell_communication()`: run CellChat, NicheNet, or LIANA communication inference
- `sn_store_cell_communication()`: persist communication results
- `sn_get_cell_communication_result()`: retrieve communication results
- `sn_run_regulatory_activity()`: infer DoRothEA TF activity or PROGENy pathway activity with decoupleR
- `sn_store_regulatory_activity()`: persist regulatory activity results
- `sn_get_regulatory_activity_result()`: retrieve regulatory activity results

## Interpretation and Reporting

- `sn_list_results()`: list stored DE, enrichment, milo, deconvolution, communication, regulatory activity, and interpretation results
- `sn_get_de_result()`: retrieve stored DE
- `sn_get_enrichment_result()`: retrieve stored enrichment
- `sn_get_interpretation_result()`: retrieve stored interpretation
- `sn_store_enrichment()`: store enrichment results for downstream reuse
- `sn_prepare_annotation_evidence()`: prepare cluster annotation evidence
- `sn_prepare_de_evidence()`: prepare DE evidence
- `sn_prepare_enrichment_evidence()`: prepare enrichment evidence
- `sn_prepare_results_evidence()`: prepare manuscript-style results evidence
- `sn_build_prompt()`: turn structured evidence into an LLM or human prompt bundle
- `sn_test_llm_provider()`: test the provider path
- `sn_make_ellmer_provider()`: create an `ellmer`-backed provider
- `sn_run_llm()`: call a provider on a message list
- `sn_interpret_annotation()`: annotation interpretation, including agentic mode
- `sn_interpret_de()`: DE interpretation
- `sn_interpret_enrichment()`: pathway interpretation
- `sn_write_results()`: manuscript-style results text
- `sn_write_figure_legend()`: figure legends
- `sn_write_presentation_summary()`: presentation-style summaries

## Visualization

- `sn_plot_dim()`: embedding plot for categorical labels
- `sn_plot_feature()`: feature plot or density-style feature map with rasterized points that preserve `pt_size`
- `sn_plot_heatmap()`: selected-gene heatmap in cell-level or group-averaged mode
- `sn_plot_violin()`: violin plot
- `sn_plot_dot()`: dot plot
- `sn_plot_boxplot()`: box plot
- `sn_plot_barplot()`: bar plot
- `sn_plot_composition()`: composition bar charts
- `sn_plot_milo()`: milo differential-abundance visualization
- `sn_list_palettes()`: list available palettes
- `sn_get_palette()`: resolve palette colors

## Utilities and Project Setup

- `sn_check_version()`: check package version
- `sn_install_shennong()`: install Shennong
- `sn_check_file()`: verify file paths
- `sn_get_codex_skill_path()`: locate packaged skills and template assets
- `sn_list_dependencies()`: list required and recommended package dependencies
- `sn_install_dependencies()`: install missing dependencies
- `sn_initialize_project()`: initialize a governed analysis project
- `sn_install_codex_skill()`: install shipped Codex skills
- `sn_set_path()`: create a directory path if needed

## Selection shortcuts

- If the task starts from counts or a path:
  use `sn_initialize_seurat_object()`
- If the task is clustering or integration:
  use `sn_run_cluster()`
- If the task is markers or contrasts:
  use `sn_find_de()`
- If the task is pathways:
  use `sn_enrich()`
- If the task is result reuse:
  use `sn_list_results()` plus `sn_get_*_result()`
- If the task is annotation or interpretation:
  use `sn_prepare_*_evidence()` and `sn_interpret_*()`
- If the task is project bootstrap:
  use `sn_initialize_project()`
