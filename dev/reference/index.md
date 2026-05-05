# Package index

## Data and Import

- [`marker_genes`](https://songqi.org/shennong/dev/reference/marker_genes.md)
  : Pan-Immune CellTypist Metadata
- [`hom_genes`](https://songqi.org/shennong/dev/reference/hom_genes.md)
  : Human-Mouse Homologous Genes Table
- [`pbmc_small`](https://songqi.org/shennong/dev/reference/pbmc_small.md)
  : Small Built-In PBMC Seurat Object
- [`pbmc_small_raw`](https://songqi.org/shennong/dev/reference/pbmc_small_raw.md)
  : Small Built-In PBMC Raw Counts Matrix
- [`shennong_gene_annotations`](https://songqi.org/shennong/dev/reference/shennong_gene_annotations.md)
  : Bundled Shennong Gene Annotations
- [`shennong_signature_catalog`](https://songqi.org/shennong/dev/reference/shennong_signature_catalog.md)
  : Bundled Shennong Signature Catalog
- [`sn_load_data()`](https://songqi.org/shennong/dev/reference/sn_load_data.md)
  : Load example datasets from Zenodo
- [`sn_download_zenodo()`](https://songqi.org/shennong/dev/reference/sn_download_zenodo.md)
  : Download files from a Zenodo record
- [`sn_list_10x_paths()`](https://songqi.org/shennong/dev/reference/sn_list_10x_paths.md)
  : List detected 10x Genomics output paths under a root folder
- [`sn_read()`](https://songqi.org/shennong/dev/reference/sn_read.md)
  [`.import.rio_bpcells()`](https://songqi.org/shennong/dev/reference/sn_read.md)
  [`.import.rio_10x()`](https://songqi.org/shennong/dev/reference/sn_read.md)
  [`.import.rio_10x_spatial()`](https://songqi.org/shennong/dev/reference/sn_read.md)
  [`.import.rio_starsolo()`](https://songqi.org/shennong/dev/reference/sn_read.md)
  [`.import.rio_h5ad()`](https://songqi.org/shennong/dev/reference/sn_read.md)
  [`.import.rio_h5()`](https://songqi.org/shennong/dev/reference/sn_read.md)
  [`.import.rio_qs()`](https://songqi.org/shennong/dev/reference/sn_read.md)
  [`.import.rio_qs2()`](https://songqi.org/shennong/dev/reference/sn_read.md)
  [`.import.rio_gmt()`](https://songqi.org/shennong/dev/reference/sn_read.md)
  : Read tabular and bioinformatics file formats
- [`sn_write()`](https://songqi.org/shennong/dev/reference/sn_write.md)
  [`.export.rio_bpcells()`](https://songqi.org/shennong/dev/reference/sn_write.md)
  [`.export.rio_h5ad()`](https://songqi.org/shennong/dev/reference/sn_write.md)
  [`.export.rio_h5()`](https://songqi.org/shennong/dev/reference/sn_write.md)
  [`.export.rio_qs()`](https://songqi.org/shennong/dev/reference/sn_write.md)
  [`.export.rio_qs2()`](https://songqi.org/shennong/dev/reference/sn_write.md)
  : Write tabular and bioinformatics file formats
- [`sn_upload_zenodo()`](https://songqi.org/shennong/dev/reference/sn_upload_zenodo.md)
  : Upload reusable data files to Zenodo
- [`sn_add_data_from_anndata()`](https://songqi.org/shennong/dev/reference/sn_add_data_from_anndata.md)
  : Add metadata and embeddings exported from AnnData

## Preprocessing and QC

- [`sn_get_species()`](https://songqi.org/shennong/dev/reference/sn_get_species.md)
  : Retrieve or infer species information
- [`sn_initialize_seurat_object()`](https://songqi.org/shennong/dev/reference/sn_initialize_seurat_object.md)
  : Initialize a Seurat object with optional QC metrics
- [`sn_standardize_gene_symbols()`](https://songqi.org/shennong/dev/reference/sn_standardize_gene_symbols.md)
  : Standardize gene symbols in a vector, count matrix, or Seurat object
- [`sn_normalize_data()`](https://songqi.org/shennong/dev/reference/sn_normalize_data.md)
  : Normalize data in a Seurat object
- [`sn_score_cell_cycle()`](https://songqi.org/shennong/dev/reference/sn_score_cell_cycle.md)
  : Score Cell Cycle Phases
- [`sn_filter_genes()`](https://songqi.org/shennong/dev/reference/sn_filter_genes.md)
  : Filter genes based on the number of expressing cells
- [`sn_filter_cells()`](https://songqi.org/shennong/dev/reference/sn_filter_cells.md)
  : Filter cells in a Seurat object based on QC metrics
- [`sn_assess_qc()`](https://songqi.org/shennong/dev/reference/sn_assess_qc.md)
  : Assess overall QC status and before/after filtering outcomes
- [`sn_find_doublets()`](https://songqi.org/shennong/dev/reference/sn_find_doublets.md)
  : Find doublets using scDblFinder
- [`sn_remove_ambient_contamination()`](https://songqi.org/shennong/dev/reference/sn_remove_ambient_contamination.md)
  : Remove ambient RNA contamination from counts

## Clustering and Integration

- [`sn_run_cluster()`](https://songqi.org/shennong/dev/reference/sn_run_cluster.md)
  : Run clustering for a single dataset or batch integration workflow
- [`sn_run_scvi()`](https://songqi.org/shennong/dev/reference/sn_run_scvi.md)
  [`sn_run_scanvi()`](https://songqi.org/shennong/dev/reference/sn_run_scvi.md)
  : Run scVI or scANVI integration through Shennong
- [`sn_transfer_labels()`](https://songqi.org/shennong/dev/reference/sn_transfer_labels.md)
  : Transfer labels from a Seurat reference to a query object
- [`sn_simulate()`](https://songqi.org/shennong/dev/reference/sn_simulate.md)
  : Simulate single-cell counts with scDesign3
- [`sn_simulate_scdesign3()`](https://songqi.org/shennong/dev/reference/sn_simulate_scdesign3.md)
  : Simulate single-cell counts with scDesign3
- [`sn_sweep_cluster_resolution()`](https://songqi.org/shennong/dev/reference/sn_sweep_cluster_resolution.md)
  : Sweep clustering resolutions and summarize cluster-quality
  diagnostics

## Diagnostics and Benchmarking

- [`sn_assess_integration()`](https://songqi.org/shennong/dev/reference/sn_assess_integration.md)
  : Assess integration quality across multiple metrics
- [`sn_detect_rare_cells()`](https://songqi.org/shennong/dev/reference/sn_detect_rare_cells.md)
  : Detect rare cells with native or optional rare-cell backends
- [`sn_calculate_lisi()`](https://songqi.org/shennong/dev/reference/sn_calculate_lisi.md)
  : Calculate LISI scores from a Seurat embedding
- [`sn_calculate_silhouette()`](https://songqi.org/shennong/dev/reference/sn_calculate_silhouette.md)
  : Calculate silhouette widths from a Seurat embedding
- [`sn_calculate_graph_connectivity()`](https://songqi.org/shennong/dev/reference/sn_calculate_graph_connectivity.md)
  : Calculate graph connectivity for a grouping label
- [`sn_calculate_pcr_batch()`](https://songqi.org/shennong/dev/reference/sn_calculate_pcr_batch.md)
  : Calculate PCR batch_by effect scores
- [`sn_calculate_variance_explained()`](https://songqi.org/shennong/dev/reference/sn_calculate_variance_explained.md)
  : Rank metadata variables by embedding variance explained
- [`sn_calculate_clustering_agreement()`](https://songqi.org/shennong/dev/reference/sn_calculate_clustering_agreement.md)
  : Calculate agreement between clusters and reference labels
- [`sn_calculate_isolated_label_score()`](https://songqi.org/shennong/dev/reference/sn_calculate_isolated_label_score.md)
  : Calculate isolated-label preservation scores
- [`sn_calculate_cluster_purity()`](https://songqi.org/shennong/dev/reference/sn_calculate_cluster_purity.md)
  : Calculate cluster purity against a reference label
- [`sn_calculate_cluster_entropy()`](https://songqi.org/shennong/dev/reference/sn_calculate_cluster_entropy.md)
  : Calculate cluster_by entropy for a categorical label
- [`sn_identify_challenging_groups()`](https://songqi.org/shennong/dev/reference/sn_identify_challenging_groups.md)
  : Identify rare or difficult-to-separate groups
- [`sn_calculate_rogue()`](https://songqi.org/shennong/dev/reference/sn_calculate_rogue.md)
  : Calculate ROGUE score for Seurat Object

## Annotation, Markers, and Pathways

- [`sn_run_celltypist()`](https://songqi.org/shennong/dev/reference/sn_run_celltypist.md)
  : Run CellTypist for automated cell type annotation
- [`sn_find_de()`](https://songqi.org/shennong/dev/reference/sn_find_de.md)
  : Run differential expression analysis on a Seurat object
- [`sn_enrich()`](https://songqi.org/shennong/dev/reference/sn_enrich.md)
  : Run gene set enrichment analysis
- [`sn_list_signatures()`](https://songqi.org/shennong/dev/reference/sn_list_signatures.md)
  : List bundled Shennong signatures
- [`sn_get_signatures()`](https://songqi.org/shennong/dev/reference/sn_get_signatures.md)
  : Retrieve bundled Shennong signature genes by category or tree path
- [`sn_add_signature()`](https://songqi.org/shennong/dev/reference/sn_add_signature.md)
  : Add a signature to the editable source registry
- [`sn_update_signature()`](https://songqi.org/shennong/dev/reference/sn_update_signature.md)
  : Update a signature in the editable source registry
- [`sn_delete_signature()`](https://songqi.org/shennong/dev/reference/sn_delete_signature.md)
  : Delete a signature from the editable source registry

## Composition and Comparative Analysis

- [`sn_calculate_composition()`](https://songqi.org/shennong/dev/reference/sn_calculate_composition.md)
  : Calculate composition proportions
- [`sn_calculate_roe()`](https://songqi.org/shennong/dev/reference/sn_calculate_roe.md)
  : Calculate observed-over-expected enrichment
- [`sn_compare_composition()`](https://songqi.org/shennong/dev/reference/sn_compare_composition.md)
  : Compare sample-level composition between groups
- [`sn_run_milo()`](https://songqi.org/shennong/dev/reference/sn_run_milo.md)
  : Run neighborhood differential abundance testing with miloR
- [`sn_store_milo()`](https://songqi.org/shennong/dev/reference/sn_store_milo.md)
  : Store a miloR differential-abundance result on a Seurat object
- [`sn_get_milo_result()`](https://songqi.org/shennong/dev/reference/sn_get_milo_result.md)
  : Retrieve a stored miloR result from a Seurat object
- [`sn_deconvolve_bulk()`](https://songqi.org/shennong/dev/reference/sn_deconvolve_bulk.md)
  : Run bulk RNA-seq deconvolution with single-cell references
- [`sn_set_cibersortx_credentials()`](https://songqi.org/shennong/dev/reference/sn_set_cibersortx_credentials.md)
  : Store local CIBERSORTx credentials for container execution
- [`sn_store_deconvolution()`](https://songqi.org/shennong/dev/reference/sn_store_deconvolution.md)
  : Store a deconvolution result on a Seurat object
- [`sn_get_deconvolution_result()`](https://songqi.org/shennong/dev/reference/sn_get_deconvolution_result.md)
  : Retrieve a stored deconvolution result from a Seurat object

## Communication and Regulatory Activity

- [`sn_run_cell_communication()`](https://songqi.org/shennong/dev/reference/sn_run_cell_communication.md)
  : Run cell-cell communication inference
- [`sn_store_cell_communication()`](https://songqi.org/shennong/dev/reference/sn_store_cell_communication.md)
  : Store a cell-cell communication result on a Seurat object
- [`sn_get_cell_communication_result()`](https://songqi.org/shennong/dev/reference/sn_get_cell_communication_result.md)
  : Retrieve a stored cell-cell communication result
- [`sn_run_regulatory_activity()`](https://songqi.org/shennong/dev/reference/sn_run_regulatory_activity.md)
  : Infer transcription-factor or pathway activity
- [`sn_store_regulatory_activity()`](https://songqi.org/shennong/dev/reference/sn_store_regulatory_activity.md)
  : Store regulatory activity results on a Seurat object
- [`sn_get_regulatory_activity_result()`](https://songqi.org/shennong/dev/reference/sn_get_regulatory_activity_result.md)
  : Retrieve stored regulatory activity results

## Interpretation and Reporting

- [`sn_list_results()`](https://songqi.org/shennong/dev/reference/sn_list_results.md)
  : List stored Shennong analysis and interpretation results on a Seurat
  object

- [`sn_get_de_result()`](https://songqi.org/shennong/dev/reference/sn_get_de_result.md)
  : Retrieve a stored DE result from a Seurat object

- [`sn_get_enrichment_result()`](https://songqi.org/shennong/dev/reference/sn_get_enrichment_result.md)
  : Retrieve a stored enrichment result from a Seurat object

- [`sn_get_interpretation_result()`](https://songqi.org/shennong/dev/reference/sn_get_interpretation_result.md)
  : Retrieve a stored interpretation result from a Seurat object

- [`sn_store_enrichment()`](https://songqi.org/shennong/dev/reference/sn_store_enrichment.md)
  : Store an enrichment result on a Seurat object

- [`sn_prepare_annotation_evidence()`](https://songqi.org/shennong/dev/reference/sn_prepare_annotation_evidence.md)
  : Prepare cluster-annotation evidence from a Seurat object

- [`sn_prepare_de_evidence()`](https://songqi.org/shennong/dev/reference/sn_prepare_de_evidence.md)
  : Prepare differential-expression evidence from a stored DE result

- [`sn_prepare_enrichment_evidence()`](https://songqi.org/shennong/dev/reference/sn_prepare_enrichment_evidence.md)
  : Prepare enrichment evidence

- [`sn_prepare_results_evidence()`](https://songqi.org/shennong/dev/reference/sn_prepare_results_evidence.md)
  : Prepare manuscript-style results evidence

- [`sn_build_prompt()`](https://songqi.org/shennong/dev/reference/sn_build_prompt.md)
  : Build an LLM prompt from structured Shennong evidence

- [`sn_test_llm_provider()`](https://songqi.org/shennong/dev/reference/sn_test_llm_provider.md)
  : Test whether an ellmer-backed LLM provider is reachable and usable

- [`sn_make_ellmer_provider()`](https://songqi.org/shennong/dev/reference/sn_make_ellmer_provider.md)
  :

  Create an ellmer-backed provider for Shennong interpretation helpers

- [`sn_run_llm()`](https://songqi.org/shennong/dev/reference/sn_run_llm.md)
  : Run an LLM provider on a prepared message list

- [`sn_interpret_annotation()`](https://songqi.org/shennong/dev/reference/sn_interpret_annotation.md)
  : Interpret cluster markers for cell-type annotation

- [`sn_interpret_de()`](https://songqi.org/shennong/dev/reference/sn_interpret_de.md)
  : Interpret a stored differential-expression result

- [`sn_interpret_enrichment()`](https://songqi.org/shennong/dev/reference/sn_interpret_enrichment.md)
  : Interpret a stored enrichment result

- [`sn_write_results()`](https://songqi.org/shennong/dev/reference/sn_write_results.md)
  : Write a manuscript-style results summary from stored analysis
  outputs

- [`sn_write_figure_legend()`](https://songqi.org/shennong/dev/reference/sn_write_figure_legend.md)
  : Write a figure legend from stored analysis outputs

- [`sn_write_presentation_summary()`](https://songqi.org/shennong/dev/reference/sn_write_presentation_summary.md)
  : Write a presentation-style summary from stored analysis outputs

## Visualization

- [`sn_plot_dim()`](https://songqi.org/shennong/dev/reference/sn_plot_dim.md)
  : Create a dimensionality reduction plot for categorical data
- [`sn_plot_feature()`](https://songqi.org/shennong/dev/reference/sn_plot_feature.md)
  : Plot feature expression in reduced dimensions
- [`sn_plot_heatmap()`](https://songqi.org/shennong/dev/reference/sn_plot_heatmap.md)
  : Plot a heatmap for selected genes or features
- [`sn_plot_violin()`](https://songqi.org/shennong/dev/reference/sn_plot_violin.md)
  : Plot a violin plot with categorical groups
- [`sn_plot_dot()`](https://songqi.org/shennong/dev/reference/sn_plot_dot.md)
  : Plot a dot plot with categorical groups
- [`sn_plot_boxplot()`](https://songqi.org/shennong/dev/reference/sn_plot_boxplot.md)
  : Create a boxplot from a data frame
- [`sn_plot_barplot()`](https://songqi.org/shennong/dev/reference/sn_plot_barplot.md)
  : Create a bar plot from a data frame
- [`sn_plot_composition()`](https://songqi.org/shennong/dev/reference/sn_plot_composition.md)
  : Plot grouped composition-style bar charts
- [`sn_plot_milo()`](https://songqi.org/shennong/dev/reference/sn_plot_milo.md)
  : Plot miloR neighborhood differential-abundance results
- [`sn_list_palettes()`](https://songqi.org/shennong/dev/reference/sn_list_palettes.md)
  : List available color palettes
- [`sn_get_palette()`](https://songqi.org/shennong/dev/reference/sn_get_palette.md)
  : Resolve a palette into explicit colors

## Utilities and Project Setup

- [`sn_check_version()`](https://songqi.org/shennong/dev/reference/sn_check_version.md)
  : Check whether Shennong is up to date
- [`sn_install_shennong()`](https://songqi.org/shennong/dev/reference/sn_install_shennong.md)
  : Install Shennong from CRAN, GitHub, or a local source
- [`sn_check_file()`](https://songqi.org/shennong/dev/reference/sn_check_file.md)
  : Check if files exist
- [`sn_get_codex_skill_path()`](https://songqi.org/shennong/dev/reference/sn_get_codex_skill_path.md)
  : Return installed Shennong Codex asset paths
- [`sn_list_dependencies()`](https://songqi.org/shennong/dev/reference/sn_list_dependencies.md)
  : List Shennong runtime and recommended R package dependencies
- [`sn_install_dependencies()`](https://songqi.org/shennong/dev/reference/sn_install_dependencies.md)
  : Install missing Shennong dependencies in one step
- [`sn_check_pixi()`](https://songqi.org/shennong/dev/reference/sn_check_pixi.md)
  [`sn_install_pixi()`](https://songqi.org/shennong/dev/reference/sn_check_pixi.md)
  [`sn_ensure_pixi()`](https://songqi.org/shennong/dev/reference/sn_check_pixi.md)
  : Check, install, and configure pixi for Shennong Python backends
- [`sn_pixi_paths()`](https://songqi.org/shennong/dev/reference/sn_pixi_paths.md)
  : Inspect Shennong pixi runtime paths
- [`sn_list_pixi_environments()`](https://songqi.org/shennong/dev/reference/sn_list_pixi_environments.md)
  : List bundled pixi environment configs
- [`sn_pixi_config_path()`](https://songqi.org/shennong/dev/reference/sn_pixi_config_path.md)
  : Locate a bundled pixi config
- [`sn_prepare_pixi_environment()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  [`sn_call_pixi_environment()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  [`sn_call_scvi()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  [`sn_call_scanvi()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  [`sn_call_scarches()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  [`sn_call_scpoli()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  [`sn_call_infercnvpy()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  [`sn_call_cellphonedb()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  [`sn_call_cell2location()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  [`sn_call_tangram()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  [`sn_call_squidpy()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  [`sn_call_spatialdata()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  [`sn_call_stlearn()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  : Prepare or call a Shennong pixi environment
- [`sn_run_scarches()`](https://songqi.org/shennong/dev/reference/sn_run_scarches.md)
  [`sn_run_scpoli()`](https://songqi.org/shennong/dev/reference/sn_run_scarches.md)
  [`sn_run_infercnvpy()`](https://songqi.org/shennong/dev/reference/sn_run_scarches.md)
  [`sn_run_cellphonedb()`](https://songqi.org/shennong/dev/reference/sn_run_scarches.md)
  [`sn_run_cell2location()`](https://songqi.org/shennong/dev/reference/sn_run_scarches.md)
  [`sn_run_tangram()`](https://songqi.org/shennong/dev/reference/sn_run_scarches.md)
  [`sn_run_squidpy()`](https://songqi.org/shennong/dev/reference/sn_run_scarches.md)
  [`sn_run_spatialdata()`](https://songqi.org/shennong/dev/reference/sn_run_scarches.md)
  [`sn_run_stlearn()`](https://songqi.org/shennong/dev/reference/sn_run_scarches.md)
  : Run a Python analysis command through a managed Shennong pixi
  environment
- [`sn_detect_accelerator()`](https://songqi.org/shennong/dev/reference/sn_detect_accelerator.md)
  : Detect local accelerator support for pixi-managed Python methods
- [`sn_configure_pixi_mirror()`](https://songqi.org/shennong/dev/reference/sn_configure_pixi_mirror.md)
  : Configure pixi mirrors for Shennong runtime environments
- [`sn_initialize_project()`](https://songqi.org/shennong/dev/reference/sn_initialize_project.md)
  : Initialize a Shennong analysis project
- [`sn_install_codex_skill()`](https://songqi.org/shennong/dev/reference/sn_install_codex_skill.md)
  : Install the bundled Shennong Codex skill for end-user agents
- [`sn_set_path()`](https://songqi.org/shennong/dev/reference/sn_set_path.md)
  : Create a directory path if needed
