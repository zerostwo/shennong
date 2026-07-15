# Package index

## Analysis Methods and Result Contract

- [`sn_list_methods()`](https://songqi.org/shennong/dev/reference/sn_list_methods.md)
  : List registered Shennong analysis methods
- [`sn_method_status()`](https://songqi.org/shennong/dev/reference/sn_method_status.md)
  : Report the status of a registered Shennong method
- [`sn_store_result()`](https://songqi.org/shennong/dev/reference/sn_store_result.md)
  : Store a Shennong analysis result on a Seurat object
- [`sn_get_result()`](https://songqi.org/shennong/dev/reference/sn_get_result.md)
  : Retrieve a stored Shennong analysis result
- [`sn_list_results()`](https://songqi.org/shennong/dev/reference/sn_list_results.md)
  : List stored Shennong analysis and interpretation results on a Seurat
  object
- [`sn_delete_result()`](https://songqi.org/shennong/dev/reference/sn_delete_result.md)
  : Delete a stored Shennong analysis result
- [`sn_validate_result()`](https://songqi.org/shennong/dev/reference/sn_validate_result.md)
  : Validate a Shennong analysis result

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
- [`sn_list_datasets()`](https://songqi.org/shennong/dev/reference/sn_list_datasets.md)
  : List datasets available through Shennong
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
- [`sn_convert_bpcells()`](https://songqi.org/shennong/dev/reference/sn_convert_bpcells.md)
  : Convert Seurat assay layers to BPCells-backed matrices
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
- [`sn_run_multimodal()`](https://songqi.org/shennong/dev/reference/sn_run_multimodal.md)
  : Run multimodal clustering through the unified clustering workflow
- [`sn_run_scvi()`](https://songqi.org/shennong/dev/reference/sn_run_scvi.md)
  [`sn_run_scanvi()`](https://songqi.org/shennong/dev/reference/sn_run_scvi.md)
  : Run scVI or scANVI integration through Shennong
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
- [`sn_transfer_labels()`](https://songqi.org/shennong/dev/reference/sn_transfer_labels.md)
  : Transfer labels from a Seurat reference to a query object
- [`sn_prepare_label_transfer_reference()`](https://songqi.org/shennong/dev/reference/sn_prepare_label_transfer_reference.md)
  : Prepare a compact label-transfer reference
- [`sn_simulate()`](https://songqi.org/shennong/dev/reference/sn_simulate.md)
  : Simulate single-cell counts with scDesign3
- [`sn_simulate_scdesign3()`](https://songqi.org/shennong/dev/reference/sn_simulate_scdesign3.md)
  : Simulate single-cell counts with scDesign3
- [`sn_sweep_cluster_resolution()`](https://songqi.org/shennong/dev/reference/sn_sweep_cluster_resolution.md)
  : Sweep clustering resolutions and summarize cluster-quality
  diagnostics

## Python and Pixi Backends

- [`sn_check_pixi()`](https://songqi.org/shennong/dev/reference/sn_check_pixi.md)
  [`sn_install_pixi()`](https://songqi.org/shennong/dev/reference/sn_check_pixi.md)
  [`sn_ensure_pixi()`](https://songqi.org/shennong/dev/reference/sn_check_pixi.md)
  : Check, install, and configure pixi for Shennong Python backends
- [`sn_prepare_pixi_environment()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  [`sn_call_pixi_environment()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  [`sn_call_scvi()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  [`sn_call_scanvi()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  [`sn_call_mmochi()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  [`sn_call_scarches()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  [`sn_call_scpoli()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  [`sn_call_infercnvpy()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  [`sn_call_trajectory()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
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

- [`sn_run_annotation()`](https://songqi.org/shennong/dev/reference/sn_run_annotation.md)
  : Run traceable cell-type annotation
- [`sn_annotation_consensus()`](https://songqi.org/shennong/dev/reference/sn_annotation_consensus.md)
  : Build consensus annotation labels from evidence
- [`sn_annotation_confidence()`](https://songqi.org/shennong/dev/reference/sn_annotation_confidence.md)
  : Calibrate confidence from annotation evidence
- [`sn_map_cell_ontology()`](https://songqi.org/shennong/dev/reference/sn_map_cell_ontology.md)
  : Map cell labels to Cell Ontology identifiers
- [`sn_review_annotation()`](https://songqi.org/shennong/dev/reference/sn_review_annotation.md)
  : Review stored annotation evidence and low-confidence labels
- [`sn_plot_annotation_confidence()`](https://songqi.org/shennong/dev/reference/sn_plot_annotation_confidence.md)
  : Plot annotation confidence
- [`sn_plot_annotation_markers()`](https://songqi.org/shennong/dev/reference/sn_plot_annotation_markers.md)
  : Plot annotation marker evidence
- [`sn_plot_annotation_confusion()`](https://songqi.org/shennong/dev/reference/sn_plot_annotation_confusion.md)
  : Plot annotation confusion against known labels
- [`sn_run_celltypist()`](https://songqi.org/shennong/dev/reference/sn_run_celltypist.md)
  : Run CellTypist for automated cell type annotation
- [`sn_find_de()`](https://songqi.org/shennong/dev/reference/sn_find_de.md)
  : Run differential expression analysis on a Seurat object
- [`sn_annotate_de_features()`](https://songqi.org/shennong/dev/reference/sn_annotate_de_features.md)
  : Annotate DE or marker genes by feature class
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
- [`sn_score_programs()`](https://songqi.org/shennong/dev/reference/sn_score_programs.md)
  : Score gene programs in cells or aggregated samples
- [`sn_test_programs()`](https://songqi.org/shennong/dev/reference/sn_test_programs.md)
  : Test program activity between conditions
- [`sn_plot_program_activity()`](https://songqi.org/shennong/dev/reference/sn_plot_program_activity.md)
  : Plot program activity distributions
- [`sn_plot_program_heatmap()`](https://songqi.org/shennong/dev/reference/sn_plot_program_heatmap.md)
  : Plot a program activity heatmap
- [`sn_discover_programs()`](https://songqi.org/shennong/dev/reference/sn_discover_programs.md)
  : Discover latent gene programs
- [`sn_plot_discovered_programs()`](https://songqi.org/shennong/dev/reference/sn_plot_discovered_programs.md)
  : Plot discovered gene programs
- [`sn_run_grn()`](https://songqi.org/shennong/dev/reference/sn_run_grn.md)
  : Infer and summarize a gene regulatory network
- [`sn_plot_regulon()`](https://songqi.org/shennong/dev/reference/sn_plot_regulon.md)
  : Plot a gene regulatory network result

## Trajectory and Dynamic Genes

- [`sn_run_trajectory()`](https://songqi.org/shennong/dev/reference/sn_run_trajectory.md)
  : Infer trajectories and test dynamic genes
- [`sn_run_velocity()`](https://songqi.org/shennong/dev/reference/sn_run_velocity.md)
  : Run RNA velocity with the managed scVelo backend
- [`sn_plot_velocity()`](https://songqi.org/shennong/dev/reference/sn_plot_velocity.md)
  : Plot RNA velocity vectors
- [`sn_run_fate()`](https://songqi.org/shennong/dev/reference/sn_run_fate.md)
  : Infer terminal states and fate probabilities with CellRank
- [`sn_plot_fate()`](https://songqi.org/shennong/dev/reference/sn_plot_fate.md)
  : Plot CellRank fate probabilities
- [`sn_plot_trajectory()`](https://songqi.org/shennong/dev/reference/sn_plot_trajectory.md)
  : Plot an inferred trajectory on its embedding
- [`sn_plot_pseudotime()`](https://songqi.org/shennong/dev/reference/sn_plot_pseudotime.md)
  : Plot pseudotime on the trajectory embedding
- [`sn_plot_lineage_probability()`](https://songqi.org/shennong/dev/reference/sn_plot_lineage_probability.md)
  : Plot lineage assignment probability
- [`sn_plot_dynamic_heatmap()`](https://songqi.org/shennong/dev/reference/sn_plot_dynamic_heatmap.md)
  : Plot fitted dynamic-gene trends as a heatmap
- [`sn_plot_gene_trend()`](https://songqi.org/shennong/dev/reference/sn_plot_gene_trend.md)
  : Plot fitted expression trends for selected genes
- [`sn_plot_branch_comparison()`](https://songqi.org/shennong/dev/reference/sn_plot_branch_comparison.md)
  : Plot branch-specific dynamic-gene evidence

## Spatial Analysis

- [`sn_run_spatial()`](https://songqi.org/shennong/dev/reference/sn_run_spatial.md)
  : Unified spatial workflow dispatcher
- [`sn_find_spatial_features()`](https://songqi.org/shennong/dev/reference/sn_find_spatial_features.md)
  : Find spatially variable features
- [`sn_find_spatial_domains()`](https://songqi.org/shennong/dev/reference/sn_find_spatial_domains.md)
  : Identify spatial domains
- [`sn_run_spatial_neighborhood()`](https://songqi.org/shennong/dev/reference/sn_run_spatial_neighborhood.md)
  : Analyze spatial neighborhoods
- [`sn_run_spatial_deconvolution()`](https://songqi.org/shennong/dev/reference/sn_run_spatial_deconvolution.md)
  : Run spatial deconvolution through cell2location
- [`sn_run_spatial_mapping()`](https://songqi.org/shennong/dev/reference/sn_run_spatial_mapping.md)
  : Map single cells to space through Tangram
- [`sn_integrate_spatial()`](https://songqi.org/shennong/dev/reference/sn_integrate_spatial.md)
  : Integrate spatial samples with an explicit backend adapter
- [`sn_run_spatial_communication()`](https://songqi.org/shennong/dev/reference/sn_run_spatial_communication.md)
  : Add spatial distance evidence to a communication result
- [`sn_plot_spatial()`](https://songqi.org/shennong/dev/reference/sn_plot_spatial.md)
  : Plot spatial coordinates colored by metadata
- [`sn_plot_spatial_feature()`](https://songqi.org/shennong/dev/reference/sn_plot_spatial_feature.md)
  : Plot expression in spatial coordinates
- [`sn_plot_spatial_domain()`](https://songqi.org/shennong/dev/reference/sn_plot_spatial_domain.md)
  : Plot spatial-domain assignments
- [`sn_plot_spatial_svg()`](https://songqi.org/shennong/dev/reference/sn_plot_spatial_svg.md)
  : Plot spatial-feature statistics
- [`sn_plot_spatial_neighborhood()`](https://songqi.org/shennong/dev/reference/sn_plot_spatial_neighborhood.md)
  : Plot spatial-neighborhood enrichment or co-occurrence
- [`sn_plot_spatial_deconvolution()`](https://songqi.org/shennong/dev/reference/sn_plot_spatial_deconvolution.md)
  : Plot spatial deconvolution proportions
- [`sn_plot_spatial_communication()`](https://songqi.org/shennong/dev/reference/sn_plot_spatial_communication.md)
  : Plot spatially constrained communication

## Composition and Comparative Analysis

- [`sn_test_abundance()`](https://songqi.org/shennong/dev/reference/sn_test_abundance.md)
  : Test differential abundance across biological samples
- [`sn_plot_abundance()`](https://songqi.org/shennong/dev/reference/sn_plot_abundance.md)
  : Plot differential-abundance effects
- [`sn_prioritize_states()`](https://songqi.org/shennong/dev/reference/sn_prioritize_states.md)
  : Prioritize phenotype-responsive or rare cell states
- [`sn_plot_state_priority()`](https://songqi.org/shennong/dev/reference/sn_plot_state_priority.md)
  : Plot cell-state priority scores
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

## Bulk Transcriptomics

- [`sn_run_bulk()`](https://songqi.org/shennong/dev/reference/sn_run_bulk.md)
  : Run a bulk transcriptomics workflow
- [`sn_assess_bulk_qc()`](https://songqi.org/shennong/dev/reference/sn_assess_bulk_qc.md)
  : Assess bulk transcriptomics sample quality
- [`sn_find_bulk_de()`](https://songqi.org/shennong/dev/reference/sn_find_bulk_de.md)
  : Find differential expression in bulk transcriptomics data
- [`sn_score_bulk_pathways()`](https://songqi.org/shennong/dev/reference/sn_score_bulk_pathways.md)
  : Score pathways in bulk expression samples
- [`sn_run_wgcna()`](https://songqi.org/shennong/dev/reference/sn_run_wgcna.md)
  : Run weighted gene co-expression network analysis
- [`sn_run_survival()`](https://songqi.org/shennong/dev/reference/sn_run_survival.md)
  : Run Cox proportional-hazards models for bulk features
- [`sn_run_clinical_association()`](https://songqi.org/shennong/dev/reference/sn_run_clinical_association.md)
  : Associate bulk features with clinical variables
- [`sn_plot_bulk_qc()`](https://songqi.org/shennong/dev/reference/sn_plot_bulk_qc.md)
  : Plot bulk sample quality metrics
- [`sn_plot_bulk_pca()`](https://songqi.org/shennong/dev/reference/sn_plot_bulk_pca.md)
  : Plot bulk sample PCA
- [`sn_plot_sample_correlation()`](https://songqi.org/shennong/dev/reference/sn_plot_sample_correlation.md)
  : Plot bulk sample correlation
- [`sn_plot_bulk_de()`](https://songqi.org/shennong/dev/reference/sn_plot_bulk_de.md)
  : Plot bulk differential expression
- [`sn_plot_wgcna()`](https://songqi.org/shennong/dev/reference/sn_plot_wgcna.md)
  : Plot WGCNA modules or trait associations
- [`sn_plot_survival()`](https://songqi.org/shennong/dev/reference/sn_plot_survival.md)
  : Plot bulk survival associations

## Communication and Regulatory Activity

- [`sn_run_cell_communication()`](https://songqi.org/shennong/dev/reference/sn_run_cell_communication.md)
  : Run cell-cell communication inference
- [`sn_store_cell_communication()`](https://songqi.org/shennong/dev/reference/sn_store_cell_communication.md)
  : Store a cell-cell communication result on a Seurat object
- [`sn_get_cell_communication_result()`](https://songqi.org/shennong/dev/reference/sn_get_cell_communication_result.md)
  : Retrieve a stored cell-cell communication result
- [`sn_plot_communication()`](https://songqi.org/shennong/dev/reference/sn_plot_communication.md)
  : Plot standardized cell-cell communication results
- [`sn_plot_ligand_target()`](https://songqi.org/shennong/dev/reference/sn_plot_ligand_target.md)
  : Plot NicheNet or MultiNicheNet ligand-target evidence
- [`sn_plot_communication_comparison()`](https://songqi.org/shennong/dev/reference/sn_plot_communication_comparison.md)
  : Plot sample-aware differential communication effects
- [`sn_run_regulatory_activity()`](https://songqi.org/shennong/dev/reference/sn_run_regulatory_activity.md)
  : Infer transcription-factor or pathway activity
- [`sn_store_regulatory_activity()`](https://songqi.org/shennong/dev/reference/sn_store_regulatory_activity.md)
  : Store regulatory activity results on a Seurat object
- [`sn_get_regulatory_activity_result()`](https://songqi.org/shennong/dev/reference/sn_get_regulatory_activity_result.md)
  : Retrieve stored regulatory activity results

## CNV, Malignancy, and Metabolism

- [`sn_run_cnv()`](https://songqi.org/shennong/dev/reference/sn_run_cnv.md)
  : Run copy-number and malignancy analysis
- [`sn_plot_cnv()`](https://songqi.org/shennong/dev/reference/sn_plot_cnv.md)
  : Plot CNV, malignancy, and subclone results
- [`sn_metabolic_signatures()`](https://songqi.org/shennong/dev/reference/sn_metabolic_signatures.md)
  : Retrieve curated core metabolic signatures
- [`sn_run_metabolism()`](https://songqi.org/shennong/dev/reference/sn_run_metabolism.md)
  : Run unified single-cell metabolic activity analysis
- [`sn_plot_metabolism()`](https://songqi.org/shennong/dev/reference/sn_plot_metabolism.md)
  : Plot metabolic pathway activity and differential results

## Interpretation and Reporting

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

- [`sn_list_figure_profiles()`](https://songqi.org/shennong/dev/reference/sn_list_figure_profiles.md)
  : List publication figure profiles
- [`sn_figure_spec()`](https://songqi.org/shennong/dev/reference/sn_figure_spec.md)
  : Inspect or calculate a Shennong figure specification
- [`sn_recommend_figure_size()`](https://songqi.org/shennong/dev/reference/sn_recommend_figure_size.md)
  : Recommend output dimensions for a figure
- [`sn_apply_figure_profile()`](https://songqi.org/shennong/dev/reference/sn_apply_figure_profile.md)
  : Apply a generic publication profile to a plot
- [`sn_validate_figure()`](https://songqi.org/shennong/dev/reference/sn_validate_figure.md)
  : Validate a Shennong figure before export
- [`sn_save_figure()`](https://songqi.org/shennong/dev/reference/sn_save_figure.md)
  : Save a publication figure with deterministic dimensions
- [`sn_export_figure()`](https://songqi.org/shennong/dev/reference/sn_export_figure.md)
  : Export a publication figure
- [`sn_export_figure_bundle()`](https://songqi.org/shennong/dev/reference/sn_export_figure_bundle.md)
  : Export a figure, source data, specification, and manifest bundle
- [`sn_plot_qc()`](https://songqi.org/shennong/dev/reference/sn_plot_qc.md)
  : Plot QC assessment summaries
- [`sn_plot_qc_thresholds()`](https://songqi.org/shennong/dev/reference/sn_plot_qc_thresholds.md)
  : Plot cell-level QC thresholds
- [`sn_plot_doublets()`](https://songqi.org/shennong/dev/reference/sn_plot_doublets.md)
  : Plot doublet classifications in an embedding
- [`sn_plot_ambient_correction()`](https://songqi.org/shennong/dev/reference/sn_plot_ambient_correction.md)
  : Plot ambient RNA correction totals
- [`sn_plot_hvg()`](https://songqi.org/shennong/dev/reference/sn_plot_hvg.md)
  : Plot highly variable feature diagnostics
- [`sn_plot_elbow()`](https://songqi.org/shennong/dev/reference/sn_plot_elbow.md)
  : Plot PCA elbow diagnostics
- [`sn_plot_cluster_tree()`](https://songqi.org/shennong/dev/reference/sn_plot_cluster_tree.md)
  : Plot cluster transitions across resolutions
- [`sn_plot_resolution_sweep()`](https://songqi.org/shennong/dev/reference/sn_plot_resolution_sweep.md)
  : Plot a clustering resolution sweep
- [`sn_plot_integration()`](https://songqi.org/shennong/dev/reference/sn_plot_integration.md)
  : Plot integration assessment metrics
- [`sn_plot_reference_projection()`](https://songqi.org/shennong/dev/reference/sn_plot_reference_projection.md)
  : Plot reference annotation projection
- [`sn_plot_de()`](https://songqi.org/shennong/dev/reference/sn_plot_de.md)
  : Plot a differential-expression result
- [`sn_plot_enrichment()`](https://songqi.org/shennong/dev/reference/sn_plot_enrichment.md)
  : Plot enrichment results
- [`sn_plot_gsea()`](https://songqi.org/shennong/dev/reference/sn_plot_gsea.md)
  : Plot a GSEA running-score curve or summary
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
  [`sn_call_mmochi()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  [`sn_call_scarches()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  [`sn_call_scpoli()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  [`sn_call_infercnvpy()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  [`sn_call_trajectory()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  [`sn_call_cellphonedb()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  [`sn_call_cell2location()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  [`sn_call_tangram()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  [`sn_call_squidpy()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  [`sn_call_spatialdata()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  [`sn_call_stlearn()`](https://songqi.org/shennong/dev/reference/sn_prepare_pixi_environment.md)
  : Prepare or call a Shennong pixi environment
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
