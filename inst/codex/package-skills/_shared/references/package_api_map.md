# Shennong Package API Map

This reference summarizes the exported Shennong surface so an agent can choose
the right `sn_*` entry point without scanning the whole package.

Core object rule:
- Current workflows are centered on Seurat objects.
- Standalone bulk workflows accept feature-by-sample matrices, lists, and
  `SummarizedExperiment` objects and return the same validated result contract.
- Prefer Shennong APIs over raw Seurat calls when the package already exposes
  the needed behavior.
- Prefer stored-result workflows over ad hoc `object@misc` access.

## Data and Import

- `sn_list_datasets()`: list sample-level datasets available through the current Shennong public Zenodo collection, plus bundled examples with `source = "all"`
- `sn_load_data()`: load packaged or remote example datasets; vectorized bundled datasets such as `c("pbmc1k", "pbmc3k")` return a merged Seurat object for filtered matrices or a named list for raw matrices. Public Shennong collection samples can be loaded by sample ID from `sn_list_datasets()` or by `study_id` plus `sample_id`; the study ZIP is cached and only the requested filtered/raw H5 or metrics file is extracted. Use `backend = "api"` for a lazy ShennongData 0.2 resource handle, with assay/layer selection in `api_args`; set `lazy = FALSE` only for explicit collection.
- `sn_download_zenodo()`: download reusable files from public Zenodo records without a token, or from restricted/private records when a token is supplied
- `sn_list_10x_paths()`: discover 10x `outs/`, filtered/raw matrices, H5, or metrics paths
- `sn_read()`: import tabular, serialized (`qs`/`qs2`), and bioinformatics formats
- `sn_write()`: export tabular, serialized (`qs`/`qs2`), and supported omics formats
- registered `rio` handlers: `.import.rio_10x()`, `.import.rio_10x_spatial()`, `.import.rio_starsolo()`, `.import.rio_gmt()`, `.import.rio_h5()`, `.import.rio_h5ad()`, `.import.rio_qs()`, `.import.rio_qs2()`, `.import.rio_bpcells()`, `.export.rio_h5()`, `.export.rio_h5ad()`, `.export.rio_qs()`, `.export.rio_qs2()`, and `.export.rio_bpcells()` are the package-level import/export hooks used by `sn_read()` and `sn_write()`
- `sn_convert_bpcells()`: write selected Seurat assay layers to BPCells matrix directories and rebind those layers in the returned object
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

- `sn_run_cluster()`: single-dataset clustering or batch integration; supports Seurat log-normalization, SCTransform, CITE-seq workflows with `modality = "cite_seq"` and `multimodal_method = "wnn"` / `"totalvi"` / `"coralysis"` / `"mmochi"`, Harmony, native Coralysis, Seurat CCA/RPCA, and pixi-managed scVI/scANVI/totalVI/MMoCHi integration. CITE-seq WNN combines RNA PCA with ADT CLR normalization/PCA, clusters on `wsnn`, and returns `wnn.umap`; CITE-seq totalVI writes RNA and ADT counts to the shared scVI-family pixi backend; CITE-seq Coralysis runs native Coralysis on the ADT protein assay; CITE-seq MMoCHi runs ADT landmark registration across `batch` or in single-sample mode when `batch = NULL`, stores the corrected protein matrix as an assay layer when supported and otherwise under `object@misc$mmochi$corrected_protein`, and clusters on a protein-derived `mmochi` reduction. SCTransform integration currently uses Harmony. Coralysis, scVI/scANVI, totalVI, and MMoCHi skip redundant Seurat PCA stages that their backends do not consume. Native Coralysis stores the trained SingleCellExperiment under `object@misc$coralysis` by default so the returned object can be used directly for label transfer; set `integration_control = list(store_sce = FALSE)` only for clustering-only runs. For `integration_method = "scvi"` or `"scanvi"`, Shennong writes selected counts/metadata under `~/.shennong/runs/`, manages the shared scVI-family pixi project under `~/.shennong/pixi/scvi/`, imports the learned latent reduction, then continues Seurat neighbors/clustering/UMAP; scANVI requires `integration_control = list(label_by = ...)`. Use `integration_control = list(accelerator = "auto", mirror = "auto")` for CUDA/CPU auto-selection and Shennong-level mirror configuration. `block_genes` can mix bundled signature queries such as `cellCycle.G2M`, `ribo`, and `mito` with custom gene symbols before internally selected HVGs are stored in log-normalization and SCTransform workflows. Rare-aware feature augmentation can combine `gini` and `local_markers`, with advanced thresholds kept in `rare_feature_control`. Use `hvg_features` to merge user-supplied marker genes into the backend feature set and, for PCA-based workflows, the final ScaleData/PCA feature set. Use `umap_control = list(...)` to tune `Seurat::RunUMAP()` arguments such as `n.neighbors`, `min.dist`, `spread`, and `reduction.name` without changing the clustering graph. Re-running on the returned object reuses matching stages by default; use `rerun_from` or `reuse = FALSE` for forced recompute. Leiden clustering auto-installs `leidenbase` unless `auto_install = FALSE`.
- `sn_run_multimodal()`: explicit CITE-seq wrapper over `sn_run_cluster()` for
  WNN, totalVI, Coralysis, or MMoCHi; it preserves the clustering return
  contract and forwards all workflow controls.
- `sn_run_scvi()` / `sn_run_scanvi()`: explicit wrappers for the corresponding `sn_run_cluster()` integration methods.
- `sn_transfer_labels()`: query-first reference `label_by` transfer wrapper. Defaults to Seurat anchors, can use `method = "coralysis"` for native Coralysis `ReferenceMapping()` when the reference stores a trained Coralysis SingleCellExperiment, and supports semi-supervised scVI-family transfer with `method = "scanvi"` or `method = "scarches"`.
- `sn_prepare_label_transfer_reference()`: create compact transfer-ready references. Coralysis output is a minimal SingleCellExperiment with trained models, PCA model, feature names, and labels; Seurat/scANVI/scArches output is a slim Seurat reference with selected assay layers and labels.
- `sn_simulate()`: method-based simulation entry point; currently supports `method = "scdesign3"` for Seurat or SingleCellExperiment inputs and returns Seurat, SingleCellExperiment, sparse counts, or the raw scDesign3 result.
- `sn_simulate_scdesign3()`: backend-specific scDesign3 wrapper when direct control of the scDesign3 design arguments is needed.
- `sn_plot_heatmap()`: focused heatmap for user-selected genes, with cell-level and group-averaged modes, optional grouping/splitting, and selected-feature scaling.

## Python Runtime Helpers

- `sn_check_pixi()`: check whether a pixi executable is available.
- `sn_install_pixi()` / `sn_ensure_pixi()`: install or ensure the standalone pixi binary when Python backends need it.
- `sn_pixi_paths()`: inspect the `~/.shennong/pixi/` layout for scVI/scANVI and other Python method families.
- `sn_list_pixi_environments()` / `sn_pixi_config_path()`: discover bundled pixi configs under `inst/pixi/`.
- `sn_prepare_pixi_environment()` / `sn_call_pixi_environment()`: materialize a bundled config into `~/.shennong/pixi/<family>/` and run commands inside it.
- `sn_call_scvi()`, `sn_call_scanvi()`, `sn_call_mmochi()`, `sn_call_scarches()`, `sn_call_scpoli()`, `sn_call_infercnvpy()`, `sn_call_trajectory()`, `sn_call_cellphonedb()`, `sn_call_cell2location()`, `sn_call_tangram()`, `sn_call_squidpy()`, `sn_call_spatialdata()`, `sn_call_stlearn()`: environment-specific command-call helpers. `scanvi` shares the `scvi` environment; `scpoli` shares the `scarches` environment.
- `sn_run_scarches(object = ...)`, `sn_run_scpoli(object = ...)`, `sn_run_infercnvpy(object = ...)`, `sn_run_cellphonedb(object = ...)`, `sn_run_cell2location(object = ...)`, `sn_run_tangram(object = ...)`, `sn_run_squidpy(object = ...)`, `sn_run_spatialdata(object = ...)`, `sn_run_stlearn(object = ...)`: object-level Python wrappers. They export Seurat input under `~/.shennong/runs/`, run family-local scripts from `inst/pixi/<family>/scripts/`, import cell-level metadata/reductions when produced, and record manifests under `object@misc`.
- Use the `sn_call_*()` helpers for direct command execution in managed Python environments. Object-level `sn_run_*()` wrappers require a Seurat object and should be used only for package workflows that export/import analysis state.
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

- `sn_run_annotation()`: marker/reference annotation mainline with consensus,
  SingleR, CellTypist, Seurat, Symphony, scmap, and scANVI backends; stores
  cell/cluster labels, confidence, hierarchy, evidence, ontology IDs, raw
  backend predictions, diagnostics, and provenance
- `sn_annotation_consensus()` / `sn_annotation_confidence()`: combine and
  calibrate long-form marker/reference evidence without LLM label overrides
- `sn_map_cell_ontology()`: map labels against the bundled versioned Cell
  Ontology snapshot or a project mapping
- `sn_review_annotation()`: inspect low-confidence cells/clusters and evidence
- `sn_plot_annotation_confidence()` / `sn_plot_annotation_markers()` /
  `sn_plot_annotation_confusion()`: result-aware annotation diagnostics
- `sn_run_celltypist()`: external CellTypist-based annotation
- `sn_find_de()`: markers, contrasts, and pseudobulk DE
- `sn_annotate_de_features()`: flag marker/DE genes that encode TFs, surface/plasma-membrane proteins, cytokines, or chemokines
- `sn_enrich()`: ORA or GSEA from vectors, tables, or stored DE
- `sn_list_signatures()`: list bundled signatures
- `sn_get_signatures()`: retrieve signatures by path or category
- `sn_add_signature()`: add a signature to the editable registry
- `sn_update_signature()`: update a signature in the editable registry
- `sn_delete_signature()`: delete a signature from the editable registry
- `sn_score_programs()`: UCell, AUCell, GSVA, ssGSEA, or mean program scoring
  with feature-coverage diagnostics, stored long-form scores, and cell metadata
- `sn_test_programs()`: sample-aware program activity comparisons; when
  `sample_by` is present, cells are aggregated before inference
- `sn_plot_program_activity()` / `sn_plot_program_heatmap()`: result-aware
  program score distributions and heatmaps
- `sn_discover_programs()`: multi-restart NMF discovery or explicit cNMF and
  Hotspot adapters; retrieve weights and activity with
  `sn_get_result(object, "program_discovery", name)`
- `sn_plot_discovered_programs()`: gene-weight, activity, and restart plots
- `sn_run_grn()`: GENIE3 inference or explicit pySCENIC/SCENIC/GRNBoost2
  adapters with unified edge, regulon, activity, and specificity tables
- `sn_plot_regulon()`: network, activity, and group-specificity plots

## Trajectory and Dynamic Genes

- `sn_run_trajectory()`: direct Slingshot or Monocle 3 inference and an
  explicit Palantir runner/result adapter, with per-cell pseudotime, lineage
  probabilities, terminal states, and optional tradeSeq dynamic/branch tests
  plus fitted trends
- `sn_plot_trajectory()` / `sn_plot_pseudotime()` /
  `sn_plot_lineage_probability()`: embedding views backed by the stored result
- `sn_plot_dynamic_heatmap()` / `sn_plot_gene_trend()` /
  `sn_plot_branch_comparison()`: tradeSeq trend and branch-test views
- `sn_run_velocity()` / `sn_plot_velocity()`: managed scVelo inference from
  spliced/unspliced layers with projected vectors, transition evidence,
  pseudotime, and confidence
- `sn_run_fate()` / `sn_plot_fate()`: CellRank GPCCA terminal states, fate
  probabilities, and optional lineage drivers from a stored velocity result

## Spatial Workflows

- `sn_run_spatial()`: dispatch QC, SVG, domain, neighborhood, deconvolution,
  mapping, integration, or communication tasks
- `sn_find_spatial_features()`: Moran's I with permutation evidence, nnSVG,
  or explicit SPARK-X adapters
- `sn_find_spatial_domains()`: optional BANKSY or explicit
  stLearn/BayesSpace/CellCharter adapters
- `sn_run_spatial_neighborhood()`: memory-bounded KNN graph, permutation
  enrichment, and distance-bin co-occurrence
- `sn_run_spatial_deconvolution()` / `sn_run_spatial_mapping()`: stable aliases
  for the existing cell2location and Tangram object workflows
- `sn_integrate_spatial()`: explicit STAligner/Harmony/custom result adapter
- `sn_run_spatial_communication()`: augment a stored communication result with
  group distance evidence and optional distance filtering
- `sn_plot_spatial*()`: result-aware coordinate, SVG, domain, neighborhood,
  deconvolution, and communication figures with fixed spatial aspect

## Bulk Transcriptomics

- `sn_run_bulk()`: dispatcher for standalone QC, DE, pathway, network, and
  survival workflows
- `sn_assess_bulk_qc()`: library size, detected features, distributions, PCA,
  sample correlation, and robust outlier evidence
- `sn_find_bulk_de()`: validated fixed/mixed design and explicit contrast with
  automatic or direct edgeR, DESeq2, limma-voom, limma, and dream backends
- `sn_score_bulk_pathways()`: mean, GSVA, or ssGSEA sample scores with gene-set
  coverage diagnostics
- `sn_run_wgcna()`: weighted co-expression modules, eigengenes, soft-power
  evidence, and sample trait associations
- `sn_run_survival()` / `sn_run_clinical_association()`: sample-level Cox and
  phenotype models using expression features or metadata scores
- `sn_plot_bulk_qc()` / `sn_plot_bulk_pca()` /
  `sn_plot_sample_correlation()` / `sn_plot_bulk_de()` / `sn_plot_wgcna()` /
  `sn_plot_survival()`: result-aware bulk figures

## Publication Figures

- `sn_list_figure_profiles()`: inspect generic screen, column, page, and slide
  constraints; journal requirements must still be checked at submission time
- `sn_figure_spec()` / `sn_recommend_figure_size()`: calculate canvas, point,
  alpha, font, line, legend, raster, layout, and pagination recommendations
  from plot/data metadata without rendering large synthetic inputs
- `sn_apply_figure_profile()`: attach profile styling/specification while
  preserving the native ggplot/patchwork class
- `sn_validate_figure()`: structured preflight checks and suggested actions
- `sn_save_figure()` / `sn_export_figure()`: deterministic PDF, SVG, TIFF, or
  PNG output independent of the interactive device
- `sn_export_figure_bundle()`: figures plus available source data, spec,
  session, checksums, validation, and JSON manifest
- `sn_plot_de()` / `sn_plot_enrichment()` / `sn_plot_gsea()`: standardized
  result-aware plots whose evidence tables can be exported in the bundle
- `sn_plot_qc*()` / `sn_plot_doublets()` /
  `sn_plot_ambient_correction()` / `sn_plot_hvg()` / `sn_plot_elbow()` /
  `sn_plot_cluster_tree()` / `sn_plot_resolution_sweep()` /
  `sn_plot_integration()` / `sn_plot_reference_projection()`: core diagnostic
  figure surface

## Differential Abundance and State Priority

- `sn_test_abundance()`: sample-level Propeller/permutation/scCODA or
  neighborhood-level Milo through one versioned result contract; scCODA and
  pertpy outputs enter through an explicit runner/result adapter, and posterior
  inclusion probabilities are never relabeled as frequentist p values
- `sn_plot_abundance()`: standardized effect view for stored abundance results
- `sn_prioritize_states()`: sample-held-out perturbation separability, explicit
  bulk-input Scissor, or RareQ discovery plus sample-level association
- `sn_plot_state_priority()`: ranked state-priority view

## CNV, Malignancy, and Metabolism

- `sn_run_cnv()`: unified inferCNVpy/CopyKAT analysis with declared normal
  references, malignancy scores/calls, subclones, chromosome evidence, sample
  summaries, CNV UMAP, and expression association
- `sn_plot_cnv()`: chromosome heatmap, CNV UMAP, malignancy distribution,
  sample summary, or CNV-expression association from a stored result
- `sn_metabolic_signatures()`: curated core metabolic pathway gene sets
- `sn_run_metabolism()`: UCell/GSVA/ssGSEA/mean pathway scoring plus
  scMetabolism or explicit scFEA/Compass adapters; condition tests aggregate to
  `sample_by` first
- `sn_plot_metabolism()`: pathway activity, sample heatmap, sample comparison,
  and differential-effect views

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
- `sn_run_cell_communication()`: run LIANA, CellChat, CellPhoneDB, NicheNet, or
  MultiNicheNet alone or as a cross-method consensus, with optional
  sample-level condition comparison
- `sn_store_cell_communication()`: persist communication results
- `sn_get_cell_communication_result()`: retrieve communication results
- `sn_plot_communication()`: bubble, heatmap, network, chord, or river view of
  standardized interactions
- `sn_plot_ligand_target()` / `sn_plot_communication_comparison()`: inspect
  ligand-target evidence and sample-level condition effects
- `sn_run_regulatory_activity()`: infer DoRothEA TF activity or PROGENy pathway activity with decoupleR
- `sn_store_regulatory_activity()`: persist regulatory activity results
- `sn_get_regulatory_activity_result()`: retrieve regulatory activity results

## Interpretation and Reporting

- `sn_list_methods()` / `sn_method_status()`: discover registered current and
  roadmap backends, their default status, runtime, dependencies, install action,
  requirements, outputs, and current availability
- `sn_store_result()` / `sn_get_result()` / `sn_delete_result()`: manage any
  versioned Shennong result type through the generic analysis-result contract
- `sn_validate_result()`: validate schema, typed containers, and provenance
- `sn_list_results()`: list registered and generic stored results, optionally
  filtered by analysis `type`
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
- `sn_plot_feature()`: feature plot or density-style feature map; rasterized points preserve `pt_size`, with a vector fallback when `ggrastr` is unavailable
- `sn_plot_heatmap()`: selected-gene heatmap in cell-level or group-averaged mode
- `sn_plot_violin()`: violin plot
- `sn_plot_dot()`: dot plot; pass a named feature list for separate free-width
  marker-group panels
- `sn_plot_boxplot()`: box plot
- `sn_plot_barplot()`: bar plot
- `sn_plot_composition()`: composition bar charts
- `sn_plot_milo()`: milo differential-abundance visualization
- `sn_list_palettes()`: list available palettes
- `sn_get_palette()`: resolve palette colors

## Utilities and Project Setup

- `sn_check_version()`: check package version against CRAN or a GitHub `source` / `ref`
- `sn_install_shennong()`: install Shennong from CRAN, GitHub `source` / `ref`, or a local `source`
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
  use `sn_list_results()` plus `sn_get_result()` or a specialized
  `sn_get_*_result()`
- If the task is annotation or interpretation:
  use `sn_prepare_*_evidence()` and `sn_interpret_*()`
- If the task is project bootstrap:
  use `sn_initialize_project()`
