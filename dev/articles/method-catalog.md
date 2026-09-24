# Choose a backend

Start with [Get
started](https://zerostwo.github.io/shennong/dev/articles/get-started.md)
for a complete small example. Use this catalog when choosing software
for a specific analysis. A listed adapter is not a claim that the
backend is installed or validated for every input; check its
requirements and recorded evidence before running it.

## Workflows and execution requirements

Shennong provides stable `sn_*` entry points over R packages,
command-line programs, and Python workflows. The table below summarizes
the current analysis surface by data-analysis module. “Managed pixi/CLI”
means Shennong prepares the input, runs the backend in a
project-independent environment, and imports the result. “Adapter” means
Shennong standardizes an existing result or a result returned by a
user-supplied runner; the external software is not silently installed or
executed.

| Analysis module | Main Shennong entry point | Supported software and methods | Execution model |
|----|----|----|----|
| Data import and storage | [`sn_read()`](https://zerostwo.github.io/shennong/dev/reference/sn_read.md), [`sn_write()`](https://zerostwo.github.io/shennong/dev/reference/sn_write.md), [`sn_initialize_seurat_object()`](https://zerostwo.github.io/shennong/dev/reference/sn_initialize_seurat_object.md), [`sn_set_layer_backend()`](https://zerostwo.github.io/shennong/dev/reference/sn_set_layer_backend.md) | 10x Genomics, STARsolo, H5/H5AD, AnnData, BPCells, qs2, GMT and rio | Native R and format adapters |
| Preprocessing and QC | [`sn_normalize_data()`](https://zerostwo.github.io/shennong/dev/reference/sn_normalize_data.md), [`sn_find_doublets()`](https://zerostwo.github.io/shennong/dev/reference/sn_find_doublets.md), [`sn_remove_ambient_contamination()`](https://zerostwo.github.io/shennong/dev/reference/sn_remove_ambient_contamination.md) | Seurat log-normalization, SCTransform/glmGamPoi, scran, scDblFinder, SoupX, decontX and HGNChelper | Native R |
| Clustering and batch integration | [`sn_run_cluster()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_cluster.md) | Seurat CCA/RPCA, Harmony, Coralysis, scVI and scANVI; Louvain, multilevel Louvain, SLM and Leiden clustering | Native R or managed pixi |
| CITE-seq and multimodal integration | [`sn_run_multimodal()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_multimodal.md) | Seurat WNN, totalVI, Coralysis and MMoCHi | Native R or managed pixi |
| Reference mapping and simulation | [`sn_transfer_labels()`](https://zerostwo.github.io/shennong/dev/reference/sn_transfer_labels.md), [`sn_simulate()`](https://zerostwo.github.io/shennong/dev/reference/sn_simulate.md) | Seurat anchors, Coralysis, scANVI, scArches, scPoli and scDesign3 | Native R or managed pixi |
| Cell-type annotation | [`sn_run_annotation()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_annotation.md) | SingleR, CellTypist, Seurat label transfer, Symphony, scmap and scANVI, with Cell Ontology mapping | Native R, CLI or managed pixi |
| Marker and differential-expression analysis | [`sn_find_de()`](https://zerostwo.github.io/shennong/dev/reference/sn_find_de.md) | Single-cell Seurat tests including Wilcoxon and COSG; pseudobulk and standalone bulk DE with DESeq2, edgeR, limma/limma-voom and dream/variancePartition | Native R; input type selects single-cell or bulk automatically; [`sn_find_bulk_de()`](https://zerostwo.github.io/shennong/dev/reference/sn_find_bulk_de.md) remains compatible |
| Enrichment and regulatory activity | [`sn_run_enrichment()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_enrichment.md), [`sn_run_regulatory_activity()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_regulatory_activity.md) | clusterProfiler ORA/GSEA, GO, KEGG, MSigDB/msigdbr, decoupleR, DoRothEA and PROGENy | Native R |
| Gene-set scoring and program discovery | [`sn_score_programs()`](https://zerostwo.github.io/shennong/dev/reference/sn_score_programs.md), [`sn_discover_programs()`](https://zerostwo.github.io/shennong/dev/reference/sn_discover_programs.md) | UCell, AUCell, GSVA, ssGSEA, sparse mean scoring, multi-restart NMF, cNMF and Hotspot | Native R; cNMF/Hotspot adapters |
| Gene-regulatory networks | [`sn_run_grn()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_grn.md) | GENIE3, pySCENIC, R SCENIC and GRNBoost2/arboreto | Native R for GENIE3; external-result adapters for SCENIC/GRNBoost2 |
| Trajectory and cell dynamics | [`sn_run_trajectory()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_trajectory.md), [`sn_run_velocity()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_velocity.md), [`sn_run_fate()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_fate.md) | Slingshot, Monocle 3, Palantir, tradeSeq, scVelo, RegVelo and CellRank | Native R or managed pixi |
| Differential abundance and state prioritization | [`sn_test_abundance()`](https://zerostwo.github.io/shennong/dev/reference/sn_test_abundance.md), [`sn_prioritize_states()`](https://zerostwo.github.io/shennong/dev/reference/sn_prioritize_states.md), [`sn_run_scissor()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_scissor.md) | Propeller/speckle, Milo/miloR, scCODA/pertpy, sample-aware permutation, Augur-inspired prioritization, Scissor and RareQ | Native R, managed runner or adapter |
| Cell-cell communication | [`sn_run_cell_communication()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_cell_communication.md) | LIANA, CellChat, CellPhoneDB, NicheNet and MultiNicheNet, including cross-method consensus | Native R or managed pixi |
| CNV and malignant-state analysis | [`sn_run_cnv()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_cnv.md) | infercnvpy and CopyKAT, with malignancy scoring, subclones and chromosome summaries | Managed pixi or native R |
| Metabolic analysis | [`sn_run_metabolism()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_metabolism.md) | UCell, GSVA, ssGSEA, mean scoring, scMetabolism, scFEA and Compass | Native R; scFEA/Compass runner-result adapters |
| Spatial transcriptomics | [`sn_run_spatial()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_spatial.md) | Moran’s I/Squidpy, nnSVG, BANKSY, stLearn, cell2location, Tangram, SPARK-X, BayesSpace, CellCharter, STAligner and Harmony | Native R, managed pixi or adapter |
| Bulk transcriptomics and clinical analysis | [`sn_run_bulk()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_bulk.md), [`sn_run_survival()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_survival.md), [`sn_run_bulk_deconvolution()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_bulk_deconvolution.md) | edgeR, DESeq2, limma/limma-voom, dream, GSVA/ssGSEA, WGCNA, adjusted Cox, Kaplan-Meier/log-rank, proportional-hazards diagnostics, BayesPrism and CIBERSORTx | Native R or local container backend |
| Integration diagnostics | [`sn_assess_integration()`](https://zerostwo.github.io/shennong/dev/reference/sn_assess_integration.md) | LISI, silhouette, graph connectivity, PCR batch effect, clustering agreement, isolated-label score, entropy, purity and ROGUE | Native R |
| Result validation and migration | [`sn_validate_result()`](https://zerostwo.github.io/shennong/dev/reference/sn_validate_result.md), [`sn_audit_results()`](https://zerostwo.github.io/shennong/dev/reference/sn_audit_results.md), [`sn_upgrade_results()`](https://zerostwo.github.io/shennong/dev/reference/sn_upgrade_results.md) | Versioned `2.0.0` analytical-result envelope, canonical primary tables, explicit legacy migration and registered artifact audit | Native R; audit is read-only and migration is explicit |
| Optional R acceleration | [`sn_check_acceleration()`](https://zerostwo.github.io/shennong/dev/reference/sn_check_acceleration.md), [`sn_with_acceleration()`](https://zerostwo.github.io/shennong/dev/reference/sn_with_acceleration.md) | Companion ShennongOpt package with guarded operation-scoped patches for Seurat RunPCA/ScaleData/FindVariableFeatures/FindNeighbors/FindClusters, SeuratObject merge/JoinLayers, scran, decontX, scDblFinder, Coralysis, UCell, LISI and Rogue | Same scientific API with optional acceleration; uncovered hot paths run plain upstream |
| Local and remote runtime observability | [`sn_enable_usage_tracking()`](https://zerostwo.github.io/shennong/dev/reference/sn_enable_usage_tracking.md), [`sn_summarize_usage()`](https://zerostwo.github.io/shennong/dev/reference/sn_summarize_usage.md) | Opt-in timing/call counts for all non-control exports, explicit mode, sanitized parameter fingerprints and consent-gated DBI delivery | SQLite outbox first; managed remote sync is explicit and fail-open |
| Publication figures and reporting | [`sn_get_figure_spec()`](https://zerostwo.github.io/shennong/dev/reference/sn_get_figure_spec.md), [`sn_export_figure_bundle()`](https://zerostwo.github.io/shennong/dev/reference/sn_export_figure_bundle.md), [`sn_write_results()`](https://zerostwo.github.io/shennong/dev/reference/sn_write_results.md) | ggplot2/patchwork, ggrastr, SVG/TIFF/PDF/PNG export, source-data bundles and optional ellmer-backed interpretation | Native R with optional LLM provider |

All methods shipped in `inst/methods/` currently have an implemented
Shennong entry point or explicit adapter. Optional R packages,
command-line programs, pixi environments, credentials, references, and
model files are still required when the selected backend depends on
them. Inspect the registry and the current machine before starting a
workflow:

``` r

# Every registered backend and whether it can run in the current session
sn_list_methods()
sn_list_methods(task = "trajectory")
sn_list_methods(available = TRUE)

# Runtime, dependency, install action, inputs and outputs for one backend
sn_get_method_status("cellrank", task = "fate")

# Install missing R dependencies or prepare a managed Python environment
sn_install_dependencies(scope = "recommended")
sn_prepare_pixi_environment("trajectory", install_environment = TRUE)
```

For shared argument names and backend control lists, see [Parameters and
results](https://zerostwo.github.io/shennong/dev/articles/parameters-and-results.md).
For result schemas, validation, and migration, see [Manage stored
results](https://zerostwo.github.io/shennong/dev/articles/analysis-results-and-methods.md).
