
# Shennong

<!-- badges: start -->

[![R-CMD-check](https://github.com/zerostwo/shennong/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/zerostwo/shennong/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/zerostwo/shennong/branch/main/graph/badge.svg)](https://app.codecov.io/gh/zerostwo/shennong?branch=main)
[![lifecycle](https://img.shields.io/badge/lifecycle-Experimental-important.svg)](https://lifecycle.r-lib.org/articles/stages.html)
<!-- badges: end -->

<!-- README.md is generated from README.Rmd. Please edit that file -->

`Shennong` is an experimental R package for single-cell, multimodal,
spatial, and bulk transcriptomics workflows. Seurat objects remain the
primary single-cell contract, while standalone bulk analyses accept
ordinary matrices and `SummarizedExperiment` objects. The package
provides reproducible entry points for preprocessing, integration,
annotation, differential testing, biological-state modeling, publication
figures, and interpretation-ready result storage.

## Choose A Research Path

The documentation is organized around scientific decisions, not isolated
function names. Start with the [research workflow
map](https://songqi.org/shennong/dev/articles/research-workflow-map.html),
then follow the path that matches the data-generating process:

| Research setting | Evidence path | Real-data narrative |
|----|----|----|
| Longitudinal or multi-sample single-cell study | QC → baseline clustering → integration diagnostics → annotation → sample-level change | Kotliarov PBMC CITE-seq |
| Tumour ecosystem and clinical translation | malignant-state evidence → programs/GRNs → pathways → bulk cohort and survival | GSE72056 + TCGA-SKCM |
| Ordered cell-state transition | topology → gene trends → velocity → terminal-state probabilities | Hermann spermatogenesis |
| Spatial tissue organization | coordinates → spatial features → domains/neighborhoods → mapping | Visium lymph node |

Every real-data article is check-safe: it executes only when a local
fixture, its optional dependencies, and `SHENNONG_RUN_VIGNETTES=true`
are all present. Absent backends remain visibly unavailable instead of
producing substitute results.

## Installation

Install the current development version from GitHub:

``` r
install.packages("remotes")
remotes::install_github("zerostwo/shennong")
```

List Shennong’s required and recommended R packages, then install the
missing ones in one step:

``` r
deps <- sn_list_dependencies()
deps

sn_install_dependencies(scope = "required")
```

## One-Command Analysis Software

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
| Data import and storage | `sn_read()`, `sn_write()`, `sn_initialize_seurat_object()`, `sn_set_layer_backend()` | 10x Genomics, STARsolo, H5/H5AD, AnnData, BPCells, qs2, GMT and rio | Native R and format adapters |
| Preprocessing and QC | `sn_normalize_data()`, `sn_find_doublets()`, `sn_remove_ambient_contamination()` | Seurat log-normalization, SCTransform/glmGamPoi, scran, scDblFinder, SoupX, decontX and HGNChelper | Native R |
| Clustering and batch integration | `sn_run_cluster()` | Seurat CCA/RPCA, Harmony, Coralysis, scVI and scANVI; Louvain, multilevel Louvain, SLM and Leiden clustering | Native R or managed pixi |
| CITE-seq and multimodal integration | `sn_run_multimodal()` | Seurat WNN, totalVI, Coralysis and MMoCHi | Native R or managed pixi |
| Reference mapping and simulation | `sn_transfer_labels()`, `sn_simulate()` | Seurat anchors, Coralysis, scANVI, scArches, scPoli and scDesign3 | Native R or managed pixi |
| Cell-type annotation | `sn_run_annotation()` | SingleR, CellTypist, Seurat label transfer, Symphony, scmap and scANVI, with Cell Ontology mapping | Native R, CLI or managed pixi |
| Marker and differential-expression analysis | `sn_find_de()` | Single-cell Seurat tests including Wilcoxon and COSG; pseudobulk and standalone bulk DE with DESeq2, edgeR, limma/limma-voom and dream/variancePartition | Native R; input type selects single-cell or bulk automatically; `sn_find_bulk_de()` remains compatible |
| Enrichment and regulatory activity | `sn_run_enrichment()`, `sn_run_regulatory_activity()` | clusterProfiler ORA/GSEA, GO, KEGG, MSigDB/msigdbr, decoupleR, DoRothEA and PROGENy | Native R |
| Gene-set scoring and program discovery | `sn_score_programs()`, `sn_discover_programs()` | UCell, AUCell, GSVA, ssGSEA, sparse mean scoring, multi-restart NMF, cNMF and Hotspot | Native R; cNMF/Hotspot adapters |
| Gene-regulatory networks | `sn_run_grn()` | GENIE3, pySCENIC, R SCENIC and GRNBoost2/arboreto | Native R for GENIE3; external-result adapters for SCENIC/GRNBoost2 |
| Trajectory and cell dynamics | `sn_run_trajectory()`, `sn_run_velocity()`, `sn_run_fate()` | Slingshot, Monocle 3, Palantir, tradeSeq, scVelo, RegVelo and CellRank | Native R or managed pixi |
| Differential abundance and state prioritization | `sn_test_abundance()`, `sn_prioritize_states()`, `sn_run_scissor()` | Propeller/speckle, Milo/miloR, scCODA/pertpy, sample-aware permutation, Augur-inspired prioritization, Scissor and RareQ | Native R, managed runner or adapter |
| Cell-cell communication | `sn_run_cell_communication()` | LIANA, CellChat, CellPhoneDB, NicheNet and MultiNicheNet, including cross-method consensus | Native R or managed pixi |
| CNV and malignant-state analysis | `sn_run_cnv()` | infercnvpy and CopyKAT, with malignancy scoring, subclones and chromosome summaries | Managed pixi or native R |
| Metabolic analysis | `sn_run_metabolism()` | UCell, GSVA, ssGSEA, mean scoring, scMetabolism, scFEA and Compass | Native R; scFEA/Compass runner-result adapters |
| Spatial transcriptomics | `sn_run_spatial()` | Moran’s I/Squidpy, nnSVG, BANKSY, stLearn, cell2location, Tangram, SPARK-X, BayesSpace, CellCharter, STAligner and Harmony | Native R, managed pixi or adapter |
| Bulk transcriptomics and clinical analysis | `sn_run_bulk()`, `sn_run_survival()`, `sn_run_bulk_deconvolution()` | edgeR, DESeq2, limma/limma-voom, dream, GSVA/ssGSEA, WGCNA, adjusted Cox, Kaplan-Meier/log-rank, proportional-hazards diagnostics, BayesPrism and CIBERSORTx | Native R or local container backend |
| Integration diagnostics | `sn_assess_integration()` | LISI, silhouette, graph connectivity, PCR batch effect, clustering agreement, isolated-label score, entropy, purity and ROGUE | Native R |
| Result validation and migration | `sn_validate_result()`, `sn_audit_results()`, `sn_upgrade_results()` | Versioned `2.0.0` analytical-result envelope, canonical primary tables, explicit legacy migration and registered artifact audit | Native R; audit is read-only and migration is explicit |
| Optional R acceleration | `sn_check_acceleration()`, `sn_with_acceleration()` | Companion ShennongOpt package with guarded operation-scoped patches for Seurat RunPCA/ScaleData/FindVariableFeatures/FindNeighbors/FindClusters, SeuratObject merge/JoinLayers, scran, decontX, scDblFinder, Coralysis, UCell, LISI and Rogue | Same scientific API with optional acceleration; uncovered hot paths run plain upstream |
| Local and remote runtime observability | `sn_enable_usage_tracking()`, `sn_summarize_usage()` | Opt-in timing/call counts for all non-control exports, explicit mode, sanitized parameter fingerprints and consent-gated DBI delivery | SQLite outbox first; managed remote sync is explicit and fail-open |
| Publication figures and reporting | `sn_get_figure_spec()`, `sn_export_figure_bundle()`, `sn_write_results()` | ggplot2/patchwork, ggrastr, SVG/TIFF/PDF/PNG export, source-data bundles and optional ellmer-backed interpretation | Native R with optional LLM provider |

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

Analytical outputs use `schema_version = "2.0.0"`, one canonical
`(analysis_type, result_id)` identity, and their principal data frame in
`tables$primary`. All 33 built-in result types have minimum semantic checks for
their identities and numerical evidence. Future schemas, malformed primary
tables, and storage-key mismatches fail closed; audit older objects before
reusing stored evidence:

``` r
audit <- sn_audit_results(object)
object <- sn_upgrade_results(object) # only after reviewing legacy/invalid rows
sn_list_results(object)
```

Optional R acceleration is provided by the companion **ShennongOpt**
package and is never activated during package loading. When an
accelerated workflow first reaches a covered hot path—Seurat RunPCA,
ScaleData, FindVariableFeatures, FindNeighbors, FindClusters, scran,
decontX, scDblFinder, Coralysis, UCell, LISI, or Rogue—Shennong checks
the installed ShennongOpt build and the upstream package version, and
activates only validated, operation-scoped patches for that single
workflow call; the pre-call patch state is restored afterward.
SeuratObject merge/JoinLayers patches are available through explicit
`sn_enable_acceleration()` requests. Uncovered hot paths always run
plain upstream, so acceleration never changes scientific contracts.

``` r
sn_check_acceleration()

pbmc <- sn_find_doublets(pbmc)

# Session-wide opt-out for automatic acceleration.
options(shennong.acceleration = FALSE)
# Environment alternative: SHENNONG_ACCELERATION_DISABLED=true.
```

The opt-outs do not deactivate a patch that you activated manually;
explicit `sn_enable_acceleration()` and `sn_with_acceleration()` calls
are intentional overrides that ignore them. Thread budgets can be set
through the `shennong.opt.threads` option or the `SHENNONG_OPT_THREADS`
environment variable. Shennong does not install ShennongOpt, its
upstream packages, or any Python environment.

Workflow timing and usage statistics are also opt-in. Supply a private
local SQLite path and explicitly distinguish development from real
production runs:

``` r
sn_enable_usage_tracking(
  "analysis-private/shennong-usage.sqlite",
  mode = "development",
  display = TRUE
)

object <- sn_run_cluster(
  object,
  batch = "sample",
  integration_method = "harmony",
  resolution = 0.6
)

sn_disable_usage_tracking()

# Which methods are common, and which consume the most total time?
sn_summarize_usage(
  "analysis-private/shennong-usage.sqlite",
  sort_by = "calls"
)
sn_summarize_usage(
  "analysis-private/shennong-usage.sqlite",
  sort_by = "total_seconds"
)
```

The recorder instruments every export except its ten usage-control APIs
(265 of the current 275 exports): analysis, plotting, get/list, storage,
IO, validation, installation, administration, rio adapters and low-level
backend calls. Its ten usage-control APIs are excluded to avoid
recursive writes. Instrumentation replaces namespace bindings only after
`sn_enable_usage_tracking()`; a function reference saved or imported
before then is not intercepted. Enable at process startup for a usage
study and confirm this boundary with
`sn_check_usage_tracking()$covers_preexisting_function_references`. It
records call/configuration ordinals, status and timing while redacting
objects, matrices, gene/signature/cell/sample/patient values, paths,
credentials, prompts and free text. Loading Shennong never creates a
usage database. Managed deployments can use an explicit
`remote_research` consent, local SQLite outbox and
`sn_flush_usage_tracking()` to deliver sanitized rows through a DBI
connection factory. See the runtime observability article for remote
security, nesting, privacy and acceleration evidence semantics.

The generated public-API parameter inventory plus clustering and
enrichment method matrices are static completeness/admission artifacts.
They prove that functions, formals, selectors, dispatch cells and
required cases are classified; they do not prove that every case has run
or passed an upstream conformance comparison.

## Agent And MCP Integration

Shennong ships installable Agent Skills plus a read-only MCP server. The
MCP surface lets an agent discover registered methods, inspect exact
installed R help, and retrieve workflow recipes; it does not execute
arbitrary R code or modify analysis files.

``` r
# Install all packaged Shennong usage skills for local agents.
sn_install_codex_skill(path = "~/.agents/skills", type = "package_skills")

# Use this command/argument pair in any stdio-capable MCP client.
sn_get_mcp_server_config()

# Equivalent direct server command:
# Rscript -e 'Shennong::sn_run_mcp_server()'
```

## Analysis Data

Dataset discovery, download, caching, and publication are owned by the
`ShennongData` package. Shennong accepts materialized matrices, file
paths, and Seurat objects; it no longer bundles analysis datasets or
exports a data distribution API. The executable website uses small,
provenance-tracked public data subsets materialized under the local
`SHENNONG_REAL_DATA_DIR` cache. Those data files are intentionally
excluded from this repository.

## Quick Start

``` r
library(Shennong)

pbmc <- qs2::qs_read(file.path(
  Sys.getenv("SHENNONG_REAL_DATA_DIR"),
  "single-cell",
  "kotliarov_pbmc.qs2"
))

pbmc <- sn_filter_cells(
  pbmc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  plot = FALSE
)
pbmc <- sn_filter_genes(pbmc, min_cells = 3, plot = FALSE)
pbmc <- sn_run_cluster(
  pbmc,
  normalization_method = "seurat",
  resolution = 0.6
)

sn_plot_dim(pbmc, group_by = "seurat_clusters", label = TRUE)
```

## Integration Example

The same public PBMC fixture can be used for batch-aware integration:

``` r
pbmc_integrated <- sn_run_cluster(
  pbmc,
  batch = "real_batch",
  normalization_method = "seurat",
  resolution = 0.6
)

sn_plot_dim(pbmc_integrated, group_by = "real_batch")
sn_plot_dim(pbmc_integrated, group_by = "seurat_clusters", label = TRUE)
```

You can summarize integration quality and surface rare or difficult
groups directly from the integrated object:

``` r
metrics <- sn_assess_integration(
  pbmc_integrated,
  batch_by = "real_batch",
  cluster_by = "seurat_clusters",
  reduction = "harmony",
  baseline_reduction = "pca"
)

metrics$summary
metrics$per_group$isolated_label_score
metrics$per_group$cluster_entropy
metrics$per_group$cluster_purity
metrics$per_group$challenging_groups
```

## Differential Expression And Enrichment

``` r
pbmc <- sn_find_de(
  pbmc,
  analysis = "markers",
  group_by = "seurat_clusters",
  layer = "data",
  store_name = "cluster_markers",
  return_object = TRUE,
  verbose = FALSE
)

pbmc <- sn_run_enrichment(
  x = pbmc,
  source_de_name = "cluster_markers",
  gene_clusters = gene ~ cluster,
  database = c("GOBP", "H"),
  species = "human",
  universe = rownames(pbmc),
  store_name = "cluster_pathways",
  pvalue_cutoff = 0.05
)
```

The same `sn_find_de()` entry point accepts a feature-by-sample matrix,
list, or `SummarizedExperiment` for standalone bulk analysis:

``` r
bulk_de <- sn_find_de(
  counts,
  metadata = sample_data,
  design = ~ batch + condition,
  contrast = c("condition", "tumor", "normal"),
  method = "auto"
)
```

## Documentation

Longer workflow articles are available in the package site and
vignettes, grouped by the decision they support:

- Build the study: [data and
  projects](https://songqi.org/shennong/dev/articles/data-io-projects.html),
  [preprocessing and
  QC](https://songqi.org/shennong/dev/articles/preprocessing-qc.html),
  [clustering and
  integration](https://songqi.org/shennong/dev/articles/clustering.html),
  and [annotation and
  pathways](https://songqi.org/shennong/dev/articles/annotation-pathways.html).
- Explain biological change:
  [composition](https://songqi.org/shennong/dev/articles/composition-analysis.html),
  [differential
  abundance](https://songqi.org/shennong/dev/articles/abundance-priority.html),
  [communication](https://songqi.org/shennong/dev/articles/communication-consensus.html),
  [trajectory](https://songqi.org/shennong/dev/articles/trajectory-dynamics.html),
  and [spatial
  analysis](https://songqi.org/shennong/dev/articles/spatial-workflows.html).
- Translate and report: [bulk
  transcriptomics](https://songqi.org/shennong/dev/articles/bulk-transcriptomics.html),
  [stored-result
  contracts](https://songqi.org/shennong/dev/articles/analysis-results-and-methods.html),
  [publication
  figures](https://songqi.org/shennong/dev/articles/publication-figures.html),
  and [runtime
  observability](https://songqi.org/shennong/dev/articles/runtime-observability.html).

## Status

`Shennong` is still experimental. The package currently prioritizes a
clean and consistent workflow surface over backward compatibility across
early versions.
