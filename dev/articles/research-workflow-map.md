# From samples to evidence: a research workflow map

For your first runnable analysis, use [Get
started](https://zerostwo.github.io/shennong/dev/articles/get-started.md).
Use this map to choose a larger study workflow after you know the core
API.

Start with the question

Shennong’s functions form an evidence chain rather than a menu of
unrelated methods. Start with the biological unit of replication,
preserve the source assay and layer, compare plausible methods, and
finish with a stored table and an inspectable figure.

observations**→** quality control**→** representation**→** biological
comparison**→** stored evidence**→** figure + provenance

Use this page to choose a route. The linked articles contain the
executable real-data analysis and its visual checkpoints.

## What the examples do—and do not—claim

The website uses locally materialized subsets of public observations.
They are not bundled in the package and the articles never download them
during `R CMD check`. A normal build shows the code and input contract.
A designated real-data build executes an article only when its fixture
and dependencies are present:

``` sh
SHENNONG_RUN_VIGNETTES=true \
SHENNONG_REAL_DATA_DIR=/path/to/pkgdown-real \
Rscript -e 'pkgdown::build_site()'
```

core uses a bounded local fixture and the default workflow. extended
additionally requires `SHENNONG_REAL_PROFILE=all` and optional
dependencies. External services, credentials, model files, or imported
backend results remain explicitly opt-in. A displayed plot is therefore
evidence from an executed chunk, not a placeholder manufactured when a
backend is absent.

## Choose a scientific narrative

### Longitudinal immune cohort

Use the real Kotliarov PBMC CITE-seq subset to ask how RNA and ADT cell
states, batch, and vaccine response vary across 20 biological samples.
Begin with
[preprocessing](https://zerostwo.github.io/shennong/dev/articles/preprocessing-qc.md),
continue through
[integration](https://zerostwo.github.io/shennong/dev/articles/clustering.md),
[annotation and
pathways](https://zerostwo.github.io/shennong/dev/articles/annotation-pathways.md),
then test
[composition](https://zerostwo.github.io/shennong/dev/articles/composition-analysis.md)
or [differential
abundance](https://zerostwo.github.io/shennong/dev/articles/abundance-priority.md).

### Tumour ecosystem and clinical translation

Use author-normalized malignant cells from GSE72056 for programs,
regulatory networks, malignancy evidence, and metabolism. Connect those
hypotheses to TCGA-SKCM expression and follow-up without relabelling TPM
as counts. Follow [program
discovery](https://zerostwo.github.io/shennong/dev/articles/program-discovery-grn.md),
[CNV and
metabolism](https://zerostwo.github.io/shennong/dev/articles/cnv-metabolism.md),
and [bulk clinical
analysis](https://zerostwo.github.io/shennong/dev/articles/bulk-transcriptomics.md).

### Cell-state dynamics

Use the Hermann spermatogenesis fixture when the question is ordered
change, not only discrete clusters. Separate trajectory topology, gene
trends, velocity, and terminal-state probabilities. Follow [trajectory
and dynamic
genes](https://zerostwo.github.io/shennong/dev/articles/trajectory-dynamics.md),
then the optional [velocity and
fate](https://zerostwo.github.io/shennong/dev/articles/velocity-fate.md)
path.

### Tissue organization

Use the real Visium lymph-node sections to keep expression and
coordinates together. Inspect native neighborhoods and spatial features
before introducing an extended backend or reference mapping. Follow the
[spatial
workflow](https://zerostwo.github.io/shennong/dev/articles/spatial-workflows.md).

## The core single-cell evidence chain

1.  **Define the observation and replicate.** Materialize the source
    once, retain sample metadata, and decide whether the scientific
    replicate is a donor, capture, section, or bulk sample.
2.  **Audit quality before filtering.** Plot distributions and
    thresholds, record flags, and retain the pre-filter evidence needed
    to explain cell and feature loss.
3.  **Choose a representation.** Cluster without integration as a
    baseline, then compare batch-aware methods against both mixing and
    biological separation metrics.
4.  **Name and test states.** Join marker evidence, reference evidence,
    signatures, and ontology terms. Test replicate-aware changes instead
    of treating cells as independent samples.
5.  **Store the result, then plot it.** Retrieve the exact result name
    from the object so tables, figures, prompts, and exported source
    data all refer to the same analysis.
6.  **Report runtime and provenance.** Record backend versions,
    parameter fingerprints, and optional acceleration separately from
    the scientific result.

### 1. Materialize once and retain sample identity

Dataset discovery and publication belong to ShennongData. Shennong
begins with a materialized matrix, file, or object and preserves the
metadata that defines the study design. See [Data input, output, and
project
setup](https://zerostwo.github.io/shennong/dev/articles/data-io-projects.md).

``` r

library(Shennong)

pbmc <- qs2::qs_read(file.path(
  Sys.getenv("SHENNONG_REAL_DATA_DIR"),
  "single-cell",
  "kotliarov_pbmc.qs2"
))

table(pbmc$real_sample, pbmc$real_response)
```

### 2. Make QC decisions visible

The filtering calls return an analysis object; plotting functions expose
the decision boundary and retained fraction. The complete real-data
sequence is in [Preprocessing and
QC](https://zerostwo.github.io/shennong/dev/articles/preprocessing-qc.md).

``` r

pbmc <- sn_filter_cells(
  pbmc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  plot = FALSE
)
pbmc <- sn_filter_genes(pbmc, min_cells = 3, plot = FALSE)

sn_plot_qc_thresholds(
  pbmc,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
)
```

### 3. Establish a baseline before integration

[`sn_run_cluster()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_cluster.md)
exposes the normalization and integration choices in one call. Compare
batch mixing and biological conservation; a UMAP alone is not an
integration test. See [Clustering and
integration](https://zerostwo.github.io/shennong/dev/articles/clustering.md),
[layer-aware
workflows](https://zerostwo.github.io/shennong/dev/articles/layered-workflows.md),
and [integration
metrics](https://zerostwo.github.io/shennong/dev/articles/metrics-diagnostics.md).

``` r

pbmc <- sn_run_cluster(
  pbmc,
  normalization_method = "seurat",
  integration_method = "harmony",
  batch_by = "real_batch",
  resolution = 0.6
)

sn_plot_dim(pbmc, reduction = "umap", group_by = "real_batch")
sn_plot_dim(pbmc, reduction = "umap", group_by = "seurat_clusters", label = TRUE)

assessment <- sn_assess_integration(
  pbmc,
  batch_by = "real_batch",
  cluster_by = "seurat_clusters",
  reduction = "harmony",
  baseline_reduction = "pca"
)
sn_plot_integration(assessment)
```

### 4. Connect states to genes and pathways

Store marker results on the object, then let downstream plots and
enrichment retrieve them by name. ORA and GSEA answer different
questions: ORA tests a selected gene set against an explicit universe,
whereas GSEA uses the complete ranked statistic. See [Markers,
signatures, pathways, and
annotation](https://zerostwo.github.io/shennong/dev/articles/annotation-pathways.md).

``` r

pbmc <- sn_find_de(
  pbmc,
  analysis = "markers",
  group_by = "seurat_clusters",
  result_id = "cluster_markers",
  return_object = TRUE,
  verbose = FALSE
)

sn_plot_dot(pbmc, features = "top_markers", result_id = "cluster_markers")

pbmc <- sn_run_enrichment(
  pbmc,
  source_de_result_id = "cluster_markers",
  gene_clusters = gene ~ cluster,
  database = c("GOBP", "H"),
  universe = rownames(pbmc),
  result_id = "cluster_pathways",
  return_object = TRUE
)

pathways <- sn_get_enrichment_result(
  pbmc,
  result_id = "cluster_pathways.GOBP"
)
sn_plot_enrichment(pathways, type = "dot")
```

### 5. Return to the sample as the replicate

Cell counts are descriptive until they are summarized within biological
samples. Compare sample-level proportions or use a
differential-abundance model with an explicit design. See [Composition
analysis](https://zerostwo.github.io/shennong/dev/articles/composition-analysis.md)
and [Differential abundance and state
priority](https://zerostwo.github.io/shennong/dev/articles/abundance-priority.md).

``` r

composition <- sn_calculate_composition(
  pbmc,
  group_by = c("real_sample", "real_response"),
  variable = "seurat_clusters"
)

sn_plot_composition(
  pbmc,
  x = real_response,
  fill = seurat_clusters,
  type = "sample_boxplot",
  sample_by = "real_sample"
)
```

### 6. Extend only when the question requires it

The common object and result contract supports several distinct
extensions:

- [cell-cell
  communication](https://zerostwo.github.io/shennong/dev/articles/communication-consensus.md)
  for sender-receiver evidence and cross-method consensus;
- [program discovery and
  GRNs](https://zerostwo.github.io/shennong/dev/articles/program-discovery-grn.md)
  for latent programs and regulatory hypotheses;
- [trajectory](https://zerostwo.github.io/shennong/dev/articles/trajectory-dynamics.md)
  and
  [velocity/fate](https://zerostwo.github.io/shennong/dev/articles/velocity-fate.md)
  for directed or ordered change;
- [spatial
  analysis](https://zerostwo.github.io/shennong/dev/articles/spatial-workflows.md)
  for coordinate-aware expression;
- [bulk
  deconvolution](https://zerostwo.github.io/shennong/dev/articles/bulk-deconvolution.md)
  and [bulk
  transcriptomics](https://zerostwo.github.io/shennong/dev/articles/bulk-transcriptomics.md)
  for cohort-level translation.

Each extension has its own admissible input. A normalized tumour matrix
is not silently treated as UMI counts, and an absent external backend
never produces a substitute result.

### 7. Retrieve, validate, and publish the same evidence

Discover a stored result before retrieving it, validate its versioned
envelope, then export the plot together with source data and a figure
specification. See [result and method
contracts](https://zerostwo.github.io/shennong/dev/articles/analysis-results-and-methods.md),
[publication
figures](https://zerostwo.github.io/shennong/dev/articles/publication-figures.md),
and [interpretation and
reporting](https://zerostwo.github.io/shennong/dev/articles/interpretation-workflows.md).

``` r

sn_list_results(pbmc)
markers <- sn_get_result(pbmc, type = "de", result_id = "cluster_markers")
sn_validate_result(markers)

figure <- sn_plot_de(markers, type = "effect")
sn_export_figure_bundle(
  figure,
  path = "figures/cluster-markers"
)
```

## Module directory

| Research decision | Principal entry points | Visual checkpoint | Detailed article |
|----|----|----|----|
| Materialize and govern inputs | [`sn_read()`](https://zerostwo.github.io/shennong/dev/reference/sn_read.md), [`sn_write()`](https://zerostwo.github.io/shennong/dev/reference/sn_write.md), [`sn_initialize_project()`](https://zerostwo.github.io/shennong/dev/reference/sn_initialize_project.md) | input and metadata audit tables | [Data and projects](https://zerostwo.github.io/shennong/dev/articles/data-io-projects.md) |
| Filter and normalize | [`sn_filter_cells()`](https://zerostwo.github.io/shennong/dev/reference/sn_filter_cells.md), [`sn_filter_genes()`](https://zerostwo.github.io/shennong/dev/reference/sn_filter_genes.md), [`sn_normalize_data()`](https://zerostwo.github.io/shennong/dev/reference/sn_normalize_data.md) | [`sn_plot_qc_thresholds()`](https://zerostwo.github.io/shennong/dev/reference/sn_plot_qc_thresholds.md), [`sn_plot_qc()`](https://zerostwo.github.io/shennong/dev/reference/sn_plot_qc.md) | [Preprocessing and QC](https://zerostwo.github.io/shennong/dev/articles/preprocessing-qc.md) |
| Cluster or integrate modalities | [`sn_run_cluster()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_cluster.md), [`sn_run_multimodal()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_multimodal.md) | [`sn_plot_dim()`](https://zerostwo.github.io/shennong/dev/reference/sn_plot_dim.md), [`sn_plot_resolution_sweep()`](https://zerostwo.github.io/shennong/dev/reference/sn_plot_resolution_sweep.md) | [Clustering](https://zerostwo.github.io/shennong/dev/articles/clustering.md) |
| Respect assay layers | layer-aware arguments across preprocessing, DE, and scoring | embedding and layer audit | [Layered workflows](https://zerostwo.github.io/shennong/dev/articles/layered-workflows.md) |
| Validate integration | [`sn_assess_integration()`](https://zerostwo.github.io/shennong/dev/reference/sn_assess_integration.md), [`sn_compare_integrations()`](https://zerostwo.github.io/shennong/dev/reference/sn_compare_integrations.md) | [`sn_plot_integration()`](https://zerostwo.github.io/shennong/dev/reference/sn_plot_integration.md) | [Metrics and diagnostics](https://zerostwo.github.io/shennong/dev/articles/metrics-diagnostics.md) |
| Annotate and test features | [`sn_run_annotation()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_annotation.md), [`sn_find_de()`](https://zerostwo.github.io/shennong/dev/reference/sn_find_de.md), [`sn_run_enrichment()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_enrichment.md) | marker dot plots, ORA dots, GSEA curves | [Annotation and pathways](https://zerostwo.github.io/shennong/dev/articles/annotation-pathways.md) |
| Score or discover programs | [`sn_score_programs()`](https://zerostwo.github.io/shennong/dev/reference/sn_score_programs.md), [`sn_discover_programs()`](https://zerostwo.github.io/shennong/dev/reference/sn_discover_programs.md), [`sn_run_grn()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_grn.md) | activity heatmaps and regulon networks | [Programs and GRNs](https://zerostwo.github.io/shennong/dev/articles/program-discovery-grn.md) |
| Compare sample composition | [`sn_calculate_composition()`](https://zerostwo.github.io/shennong/dev/reference/sn_calculate_composition.md), [`sn_test_abundance()`](https://zerostwo.github.io/shennong/dev/reference/sn_test_abundance.md) | composition bars and abundance effects | [Composition](https://zerostwo.github.io/shennong/dev/articles/composition-analysis.md), [abundance](https://zerostwo.github.io/shennong/dev/articles/abundance-priority.md) |
| Model communication | [`sn_run_cell_communication()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_cell_communication.md) | bubble, heatmap, network, ligand-target plot | [Communication](https://zerostwo.github.io/shennong/dev/articles/communication-consensus.md) |
| Follow dynamic states | [`sn_run_trajectory()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_trajectory.md), [`sn_run_velocity()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_velocity.md), [`sn_run_fate()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_fate.md) | topology, trends, velocity, fate | [Trajectory](https://zerostwo.github.io/shennong/dev/articles/trajectory-dynamics.md), [velocity](https://zerostwo.github.io/shennong/dev/articles/velocity-fate.md) |
| Explain tumour states | [`sn_run_cnv()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_cnv.md), [`sn_run_metabolism()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_metabolism.md) | CNV heatmap and pathway activity | [CNV and metabolism](https://zerostwo.github.io/shennong/dev/articles/cnv-metabolism.md) |
| Keep tissue coordinates | [`sn_run_spatial()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_spatial.md) and spatial sub-workflows | features, domains, neighborhoods, mapping | [Spatial workflows](https://zerostwo.github.io/shennong/dev/articles/spatial-workflows.md) |
| Translate to cohorts | [`sn_run_bulk()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_bulk.md), [`sn_find_de()`](https://zerostwo.github.io/shennong/dev/reference/sn_find_de.md), [`sn_run_survival()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_survival.md), [`sn_run_bulk_deconvolution()`](https://zerostwo.github.io/shennong/dev/reference/sn_run_bulk_deconvolution.md) | QC, PCA, DE, pathways, survival, fractions | [Bulk](https://zerostwo.github.io/shennong/dev/articles/bulk-transcriptomics.md), [deconvolution](https://zerostwo.github.io/shennong/dev/articles/bulk-deconvolution.md) |
| Build publication evidence | result retrieval, figure specification, interpretation helpers | exported figure and source-data bundle | [Visualization](https://zerostwo.github.io/shennong/dev/articles/visualization.md), [publication](https://zerostwo.github.io/shennong/dev/articles/publication-figures.md), [interpretation](https://zerostwo.github.io/shennong/dev/articles/interpretation-workflows.md) |
| Profile the workflow | usage tracking and acceleration controls | timing summaries | [Runtime observability](https://zerostwo.github.io/shennong/dev/articles/runtime-observability.md) |

## Read acceleration evidence as four separate gates

An active patch is not proof that a Shennong call used a fast kernel.
The full sequence is:

### 1. Activation

The replacement is installed in a namespace for a bounded scope. This
reports availability, not use.

### 2. Target intersection

The workflow must actually call the replaced upstream symbol. A loaded
patch for an unused backend has no effect.

### 3. Guard accepted

The real call shape, input class, operation, version, and parameters
must match the validated fast-path guard. Otherwise the replacement must
fall back.

### 4. Fast-path hit

Runtime evidence must show that the guarded kernel executed. Shennong
does not infer this from activation alone.

This distinction prevents scientifically unsafe false hits. The pinned
JoinLayers and Coralysis replacements can be active while Shennong’s
real call shape fails their benchmark-specific guards and falls back. An
older broad Seurat `FindIntegrationAnchors` replacement lacked a
reduction guard, so an RPCA request could incorrectly enter a
CCA-oriented path. The tradeSeq patch globally replaced `mgcv` behavior
without an exact dependency guard. Shennong therefore uses strict
source/version checks, operation-specific scopes, and temporary
`with_disabled` boundaries around incompatible calls. Those patches
remain explicit-only until a three-arm contract proves the actual
Shennong workflow. See [Analysis results and
methods](https://zerostwo.github.io/shennong/dev/articles/analysis-results-and-methods.md)
for the measured intersections and quarantine rationale.

## Find the next article from an output

If the output is a table, call
[`sn_list_results()`](https://zerostwo.github.io/shennong/dev/reference/sn_list_results.md)
and retrieve its exact stored name before interpretation. If it is a
figure, keep the source table beside the export. If it is only an
embedding, continue to diagnostics or biological testing before treating
it as evidence. Runtime statistics belong in a separate opt-in database
and never replace the result stored on the analysis object.
