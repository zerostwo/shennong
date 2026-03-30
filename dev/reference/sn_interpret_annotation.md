# Interpret cluster markers for cell-type annotation

Interpret cluster markers for cell-type annotation

## Usage

``` r
sn_interpret_annotation(
  object,
  de_name = NULL,
  cluster_col = "seurat_clusters",
  n_markers = 10,
  marker_selection = c("specific", "top"),
  enrichment_name = NULL,
  n_terms = 5,
  enrichment_selection = c("specific", "top"),
  include_qc = TRUE,
  reduction = "umap",
  n_neighbor_clusters = 3,
  background = NULL,
  annotation_mode = c("single_pass", "agentic"),
  output_format = c("llm", "human"),
  provider = NULL,
  model = NULL,
  reasoning_effort = NULL,
  include_json_schema = TRUE,
  apply_metadata = TRUE,
  metadata_prefix = "sn_annotation",
  metadata_fields = c("primary_label", "broad_label", "confidence", "status",
    "risk_flags"),
  label_candidates = NULL,
  label_style = c("title", "snake", "asis"),
  return_prompt = FALSE,
  store_name = "default",
  return_object = TRUE,
  show_progress = interactive(),
  ...
)
```

## Arguments

- object:

  A `Seurat` object.

- de_name:

  Optional stored marker-result name. When omitted, Shennong prefers
  `"default"`, then a single available result, and otherwise the most
  recent marker result.

- cluster_col:

  Metadata column containing cluster labels.

- n_markers:

  Number of top markers per cluster.

- marker_selection:

  How to choose marker genes for annotation evidence: `"specific"`
  prefers genes that are relatively unique to one cluster, while `"top"`
  keeps the raw top-ranked genes.

- enrichment_name:

  Optional stored enrichment result used to add cluster-level functional
  evidence.

- n_terms:

  Number of enrichment terms per cluster when `enrichment_name` is
  supplied.

- enrichment_selection:

  How to choose pathway/function terms for annotation evidence:
  `"specific"` prefers terms concentrated in fewer clusters, while
  `"top"` keeps the raw top-ranked terms.

- include_qc:

  Logical; whether to include cluster-level QC summaries in the evidence
  bundle.

- reduction:

  Optional dimensional reduction name used to summarize cluster
  neighborhood geometry, for example `"umap"`. Use `NULL` to disable
  geometry evidence.

- n_neighbor_clusters:

  Number of nearest clusters to report from the reduction centroid
  distances.

- background:

  Optional study-specific background information to provide additional
  interpretation context.

- annotation_mode:

  Annotation workflow mode. `"single_pass"` sends one compact prompt.
  `"agentic"` runs a broad-pass annotation followed by a focused
  refinement pass on ambiguous or ILC-relevant clusters, then reconciles
  the result against canonical lineage hints.

- output_format:

  One of `"llm"` for a model-ready prompt bundle or `"human"` for a
  human-readable summary.

- provider:

  Optional model provider function. When left `NULL`, Shennong will try
  to construct an ellmer-backed provider from `OPENAI_API_KEY` plus
  optional `OPENAI_BASE_URL` and `OPENAI_MODEL` environment variables.

- model:

  Optional model identifier.

- reasoning_effort:

  Optional reasoning effort forwarded to compatible GPT-5
  chat-completions providers, for example `"minimal"`, `"low"`,
  `"medium"`, `"high"`, or `"xhigh"` when the upstream endpoint supports
  it.

- include_json_schema:

  Logical; whether to request structured JSON output from the provider.
  Defaults to `TRUE` for annotation workflows.

- apply_metadata:

  Logical; if `TRUE` and a structured annotation response is returned,
  map the cluster labels back onto each cell in the Seurat metadata.

- metadata_prefix:

  Prefix used for metadata columns written back to the Seurat object
  when `apply_metadata = TRUE`.

- metadata_fields:

  Annotation fields to write back into Seurat metadata when
  `apply_metadata = TRUE`. Defaults to core fields only:
  `primary_label`, `broad_label`, `confidence`, `status`, and
  `risk_flags`. Detailed evidence remains in the stored interpretation
  table under `object@misc`.

- label_candidates:

  Optional vector of candidate cell-type labels or broad lineages that
  should be treated as annotation priors for sorted or enriched
  datasets. These are injected into the prompt as constraints, but the
  model may still return a broader label when evidence is weak.

- label_style:

  Naming style used to normalize returned cell-type labels. One of
  `"title"`, `"snake"`, or `"asis"`.

- return_prompt:

  If `TRUE`, return the prompt bundle without calling the provider.

- store_name:

  Name used under `object@misc$interpretation_results`.

- return_object:

  If `TRUE`, return the updated Seurat object.

- show_progress:

  Logical; if `TRUE`, emit step-wise progress logs and, when cli is
  available, a console progress bar while waiting for the LLM response.

- ...:

  Additional arguments forwarded to `provider`.

## Value

A prompt bundle, response, or updated `Seurat` object.

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
    store_name = "celltype_markers", return_object = TRUE, verbose = FALSE
  )
  prompt <- sn_interpret_annotation(
    obj,
    de_name = "celltype_markers",
    cluster_col = "cell_type",
    return_prompt = TRUE
  )
  prompt$task
}
#> INFO [2026-03-30 22:14:06] Initializing Seurat object for project: Shennong.
#> INFO [2026-03-30 22:14:06] Running QC metrics for human.
#> INFO [2026-03-30 22:14:06] Seurat object initialization complete.
#> INFO [2026-03-30 22:14:07] [sn_interpret_annotation] Starting interpretation workflow.
#> INFO [2026-03-30 22:14:07] [sn_interpret_annotation] Step 1/5: Preparing annotation evidence (elapsed 0.0s).
#> As of Seurat v5, we recommend using AggregateExpression to perform pseudo-bulk analysis.
#> This message is displayed once per session.
#> INFO [2026-03-30 22:14:07] [sn_interpret_annotation] Step 2/5: Building annotation prompt (elapsed 0.1s).
#> INFO [2026-03-30 22:14:07] [sn_interpret_annotation] Prompt prepared (total elapsed 0.1s).
#> [1] "annotation"
```
