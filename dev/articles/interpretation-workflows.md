# Interpretation and reporting

Shennong’s interpretation layer is intentionally downstream of analysis.
It does not replace markers, enrichment, QC, or clustering; it packages
those results into structured evidence, then builds prompts or calls a
provider only when the user asks.

This article uses a deterministic subset of the public Kotliarov PBMC
CITE-seq cohort as the evidence source. The local fixture retains real
donor, batch, response, RNA, and ADT information. Provider calls remain
opt-in because they require credentials and create external state.

## Create analysis evidence first

``` r

library(Shennong)
library(Seurat)
library(dplyr)

pbmc <- qs2::qs_read(pbmc_fixture)

tibble::tibble(
  source = "Kotliarov PBMC CITE-seq",
  cells = ncol(pbmc),
  genes = nrow(pbmc),
  samples = dplyr::n_distinct(pbmc$real_sample),
  batches = dplyr::n_distinct(pbmc$real_batch),
  assays = paste(SeuratObject::Assays(pbmc), collapse = ", ")
)

pbmc <- sn_run_cluster(
  object = pbmc,
  normalization_method = "seurat",
  nfeatures = 1500,
  dims = 1:15,
  resolution = 0.6,
  species = "human",
  verbose = FALSE
)

pbmc <- sn_find_de(
  object = pbmc,
  analysis = "markers",
  group_by = "seurat_clusters",
  layer = "data",
  result_id = "cluster_markers",
  return_object = TRUE,
  verbose = FALSE
)
```

The interpretation functions expect stored results. That keeps prompts
reproducible: a report can say which DE result, enrichment result, and
cluster column were used.

``` r

sn_list_results(pbmc)
```

## Prepare evidence before building prompts

Evidence helpers return ordinary R lists and tables. Inspect them before
using an LLM; this is the step that keeps interpretation transparent.

``` r

annotation_evidence <- sn_prepare_annotation_evidence(
  object = pbmc,
  de_result_id = "cluster_markers",
  cluster_by = "seurat_clusters",
  n_markers = 8,
  include_qc = TRUE,
  reduction = "umap"
)

de_evidence <- sn_prepare_de_evidence(
  object = pbmc,
  de_result_id = "cluster_markers",
  n_genes = 10
)

results_evidence <- sn_prepare_results_evidence(
  object = pbmc,
  cluster_de_result_id = "cluster_markers",
  cluster_by = "seurat_clusters",
  n_markers = 5
)

names(annotation_evidence)
head(annotation_evidence$cluster_summary)
```

If enrichment has been run, add pathway evidence the same way.

``` r

pbmc <- sn_run_enrichment(
  x = pbmc,
  source_de_result_id = "cluster_markers",
  species = "human",
  database = "GOBP",
  result_id = "cluster_gobp",
  return_object = TRUE
)

enrichment_evidence <- sn_prepare_enrichment_evidence(
  object = pbmc,
  enrichment_result_id = "cluster_gobp",
  n_terms = 8
)
```

## Build a prompt without calling a model

[`sn_build_prompt()`](https://zerostwo.github.io/shennong/dev/reference/sn_build_prompt.md)
is useful even when you do not plan to call a model from R. It gives you
a model-ready or human-readable representation of the exact evidence
bundle.

``` r

annotation_prompt <- sn_build_prompt(
  evidence = annotation_evidence,
  task = "annotation",
  audience = "scientist",
  language = "en",
  background = "Kotliarov PBMC CITE-seq response cohort clustered with Shennong.",
  output_format = "llm",
  include_json_schema = TRUE
)

names(annotation_prompt)
annotation_prompt$task
```

The high-level helper can return the same prompt directly.

``` r

annotation_prompt <- sn_interpret_annotation(
  object = pbmc,
  de_result_id = "cluster_markers",
  cluster_by = "seurat_clusters",
  return_prompt = TRUE,
  output_format = "llm",
  background = "Annotate Kotliarov PBMC clusters conservatively."
)

annotation_prompt$task
```

## Provider calls are opt-in

Shennong uses provider functions at the package boundary. You can create
an `ellmer` provider from environment variables, test it, then pass it
to the interpretation helpers.

``` r

provider <- sn_make_ellmer_provider(
  model = Sys.getenv("OPENAI_MODEL", "gpt-5.1"),
  reasoning_effort = "medium"
)

sn_test_llm_provider(provider = provider)

response <- sn_run_llm(
  messages = annotation_prompt$messages,
  provider = provider,
  model = Sys.getenv("OPENAI_MODEL", "gpt-5.1")
)
```

## Interpret specific evidence types

The DE and enrichment helpers use the same pattern: return a prompt by
default for review, or pass a provider to store a response.

``` r

de_prompt <- sn_interpret_de(
  object = pbmc,
  de_result_id = "cluster_markers",
  n_genes = 10,
  return_prompt = TRUE,
  background = "Summarize Kotliarov PBMC cluster marker programs."
)

de_prompt$task
```

``` r

enrichment_prompt <- sn_interpret_enrichment(
  object = pbmc,
  enrichment_result_id = "cluster_gobp",
  n_terms = 8,
  return_prompt = TRUE,
  background = "Summarize pathway programs by cluster."
)
```

## Write results, legends, and presentation summaries

The reporting helpers assemble multiple evidence sources into
audience-specific outputs. They can return prompts for manual review or
store generated responses in the canonical Shennong result registry.

``` r

results_prompt <- sn_write_results(
  object = pbmc,
  cluster_de_result_id = "cluster_markers",
  cluster_by = "seurat_clusters",
  return_prompt = TRUE,
  background = "Write a concise Results paragraph for Kotliarov PBMC clustering."
)

legend_prompt <- sn_write_figure_legend(
  object = pbmc,
  cluster_de_result_id = "cluster_markers",
  cluster_by = "seurat_clusters",
  return_prompt = TRUE,
  background = "Describe a UMAP and marker dot plot."
)

summary_prompt <- sn_write_presentation_summary(
  object = pbmc,
  cluster_de_result_id = "cluster_markers",
  cluster_by = "seurat_clusters",
  return_prompt = TRUE,
  background = "Prepare a short collaborator-facing Kotliarov PBMC summary."
)

c(results_prompt$task, legend_prompt$task, summary_prompt$task)
```

When a provider is used, retrieve stored outputs by name.

``` r

pbmc <- sn_interpret_de(
  object = pbmc,
  de_result_id = "cluster_markers",
  provider = provider,
  result_id = "cluster_marker_summary",
  return_object = TRUE
)

sn_get_interpretation_result(
  pbmc,
  result_id = "cluster_marker_summary"
)
```

The practical workflow is: compute evidence, inspect evidence, build a
prompt, then call a provider only when the prompt is good enough to
send.
