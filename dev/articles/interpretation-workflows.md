# Interpretation and reporting

Shennong’s interpretation layer is intentionally downstream of analysis.
It does not replace markers, enrichment, QC, or clustering; it packages
those results into structured evidence, then builds prompts or calls a
provider only when the user asks.

This article uses PBMC3k markers as the evidence source. Provider calls
are shown but not executed by default.

## Create analysis evidence first

``` r

library(Shennong)
library(Seurat)
library(dplyr)

pbmc <- sn_load_data("pbmc3k")
#> INFO [2026-05-05 21:29:48] Initializing Seurat object for project: pbmc3k.
#> INFO [2026-05-05 21:29:48] Running QC metrics for human.
#> INFO [2026-05-05 21:29:49] Seurat object initialization complete.

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
  store_name = "cluster_markers",
  return_object = TRUE,
  verbose = FALSE
)
```

The interpretation functions expect stored results. That keeps prompts
reproducible: a report can say which DE result, enrichment result, and
cluster column were used.

``` r

sn_list_results(pbmc)
#> # A tibble: 1 × 8
#>   collection type  name        analysis method created_at n_rows source
#>   <chr>      <chr> <chr>       <chr>    <chr>  <chr>       <int> <chr> 
#> 1 de_results de    cluster_ma… markers  wilcox 2026-05-0…   5049 NA
```

## Prepare evidence before building prompts

Evidence helpers return ordinary R lists and tables. Inspect them before
using an LLM; this is the step that keeps interpretation transparent.

``` r

annotation_evidence <- sn_prepare_annotation_evidence(
  object = pbmc,
  de_name = "cluster_markers",
  cluster_by = "seurat_clusters",
  n_markers = 8,
  include_qc = TRUE,
  reduction = "umap"
)

de_evidence <- sn_prepare_de_evidence(
  object = pbmc,
  de_name = "cluster_markers",
  n_genes = 10
)

results_evidence <- sn_prepare_results_evidence(
  object = pbmc,
  cluster_de_name = "cluster_markers",
  cluster_by = "seurat_clusters",
  n_markers = 5
)

names(annotation_evidence)
#>  [1] "task"                      "cluster_col"              
#>  [3] "source_de_name"            "source_enrichment_name"   
#>  [5] "analysis_method"           "species"                  
#>  [7] "marker_selection"          "enrichment_selection"     
#>  [9] "geometry_reduction"        "cluster_summary"          
#> [11] "top_marker_table"          "enrichment_summary"       
#> [13] "qc_summary"                "lineage_hints"            
#> [15] "canonical_marker_snapshot" "geometry_summary"         
#> [17] "caveats"
head(annotation_evidence$cluster_summary)
#> # A tibble: 6 × 15
#>   cluster n_cells fraction sample_distribution top_markers             
#>   <chr>     <int>    <dbl> <chr>               <chr>                   
#> 1 0           590   0.214  pbmc3k (590)        AQP3, CD40LG, LMNA, TRA…
#> 2 1           563   0.205  pbmc3k (563)        CCR7, APBA2, TRABD2A, T…
#> 3 2           479   0.174  pbmc3k (479)        FOLR3, S100A12, S100A8,…
#> 4 3           360   0.131  pbmc3k (360)        IGLC3, FCRLA, IGHD, TCL…
#> 5 4           345   0.125  pbmc3k (345)        GZMK, TRGC2, CD8A, KLRG…
#> 6 5           171   0.0621 pbmc3k (171)        CKB, CDKN1C, RDUR, MS4A…
#> # ℹ 10 more variables: median_nFeature_RNA <dbl>,
#> #   median_nCount_RNA <dbl>, median_percent.mt <dbl>,
#> #   median_percent.ribo <dbl>, median_percent.hb <dbl>,
#> #   doublet_fraction <dbl>, heuristic_hint <chr>,
#> #   heuristic_rationale <chr>, heuristic_top_signatures <chr>,
#> #   nearest_clusters <chr>
```

If enrichment has been run, add pathway evidence the same way.

``` r

pbmc <- sn_enrich(
  x = pbmc,
  source_de_name = "cluster_markers",
  species = "human",
  database = "GOBP",
  store_name = "cluster_gobp",
  return_object = TRUE
)

enrichment_evidence <- sn_prepare_enrichment_evidence(
  object = pbmc,
  enrichment_name = "cluster_gobp",
  n_terms = 8
)
```

## Build a prompt without calling a model

[`sn_build_prompt()`](https://songqi.org/shennong/dev/reference/sn_build_prompt.md)
is useful even when you do not plan to call a model from R. It gives you
a model-ready or human-readable representation of the exact evidence
bundle.

``` r

annotation_prompt <- sn_build_prompt(
  evidence = annotation_evidence,
  task = "annotation",
  audience = "scientist",
  language = "en",
  background = "PBMC3k single-cell RNA-seq clustered with Shennong.",
  output_format = "llm",
  include_json_schema = TRUE
)

names(annotation_prompt)
#> [1] "output_format" "task"          "system"        "user"         
#> [5] "messages"      "evidence"
annotation_prompt$task
#> [1] "annotation"
```

The high-level helper can return the same prompt directly.

``` r

annotation_prompt <- sn_interpret_annotation(
  object = pbmc,
  de_name = "cluster_markers",
  cluster_by = "seurat_clusters",
  return_prompt = TRUE,
  output_format = "llm",
  background = "Annotate PBMC3k clusters conservatively."
)
#> INFO [2026-05-05 21:30:53] [sn_interpret_annotation] Starting interpretation workflow.
#> INFO [2026-05-05 21:30:53] [sn_interpret_annotation] Step 1/5: Preparing annotation evidence (elapsed 0.0s).
#> INFO [2026-05-05 21:30:54] [sn_interpret_annotation] Step 2/5: Building annotation prompt (elapsed 0.3s).
#> INFO [2026-05-05 21:30:54] [sn_interpret_annotation] Prompt prepared (total elapsed 0.4s).

annotation_prompt$task
#> [1] "annotation"
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
  de_name = "cluster_markers",
  n_genes = 10,
  return_prompt = TRUE,
  background = "Summarize PBMC3k cluster marker programs."
)
#> INFO [2026-05-05 21:30:54] [sn_interpret_de] Starting interpretation workflow.
#> INFO [2026-05-05 21:30:54] [sn_interpret_de] Step 1/4: Preparing DE evidence (elapsed 0.0s).
#> INFO [2026-05-05 21:30:54] [sn_interpret_de] Step 2/4: Building interpretation prompt (elapsed 0.1s).
#> INFO [2026-05-05 21:30:54] [sn_interpret_de] Prompt prepared (total elapsed 0.1s).

de_prompt$task
#> [1] "de"
```

``` r

enrichment_prompt <- sn_interpret_enrichment(
  object = pbmc,
  enrichment_name = "cluster_gobp",
  n_terms = 8,
  return_prompt = TRUE,
  background = "Summarize pathway programs by cluster."
)
```

## Write results, legends, and presentation summaries

The reporting helpers assemble multiple evidence sources into
audience-specific outputs. They can return prompts for manual review or
store generated responses under `object@misc$interpretation_results`.

``` r

results_prompt <- sn_write_results(
  object = pbmc,
  cluster_de_name = "cluster_markers",
  cluster_by = "seurat_clusters",
  return_prompt = TRUE,
  background = "Write a concise Results paragraph for PBMC3k clustering."
)
#> INFO [2026-05-05 21:30:54] [sn_write_results] Starting interpretation workflow.
#> INFO [2026-05-05 21:30:54] [sn_write_results] Step 1/4: Preparing results evidence (elapsed 0.0s).
#> INFO [2026-05-05 21:30:55] [sn_write_results] Step 2/4: Building writing prompt (elapsed 0.3s).
#> INFO [2026-05-05 21:30:55] [sn_write_results] Prompt prepared (total elapsed 0.3s).

legend_prompt <- sn_write_figure_legend(
  object = pbmc,
  cluster_de_name = "cluster_markers",
  cluster_by = "seurat_clusters",
  return_prompt = TRUE,
  background = "Describe a UMAP and marker dot plot."
)
#> INFO [2026-05-05 21:30:55] [sn_write_figure_legend] Starting interpretation workflow.
#> INFO [2026-05-05 21:30:55] [sn_write_figure_legend] Step 1/4: Preparing legend evidence (elapsed 0.0s).
#> INFO [2026-05-05 21:30:55] [sn_write_figure_legend] Step 2/4: Building legend prompt (elapsed 0.3s).
#> INFO [2026-05-05 21:30:55] [sn_write_figure_legend] Prompt prepared (total elapsed 0.3s).

summary_prompt <- sn_write_presentation_summary(
  object = pbmc,
  cluster_de_name = "cluster_markers",
  cluster_by = "seurat_clusters",
  return_prompt = TRUE,
  background = "Prepare a short collaborator-facing PBMC3k summary."
)
#> INFO [2026-05-05 21:30:55] [sn_write_presentation_summary] Starting interpretation workflow.
#> INFO [2026-05-05 21:30:55] [sn_write_presentation_summary] Step 1/4: Preparing presentation evidence (elapsed 0.0s).
#> INFO [2026-05-05 21:30:55] [sn_write_presentation_summary] Step 2/4: Building presentation prompt (elapsed 0.3s).
#> INFO [2026-05-05 21:30:55] [sn_write_presentation_summary] Prompt prepared (total elapsed 0.3s).

c(results_prompt$task, legend_prompt$task, summary_prompt$task)
#> [1] "results"              "figure_legend"       
#> [3] "presentation_summary"
```

When a provider is used, retrieve stored outputs by name.

``` r

pbmc <- sn_interpret_de(
  object = pbmc,
  de_name = "cluster_markers",
  provider = provider,
  store_name = "cluster_marker_summary",
  return_object = TRUE
)

sn_get_interpretation_result(
  pbmc,
  interpretation_name = "cluster_marker_summary"
)
```

The practical workflow is: compute evidence, inspect evidence, build a
prompt, then call a provider only when the prompt is good enough to
send.
