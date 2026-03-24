# Interpretation workflow

This article shows how the new interpretation layer fits on top of the
current Shennong analysis workflow. The key idea is that analysis stays
deterministic, while interpretation reuses stored DE and enrichment
results to build evidence bundles, prompts, and optional LLM-style
summaries.

The workflow below demonstrates:

- marker discovery with
  [`sn_find_de()`](https://songqi.org/shennong/reference/sn_find_de.md)
- enrichment storage with
  [`sn_enrich()`](https://songqi.org/shennong/reference/sn_enrich.md) /
  [`sn_store_enrichment()`](https://songqi.org/shennong/reference/sn_store_enrichment.md)
- result discovery with
  [`sn_list_results()`](https://songqi.org/shennong/reference/sn_list_results.md)
- direct retrieval with
  [`sn_get_de_result()`](https://songqi.org/shennong/reference/sn_get_de_result.md)
  and
  [`sn_get_enrichment_result()`](https://songqi.org/shennong/reference/sn_get_enrichment_result.md)
- evidence preparation for annotation, DE, and writing tasks
- prompt generation with
  [`sn_build_prompt()`](https://songqi.org/shennong/reference/sn_build_prompt.md)
- provider-based interpretation with
  [`sn_interpret_annotation()`](https://songqi.org/shennong/reference/sn_interpret_annotation.md)

``` r
library(Shennong)
library(dplyr)
library(knitr)
library(Seurat)
```

## Marker discovery remains the analysis layer

Interpretation starts from stored DE results rather than recomputing
statistics. Here,
[`sn_find_de()`](https://songqi.org/shennong/reference/sn_find_de.md)
has already saved cluster markers into the Seurat object, and we
retrieve them through
[`sn_get_de_result()`](https://songqi.org/shennong/reference/sn_get_de_result.md).

``` r
knitr::kable(top_cluster_markers |> dplyr::select(cluster, gene, avg_log2FC, p_val_adj))
```

| cluster | gene            | avg_log2FC | p_val_adj |
|:--------|:----------------|-----------:|----------:|
| 0       | BPI             |   5.812083 |         0 |
| 0       | ENSG00000289381 |   5.667320 |         0 |
| 0       | MTARC1          |   5.624992 |         0 |

## Store enrichment for later interpretation

The interpretation layer expects enrichment results to be stored in the
Seurat object.
[`sn_enrich()`](https://songqi.org/shennong/reference/sn_enrich.md) can
now store them directly when you pass `object =`.

``` r
knitr::kable(top_pathways)
```

| ID           | Description                                            |       NES |  p.adjust |
|:-------------|:-------------------------------------------------------|----------:|----------:|
| <GO:0042776> | proton motive force-driven mitochondrial ATP synthesis | -2.515028 | 0.0002216 |
| <GO:0006120> | mitochondrial electron transport, NADH to ubiquinone   | -2.099828 | 0.0427655 |
| <GO:0032543> | mitochondrial translation                              | -1.960416 | 0.0651602 |
| <GO:0098581> | detection of external biotic stimulus                  |  1.921567 | 0.0015675 |
| <GO:0009595> | detection of biotic stimulus                           |  1.912715 | 0.0031257 |

## Discover what is currently stored on the object

Use
[`sn_list_results()`](https://songqi.org/shennong/reference/sn_list_results.md)
to see which DE, enrichment, and interpretation artifacts are available
for reuse.

``` r
knitr::kable(result_index)
```

| collection             | type           | name                    | analysis | method | created_at              | n_rows | source          |
|:-----------------------|:---------------|:------------------------|:---------|:-------|:------------------------|-------:|:----------------|
| de_results             | de             | cluster_markers         | markers  | wilcox | 2026-03-24 17:53:20 UTC |  17203 | NA              |
| enrichment_results     | enrichment     | cluster0_gsea           | gsea     | NA     | 2026-03-24 17:53:54 UTC |   2646 | cluster_markers |
| interpretation_results | interpretation | cluster_annotation_note | NA       | NA     | 2026-03-24 17:53:55 UTC |      0 | NA              |

## Prepare structured evidence

`sn_prepare_*_evidence()` turns stored outputs into compact, reusable
evidence bundles. These are the objects that can later be turned into
prompts or passed into writing helpers.

``` r
knitr::kable(annotation_evidence$cluster_summary)
```

| cluster | n_cells |  fraction | top_markers                                                               |
|:--------|--------:|----------:|:--------------------------------------------------------------------------|
| 0       |     229 | 0.1884774 | BPI, ENSG00000289381, MTARC1, MCEMP1, VNN3P                               |
| 1       |     172 | 0.1415638 | IATPR, TTC39C-AS1, PI16, TNFRSF4, CRIP2                                   |
| 2       |     149 | 0.1226337 | EDAR, ANKRD55, ADTRP, TSHZ2, ENSG00000271774                              |
| 3       |     144 | 0.1185185 | TCL1A, ENSG00000257275, SCN3A, ENSG00000224610, KCNG1                     |
| 4       |     133 | 0.1094650 | LINC02446, CD8B, ENSG00000310107, CRTAM, GZMH                             |
| 5       |     125 | 0.1028807 | CLEC10A, FCER1A, CYP2S1, DTNA, ZBTB46                                     |
| 6       |      75 | 0.0617284 | SLC4A10, ENSG00000228033, IL23R, ADAM12, LINC01644                        |
| 7       |      61 | 0.0502058 | IGHA1, IGHG2, IGHG3, IGHG1, SSPN                                          |
| 8       |      54 | 0.0444444 | ENSG00000288970, ENSG00000291157, MYOM2, ENSG00000288782, ENSG00000294329 |
| 9       |      32 | 0.0263374 | MT-CYB, MT-CO3, MT-ATP6, MT-CO2, MT-ND3                                   |
| 10      |      26 | 0.0213992 | LYPD2, ENSG00000301038, UICLM, PPP1R17, LINC02345                         |
| 11      |      15 | 0.0123457 | CLDN5, PF4V1, CMTM5, PDZK1IP1, PDGFA-DT                                   |

``` r
knitr::kable(enrichment_evidence$top_terms)
```

| ID           | Description                                            |       NES |  p.adjust |
|:-------------|:-------------------------------------------------------|----------:|----------:|
| <GO:0042776> | proton motive force-driven mitochondrial ATP synthesis | -2.515028 | 0.0002216 |
| <GO:0006120> | mitochondrial electron transport, NADH to ubiquinone   | -2.099828 | 0.0427655 |
| <GO:0032543> | mitochondrial translation                              | -1.960416 | 0.0651602 |
| <GO:0098581> | detection of external biotic stimulus                  |  1.921567 | 0.0015675 |
| <GO:0009595> | detection of biotic stimulus                           |  1.912715 | 0.0031257 |

## Build prompts without binding to a specific provider

The first-class package object for the LLM layer is a prompt bundle, not
a network call. This keeps analysis and interpretation cleanly
separated.

``` r
cat(substr(annotation_prompt$user, 1, 900))
#> Task: annotation.
#> 
#> Audience: scientist.
#> 
#> Language: en.
#> 
#> Target style: concise annotation note.
#> 
#> Task instructions: Interpret cluster-level marker evidence to support cell type annotation. For each cluster, identify the most plausible cell type or state. Explicitly cite the top markers that support the label and mention conflicting markers when relevant. Flag ambiguous clusters and suggest additional markers or orthogonal checks that would improve confidence. Return one concise cluster-by-cluster annotation table followed by a short narrative summary.
#> 
#> 
#> 
#> 
#> 
#> Evidence:
#> 
#> task:
#> annotation
#> 
#> cluster_col:
#> seurat_clusters
#> 
#> source_de_name:
#> cluster_markers
#> 
#> analysis_method:
#> wilcox
#> 
#> species:
#> human
#> 
#> cluster_summary:
#> # A tibble: 8 × 4
#>   cluster n_cells fraction top_markers                                          
#>   <fct>     <int>    
```

## Generate manuscript-style prompts

High-level helpers such as
[`sn_write_results()`](https://songqi.org/shennong/reference/sn_write_results.md)
can return prompts directly.

## Initialize a Codex-ready analysis project

For longer-running projects, it helps to initialize a governed analysis
repository rather than improvising the folder structure. The packaged
project template used by
[`sn_initialize_codex_project()`](https://songqi.org/shennong/reference/sn_initialize_codex_project.md)
creates `AGENTS.md`, `memory/`, `docs/standards/`, `skills/`, `config/`,
`data/`, `scripts/`, `notebooks/`, `runs/`, and `results/`.

``` r
codex_project <- file.path(tempdir(), "shennong-codex-project")

created <- sn_initialize_codex_project(
  path = codex_project,
  project_name = "PBMC interpretation pilot",
  objective = "Build a reproducible PBMC interpretation workflow with Shennong.",
  overwrite = TRUE
)

knitr::kable(
  tibble::tibble(
    file = names(created)[-1],
    path = unname(unlist(created[-1]))
  )
)
```

| file             | path                                                                                       |
|:-----------------|:-------------------------------------------------------------------------------------------|
| readme           | /tmp/RtmpmxFJ1X/shennong-codex-project/README.md                                           |
| config           | /tmp/RtmpmxFJ1X/shennong-codex-project/config                                              |
| config_default   | /tmp/RtmpmxFJ1X/shennong-codex-project/config/default.yaml                                 |
| data             | /tmp/RtmpmxFJ1X/shennong-codex-project/data                                                |
| data_raw         | /tmp/RtmpmxFJ1X/shennong-codex-project/data/raw                                            |
| data_processed   | /tmp/RtmpmxFJ1X/shennong-codex-project/data/processed                                      |
| data_metadata    | /tmp/RtmpmxFJ1X/shennong-codex-project/data/metadata                                       |
| scripts          | /tmp/RtmpmxFJ1X/shennong-codex-project/scripts                                             |
| notebooks        | /tmp/RtmpmxFJ1X/shennong-codex-project/notebooks                                           |
| runs             | /tmp/RtmpmxFJ1X/shennong-codex-project/runs                                                |
| results          | /tmp/RtmpmxFJ1X/shennong-codex-project/results                                             |
| results_figures  | /tmp/RtmpmxFJ1X/shennong-codex-project/results/figures                                     |
| results_tables   | /tmp/RtmpmxFJ1X/shennong-codex-project/results/tables                                      |
| results_reports  | /tmp/RtmpmxFJ1X/shennong-codex-project/results/reports                                     |
| agents_md        | /tmp/RtmpmxFJ1X/shennong-codex-project/AGENTS.md                                           |
| agents           | /tmp/RtmpmxFJ1X/shennong-codex-project/AGENTS.md                                           |
| memory           | /tmp/RtmpmxFJ1X/shennong-codex-project/memory                                              |
| memory_decisions | /tmp/RtmpmxFJ1X/shennong-codex-project/memory/Decisions.md                                 |
| decisions        | /tmp/RtmpmxFJ1X/shennong-codex-project/memory/Decisions.md                                 |
| memory_plan      | /tmp/RtmpmxFJ1X/shennong-codex-project/memory/Plan.md                                      |
| plan             | /tmp/RtmpmxFJ1X/shennong-codex-project/memory/Plan.md                                      |
| memory_prompt    | /tmp/RtmpmxFJ1X/shennong-codex-project/memory/Prompt.md                                    |
| prompt           | /tmp/RtmpmxFJ1X/shennong-codex-project/memory/Prompt.md                                    |
| memory_status    | /tmp/RtmpmxFJ1X/shennong-codex-project/memory/Status.md                                    |
| status           | /tmp/RtmpmxFJ1X/shennong-codex-project/memory/Status.md                                    |
| standards        | /tmp/RtmpmxFJ1X/shennong-codex-project/docs/standards                                      |
| conventions      | /tmp/RtmpmxFJ1X/shennong-codex-project/docs/standards/BioinformaticsAnalysisConventions.md |
| skills           | /tmp/RtmpmxFJ1X/shennong-codex-project/skills                                              |

The initialized project keeps durable operating rules in `AGENTS.md`,
project state in `memory/`, enforceable directory and naming rules in
`docs/standards/BioinformaticsAnalysisConventions.md`, and reusable
governance procedures in `skills/`.

``` r
cat(substr(results_prompt$user, 1, 900))
#> Task: results.
#> 
#> Audience: scientist.
#> 
#> Language: en.
#> 
#> Target style: manuscript-style Results section.
#> 
#> Task instructions: Write a manuscript-style Results subsection using only the supplied evidence. Keep the tone formal, precise, and evidence-based. Integrate cluster, DE, and enrichment findings into a coherent paragraph sequence instead of bullet fragments. Do not claim validation beyond the provided evidence.
#> 
#> 
#> 
#> 
#> 
#> Evidence:
#> 
#> task:
#> results
#> 
#> dataset:
#> n_cells:
#> 1215
#> 
#> n_features:
#> 54872
#> 
#> cluster_col:
#> seurat_clusters
#> 
#> clusters:
#> 12
#> 
#> cluster_summary:
#> # A tibble: 8 × 3
#>   cluster n_cells fraction
#>   <fct>     <int>    <dbl>
#> 1 0           229   0.188 
#> 2 1           172   0.142 
#> 3 2           149   0.123 
#> 4 3           144   0.119 
#> 5
```

## Plug in a provider when needed

If you provide a function that accepts `messages` and returns text, the
same high-level helpers can store the generated interpretation back into
the Seurat object.

``` r
annotation_response
#> [1] "Mock interpretation generated from 2 messages using demo-model"
```

## Interpretation results are stored alongside analysis results

``` r
names(pbmc_interpreted@misc$interpretation_results)
#> [1] "cluster_annotation_note"
```
