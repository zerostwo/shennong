# Interpretation and reporting workflow

This article shows how the new interpretation layer fits on top of the
current Shennong analysis workflow. The key idea is that analysis stays
deterministic, while interpretation reuses stored DE and enrichment
results to build evidence bundles, prompts, and optional LLM-style
summaries.

The workflow below demonstrates:

- marker discovery with
  [`sn_find_de()`](https://songqi.org/shennong/dev/reference/sn_find_de.md)
- enrichment storage with
  [`sn_enrich()`](https://songqi.org/shennong/dev/reference/sn_enrich.md)
  /
  [`sn_store_enrichment()`](https://songqi.org/shennong/dev/reference/sn_store_enrichment.md)
- result discovery with
  [`sn_list_results()`](https://songqi.org/shennong/dev/reference/sn_list_results.md)
- direct retrieval with
  [`sn_get_de_result()`](https://songqi.org/shennong/dev/reference/sn_get_de_result.md)
  and
  [`sn_get_enrichment_result()`](https://songqi.org/shennong/dev/reference/sn_get_enrichment_result.md)
- evidence preparation for annotation, DE, and writing tasks
- prompt generation with
  [`sn_build_prompt()`](https://songqi.org/shennong/dev/reference/sn_build_prompt.md)
- provider-based interpretation with
  [`sn_interpret_annotation()`](https://songqi.org/shennong/dev/reference/sn_interpret_annotation.md)

When a Seurat object already contains stored marker results,
[`sn_interpret_annotation()`](https://songqi.org/shennong/dev/reference/sn_interpret_annotation.md)
can now resolve `de_name` automatically. It prefers a stored `default`
result, then a single available result, and otherwise the most recent
marker result. Passing `de_name` explicitly is still recommended when
your object carries multiple marker tables that represent different
biological questions. For annotation tasks, Shennong now includes the
full cluster evidence table in the model-facing prompt instead of
truncating the prompt to the first eight rows, so larger clustering
results are annotated more completely. It can also prefer
cluster-specific markers/pathways over generic top-ranked features,
summarize nearest-cluster geometry from reductions such as UMAP, and
explicitly pass `reasoning_effort` through compatible GPT-5
chat-completions providers. When the provider is built with
[`sn_make_ellmer_provider()`](https://songqi.org/shennong/dev/reference/sn_make_ellmer_provider.md),
annotation now uses `ellmer::chat_structured()` for the final label
table and can insert a tool-assisted focused-comparison step in
`annotation_mode = "agentic"` before the final refinement pass.

``` r
library(Shennong)
library(dplyr)
library(knitr)
library(Seurat)
```

## Marker discovery remains the analysis layer

Interpretation starts from stored DE results rather than recomputing
statistics. Here,
[`sn_find_de()`](https://songqi.org/shennong/dev/reference/sn_find_de.md)
has already saved cluster markers into the Seurat object, and we
retrieve them through
[`sn_get_de_result()`](https://songqi.org/shennong/dev/reference/sn_get_de_result.md).

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
[`sn_enrich()`](https://songqi.org/shennong/dev/reference/sn_enrich.md)
can now store them directly when you pass a Seurat object through `x` or
`object =`, and its `gene_clusters` formula now covers both grouped ORA
inputs such as `gene ~ cluster` and ranked GSEA inputs such as

When
[`sn_interpret_annotation()`](https://songqi.org/shennong/dev/reference/sn_interpret_annotation.md)
receives a structured JSON response, it also stores the parsed
cluster-level annotation table and can map the labels back onto each
cell’s metadata. In the example above, that means columns such as
`shennong_celltype_label`, `shennong_celltype_confidence`, and
`shennong_celltype_risk_flags` are available directly in
`pbmc_interpreted[[]]`. `gene ~ avg_log2FC`.

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
[`sn_list_results()`](https://songqi.org/shennong/dev/reference/sn_list_results.md)
to see which DE, enrichment, and interpretation artifacts are available
for reuse.

``` r
knitr::kable(result_index)
```

| collection             | type           | name                    | analysis | method | created_at              | n_rows | source          |
|:-----------------------|:---------------|:------------------------|:---------|:-------|:------------------------|-------:|:----------------|
| de_results             | de             | cluster_markers         | markers  | wilcox | 2026-03-30 08:52:12 UTC |  17203 | NA              |
| enrichment_results     | enrichment     | cluster0_gsea           | gsea     | NA     | 2026-03-30 08:52:45 UTC |   2646 | cluster_markers |
| interpretation_results | interpretation | cluster_annotation_note | NA       | NA     | 2026-03-30 08:52:47 UTC |      0 | NA              |

## Prepare structured evidence

`sn_prepare_*_evidence()` turns stored outputs into compact, reusable
evidence bundles. These are the objects that can later be turned into
prompts or passed into writing helpers. For annotation workflows, this
now includes a `canonical_marker_snapshot` table that can be used to
compare lineage-defining markers across clusters during focused
refinement.

``` r
knitr::kable(annotation_evidence$cluster_summary)
```

| cluster | n_cells |  fraction | top_markers                                                               | median_nFeature_RNA | median_nCount_RNA | median_percent.mt | median_percent.ribo | median_percent.hb | doublet_fraction | heuristic_hint               | heuristic_rationale                                                                                                                        | heuristic_top_signatures                               | nearest_clusters                                           |
|:--------|--------:|----------:|:--------------------------------------------------------------------------|--------------------:|------------------:|------------------:|--------------------:|------------------:|-----------------:|:-----------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------|:-------------------------------------------------------|:-----------------------------------------------------------|
| 0       |     229 | 0.1884774 | BPI, ENSG00000289381, MTARC1, MCEMP1, VNN3P                               |              3957.0 |           13263.0 |          6.889921 |          10.4243822 |         0.0000000 |              NaN | B-cell-like                  | Canonical B-cell markers dominate (BPI, ENSG00000289381, MTARC1, MCEMP1, VNN3P)                                                            | b_cell=4.18; dendritic_apc=3.32; nk_cytotoxic=1.23     | 5 (distance=4.11); 11 (distance=4.44); 10 (distance=7.14)  |
| 1       |     172 | 0.1415638 | IATPR, TTC39C-AS1, PI16, TNFRSF4, CRIP2                                   |              3473.0 |           10690.0 |          5.934298 |          25.0039520 |         0.0000000 |              NaN | T-cell-like                  | Canonical T-cell markers dominate (IATPR, TTC39C-AS1, PI16, TNFRSF4, CRIP2)                                                                | t_cell=4.36; kit_ilc=2.78; ilc2=2.52                   | 2 (distance=3.07); 6 (distance=4.04); 9 (distance=5.46)    |
| 2       |     149 | 0.1226337 | EDAR, ANKRD55, ADTRP, MAL-AS1, SARDH                                      |              3120.0 |            9485.0 |          6.028954 |          28.9157990 |         0.0000000 |              NaN | T-cell-like                  | Canonical T-cell markers dominate (EDAR, ANKRD55, ADTRP, MAL-AS1, SARDH)                                                                   | t_cell=3.46; kit_ilc=1.26; ilc2=1.12                   | 1 (distance=3.07); 9 (distance=4.17); 6 (distance=6.43)    |
| 3       |     144 | 0.1185185 | TCL1A, ENSG00000257275, SCN3A, ENSG00000224610, KCNG1                     |              2977.0 |            8427.5 |          7.027844 |          19.8495396 |         0.0000000 |              NaN | B-cell-like                  | Canonical B-cell markers dominate (TCL1A, ENSG00000257275, SCN3A, ENSG00000224610, KCNG1)                                                  | b_cell=21.26; dendritic_apc=4.52; t_cell=1.05          | 7 (distance=1.88); 11 (distance=14.72); 1 (distance=14.76) |
| 4       |     133 | 0.1094650 | LINC02446, CD8B, ENSG00000310107, CRTAM, ENSG00000303688                  |              3169.0 |            8875.0 |          7.271616 |          19.3367585 |         0.0000000 |              NaN | T/NK mixed lymphoid-like     | Both T-cell receptor and cytotoxic/NK programs are substantial (LINC02446, CD8B, ENSG00000310107, CRTAM, ENSG00000303688)                  | t_cell=3.11; nk_cytotoxic=2.74; kit_ilc=1.72           | 6 (distance=3.83); 1 (distance=6.56); 2 (distance=9.52)    |
| 5       |     125 | 0.1028807 | CLEC10A, FCER1A, CYP2S1, DTNA, ZBTB46                                     |              5240.0 |           21033.0 |          6.195413 |          12.0623713 |         0.0035853 |              NaN | B-cell-like                  | Canonical B-cell markers dominate (CLEC10A, FCER1A, CYP2S1, DTNA, ZBTB46)                                                                  | b_cell=16.58; dendritic_apc=8.62; nk_cytotoxic=1.18    | 10 (distance=3.07); 0 (distance=4.11); 11 (distance=6.86)  |
| 6       |      75 | 0.0617284 | SLC4A10, ENSG00000228033, IL23R, ADAM12, LINC01644                        |              3248.0 |            9372.0 |          8.657414 |          21.5659628 |         0.0000000 |              NaN | T/NK mixed lymphoid-like     | Both T-cell receptor and cytotoxic/NK programs are substantial (SLC4A10, ENSG00000228033, IL23R, ADAM12, LINC01644)                        | kit_ilc=5.29; ilc2=4.77; t_cell=3.50                   | 4 (distance=3.83); 1 (distance=4.04); 9 (distance=6.16)    |
| 7       |      61 | 0.0502058 | IGHA1, IGHG2, IGHG3, IGHG1, SSPN                                          |              3745.0 |           12554.0 |          6.794110 |          23.4140484 |         0.0000000 |              NaN | B-cell-like                  | Canonical B-cell markers dominate (IGHA1, IGHG2, IGHG3, IGHG1, SSPN)                                                                       | b_cell=18.57; dendritic_apc=4.22; t_cell=1.23          | 3 (distance=1.88); 11 (distance=15.42); 1 (distance=16.64) |
| 8       |      54 | 0.0444444 | ENSG00000288970, ENSG00000291157, MYOM2, ENSG00000288782, ENSG00000294329 |              3441.0 |            8378.0 |          6.311986 |          10.1437599 |         0.0000000 |              NaN | T/NK mixed lymphoid-like     | Both T-cell receptor and cytotoxic/NK programs are substantial (ENSG00000288970, ENSG00000291157, MYOM2, ENSG00000288782, ENSG00000294329) | nk_cytotoxic=20.13; kit_ilc=1.91; ilc2=1.73            | 6 (distance=9.19); 4 (distance=10.47); 9 (distance=10.91)  |
| 9       |      32 | 0.0263374 | MT-CYB, MT-CO3, MT-ATP6, MT-CO2, MT-ND3                                   |                42.5 |             872.0 |         96.021924 |           0.4498247 |         0.0000000 |              NaN | Erythroid contamination-like | Hemoglobin/erythroid markers dominate (MT-CYB, MT-CO3, MT-ATP6, MT-CO2, MT-ND3)                                                            | erythroid_contam=27.75; ilc2=0.46; kit_ilc=0.44        | 2 (distance=4.17); 1 (distance=5.46); 6 (distance=6.16)    |
| 10      |      26 | 0.0213992 | LYPD2, ENSG00000301038, UICLM, PPP1R17, LINC02345                         |              4959.5 |           19769.0 |          6.853679 |           9.9982262 |         0.0000000 |              NaN | ILC2-like                    | Canonical ILC2 program is strongest (LYPD2, ENSG00000301038, UICLM, PPP1R17, LINC02345)                                                    | b_cell=7.62; dendritic_apc=4.16; nk_cytotoxic=3.25     | 5 (distance=3.07); 0 (distance=7.14); 11 (distance=9.12)   |
| 11      |      15 | 0.0123457 | CLDN5, PF4V1, CMTM5, PDZK1IP1, PDGFA-DT                                   |               473.0 |            1011.0 |          8.003680 |           1.0396975 |         0.0636335 |              NaN | Erythroid contamination-like | Hemoglobin/erythroid markers dominate (CLDN5, PF4V1, CMTM5, PDZK1IP1, PDGFA-DT)                                                            | erythroid_contam=1.21; dendritic_apc=1.05; t_cell=0.20 | 0 (distance=4.44); 5 (distance=6.86); 10 (distance=9.12)   |

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
#> # Interpretation Request
#> 
#> ## Task Metadata
#> 
#> - Task: annotation
#> - Audience: scientist
#> - Language: en
#> - Target style: concise annotation note
#> 
#> ## Task Instructions
#> 
#> Interpret cluster-level evidence to support cell type annotation. For each cluster, integrate marker genes, pathway/function evidence, canonical lineage heuristic hints, reference predictions, QC-related caveats, and cluster-neighborhood geometry when available. You must return one annotation record for every cluster present in the evidence and must not omit clusters. For each cluster, identify the most plausible cell type or state, plus any strong alternative labels. Compare clusters against each other when that helps distinguish related states or subtypes; do not reason about each cluster in isolation. Explicitly cite the top markers and functional terms that support the label and mention conflicting evidence when relevant. P
```

## Generate manuscript-style prompts

High-level helpers such as
[`sn_write_results()`](https://songqi.org/shennong/dev/reference/sn_write_results.md)
can return prompts directly.

## Initialize a Codex-ready analysis project

For longer-running projects, it helps to initialize a governed analysis
repository rather than improvising the folder structure. The packaged
project template used by `sn_initialize_project(with_agent = TRUE)`
creates `AGENTS.md`, `.gitignore`, a project `.Rproj`, `memory/`,
`docs/standards/`, `skills/`, `config/`, `data/`, `scripts/`,
`notebooks/`, `runs/`, and `results/`.

``` r
codex_project <- file.path(tempdir(), "shennong-codex-project")

created <- sn_initialize_project(
  path = codex_project,
  project_name = "PBMC interpretation pilot",
  objective = "Build a reproducible PBMC interpretation workflow with Shennong.",
  with_agent = TRUE,
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
| readme           | /tmp/RtmpNogFQJ/shennong-codex-project/README.md                                           |
| gitignore        | /tmp/RtmpNogFQJ/shennong-codex-project/.gitignore                                          |
| rproj            | /tmp/RtmpNogFQJ/shennong-codex-project/shennong-codex-project.Rproj                        |
| config           | /tmp/RtmpNogFQJ/shennong-codex-project/config                                              |
| config_default   | /tmp/RtmpNogFQJ/shennong-codex-project/config/default.yaml                                 |
| data             | /tmp/RtmpNogFQJ/shennong-codex-project/data                                                |
| data_raw         | /tmp/RtmpNogFQJ/shennong-codex-project/data/raw                                            |
| data_processed   | /tmp/RtmpNogFQJ/shennong-codex-project/data/processed                                      |
| data_metadata    | /tmp/RtmpNogFQJ/shennong-codex-project/data/metadata                                       |
| scripts          | /tmp/RtmpNogFQJ/shennong-codex-project/scripts                                             |
| notebooks        | /tmp/RtmpNogFQJ/shennong-codex-project/notebooks                                           |
| runs             | /tmp/RtmpNogFQJ/shennong-codex-project/runs                                                |
| results          | /tmp/RtmpNogFQJ/shennong-codex-project/results                                             |
| results_figures  | /tmp/RtmpNogFQJ/shennong-codex-project/results/figures                                     |
| results_tables   | /tmp/RtmpNogFQJ/shennong-codex-project/results/tables                                      |
| results_reports  | /tmp/RtmpNogFQJ/shennong-codex-project/results/reports                                     |
| agents_md        | /tmp/RtmpNogFQJ/shennong-codex-project/AGENTS.md                                           |
| agents           | /tmp/RtmpNogFQJ/shennong-codex-project/AGENTS.md                                           |
| memory           | /tmp/RtmpNogFQJ/shennong-codex-project/memory                                              |
| memory_decisions | /tmp/RtmpNogFQJ/shennong-codex-project/memory/Decisions.md                                 |
| decisions        | /tmp/RtmpNogFQJ/shennong-codex-project/memory/Decisions.md                                 |
| memory_plan      | /tmp/RtmpNogFQJ/shennong-codex-project/memory/Plan.md                                      |
| plan             | /tmp/RtmpNogFQJ/shennong-codex-project/memory/Plan.md                                      |
| memory_prompt    | /tmp/RtmpNogFQJ/shennong-codex-project/memory/Prompt.md                                    |
| prompt           | /tmp/RtmpNogFQJ/shennong-codex-project/memory/Prompt.md                                    |
| memory_status    | /tmp/RtmpNogFQJ/shennong-codex-project/memory/Status.md                                    |
| status           | /tmp/RtmpNogFQJ/shennong-codex-project/memory/Status.md                                    |
| standards        | /tmp/RtmpNogFQJ/shennong-codex-project/docs/standards                                      |
| conventions      | /tmp/RtmpNogFQJ/shennong-codex-project/docs/standards/BioinformaticsAnalysisConventions.md |
| skills           | /tmp/RtmpNogFQJ/shennong-codex-project/skills                                              |

The initialized project keeps durable operating rules in `AGENTS.md`,
project state in `memory/`, enforceable directory and naming rules in
`docs/standards/BioinformaticsAnalysisConventions.md`, and reusable
governance procedures in `skills/`.

``` r
cat(substr(results_prompt$user, 1, 900))
```

    #> # Interpretation Request
    #> 
    #> ## Task Metadata
    #> 
    #> - Task: results
    #> - Audience: scientist
    #> - Language: en
    #> - Target style: manuscript-style Results section
    #> 
    #> ## Task Instructions
    #> 
    #> Write a manuscript-style Results subsection using only the supplied evidence. Keep the tone formal, precise, and evidence-based. Integrate cluster, DE, and enrichment findings into a coherent paragraph sequence instead of bullet fragments. Do not claim validation beyond the provided evidence.
    #> 
    #> ## Evidence
    #> 
    #> ### task
    #> 
    #> ```text
    #> results
    #> ```
    #> 
    #> ### dataset
    #> 
    #> ### n_cells
    #> 
    #> ```text
    #> 1215
    #> ```
    #> 
    #> ### n_features
    #> 
    #> ```text
    #> 54872
    #> ```
    #> 
    #> ### cluster_col
    #> 
    #> ```text
    #> seurat_clusters
    #> ```
    #> 
    #> ### clusters
    #> 
    #> ```text
    #> 12
    #> ```
    #> 
    #> ### cluster_summary
    #> 
    #> 
    #> 
    #> |cluster | n_cells|  fraction|
    #> |:-------|-------:|---------:|
    #> |0       |     229| 0.1884774|
    #> |1       |     172| 0.1415638|
    #> |2       |     149| 0.1226337|
    #> |3       |     144| 0.1185185|
    #> |4       |     133| 0.109465

## Plug in a provider when needed

If you provide a function that accepts `messages` and returns text, the
same high-level helpers can store the generated interpretation back into
the Seurat object.

For live model calls, the recommended path is now `ellmer`. You can
either pass an explicit provider returned by
[`sn_make_ellmer_provider()`](https://songqi.org/shennong/dev/reference/sn_make_ellmer_provider.md)
or rely on environment variables such as `OPENAI_API_KEY`,
`OPENAI_BASE_URL`, and `OPENAI_MODEL` so the high-level helpers can
construct an `ellmer` provider directly.

``` r
annotation_response
#> [1] "{\"cluster_annotations\":[{\"cluster\":\"0\",\"primary_label\":\"T Cell\",\"broad_label\":\"Lymphocyte\",\"confidence\":\"high\",\"status\":\"confident\",\"alternatives\":\"Activated T Cell\",\"supporting_markers\":[\"CD3D\",\"TRAC\"],\"supporting_functions\":\"T cell activation\",\"risk_flags\":\"transitional_state\",\"note\":\"Canonical T-cell markers dominate this cluster.\",\"recommended_checks\":\"Inspect IL7R and LTB\"}],\"narrative_summary\":\"Mock structured annotation generated for the vignette.\"}"
```

## Interpretation results are stored alongside analysis results

``` r
names(pbmc_interpreted@misc$interpretation_results)
#> [1] "cluster_annotation_note"
```
